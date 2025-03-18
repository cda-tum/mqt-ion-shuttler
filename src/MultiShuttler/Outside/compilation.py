import matplotlib.pyplot as plt
import networkx as nx
from qiskit.dagcircuit import DAGDependency
from qiskit import QuantumCircuit

from qiskit.converters import circuit_to_dagdependency
from qiskit.transpiler.passes import RemoveBarriers, RemoveFinalMeasurements

import math
import re
from pathlib import Path
from scheduling import pick_pz_for_2_q_gate


from graph_utils import GraphCreator, PZCreator, ProcessingZone, create_idc_dictionary
from Cycles import create_starting_config, find_path_edge_to_edge
from scheduling import get_ions
import math
import networkx as nx
import numpy as np
from datetime import datetime
from plotting import plot_state
from graph_utils import get_idx_from_idc



def is_qasm_file(filename):
    # Check if the file has a .qasm extension
    if not filename.endswith(".qasm"):
        return False

    
    file_path = Path(filename)
    with Path.open(file_path) as file:
        # Read the first line of the file (7th line, specific to MQT Bench)
        for _f in range(7):
            first_line = file.readline()
        # Check if the first line contains the OPENQASM identifier
        return "OPENQASM" in first_line



def extract_qubits_from_gate(gate_line):
    """Extract qubit numbers from a gate operation line."""
    # Regular expression to match qubits (assuming they are in the format q[<number>])
    pattern = re.compile(r"q\[(\d+)\]")
    matches = pattern.findall(gate_line)

    # Convert matched qubit numbers to integers
    return [int(match) for match in matches]


def parse_qasm(filename):
    """Parse a QASM file and return qubits used for each gate
    preserving their order."""
    gates_and_qubits = []
    # if filename is str
    if not isinstance(filename, Path):
        filename = Path(filename)
    with Path.open(filename) as file:
        for _line in file:
            line = _line.strip()

            # Check if line represents a gate operation
            if not line.startswith(
                ("OPENQASM", "include", "qreg", "creg", "gate", "barrier", "measure")
            ):
                qubits = extract_qubits_from_gate(line)
                if qubits:
                    gates_and_qubits.append(tuple(qubits))
    return gates_and_qubits


def compile(filename):
    """Compile a QASM file and return the compiled sequence of qubits."""
    # Check if the file is a valid QASM file
    if not is_qasm_file(filename):
        raise ValueError("Invalid QASM file format")

    # Parse the QASM file to extract the qubits used for each gate
    parse_qasm(filename)
    sequence = parse_qasm(filename)

    return sequence


def get_front_layer(dag):
    """Get the front layer of the DAG."""
    front_layer = []
    for node in dag.get_nodes():
        # If a node has no predecessors, it's in the front layer
        if not dag.direct_predecessors(node.node_id):
            front_layer.append(node)
    return front_layer


def remove_node(dag, node):
    """Execute a node and update the DAG (remove the node and its edges)."""
    # if dag.direct_successors(node.node_id):
    #    for successor in dag.direct_successors(node.node_id):
    #        dag._multi_graph.remove_edge(node.node_id, successor)
    dag._multi_graph.remove_node(node.node_id)


def find_best_gate(front_layer, dist_map):
    """Find the best gate to execute based on distance."""
    min_gate_cost = math.inf
    for _i, gate_node in enumerate(front_layer):
        qubit_indices = gate_node.qindices
        gate_cost = max([dist_map[qs] for qs in qubit_indices])
        # if both ions of 2-qubit gate are in pz execute 2-qubit gate
        if len(qubit_indices) == 2 and gate_cost == 0:
            return gate_node
        if gate_cost < min_gate_cost:
            min_gate_cost = gate_cost
            best_gate = gate_node
    return best_gate


def manual_copy_dag(dag):
    new_dag = DAGDependency()

    # Recreate quantum registers in the new DAG
    for qreg in dag.qregs.values():
        new_dag.add_qreg(qreg)

    # Iterate over all operation nodes in the original DAG and copy them
    for node in dag.get_nodes():
        new_dag.add_op_node(node.op, node.qargs, node.cargs)

    return new_dag


def create_initial_sequence(distance_map, filename, compilation):
    # assert file is a qasm file
    assert is_qasm_file(filename), "The file is not a valid QASM file."

    # generate sequence
    if compilation is False:
        seq = parse_qasm(filename)
        next_node = None
        dag_dep = None
    else:
        qc = QuantumCircuit.from_qasm_file(filename)
        # Remove barriers
        qc = RemoveBarriers()(qc)
        # Remove measurement operations
        qc = RemoveFinalMeasurements()(qc)

        dag_dep = circuit_to_dagdependency(qc)
        working_dag = manual_copy_dag(dag_dep)
        seq = []
        first_flag = True
        while True:
            first_gates = get_front_layer(working_dag)
            if not first_gates:
                break
            first_gate_to_excute = find_best_gate(first_gates, distance_map)
            if first_flag == True:
                next_node = first_gate_to_excute
            first_flag = False
            remove_node(working_dag, first_gate_to_excute)
            seq.append(tuple(first_gate_to_excute.qindices))
        

    flat_seq = [item for sublist in seq for item in sublist]

    return seq, flat_seq, dag_dep, next_node

def update_sequence(dag, dist_map, sequence, max_number_of_front_gates=4):
    """Get the sequence of gates from the DAG.
    Creates a new DAG and removes all gates from it while creating the sequence."""
    working_dag = manual_copy_dag(dag)
    ordered_sequence = []
    i = 0
    while i < max_number_of_front_gates:
        first_gates = get_front_layer(working_dag)
        if not first_gates:
            break
        first_gate_to_excute = find_best_gate(first_gates, dist_map)
        if i == 0:
            first_node = first_gate_to_excute
        i += 1
        remove_node(working_dag, first_gate_to_excute)
        print(tuple(first_gate_to_excute.qindices))
        sequence.remove(tuple(first_gate_to_excute.qindices))
        ordered_sequence.append(tuple(first_gate_to_excute.qindices))
    sequence = ordered_sequence + sequence
    return sequence, first_node


def get_front_layer_non_destructive(dag, virtually_processed_nodes):
    """Get the front layer of the DAG without modifying it."""
    front_layer = []
    
    for node in dag.get_nodes():
        # Skip nodes we've already processed
        if node.node_id in virtually_processed_nodes:
            continue
            
        # Check if all predecessors have been processed
        predecessors = dag.direct_predecessors(node.node_id)
        if not predecessors or all(pred in virtually_processed_nodes for pred in predecessors):
            front_layer.append(node)
            
    return front_layer

def remove_node_by_ions(dag, ions):
    """
    Remove a node from the DAG based on the qubits it operates on.
    
    Args:
        dag: The directed acyclic graph
        qubits: A tuple or list containing the qubits the gate operates on
    
    Returns:
        bool: True if a node was found and removed, False otherwise
    """
    target_node = None
    
    # Search for the node with matching qubits
    for node in dag.get_nodes():
        if tuple(node.qindices) == tuple(ions):
            target_node = node
            break
    
    # If found, remove the node
    if target_node:
        remove_node(dag, target_node)
        return True
    
    return False

#def update_dag_and_sequence(dag, sequence, processed_gates):
  #  for gate in processed_gates:
    #    remove_node_by_ions(dag, gate)
      #  sequence.remove(gate)

def map_gates_to_pzs(graph, front_layer_gates): # new TODO need to swap keys and values? Need a real gate info map? node to pz? then take indices later? node should be unique
    # create list of all gates at each processing zone
    front_layer_info = {pz.name: [] for pz in graph.pzs}
    for seq_elem in front_layer_gates:
        if len(seq_elem) == 1:
            elem = seq_elem[0]
            pz = graph.map_to_pz[elem]
            if front_layer_info[pz] == []:
                front_layer_info[pz].append(elem)
        elif len(seq_elem) == 2:
            # only pick processing zone for 2-qubit gate if not already locked -> then lock it
            if seq_elem not in graph.locked_gates:
                pz = pick_pz_for_2_q_gate(graph, seq_elem[0], seq_elem[1])
                graph.locked_gates[seq_elem] = pz
            else:
                pz = graph.locked_gates[seq_elem]
            if front_layer_info[pz] == []:
                front_layer_info[pz].append(seq_elem[0])
                front_layer_info[pz].append(seq_elem[1])

    return front_layer_info


def remove_processed_gates_from_sequence_non_destructive(graph, dag, dist_map, sequence, max_number_of_front_gates=4):
    """Get the sequence of gates from the DAG non-destructively."""

    ordered_sequence = []
    processed_nodes = set()  # Track nodes we've "virtually removed"

    for pz in range(graph.pzs):
        # Get front layer excluding already processed nodes
        front_layer = get_front_layer_non_destructive(dag, processed_nodes)
        
        if not front_layer:
            break

        first_gates = [node.indices for node in front_layer]
            
        gates_pz_map = map_gates_to_pzs(graph, first_gates)
        first_gates[pz] = find_best_gate(first_gates, dist_map) # new TODO hier first gates und oben, oben aber front layer, front layer in gate_info_list verwenden?
            
        # "Virtually remove" the node by adding it to processed_nodes
        processed_nodes.add(first_gates[pz].node_id)
        
        #print(tuple(first_gate_to_execute.qindices))
        sequence.remove(tuple(first_gates[pz].qindices))
        ordered_sequence.append(tuple(first_gates[pz].qindices))
    
    sequence = ordered_sequence + sequence
    gate_info_list = None #TODO hier für jede pz first gate speichern
    return sequence

def get_all_first_gates_and_update_sequence(graph, dag, dist_map, max_rounds=5):
    """Get the first gates from the DAG for each pz."""
    gate_info_list = {}
    ordered_sequence = []
    processed_nodes = set()  # Track nodes we've "virtually removed"
    
    for round in max_rounds:
        first_gates = {pz: [] for pz in graph.pzs}
        # Get front layer excluding already processed nodes
        first_gates = get_front_layer_non_destructive(dag, processed_nodes)
        
        for pz in range(graph.pzs):
            first_gates[pz].append(find_best_gate(first_gates, dist_map))

            # save first gates of each pz as gate_info_list
            if round == 0:
                gate_info_list[pz.name] = [first_gates[pz][0]]
            
            # "Virtually remove" the node by adding it to processed_nodes (used in get_front_layer)
            processed_nodes.add(first_gates[pz].node_id)

            ordered_sequence.append(tuple(first_gates[pz].qindices))
            sequence.remove(tuple(first_gates[pz].qindices))
    
    sequence = ordered_sequence + sequence
    
    return sequence


if __name__ == "__main__":
    """
    num_ion_chains = 3
    distance_map = {0: 10, 1: 33, 2: 5}#, 3: 11, 4: 2, 5: 3}
    partition =  {'pz1': [2], 'pz2': [1], 'pz3': [0]}#{'pz1': [2, 4], 'pz2': [5, 1], 'pz3': [3], 'pz4': [0]}
    filename = f"QASM_files/qft_no_swaps_nativegates_quantinuum_tket/qft_no_swaps_nativegates_quantinuum_tket_{num_ion_chains}.qasm"
    seq, flat_seq, dag_dep, next_node_initial = create_initial_sequence(
                distance_map, filename, compilation=True
            )
    print('initial sequence: \n', seq)

    processed_gates = [(2, 1), (2,), (0,)]
    update_dag_and_sequence(dag_dep, seq, processed_gates=processed_gates)
    print('updated sequence: \n', seq)

    seq = remove_processed_gates_from_sequence_non_destructive(None, dag_dep, distance_map, seq, max_number_of_front_gates=40)


    
    while True:
        processed_gates = [tuple(node.qindices) for node in dag_dep.get_nodes()][:3]
        update_dag_and_sequence(dag_dep, seq, processed_gates=processed_gates)
        seq = restructure_sequence_non_destructive(dag_dep, distance_map, seq, max_number_of_front_gates=40)

        print(len([node for node in dag_dep.get_nodes()]))

        if len([node for node in dag_dep.get_nodes()]) == 0:
            break
    """
    
    #dag_dep.draw(filename='dags.png')


    m, n, v, h = 3, 3, 1, 1
    number_of_pz = 3
    failing_junctions = 0
    seed = 0
    plot = False
    save = False

    height = -1.5 # height > -1 (position of pz outside of mz) -> all edges can be ordered by sum -> if height -1 -> edge may be smaller in edge_idc (important for get_idx_from_idc())

    exit1 = (float((m-1)*v), float((n-1)*h))
    entry1 = (float((m-1)*v), float(0))
    processing_zone1 = (float((m-1)*v-height), float((n-1)*h/2))

    exit2 = (0.0, 0.0)
    entry2 = (0.0, float((n-1)*h))
    processing_zone2 = (float(height), float((n-1)*h/2))

    exit3 = (float((m-1)*v), float(0))
    entry3 = (float(0), float(0))
    processing_zone3 = (float((m-1)*v/2), float(height))

    exit4 = (float(0), float((n-1)*h))
    entry4 = (float((m-1)*v), float((n-1)*h))
    processing_zone4 = (float((m-1)*v/2), float((n-1)*h-height))

    pz1 = ProcessingZone("pz1", [exit1, entry1, processing_zone1])
    pz2 = ProcessingZone("pz2", [exit2, entry2, processing_zone2])
    pz3 = ProcessingZone("pz3", [exit3, entry3, processing_zone3])
    pz4 = ProcessingZone("pz4", [exit4, entry4, processing_zone4])
    pzs = [pz1, pz2, pz3, pz4][0:number_of_pz]

    basegraph_creator = GraphCreator(m, n, v, h, failing_junctions, pzs)
    MZ_graph = basegraph_creator.get_graph()        
    pzgraph_creator = PZCreator(m, n, v, h, failing_junctions, pzs)
    G = pzgraph_creator.get_graph()
    G.mz_graph = MZ_graph
    
    G.seed = seed

    G.idc_dict = create_idc_dictionary(G)
    G.pzs = pzs
    G.parking_edges_idxs = []
    G.pzs_name_map = {}
    G.edge_to_pz_map = {}
    for pz in G.pzs:
        G.parking_edges_idxs.append(get_idx_from_idc(G.idc_dict, pz.parking_edge))
        G.pzs_name_map[pz.name] = pz
        for edge_idx in pz.path_to_pz_idxs:
            G.edge_to_pz_map[edge_idx] = pz
        G.edge_to_pz_map[get_idx_from_idc(G.idc_dict, pz.parking_edge)] = pz
        for edge_idx in pz.path_from_pz_idxs:
            G.edge_to_pz_map[edge_idx] = pz
    print(f"parking_edges_idxs: {G.parking_edges_idxs}")
    
    G.max_num_parking = 2
    for pz in G.pzs:
        pz.max_num_parking = G.max_num_parking # if changed here (meaning pzs can hold different amounts of ions), also change in shuttle.py (check_duplicates) and check for further updates to max_num_parking

    G.plot = plot
    G.save = save
    G.arch = str([m, n, v, h])

    number_of_mz_edges = len(MZ_graph.edges())
    number_of_chains = math.ceil(.5*len(MZ_graph.edges()))
    

    # plot for paper
    # plot_state(
    #     G, (None, None), plot_ions=True, show_plot=plot, save_plot=save
    # )

    print(f"Number of chains: {number_of_chains}")
    #algorithm = "qft_no_swaps_nativegates_quantinuum_tket"
    algorithm = "full_register_access"
    qasm_file_path = (
        #f"../../../QASM_files/{algorithm}/{algorithm}_{number_of_chains}.qasm"
        f"QASM_files/{algorithm}/{algorithm}_{number_of_chains}.qasm"
    )

    edges = list(G.edges())

    create_starting_config(G, number_of_chains, seed=seed)
    G.idc_dict = create_idc_dictionary(G)
    G.state = get_ions(G)

    best_gates = True
    dist_map = {0: 10, 1: 33, 2: 5, 3: 11, 4: 2, 5: 3}
    if best_gates:
        sequence, flat_sequence, dag, next_node = create_initial_sequence(dist_map, qasm_file_path, compilation=best_gates)
    else:
        sequence = compile(qasm_file_path)
    G.sequence = sequence
    print(len(G.sequence), "len seq")

    # if there is a real tuple in sequence (2-qbuit gate) use partitioning
    # if any(len(i) > 1 for i in sequence):
    #     part = get_partition(qasm_file_path, len(G.pzs))
    #     partition = {pz.name: part[i] for i, pz in enumerate(G.pzs)}
    #     num_pzs = len(G.pzs)
    # else:
    if True:
        # else place them in the closest processing zone (equally distributed)
        # TODO double check
        partition = {pz.name: [] for pz in G.pzs}
        # Assign each ion to the closest processing zone
        for ion, position in G.state.items():
            closest_pz = None
            min_distance = float("inf")
            for pz in G.pzs:
                distance = len(find_path_edge_to_edge(G, position, pz.parking_edge))
                if distance < min_distance:
                    min_distance = distance
                    closest_pz = pz
            partition[closest_pz.name].append(ion)
        # Balance the ions among the processing zones
        all_ions = [ion for ions in partition.values() for ion in ions]
        all_ions.sort(key=lambda ion: G.state[ion])
        num_pzs = len(G.pzs)
        ions_per_pz = len(all_ions) // num_pzs
        for i, pz in enumerate(G.pzs):
            start_index = i * ions_per_pz
            end_index = start_index + ions_per_pz
            partition[pz.name] = all_ions[start_index:end_index]
        # Distribute any remaining ions
        remaining_ions = all_ions[num_pzs * ions_per_pz :]
        for i, ion in enumerate(remaining_ions):
            partition[G.pzs[i % num_pzs].name].append(ion)

    print('partition: ', partition)

    # Create a reverse mapping from element to partition name
    map_to_pz = {
        element: pz
        for pz, elements in partition.items()
        for element in elements
    }
    G.map_to_pz = map_to_pz

    # Ensure all elements are in one of the partitions
    all_partition_elements = []
    for elements in partition.values():
        all_partition_elements.extend(elements)
    unique_sequence = []
    for seq_elem in G.sequence:
        for elem in seq_elem:
            if elem not in unique_sequence:
                unique_sequence.append(elem)
    assert all(element in all_partition_elements for element in unique_sequence)

    if len(G.pzs) > 1:
        # and no element is in both partitions
        pz_sets = {pz: set(elements) for pz, elements in partition.items()}
        common_elements = set.intersection(*pz_sets.values())
        assert (
            not common_elements
        ), f"{common_elements} are overlapping in partitions"

    front_layer = get_front_layer_non_destructive(dag, virtually_processed_nodes=[])
    front_layer_gates = [node.qindices for node in front_layer]
    print(front_layer)
    Gate_info_map_new = pz_gates_map(G, front_layer_gates=front_layer_gates)
    print(Gate_info_map_new)
