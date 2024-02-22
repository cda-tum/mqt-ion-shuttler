import argparse
import json
import pathlib
import contextlib
import random
import networkx as nx
import numpy as np
from qiskit import QuantumCircuit
from qiskit.converters import circuit_to_dagdependency

from compilation import is_qasm_file, manual_copy_dag, remove_node, update_sequence
from Cycles import GraphCreator, MemoryZone, get_idx_from_idc

def create_starting_config(n_of_chains, graph, seed=None):
    if seed is not None:
        random.seed(seed)
        random_starting_traps = random.sample(range(n_of_traps), (n_of_chains))
        starting_traps = []
        traps = [edges for edges in graph.edges() if graph.get_edge_data(edges[0], edges[1])["edge_type"] == "trap"]
        for trap in random_starting_traps:
            starting_traps.append(traps[trap])
    else:
        starting_traps = [
            edges for edges in graph.edges() if graph.get_edge_data(edges[0], edges[1])["edge_type"] == "trap"
        ][: n_of_chains]
    number_of_registers = len(starting_traps)

    # place ions onto traps (ion0 on starting_trap0)
    ion_chains = {}
    for ion, idc in enumerate(starting_traps):
        ion_chains[ion] = idc

    return ion_chains, number_of_registers


def preprocess(memorygrid, sequence):
    need_rotate = [False] * len(sequence)
    while sum(need_rotate) < len(sequence):
        for i, rotate_chain in enumerate(sequence):
            edge_idc = memorygrid.ion_chains[rotate_chain]
            next_edge = memorygrid.find_next_edge(edge_idc)
            state_edges_idx = memorygrid.get_state_idxs()

            if (
                memorygrid.have_common_junction_node(edge_idc, next_edge) is False
                and get_idx_from_idc(memorygrid.idc_dict, next_edge) not in state_edges_idx
            ):
                memorygrid.ion_chains[rotate_chain] = next_edge
            else:
                need_rotate[i] = True
    return memorygrid


def create_move_list(memorygrid, sequence, max_length=10):
    """
    max_length: max length of move_list (if sequence is longer than max_length, only first max_length elements are considered)
    """
    # unique sequence is sequence without repeating elements (for move_list and 2-qubit gates)
    unique_sequence = []
    for seq_elem in sequence:
        if seq_elem not in unique_sequence:
            unique_sequence.append(seq_elem)
            if len(unique_sequence) == max_length:
                break

    path_length_sequence = {}
    move_list = []
    for i, rotate_chain in enumerate(unique_sequence):
        edge_idc = memorygrid.ion_chains[rotate_chain]
        # TODO shortest path here maybe not optimal?
        path_to_go = nx.shortest_path(
            memorygrid.graph,
            edge_idc[0],
            memorygrid.graph_creator.processing_zone,
            lambda _, __, edge_attr_dict: (edge_attr_dict["edge_type"] == "first_entry_connection") * 1e8 + 1,
        )
        path_length_sequence[rotate_chain] = len(path_to_go)

        if i == 0 or sum(
            np.array([path_length_sequence[rotate_chain]] * len(move_list))
            > np.array([path_length_sequence[chain] for chain in move_list])
        ) == len(move_list):
            move_list.append(rotate_chain)

    # add exit edges (needed in rare cases, when chain was moved into exit but dag dependency changed right after that -> chain is in exit but not in move sequence)
    for exit_connection_idc in iontrap.graph_creator.path_to_pz:
        ion = iontrap.find_chain_in_edge(exit_connection_idc)
        if ion is not None and ion not in move_list:
            move_list.insert(0, ion)

    # get chains in all entry edges and place in front
    # chain in entry must move out
    chains_in_entry_connections = []
    for ion, chain_edge_idx in enumerate(memorygrid.get_state_idxs()):
        if chain_edge_idx in memorygrid.graph_creator.path_from_pz_idxs:
            if chain_edge_idx == get_idx_from_idc(memorygrid.idc_dict, memorygrid.graph_creator.entry_edge):
                # place chain in entry at the end of move_list -> so later looping over list leads to chain in entry being first
                chains_in_entry_connections.append(ion)
            else:
                chains_in_entry_connections.insert(0, ion)

    if len(chains_in_entry_connections) > 0:
        for ion in chains_in_entry_connections:
            with contextlib.suppress(Exception):
                move_list.remove(ion)
            move_list = [ion, *move_list]

    return move_list


def create_initial_sequence(distance_map, filename):

    with open(filename) as file:
        first_line = file.readline()
        print(first_line)
    # assert file is a qasm file
    assert is_qasm_file(filename), "The file is not a valid QASM file."

    # generate sequence

    qc = QuantumCircuit.from_qasm_file(filename)
    dag_dep = circuit_to_dagdependency(qc)

    gate_ids, next_node = update_sequence(dag_dep, distance_map)
    seq = [tuple(gate) for gate in gate_ids]
    flat_seq = [item for sublist in seq for item in sublist]

    return seq, flat_seq, dag_dep, next_node


def run_simulation(iontrap, max_timesteps, seq, flat_seq, dag_dep, next_node, max_length):
    time_in_pz_counter = 0
    next_gate_is_two_qubit_gate = len(seq[0]) == 2
    gate_execution_finished = True

    timestep = 0
    while timestep < max_timesteps:
        rotate_entry = False

        # update state_idxs
        iontrap.get_state_idxs()

        ######### PREPROCESS #########
        iontrap = preprocess(iontrap, flat_seq)

        iontrap.update_distance_map()

        ######### CREATE MOVE SEQUENCE #########
        move_list = create_move_list(iontrap, flat_seq, max_length=max_length)

        ######### CREATE CIRCLES #########
        ### create circles for all chains in move_list (dictionary with chain as key and circle_idcs as value)
        chain_to_park = iontrap.find_chain_in_edge(iontrap.graph_creator.path_to_pz[-1])
        if iontrap.count_chains_in_parking() < iontrap.max_num_parking or gate_execution_finished:
            parking_open = True
        else:
            parking_open = False

        all_circles = {}
        stop_exit_edges = []
        # need to find next_edges before for bfs search of "out of entry move"
        next_edges = {}
        for rotate_chain in move_list:
            edge_idc = iontrap.ion_chains[rotate_chain]
            # if chain is needed again (is present in rest of sequence) -> move (only out of entry) towards exit instead of top left
            towards = "exit" if rotate_chain in flat_seq[1:] else (0, 0)
            next_edges[rotate_chain] = iontrap.find_next_edge(edge_idc, towards=towards)

        for rotate_chain in move_list:
            edge_idc = iontrap.ion_chains[rotate_chain]
            next_edge = next_edges[rotate_chain]

            # make edge_idc and next_edge consistent
            edge_idc, next_edge = iontrap.find_ordered_edges(edge_idc, next_edge)

            # moves in pz
            if get_idx_from_idc(iontrap.idc_dict, next_edge) in [
                *iontrap.graph_creator.path_to_pz_idxs,
                get_idx_from_idc(iontrap.idc_dict, iontrap.graph_creator.parking_edge),
            ]:
                all_circles[rotate_chain] = [edge_idc, next_edge]
                # block moves to pz if parking is full
                if (
                    get_idx_from_idc(iontrap.idc_dict, next_edge)
                    in [
                        *iontrap.graph_creator.path_to_pz_idxs,
                        get_idx_from_idc(iontrap.idc_dict, iontrap.graph_creator.parking_edge),
                    ]
                    and parking_open is False
                ) and (get_idx_from_idc(iontrap.idc_dict, next_edge) in stop_exit_edges or stop_exit_edges == []):
                    all_circles[rotate_chain] = [edge_idc, edge_idc]
                    # needed later for blocking moves to parking
                    stop_exit_edges.append(get_idx_from_idc(iontrap.idc_dict, edge_idc))

            # moves without circle
            # also if chain is moving out of entry connections (entry is handled in create_outer_circle)
            elif (
                not iontrap.check_if_edge_is_filled(next_edge)
                or get_idx_from_idc(iontrap.idc_dict, edge_idc) in iontrap.graph_creator.path_from_pz_idxs[:-1]
            ):
                all_circles[rotate_chain] = [edge_idc, next_edge]

            # moves with circle
            else:
                # create circle (deleted in create_outer_circle: in parking circle is a "stop move")
                all_circles[rotate_chain] = iontrap.create_outer_circle(edge_idc, next_edge, next_edges.values())

        # move chain out of parking edge if needed
        chains_in_parking = iontrap.find_chains_in_parking()
        # if pz full and no chain is moving out (not in state_idxs entry edge) but chain is moving in
        if (
            iontrap.count_chains_in_parking() >= iontrap.max_num_parking
            and gate_execution_finished
            and chain_to_park is not None
        ):
            # find least important chain in parking edge
            chain_to_move_out_of_pz = iontrap.find_least_import_chain_in_parking(
                flat_seq, [*chains_in_parking, chain_to_park]
            )
            if chain_to_move_out_of_pz != chain_to_park:
                # move it to entry
                rotate_entry = True
                # change its path/circle to a stop move
                all_circles[chain_to_move_out_of_pz] = [
                    iontrap.graph_creator.path_from_pz[0],
                    iontrap.graph_creator.path_from_pz[0],
                ]

        ######### FIND CIRCLES THAT CAN MOVE #########
        # find circles that can move while first seq ion is moving
        nonfree_circles, free_circle_combs = iontrap.find_nonfree_and_free_circle_idxs(all_circles)
        free_circle_seq_idxs = [move_list[0]]
        for seq_circ in move_list[1:]:
            nonfree = False
            for mov_circ in free_circle_seq_idxs:
                if (seq_circ, mov_circ) in nonfree_circles or (mov_circ, seq_circ) in nonfree_circles:
                    nonfree = True
                    break
            if nonfree is False:
                free_circle_seq_idxs.append(seq_circ)

        ######### ROTATE CIRCLES #########
        # need circles given in idxs for rotate function
        free_circle_idxs = {}
        for seq_idx in free_circle_seq_idxs:
            free_circle_idxs[seq_idx] = [
                get_idx_from_idc(iontrap.idc_dict, edge_idc) for edge_idc in all_circles[seq_idx]
            ]
            # rotate chains
            iontrap.rotate(free_circle_idxs[seq_idx])
            if rotate_entry:
                iontrap.ion_chains[chain_to_move_out_of_pz] = iontrap.graph_creator.path_from_pz[0]

        ######### UPDATE SEQUENCE / PROCESS GATE #########
        gate = seq[0]
        chains_in_parking = iontrap.find_chains_in_parking()
        if sum((gate_element in chains_in_parking) for gate_element in gate) == len(gate):
            gate_execution_finished = False
            time_in_pz_counter += 1

            print(f"\ntime step: {timestep}, gate {seq[0]} is executed,")
            time_gate = time_2qubit_gate if next_gate_is_two_qubit_gate is True else time_1qubit_gate

            if time_in_pz_counter == time_gate:
                ######### END IF SEQUENCE IS FINISHED #########
                if len(seq) == 1:
                    print("\nFull Sequence executed in %s time steps" % timestep)
                    break

                for _ in gate:
                    flat_seq.pop(0)
                time_in_pz_counter = 0
                gate_execution_finished = True

                # update dag
                remove_node(dag_dep, next_node)
                dag_dep = manual_copy_dag(dag_dep)

                gate_ids, next_node = update_sequence(dag_dep, iontrap.distance_map)
                seq = [tuple(gate) for gate in gate_ids]
                flat_seq = [item for sublist in seq for item in sublist]
                next_gate_is_two_qubit_gate = len(seq[0]) == 2

        ######### SETUP NEW TIME STEP #########
        timestep += 1


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("config_file", help="path to json config file")
    #parser.add_argument("--plot", action="store_true", help="plot grid")
    args = parser.parse_args()

    with pathlib.Path(args.config_file).open("r") as f:
        config = json.load(f)
    arch = config["arch"]
    max_timesteps = config["max_timesteps"]
    num_ion_chains = config["num_ion_chains"]
    filename = config["qu_alg"]
    
    seed = 0
    m, n, v, h = arch
    # create dummy graph
    graph = GraphCreator(m, n, v, h).get_graph()
    n_of_traps = len([trap for trap in graph.edges() if graph.get_edge_data(trap[0], trap[1])["edge_type"] == "trap"])
    ion_chains, number_of_registers = create_starting_config(num_ion_chains, graph, seed=seed)


    print(f"arch: {arch}, seed: {seed}, registers: {number_of_registers}\n")
    max_timesteps = 50000

    time_2qubit_gate = 3
    time_1qubit_gate = 1
    max_chains_in_parking = 3

    iontrap = MemoryZone(
        m,
        n,
        v,
        h,
        ion_chains,
        max_timesteps,
        max_chains_in_parking,
        time_2qubit_gate=time_2qubit_gate,
        time_1qubit_gate=time_1qubit_gate,
    )

    iontrap.update_distance_map()
    distance_map = iontrap.distance_map

    seq, flat_seq, dag_dep, next_node = create_initial_sequence(distance_map, filename)
    run_simulation(iontrap, max_timesteps, seq, flat_seq, dag_dep, next_node, max_length=10)
