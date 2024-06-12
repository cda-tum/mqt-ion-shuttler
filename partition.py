from qiskit import QuantumCircuit
from qiskit.transpiler.passes import RemoveBarriers, RemoveFinalMeasurements
import networkx as nx
from networkx.algorithms.community import kernighan_lin_bisection
import matplotlib.pyplot as plt

def read_qasm_file(file_path):
    circuit = QuantumCircuit.from_qasm_file(file_path)
    # Remove barriers
    circuit = RemoveBarriers()(circuit)
    # Remove measurement operations
    circuit = RemoveFinalMeasurements()(circuit)
    return circuit

def construct_interaction_graph(circuit):
    graph = nx.Graph()
    qubits = circuit.qubits
    
    for qubit in qubits:
        graph.add_node(qubit._index)
    
    for gate in circuit.data:
        if len(gate[1]) == 2:  # Check if it is a 2-qubit gate
            q0 = gate[1][0]._index
            q1 = gate[1][1]._index
            if graph.has_edge(q0, q1):
                graph[q0][q1]['weight'] += 1
            else:
                graph.add_edge(q0, q1, weight=1)
        elif len(gate[1]) > 2:
            raise ValueError("Circuit contains gates with more than 2 qubits")
        
    # plot graph
    nx.draw(graph, with_labels=True)
    nx.draw_networkx_edge_labels(graph, pos=nx.circular_layout(graph), edge_labels={(u, v): d['weight'] for u, v, d in graph.edges(data=True)})
    plt.show()

    return graph

def partition_graph(graph, n):
    partitions = []
    subgraphs = [graph]
    
    while len(partitions) < n:
        new_subgraphs = []
        for subgraph in subgraphs:
            if len(partitions) + len(new_subgraphs) < n:
                part1, part2 = kernighan_lin_bisection(subgraph)
                new_subgraphs.append(graph.subgraph(part1))
                new_subgraphs.append(graph.subgraph(part2))
            else:
                new_subgraphs.append(subgraph)
        subgraphs = new_subgraphs
        partitions = subgraphs
    
    return partitions

def partition(qasm_file_path, n):
    circuit = read_qasm_file(qasm_file_path)
    interaction_graph = construct_interaction_graph(circuit)
    partitions = partition_graph(interaction_graph, n)
    
    for i, partition in enumerate(partitions):
        print(f"Partition {i + 1}: {list(partition.nodes)}")

if __name__ == '__main__':
    # Example usage
    qasm_file_path = 'QASM_files/qft_no_swaps/QFT_no_swaps_nativegates_quantinuum_tket_6.qasm'
    n = 3
    partition(qasm_file_path, n)

