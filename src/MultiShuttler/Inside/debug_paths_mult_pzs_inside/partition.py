from qiskit import QuantumCircuit
from qiskit.transpiler.passes import RemoveBarriers, RemoveFinalMeasurements
import networkx as nx
from networkx.algorithms.community import kernighan_lin_bisection


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
                graph[q0][q1]["weight"] += 1
            else:
                graph.add_edge(q0, q1, weight=1)
        elif len(gate[1]) > 2:
            raise ValueError("Circuit contains gates with more than 2 qubits")

    # # plot graph
    # nx.draw(graph, with_labels=True)
    # nx.draw_networkx_edge_labels(
    #     graph, pos=nx.circular_layout(graph),
    #     edge_labels={(u, v): d['weight'] for u, v, d in graph.edges(data=True)}
    #     )
    # plt.show()

    return graph


# def partition_graph(graph, n):
#     #partitions = []
#     subgraphs = [graph]
#     assert n > 1, "Number of partitions must be greater than 1"
#     assert n <= len(graph.nodes), "Number of partitions
# must be less than or equal to the number of nodes in the graph"
#     #assert n != len(graph.nodes), "Number of partitions must
# be less than the number of nodes in the graph TODO: handle this case?"
#     while len(subgraphs) < n:
#         print('new iteration')
#         new_subgraphs = []
#         for subgraph in subgraphs:
#             #print(f"Start subgraph nodes: {list(subgraph.nodes)}",
# len(subgraph)+len(new_subgraphs),
# [list(new_subg.nodes) for new_subg in new_subgraphs])
#             print('len(subgraph):', len(subgraph))
#             print('len(new_subgraphs):', len(new_subgraphs),'\n')
#             if len(subgraphs) + len(new_subgraphs) - 1 < n:
#                 # handle cases where the subgraph has less
# than 2 nodes (happens if n = # of nodes)
#                 if len(subgraph.nodes()) < 2:
#                 #     #print(f"Subgraph nodes: {list(subgraph.nodes)}")
#                 #     #print(f"Subgraph edges: {list(subgraph.edges)}")
#                     new_subgraphs.append(subgraph)
#                     continue
#                 part1, part2 = kernighan_lin_bisection(subgraph)
#                 #print(f"Part1: {part1}, Part2: {part2}")
#                 #print(f"Subgraph nodes: {list(subgraph.nodes)}")
#                 #print(f"Subgraph edges: {list(subgraph.edges)}")
#                 new_subgraphs.append(graph.subgraph(part1).copy())
#                 new_subgraphs.append(graph.subgraph(part2).copy())
#             else:
#                 new_subgraphs.append(subgraph)
#         subgraphs = new_subgraphs
#         #partitions = subgraphs

#     return subgraphs#partitions


def partition_graph_new(graph, n):
    if n == 1:
        return [graph]

    assert n <= len(
        graph.nodes
    ), f"Number of partitions must be less or equal to the number of nodes {len(graph.nodes), graph.nodes, n}"
    partitions = [graph.copy()]
    new_partitions = []
    while len(partitions) < n:
        len_partitions = len(partitions)
        for i, partition in enumerate(partitions):
            if len_partitions + len(new_partitions) < n:
                if len(partition) < 2:
                    new_partitions.append(partition)
                    len_partitions -= 1
                    continue

                part1, part2 = kernighan_lin_bisection(partition)
                new_partitions.append(graph.subgraph(part1).copy())
                new_partitions.append(graph.subgraph(part2).copy())
            else:
                new_partitions.append(partition)
        partitions = new_partitions
        new_partitions = []
    return partitions


def get_partition(qasm_file_path, n):
    circuit = read_qasm_file(qasm_file_path)
    interaction_graph = construct_interaction_graph(circuit)
    partition_graphs = partition_graph_new(interaction_graph, n)

    partition = []
    for graph in partition_graphs:
        partition.append(list(graph.nodes))

    return partition


if __name__ == "__main__":
    # Example usage
    qasm_file_path = (
        "QASM_files/full_register_access/full_register_access_2.qasm"
        # "QASM_files/QFT_no_swaps/qft_no_swaps_nativegates_quantinuum_tket_36.qasm"
    )
    n = 5
    print(get_partition(qasm_file_path, n))
