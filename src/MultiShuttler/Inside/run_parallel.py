from graph_utils import GraphCreator, create_idc_dictionary
from Cycles import create_starting_config, find_path_edge_to_edge
from scheduling import ProcessingZone, get_ion_chains
from shuttle import main
from partition import get_partition
from compilation import compile
import math
import networkx as nx
import numpy as np
from datetime import datetime
import sys

plot = False
save = False

number_of_pzs_list = [2, 3, 4, 5, 6, 7, 8, 9, 10]
# archs = [
# [3, 3, 2, 2],
# [3, 3, 1, 1],
# [3, 3, 2, 2],
# [3, 3, 3, 3],
# [4, 4, 1, 1],
# [4, 4, 2, 2],
# [4, 4, 3, 3],
# [5, 5, 1, 1],
# [5, 5, 2, 2],
# [5, 5, 3, 3],
# [6, 6, 1, 1],
# [6, 6, 2, 2],
# [6, 6, 3, 3],
# [7, 7, 1, 1],
# [7, 7, 2, 2],

# [7, 7, 3, 3],
# [8, 8, 1, 1],
# [8, 8, 2, 2],
# [8, 8, 3, 3],
# [9, 9, 1, 1],
# [9, 9, 2, 2],
# [9, 9, 3, 3],
# [10, 10, 1, 1],
# [10, 10, 2, 2],
# [10, 10, 3, 3],
# ]

seeds = [0, 1, 2, 3, 4]
time = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")

number_of_pzs = int(sys.argv[1])
m, n, ion_chain_size_vertical, ion_chain_size_horizontal = [
    int(sys.argv[2]),
    int(sys.argv[3]),
    int(sys.argv[4]),
    int(sys.argv[5]),
]

# for m, n, ion_chain_size_vertical, ion_chain_size_horizontal in archs:
timesteps_array = []
cpu_time_array = []

for seed in seeds:
    # for number_of_pzs in number_of_pzs_list:
    start_time = datetime.now()

    graph_creator = GraphCreator(
        m, n, ion_chain_size_vertical, ion_chain_size_horizontal
    )
    G = graph_creator.get_graph()
    G.plot = plot
    G.save = save
    G.arch = str([m, n, ion_chain_size_vertical, ion_chain_size_horizontal])

    number_of_chains = math.ceil(len(G.edges()) / 2)
    print(f"Number of chains: {number_of_chains}")
    algorithm = "full_register_access"
    qasm_file_path = f"QASM_files/{algorithm}/{algorithm}_{number_of_chains}.qasm"

    # edges = list(G.edges())
    # # Select the middle edge
    # middle_index = math.ceil(len(edges) / 2)
    # middle_edge = edges[middle_index]
    # pz1 = ProcessingZone("pz1", ((0, 0), (1, 0)))
    # pz2 = ProcessingZone(
    #     "pz2",
    #     (middle_edge),
    # )
    # pz3 = ProcessingZone(
    #     "pz3", ((max(G.nodes)[0], max(G.nodes)[1] - 1),
    # (max(G.nodes)[0], max(G.nodes)[1]))
    # )

    def add_processing_zones(graph, num_zones):
        edges = list(graph.edges)
        if len(edges) < num_zones:
            raise ValueError(
                "Number of processing zones exceeds number of available edges."
            )

        # Select edges evenly spaced
        indices = np.linspace(0, len(edges) - 1, num_zones, dtype=int)
        selected_edges = [edges[i] for i in indices]

        processing_zones = []
        for idx, edge in enumerate(selected_edges):
            pz_name = f"pz{idx + 1}"
            processing_zone = ProcessingZone(pz_name, edge)
            processing_zones.append(processing_zone)
            # Optional: Set edge attributes if necessary
            nx.set_edge_attributes(graph, {edge: "processing"}, "edge_type")

        return processing_zones

    G.pzs = add_processing_zones(G, number_of_pzs)

    # G.pzs = [pz1, pz2, pz3]

    create_starting_config(G, number_of_chains, seed=seed)
    G.idc_dict = create_idc_dictionary(G)
    G.state = get_ion_chains(G)

    sequence = compile(qasm_file_path)
    G.sequence = sequence
    print(sequence)

    # if there is a real tuple in sequence (2-qbuit gate) use partitioning
    if any(len(i) > 1 for i in sequence):
        part = get_partition(qasm_file_path, len(G.pzs))
        partition = {pz.name: part[i] for i, pz in enumerate(G.pzs)}
        num_pzs = len(G.pzs)
    else:
        # else place them in the closest processing zone (equally distributed)
        # TODO double check
        partition = {pz.name: [] for pz in G.pzs}
        # Assign each ion to the closest processing zone
        for ion, position in G.state.items():
            closest_pz = None
            min_distance = float("inf")
            for pz in G.pzs:
                distance = len(find_path_edge_to_edge(G, position, pz.edge_idc))
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

    print(partition)

    # Create a reverse mapping from element to partition name
    map_to_pz = {
        element: pz for pz, elements in partition.items() for element in elements
    }
    G.map_to_pz = map_to_pz

    # Ensure all elements are in one of the partitions
    all_partition_elements = []
    for elements in partition.values():
        all_partition_elements.extend(elements)
    unique_sequence = []
    for seq_elem in sequence:
        for elem in seq_elem:
            if elem not in unique_sequence:
                unique_sequence.append(elem)
    assert all(element in all_partition_elements for element in unique_sequence)

    if len(G.pzs) > 1:
        # and no element is in both partitions
        pz_sets = {pz: set(elements) for pz, elements in partition.items()}
        common_elements = set.intersection(*pz_sets.values())
        assert not common_elements, f"{common_elements} are overlapping in partitions"

    timesteps = main(G, sequence, partition)
    end_time = datetime.now()
    cpu_time = end_time - start_time

    timesteps_array.append(timesteps)
    cpu_time_array.append(cpu_time)

    # save timesteps in a file
    with open(f"benchmarks/{time}_{num_pzs}_{algorithm}.txt", "a") as f:
        f.write(
            f"{m, n, ion_chain_size_vertical, ion_chain_size_horizontal}, #pzs: {num_pzs}, ts: {timesteps}, cpu_time: {cpu_time}, seed: {seed}\n"
        )

# calculate averages
timesteps_array = np.array(timesteps_array)
cpu_time_array = np.array(cpu_time_array)
timesteps_average = np.mean(timesteps_array)
cpu_time_average = np.mean(cpu_time_array)

# save averages
with open(f"benchmarks/{time}_{num_pzs}_{algorithm}.txt", "a") as f:
    f.write(
        f"{m, n, ion_chain_size_vertical, ion_chain_size_horizontal}, #pzs: {num_pzs}, average_ts: {timesteps_average}, average_cpu_time: {cpu_time_average}\n"
    )
