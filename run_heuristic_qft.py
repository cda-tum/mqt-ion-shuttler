import math
import time
from pathlib import Path

import numpy as np

from Cycles import GraphCreator, MemoryZone
from scheduling import create_initial_sequence, create_starting_config, run_simulation

# parser = argparse.ArgumentParser()
# parser.add_argument("config_file", help="path to json config file")
# # parser.add_argument("--plot", action="store_true", help="plot grid")
# args = parser.parse_args()

# with pathlib.Path(args.config_file).open("r") as f:
#     config = json.load(f)
# arch = config["arch"]
# max_timesteps = config["max_timesteps"]
# num_ion_chains = config["num_ion_chains"]
# filename = config["qu_alg"]

# archs = [[12, 12, 3, 3]]  # [4, 4, 2, 2], [7, 7, 1, 1], [10, 10, 1, 1], [3, 3, 5, 5]]
# arch = [10, 10, 2, 2]

# archs =  [[2, 2, 1, 100]]#, [2, 2, 1, 11], [2, 2, 1, 19], [2, 2, 1, 29], [2, 2, 1, 100],
archs = [ 
    [6, 2, 1, 1],
    [4, 4, 1, 1],
    [5, 5, 1, 1],
    [6, 6, 1, 1],
    # [4, 4, 2, 2],
    # [5, 5, 2, 2],
    # [6, 6, 2, 2],
    # [7, 7, 2, 2],
    # [8, 8, 2, 2],
    # [9, 9, 2, 2],
    # [10, 10, 2, 2],
    # [11, 11, 2, 2],
    # [12, 12, 2, 2],
] 
seeds = [0, 1, 2, 3, 4]  # [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]

compilation = True
#percs = [0.1, 0.25, 0.5, 0.75, 1]
pz = 'outer'#, 'mid']

for arch in archs:
    timestep_arr = []
    cpu_time_arr = []
    start_time = time.time()
    for seed in seeds:
        m, n, v, h = arch
        # create dummy graph
        graph = GraphCreator(m, n, v, h, pz).get_graph()
        n_of_traps = len(
            [trap for trap in graph.edges() if graph.get_edge_data(trap[0], trap[1])["edge_type"] == "trap"]
        )
        num_ion_chains = math.ceil(n_of_traps / 2)
        try:
            ion_chains, number_of_registers = create_starting_config(num_ion_chains, graph, seed=seed)
        except:
            continue
        #filename = "QASM_files/GHZ/ghz_nativegates_quantinuum_tket_%s.qasm" % num_ion_chains #qft_%squbits.qasm" % num_ion_chains
        #filename = "QASM_files/full_register_access/full_register_access_%s.qasm" % num_ion_chains
        filename = "QASM_files/Qft_no_swaps/qft_no_swaps_nativegates_quantinuum_tket_%s.qasm" % num_ion_chains
        max_timesteps = 100000000

        print(f"arch: {arch}, seed: {seed}, registers: {number_of_registers}\n")

        time_2qubit_gate = 3
        time_1qubit_gate = 1
        max_chains_in_parking = 3

        memorygrid = MemoryZone(
            m,
            n,
            v,
            h,
            ion_chains,
            max_timesteps,
            max_chains_in_parking,
            pz,
            time_2qubit_gate=time_2qubit_gate,
            time_1qubit_gate=time_1qubit_gate,
        )

        memorygrid.update_distance_map()
        seq, flat_seq, dag_dep, next_node_initial = create_initial_sequence(
            memorygrid.distance_map, filename, compilation=compilation
        )
        seq_length = len(seq)
        timestep = run_simulation(
            memorygrid, max_timesteps, seq, flat_seq, dag_dep, next_node_initial, max_length=10, show_plot=False
        )
        timestep_arr.append(timestep)
        cpu_time = time.time() - start_time
        cpu_time_arr.append(cpu_time)

        # Create a Path object for the file
        #file_path = Path("VS_SWAP_ghz_results.txt")
        file_path = Path("qft.txt")

        # with file_path.open("a") as file:
        #     line = f"& {arch[0]} {arch[1]} {arch[2]} {arch[3]} {cpu_time} {timestep} \\\\"
        #     file.write(f"{pz}\n" + line + "\n")

    timestep_mean = np.mean(timestep_arr)
    timestep_var = np.var(timestep_arr)
    cpu_time_mean = np.mean(cpu_time_arr)
    # results[j] = timestep_mean
    # cpu_time_results[j] = cpu_time_mean

    try:
        # Open the file using the Path object
        with file_path.open("a") as file:
            line = f"& {arch[0]} {arch[1]} {arch[2]} {arch[3]} & {number_of_registers}/{n_of_traps} & {seq_length} & {timestep_mean} & & {cpu_time_mean} {'s'} {' compilation='}{compilation} \\\\"
            file.write(f"array ts: {timestep_arr}\n" + line + "\n\n")
    except:
        continue