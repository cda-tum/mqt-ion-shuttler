import math
import time
from pathlib import Path
import numpy as np
import networkx as nx
from datetime import datetime
from Cycles_new import BaseGraphCreator, PZGraphCreator, MemoryZone
from scheduling import create_initial_sequence, create_starting_config, run_simulation

def run_simulation_for_architecture(arch, seeds, pz, max_timesteps, failing_junctions, compilation=True):
    """
    Runs simulations for the given architecture and seeds, logs the results.

    Args:
        arch (list): Architecture parameters.
        seeds (list): List of seed values.
        pz (str): Position of Processing zone.
        max_timesteps (int): Maximum timesteps.
        compilation (bool): Compilation flag (Gate Selection Step).

    Returns:
        tuple: (timestep_arr, cpu_time_arr, number_of_registers, n_of_traps, seq_length)
    """

    timestep_arr = []
    cpu_time_arr = []
    start_time = time.time()
    flag_seed = 0
    for seed in seeds:
        m, n, v, h = arch
        basegraph_creator = BaseGraphCreator(m, n, v, h, pz, failing_junctions)
        MZ_graph = basegraph_creator.get_graph()        
        pzgraph_creator = PZGraphCreator(m, n, v, h, pz, failing_junctions)
        graph = pzgraph_creator.get_graph()

        try:
            for node in MZ_graph.nodes():
                nx.shortest_path(MZ_graph, node, pzgraph_creator.exit)
        except:
            print('skipped architecture for seed', seed)
            continue
        
        n_of_traps = len([trap for trap in graph.edges() if graph.get_edge_data(trap[0], trap[1])["edge_type"] == "trap"])
        # num_ion_chains = math.ceil(n_of_traps / 2)
        num_ion_chains = math.ceil(((arch[0]-1)*arch[1]*arch[2] + (arch[1]-1)*arch[0]*arch[3]) / 2)

        
        try:
            ion_chains, number_of_registers = create_starting_config(num_ion_chains, graph, seed=seed)
        except:
            continue
        print(f"ion chains: {ion_chains}, number of registers: {number_of_registers}")
        filename = f"QASM_files/full_register_access/full_register_access_{num_ion_chains}.qasm"
        #filename = f"QASM_files/QFT_no_swaps/qft_no_swaps_nativegates_quantinuum_tket_{num_ion_chains}.qasm"
        print(f"arch: {arch}, seed: {seed}, registers: {number_of_registers}\n")

        time_2qubit_gate = 1
        time_1qubit_gate = 1
        max_chains_in_parking = 3

        memorygrid = MemoryZone(
            pzgraph_creator, MZ_graph, ion_chains, max_timesteps, max_chains_in_parking,
            time_2qubit_gate=time_2qubit_gate, time_1qubit_gate=time_1qubit_gate
        )

        memorygrid.update_distance_map()
        seq, flat_seq, dag_dep, next_node_initial = create_initial_sequence(
            memorygrid.distance_map, filename, compilation=compilation
        )
        seq_length = len(seq)
        #print(f"seq: {seq}")
        timestep = run_simulation(
            memorygrid, max_timesteps, seq, flat_seq, dag_dep, next_node_initial, max_length=10
        )
        timestep_arr.append(timestep)
        cpu_time = time.time() - start_time
        cpu_time_arr.append(cpu_time)

        flag_seed += 1
        if flag_seed == 5:
            break

    return timestep_arr, cpu_time_arr, number_of_registers, n_of_traps, seq_length

def log_results(arch, timestep_arr, cpu_time_arr, number_of_registers, n_of_traps, seq_length, run_folder, compilation=True):
    """
    Logs the results of the simulation to a file.

    Args:
        arch (list): Architecture parameters.
        timestep_arr (list): List of timesteps.
        cpu_time_arr (list): List of CPU times.
        number_of_registers (int): Number of registers.
        n_of_traps (int): Number of traps.
        seq_length (int): Sequence length.
        compilation (bool): Compilation flag (Gate Selection Step).
    """
    timestep_mean = np.mean(timestep_arr)
    timestep_var = np.var(timestep_arr)
    cpu_time_mean = np.mean(cpu_time_arr)
    print(cpu_time_mean)
    print(f"timestep mean: {timestep_mean}, timestep var: {timestep_var}, cpu time mean: {cpu_time_mean}")
    
    file_path = Path(run_folder) / f"paths_junctions_0jct.txt"
        #Path("paths_junctions_7jct.txt")
    try:
        with file_path.open("a") as file:
            line = (
                f"& {arch[0]} {arch[1]} {arch[2]} {arch[3]} & {number_of_registers}/{n_of_traps} & {seq_length} "
                f"& {timestep_mean} & {cpu_time_mean} s & Gate Selection={compilation} \\\\"
            )
            file.write(f"array ts: {timestep_arr}\n" + line + "\n\n")
    except:
        pass

def main():
    archs = [
        # [2, 2, 1, 5],
        # [2, 2, 1, 11],
        # [2, 2, 1, 29],
        # [2, 2, 1, 39],
        
        [3, 3, 1, 1],
        #[5, 5, 1, 1],
        #[6, 6, 1, 1],
        #[10, 10, 1, 1],
        
        #[4, 4, 2, 2],
    ]
    # 6 6 1 1 seed = 3 failing junctions = 7 does not work 
    # -> only one way to exit and two ways from entry 
    # -> infinite loop if needed ion is too far away and is blocked by ions coming from entry
    seeds = range(0, 99)
    pz = 'outer'
    max_timesteps = 1_000
    compilation = False
    failing_junctions = 0

    # Create a folder for each run with a timestamp
    run_folder = Path(f'results/run_{datetime.now().strftime("%Y%m%d")}')
    run_folder.mkdir(parents=True, exist_ok=True)

    for arch in archs:
        #try:
            timestep_arr, cpu_time_arr, number_of_registers, n_of_traps, seq_length = run_simulation_for_architecture(
            arch, seeds, pz, max_timesteps, failing_junctions, compilation=compilation
            )
            log_results(arch, timestep_arr, cpu_time_arr, number_of_registers, n_of_traps, seq_length, run_folder, compilation=compilation)
        #except:
        #    print('skipped all seeds for architecture', arch)
        #    continue

if __name__ == "__main__":
    main()
