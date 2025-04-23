import argparse
import json
import pathlib
import math
import networkx as nx
import numpy as np
from datetime import datetime
import time # Use time for CPU timing if needed

# Import your project modules
from graph_utils import GraphCreator, PZCreator, ProcessingZone, create_idc_dictionary, get_idx_from_idc
from Cycles import create_starting_config, get_ions, get_state_idxs # Removed find_path_edge_to_edge as it's in graph_utils/scheduling now
from scheduling import get_partitioned_priority_queues, create_priority_queue # Import necessary functions
from shuttle import main as run_shuttle_main # Rename imported main
from compilation import create_dag, create_updated_sequence_destructive, get_all_first_gates_and_update_sequence_non_destructive, map_front_gates_to_pzs, create_dist_dict, update_distance_map, create_initial_sequence
# from get_baseline import minimal_schedule_qiskit_dagdependency # Assuming you might want this later
from partition import get_partition
from plotting import plot_state

# --- Argument Parsing ---
parser = argparse.ArgumentParser(description="Run MQT IonShuttler")
parser.add_argument("config_file", help="Path to the JSON configuration file")
# Add any other command-line args if needed, e.g., overriding plot/save
# parser.add_argument("--plot", action="store_true", help="Show plots during execution")
# parser.add_argument("--save", action="store_true", help="Save plots to 'runs' directory")
args = parser.parse_args()

# --- Load Configuration ---
try:
    with pathlib.Path(args.config_file).open("r") as f:
        config = json.load(f)
except FileNotFoundError:
    print(f"Error: Configuration file not found at {args.config_file}")
    exit(1)
except json.JSONDecodeError:
    print(f"Error: Could not parse JSON file {args.config_file}")
    exit(1)

# --- Extract Parameters from Config ---
arch = config.get("arch")
num_pzs_config = config.get("num_pzs", 1)
seed = config.get("seed", 0) # Default seed to 0 if not provided
algorithm_name = config.get("algorithm_name")
num_ions_config = config.get("num_ions")
use_compilation = config.get("use_compilation", True)
use_paths = config.get("use_paths", False)
max_timesteps = config.get("max_timesteps", 100000) # Provide a default max
plot_flag = config.get("plot", False)
save_flag = config.get("save", False)
failing_junctions = config.get("failing_junctions", 0)
# Define base path for QASM files if needed
qasm_base_dir = config.get("qasm_base_dir", "../../../QASM_files") # Example, adjust as needed

# --- Validate Config ---
if not all([arch, algorithm_name, num_ions_config]):
    print("Error: Missing required parameters in config file (arch, algorithm_name, num_ions)")
    exit(1)
if not isinstance(arch, list) or len(arch) != 4:
    print("Error: 'arch' must be a list of 4 integers [m, n, v, h]")
    exit(1)
# Add more validation as needed...

# --- Setup ---
start_time = datetime.now()
cycle_or_paths_str = "Paths" if use_paths else "Cycles"
m, n, v, h = arch

# --- PZ Definitions (Could also be part of JSON if more complex) ---
height = -4.5 # Consider making this configurable if it varies
pz_definitions = {
    "pz1": ProcessingZone("pz1", [(float((m-1)*v), float((n-1)*h)), (float((m-1)*v), float(0)), (float((m-1)*v-height), float((n-1)*h/2))]),
    "pz2": ProcessingZone("pz2", [(0.0, 0.0), (0.0, float((n-1)*h)), (float(height), float((n-1)*h/2))]),
    "pz3": ProcessingZone("pz3", [(float((m-1)*v), float(0)), (float(0), float(0)), (float((m-1)*v/2), float(height))]),
    "pz4": ProcessingZone("pz4", [(float(0), float((n-1)*h)), (float((m-1)*v), float((n-1)*h)), (float((m-1)*v/2), float((n-1)*h-height))])
}
available_pz_names = list(pz_definitions.keys())
pzs_to_use = [pz_definitions[name] for name in available_pz_names[:num_pzs_config]]

if not pzs_to_use:
    print(f"Error: num_pzs ({num_pzs_config}) is invalid or results in no PZs selected.")
    exit(1)

print(f"Using {len(pzs_to_use)} PZs: {[pz.name for pz in pzs_to_use]}")
print(f"Architecture: {arch}, Seed: {seed}")
print(f"Algorithm: {algorithm_name}, ions: {num_ions_config}")
print(f"Compilation: {use_compilation}, Conflict Resolution: {cycle_or_paths_str}")

# --- Graph Creation ---
basegraph_creator = GraphCreator(m, n, v, h, failing_junctions, pzs_to_use)
MZ_graph = basegraph_creator.get_graph()
pzgraph_creator = PZCreator(m, n, v, h, failing_junctions, pzs_to_use)
G = pzgraph_creator.get_graph()
G.mz_graph = MZ_graph # Attach MZ graph for BFS lookups if needed by Cycles/Paths

G.seed = seed
G.idc_dict = create_idc_dictionary(G)
G.pzs = pzs_to_use # List of ProcessingZone objects
G.parking_edges_idxs = []
G.pzs_name_map = {} # Map from pz name to pz object
G.edge_to_pz_map = {} # Map from edge index to owning pz object (for non-MZ edges)
for pz in G.pzs:
    if not hasattr(pz, 'parking_edge'): # Ensure PZCreator added this
         print(f"Error: PZ {pz.name} seems malformed (missing parking_edge).")
         exit(1)
    parking_idx = get_idx_from_idc(G.idc_dict, pz.parking_edge)
    G.parking_edges_idxs.append(parking_idx)
    G.pzs_name_map[pz.name] = pz
    # Populate edge_to_pz_map for edges belonging *only* to this PZ's structure
    for edge_idx in pz.pz_edges_idx:
        G.edge_to_pz_map[edge_idx] = pz
print(f"Parking Edges Idxs: {G.parking_edges_idxs}")

G.max_num_parking = 2 # Make this configurable?
for pz in G.pzs:
    pz.max_num_parking = G.max_num_parking

G.plot = plot_flag
G.save = save_flag
G.arch = str(arch) # For plotting/logging

number_of_mz_edges = len(MZ_graph.edges())
# Ensure num_ions is determined correctly (either fixed or based on graph size)
if isinstance(num_ions_config, str) and "ceil" in num_ions_config:
    # Handle formulas like "ceil(0.5*num_edges)" - requires careful parsing
    # For simplicity, let's assume num_ions is given as an int or calculated based on MZ size
    try:
         factor = float(num_ions_config.split('*')[0].split('(')[1])
         number_of_ions = math.ceil(factor * number_of_mz_edges)
    except:
         print("Warning: Could not parse num_ions formula, defaulting to config value if integer.")
         number_of_ions = num_ions_config if isinstance(num_ions_config, int) else number_of_mz_edges # Fallback
else:
    number_of_ions = int(num_ions_config)


print(f"Number of ions: {number_of_ions}")

qasm_file_path = pathlib.Path(qasm_base_dir) / algorithm_name / f"{algorithm_name}_{number_of_ions}.qasm"

if not qasm_file_path.is_file():
    print(f"Error: QASM file not found at {qasm_file_path}")
    exit(1)

# --- Initial State & Sequence ---
create_starting_config(G, number_of_ions, seed=seed)
G.state = get_ions(G) # Get initial state {ion: edge_idc}

G.sequence = create_initial_sequence(qasm_file_path)
seq_length = len(G.sequence)
print(f'Initial sequence length: {seq_length}')

# --- Partitioning ---
partitioning = True # Make configurable?
if partitioning:
    part = get_partition(qasm_file_path, len(G.pzs))
    # Ensure partition list length matches num_pzs
    if len(part) != len(G.pzs):
         print(f"Warning: Partitioning returned {len(part)} parts, but expected {len(G.pzs)}. Adjusting...")
         # Simple fix: assign remaining qubits to the last partition, or distribute evenly.
         # This might need a more sophisticated balancing strategy.
         if len(part) < len(G.pzs):
             print("Error: Partitioning failed to produce enough parts.")
             # Handle error appropriately, maybe fall back to non-partitioned approach or exit.
             exit(1)
         else: # More parts than PZs, merge extra parts into the last ones
             part = part[:len(G.pzs)-1] + [qubit for sublist in part[len(G.pzs)-1:] for qubit in sublist]

    partition = {pz.name: part[i] for i, pz in enumerate(G.pzs)}
    print(f'Partition: {partition}')
else:
    # Fallback: Assign ions to closest PZ (example logic)
    print("Partitioning disabled. Assigning ions to closest PZ (basic method).")
    partition = {pz.name: [] for pz in G.pzs}
    # ... (implement closest PZ assignment logic as in your original file) ...
    # Make sure this logic correctly assigns *all* ions involved in the sequence.

# Create reverse map and validate partition
map_to_pz = {}
all_partition_elements = []
for pz_name, elements in partition.items():
    all_partition_elements.extend(elements)
    for element in elements:
        if element in map_to_pz:
            print(f"Warning: Qubit {element} assigned to multiple partitions ({map_to_pz[element]}, {pz_name}). Check partitioning logic.")
            # Decide on a conflict resolution strategy if needed
        map_to_pz[element] = pz_name
G.map_to_pz = map_to_pz

# Validation
unique_sequence_qubits = set(item for sublist in G.sequence for item in sublist)
missing_qubits = unique_sequence_qubits - set(all_partition_elements)
if missing_qubits:
    print(f"Error: Qubits {missing_qubits} from sequence are not in any partition.")
    # This indicates a problem with partitioning or qubit indexing.
    exit(1)
# Check for overlaps if needed (already done within map_to_pz creation loop)


# --- Compilation Setup (if enabled) ---
dag = None
dag_copy = None # Store original DAG if needed
if use_compilation:
    try:
        dag = create_dag(qasm_file_path)
        dag_copy = dag.copy() # Keep a copy of the original DAG if needed later
        # Initial DAG-based sequence update
        G.dist_dict = create_dist_dict(G)
        state_idxs = get_state_idxs(G) # {ion: edge_idx}
        G.dist_map = update_distance_map(G, state_idxs) # {ion: {pz_name: dist}}
        # This sequence update is destructive to the DAG passed in
        # sequence, flat_sequence, _ = create_updated_sequence_destructive(G, qasm_file_path, dag, compilation=use_compilation)
        # G.sequence = sequence # Overwrite initial sequence
        print("Compilation enabled, using DAG-based scheduling.")
    except Exception as e:
        print(f"Error during DAG creation or initial sequence update: {e}")
        print("Falling back to non-compiled sequence.")
        use_compilation = False # Disable compilation if setup fails
        dag = None
        G.sequence = create_initial_sequence(qasm_file_path) # Revert to basic sequence
else:
    print("Compilation disabled, using static QASM sequence.")

# --- Run Simulation ---
print("\nStarting shuttling simulation...")
# Initialize PZ states
for pz in G.pzs:
    pz.getting_processed = [] # Track nodes being processed by this PZ

# Run the main shuttling logic
final_timesteps = run_shuttle_main(G, partition, dag, cycle_or_paths_str, compilation=use_compilation)

# --- Results ---
end_time = datetime.now()
cpu_time = end_time - start_time

print(f"\nSimulation finished in {final_timesteps} timesteps.")
print(f"Total CPU time: {cpu_time}")

# # --- Benchmarking Output (Optional) ---
# bench_filename = f"benchmarks/{start_time.strftime('%Y%m%d_%H%M%S')}_{algorithm_name}.txt"
# pathlib.Path("benchmarks").mkdir(exist_ok=True)
# benchmark_output = (
#     f"{arch}, ions{number_of_ions}/pos{number_of_mz_edges}: {number_of_ions/number_of_mz_edges if number_of_mz_edges > 0 else 0:.2f}, "
#     f"#pzs: {len(pzs_to_use)}, ts: {final_timesteps}, cpu_time: {cpu_time.total_seconds():.2f}, "
#     f"gates: {seq_length}, baseline: {None}, compilation: {use_compilation}, paths: {use_paths}, "
#     f"seed: {seed}, failing_jcts: {failing_junctions}\n"
# )
# try:
#     with open(bench_filename, "a") as f:
#         f.write(benchmark_output)
#     print(f"Benchmark results appended to {bench_filename}")
# except Exception as e:
#     print(f"Warning: Could not write benchmark file: {e}")