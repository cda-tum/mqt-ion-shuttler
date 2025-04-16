from compilation import create_dag
from qiskit.dagcircuit import DAGDependency

def minimal_schedule_qiskit_dagdependency(dag, partition):
    print(partition)
    qubit_zone = {q: zone for zone, qubits in partition.items() for q in qubits}
    zone_available_time = {zone: 0 for zone in partition}
    gate_completion_time = {}

    # Initialize every node to zero
    for node in dag.topological_nodes():
        gate_completion_time[node.node_id] = 0

    # Now schedule only real gates
    for node in dag.topological_nodes():
        qubits = [q.index for q in node.qargs]
        zone = qubit_zone[qubits[0]]
        duration = 3 if len(qubits) == 2 else 1
        earliest_start = zone_available_time[zone]

        pred_ids = dag.direct_predecessors(node.node_id)
        dependency_times = [gate_completion_time[pid] for pid in pred_ids]

        if dependency_times:
            earliest_start = max(earliest_start, max(dependency_times))

        completion_time = earliest_start + duration
        zone_available_time[zone] = completion_time
        gate_completion_time[node.node_id] = completion_time

    return max(gate_completion_time.values())






# qubit_zone_partition = {'pz1': [11, 2, 5, 6, 3, 0, 8, 9, 4, 10, 7, 1]}
# algorithm = "qft_no_swaps_nativegates_quantinuum_tket"
# #algorithm = "full_register_access"
# number_of_chains = 12
# qasm_file_path = (
#     #f"../../../QASM_files/{algorithm}/{algorithm}_{number_of_chains}.qasm"
#     f"QASM_files/{algorithm}/{algorithm}_{number_of_chains}.qasm"
# )
# dag = create_dag(qasm_file_path)

# print(minimal_schedule_qiskit_dagdependency(dag, qubit_zone_partition))
# # Output: 5
