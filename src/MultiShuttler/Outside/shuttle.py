from Cycles import get_ion_chains, get_ions_in_pz_and_connections, get_ions_in_exit_connections, get_ions_in_parking, get_state_idxs  
from graph_utils import get_idx_from_idc, get_idc_from_idx
from scheduling import (
    create_move_list,
    create_priority_queue,
    get_partitioned_priority_queues,
    find_movable_cycles,
    rotate_free_cycles,
    create_cycles_for_moves,
    preprocess,
    create_gate_info_list,
    update_entry_and_exit_cycles,
    find_out_of_entry_moves
)
from plotting import plot_state
import os
from datetime import datetime
from collections import Counter

def check_duplicates(graph):
    edge_idxs_occupied = []
    for edge_idc in graph.state.values():
        edge_idxs_occupied.append(get_idx_from_idc(graph.idc_dict, edge_idc))
    # Count occurrences of each integer
    counts = Counter(edge_idxs_occupied)

    for idx, count in counts.items():
        edge_idc = get_idc_from_idx(graph.idc_dict, idx)
        if graph.get_edge_data(edge_idc[0], edge_idc[1])["edge_type"] != "parking_edge" and count > 1:
            message = f"More than one chain in edge {edge_idc}!"
            raise AssertionError(message)
        
        if graph.get_edge_data(edge_idc[0], edge_idc[1])["edge_type"] == "parking_edge" and count > graph.max_num_parking:
            message = (
                f"More than {graph.max_num_parking} chains in parking edge {edge_idc}!"
            )
            raise AssertionError(message)

def find_pz_order(graph, gate_info_list):
    # find next processing zone that will execute a gate
    pz_order = []
    for gate in graph.sequence:
        if len(gate) == 1:
            ion = gate[0]
            for pz in graph.pzs:
                if ion in gate_info_list[pz.name]:
                    pz_order.append(pz.name)
                    break
        elif len(gate) == 2:
            ion1, ion2 = gate
            for pz in graph.pzs:
                if ion1 in gate_info_list[pz.name] and ion2 in gate_info_list[pz.name]:
                    pz_order.append(pz.name)
                    break
    return pz_order


def shuttle(graph, priority_queue, partition, timestep, cycle_or_paths, unique_folder):
    gate_info_list = create_gate_info_list(graph)
    print(f"Gate info list: {gate_info_list}")

    pz_executing_gate_order = find_pz_order(graph, gate_info_list)
    print(f"Next processing zone executing gate: {pz_executing_gate_order}")

    # new: stop moves (ions that are already in the correct processing zone for a two-qubit gate)
    graph.stop_moves = []

    # swapping from Inside Shuttler
    # ion1_needed_in_pz = None
    # # "swap" ions in the same processing zone if only one is needed
    # for pz in graph.pzs:
    #     ions_at_pz = graph[pz.parking_edge[0]][pz.parking_edge[1]]["ions"]
    #     if len(ions_at_pz) == 2:
    #         ion1, ion2 = ions_at_pz
    #         if ion2 not in gate_info_list[pz.name]:
    #             graph[pz.parking_edge[0]][pz.parking_edge[1]]["ions"].remove(ion2)
    #             graph[pz.parking_edge[0]][pz.parking_edge[1]]["ions"].insert(0, ion2)
    #             print("swapped ion2, now: %s" % ions_at_pz)

    #         # find the next processing zone that will execute a gate on ion1
    #         # in case it is needed elsewhere
    #         # solution to bug where ion1 was needed later in the current pz
    #         # but also needed elsewhere
    #         # so was never rotated because other ions where rotating into that pz
    #         # and swapped with ion1 (ion1 would never be rotated then)
    #         for pz_name in pz_executing_gate_order:
    #             if ion1 in gate_info_list[pz_name]:
    #                 ion1_needed_in_pz = pz_name
    #                 break

    #         if ion1_needed_in_pz is not None:
    #             # TODO ion1 swap also necessary?
    #             if ion1_needed_in_pz != pz.name:
    #                 # ion1 not in gate_info_list[pz.name]:
    #                 graph[pz.parking_edge[0]][pz.parking_edge[1]]["ions"].remove(ion1)
    #                 graph[pz.parking_edge[0]][pz.parking_edge[1]]["ions"].insert(0, ion1)
    #                 print("swapped back ion1, now %s" % ions_at_pz)

    #         # new: this could maybe be used to time the gates
    #         # (bug where two ions were randomly in the correct pz already while only
    #         # a 1-qubit gate was executed on one of them and next gate was
    #         # the 2-qubit gate on them -> 3rd ion was moved into them)
    #         if ion1 in gate_info_list[pz.name] and ion2 in gate_info_list[pz.name]:
    #             graph.stop_moves.append(ion1)
    #             graph.stop_moves.append(ion2)

    preprocess(graph, priority_queue)

    # Update ion chains after preprocess
    graph.state = get_ion_chains(graph)
    check_duplicates(graph)
    print(f"priority queue: {priority_queue}")
    part_prio_queues = get_partitioned_priority_queues(graph, priority_queue, partition)

    all_cycles = {}
    # Iterate over all processing zones
    # create move list for each pz -> needed to get all cycles
    # priority queue later picks the cycles to rotate
    in_and_into_exit_moves = {}
    for pz in graph.pzs:
        prio_queue = part_prio_queues[pz.name]
        move_list = create_move_list(graph, prio_queue, pz)
        print(f"Priority queue for {pz.name}: {prio_queue}")
        print(f"Move list for {pz.name}: {move_list}")
        gate_execution_finished = True
        cycles, in_and_into_exit_moves = create_cycles_for_moves(graph, move_list, prio_queue, gate_execution_finished, cycle_or_paths, pz, other_next_edges=None)
        # add cycles to all_cycles
        all_cycles = {**all_cycles, **cycles}
        in_and_into_exit_moves[pz.name] = in_and_into_exit_moves

    # print(f"All cycles: {all_cycles}")
    out_of_entry_moves = find_out_of_entry_moves(graph, other_next_edges=[])
    
    for pz in graph.pzs:
        prio_queue = part_prio_queues[pz.name]
        out_of_entry_moves_pz = out_of_entry_moves[pz] if pz in out_of_entry_moves else None
        update_entry_and_exit_cycles(graph, pz, all_cycles, in_and_into_exit_moves, out_of_entry_moves_pz, gate_execution_finished, prio_queue)

    # now general priority queue picks cycles to rotate
    chains_to_rotate = find_movable_cycles(graph, all_cycles, priority_queue, cycle_or_paths)
    print(f"Chains to rotate: {chains_to_rotate}")

    rotate_free_cycles(graph, all_cycles, chains_to_rotate)

    # new: postprocess?
    # -> ions can already move into processing zone if they pass a junction
    # Update ion chains after rotate
    graph.state = get_ion_chains(graph)
    preprocess(graph, priority_queue)

    labels = ("timestep %s" % timestep, "remaining sequence: %s" % graph.sequence)

    if timestep >= 0:
        plot_state(
            graph,
            labels,
            plot_ions=True,
            show_plot=graph.plot,
            save_plot=graph.save,
            plot_cycle=False,
            plot_pzs=False,
            filename=os.path.join(
                unique_folder, "%s_timestep_%s.png" % (graph.arch, timestep)
            ),
        )


def main(graph, sequence, partition, cycle_or_paths):
    timestep = 0
    max_timesteps = 1e6
    graph.state = get_ion_chains(graph)

    unique_folder = os.path.join("runs", datetime.now().strftime("%Y%m%d_%H%M%S"))
    if graph.save is True:
        os.makedirs(unique_folder, exist_ok=True)

    plot_state(
        graph,
        labels=("Initial state", None),
        plot_ions=True,
        show_plot=graph.plot,
        save_plot=graph.save,
        plot_cycle=False,
        plot_pzs=True,
        filename=os.path.join(
            unique_folder, "%s_timestep_%s.png" % (graph.arch, timestep)
        ),
    )

    graph.in_process = []
    graph.locked_gates = {}
    while timestep < max_timesteps:    
        print(f"\nStarting timestep {timestep}")
        print(f"Remaining sequence: {graph.sequence}")

        for pz in graph.pzs:
            pz.rotate_entry = False

        # priority queue is dict with ions as keys and pz as values
        # (for 2-qubit gates pz may not match the pz of the individual ion)
        priority_queue, next_gate_at_pz = create_priority_queue(graph, sequence)

        # check if ions are already in processing zone ->
        # important for 2-qubit gates
        # -> leave ion in processing zone if needed in a 2-qubit gate
        for i in range(min(len(graph.pzs), len(sequence))):
            # only continue if previous ion was processed
            gate = sequence[i]

            if len(gate) == 2:
                ion1, ion2 = gate
                for pz in graph.pzs:
                    state1 = graph.state[ion1]
                    state2 = graph.state[ion2]
                    # append ion to in_process if it is in the correct processing zone
                    if (
                        state1 == pz.parking_edge
                        and ion1 in next_gate_at_pz[pz.name]
                        and ion2 in next_gate_at_pz[pz.name]
                    ):
                        graph.in_process.append(ion1)
                        # print(f"Added ion {ion1} to in_process")
                    if (
                        state2 == pz.parking_edge
                        and ion1 in next_gate_at_pz[pz.name]
                        and ion2 in next_gate_at_pz[pz.name]
                    ):
                        graph.in_process.append(ion2)
                        # print(f"Added ion {ion2} to in_process")

        # print('in process before shuttling:', graph.in_process)

        # shuttle one timestep
        shuttle(graph, priority_queue, partition, timestep, cycle_or_paths, unique_folder)

        # reset ions in process
        graph.in_process = []

        # Check the state of each ion in the sequence
        graph.state = get_ion_chains(graph)
        processed_ions = []
        previous_ion_processed = True
        pzs = graph.pzs.copy()
        # go through the first gates in the sequence (as many as pzs or sequence length)
        # for now, gates are processed in order
        # (can only be processed in parallel if previous gates are processed)
        for i in range(min(len(graph.pzs), len(sequence))):
            # only continue if previous ion was processed
            if not previous_ion_processed:
                break
            gate = sequence[i]
            ion_processed = False
            # wenn auf weg zu pz in anderer pz -> wird processed?
            # Problem nur für 2-qubit gate? -> TODO fix
            for pz in pzs:
                if len(gate) == 1:
                    ion = gate[0]
                    if get_idx_from_idc(graph.idc_dict, graph.state[ion]) == get_idx_from_idc(graph.idc_dict, pz.parking_edge):
                        print(f"Ion {ion} at Processing Zone {pz.name}")
                        processed_ions.insert(0, (ion,))
                        ion_processed = True
                        # remove the processing zone from the list
                        # (it can only process one ion)
                        pzs.remove(pz)
                        # graph.in_process.append(ion)
                        break
                elif len(gate) == 2:
                    ion1, ion2 = gate
                    state1 = graph.state[ion1]
                    state2 = graph.state[ion2]

                    # The following is now done at the beginning of next timestep
                    # (otherwise would do it double
                    # -> would leave the ones from last time step in in_process
                    #  -> would not move even though
                    # they are move out of pz by preprocessing)
                    # append ion to in_process if it is in the correct processing zone
                    # if state1 == pz.parking_edge and ion1 in next_gate_at_pz[pz.name]
                    # and ion2 in next_gate_at_pz[pz.name]:
                    #     graph.in_process.append(ion1)
                    # if state2 == pz.parking_edge and ion1 in next_gate_at_pz[pz.name]
                    # and ion2 in next_gate_at_pz[pz.name]: # also 1 qubit gate?
                    #     graph.in_process.append(ion2)

                    # if both ions are in the processing zone, process the gate
                    if get_idx_from_idc(graph.idc_dict, state1) == get_idx_from_idc(graph.idc_dict, pz.parking_edge) and get_idx_from_idc(graph.idc_dict, state2) == get_idx_from_idc(graph.idc_dict, pz.parking_edge):
                        print(f"Ions {ion1} and {ion2} at Processing Zone {pz.name}")
                        processed_ions.insert(0, (ion1, ion2))
                        ion_processed = True
                        # remove the processing zone from the list
                        # (it can only process one gate)
                        pzs.remove(pz)

                        # remove the locked pz of the processed two-qubit gate
                        if gate in graph.locked_gates:
                            if graph.locked_gates[gate] == pz.name:
                                graph.locked_gates.pop(gate)
                        break
                else:
                    raise ValueError("Invalid gate format")
            previous_ion_processed = ion_processed

        # Remove processed ions from the sequence
        for gate in processed_ions:
            sequence.remove(gate)

        if len(sequence) == 0:
            print(f"All ions have arrived at their destination in {timestep} timesteps")
            break

        timestep += 1

    return timestep


# TODO entern pzs nicht gleichzeitig
# pz executed gates nicht richtig (ions in pzs aber werden nicht executed)