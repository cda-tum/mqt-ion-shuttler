import logging
import numpy as np
from more_itertools import pairwise, distinct_combinations
import contextlib
import networkx as nx
from Cycles import (
    get_ions,
    get_edge_state,
    find_path_edge_to_edge,
    find_next_edge,
    find_ordered_edges,
    check_if_edge_is_filled,
    create_cycle,
    have_common_junction_node,
    get_ions_in_pz_and_connections,
    get_ions_in_exit_connections,
    get_ions_in_parking,
    find_ions_in_parking,
    find_ion_in_edge,
    get_state_idxs,
    find_least_import_ion_in_parking,
    find_conflict_cycle_idxs,
    find_path_node_to_edge
)
from Paths import create_path_via_bfs_directional, find_nonfree_paths
from graph_utils import get_idx_from_idc, get_idc_from_idx
from collections import OrderedDict, defaultdict

def preprocess(graph, priority_queue):
    need_rotate = [False] * len(priority_queue)
    while sum(need_rotate) < len(priority_queue):
        for i, rotate_ion in enumerate(priority_queue):
            pz_name = priority_queue[rotate_ion]
            pz_parking_edge = [pz.parking_edge for pz in graph.pzs if pz.name == pz_name][0]
            edge_idc = graph.state[rotate_ion]
            next_edge = find_next_edge(graph, edge_idc, pz_parking_edge)
            next_edge = tuple(sorted((next_edge[0], next_edge[1]), key=sum))

            state_edges_idc = get_edge_state(graph)

            # if next_edge is free, not a stop move (same edge) and not at junction node
            # move the ion to the next edge
            if (
                have_common_junction_node(graph, edge_idc, next_edge) is False
                and state_edges_idc[next_edge] == []
                and edge_idc != next_edge
            ):
                graph.edges[edge_idc]["ions"].remove(rotate_ion)
                graph.edges[next_edge]["ions"].append(rotate_ion)
                print(f"preprocess move {rotate_ion} from {edge_idc} to {next_edge}")
            else:
                need_rotate[i] = True


# def create_general_priority_queue(graph, sequence, max_length=10):
#     logging.debug("Entered create_priority_queue function")
#     # TODO create real flat unique sequence?
#     unique_sequence = []
#     for seq_elem in sequence:
#         for elem in seq_elem:
#             if elem not in unique_sequence:
#                 unique_sequence.append(elem)
#                 if len(unique_sequence) == max_length:
#                     break
#     logging.debug(f"Priority queue created: {unique_sequence}")
#     return unique_sequence


def get_edge_idc_by_pz_name(graph, pz_name):
    for pz in graph.pzs:
        if pz.name == pz_name:
            return pz.parking_edge
    raise ValueError(f"Processing zone with name {pz_name} not found.")


def pick_pz_for_2_q_gate(graph, ion0, ion1):
    # TODO implement a better way to pick the processing zone for 2-qubit gates
    # pick the processing zone that both ions are closest to (so sum of distances is minimal)
    min_distance = float("inf")
    for pz_name in [graph.map_to_pz[ion0], graph.map_to_pz[ion1]]:
        pz_parking_edge = get_edge_idc_by_pz_name(graph, pz_name)
        distance = len(
            find_path_edge_to_edge(graph, graph.state[ion0], pz_parking_edge)
        ) + len(find_path_edge_to_edge(graph, graph.state[ion1], pz_parking_edge))
        if distance < min_distance:
            min_distance = distance
            closest_pz = pz_name
    return closest_pz


def create_priority_queue(graph, pz_executing_gate_order, max_length=10):
    """
    Create a priority queue based on a given graph and sequence of gates.
    Also creates a dictionary of the next gate of each processing zone.

    Args:
        graph (Graph): The graph representing the QCCD architecture.
        sequence (list): The sequence of gates.
        max_length (int, optional): The maximum length of the priority queue.
        Defaults to 10.

    Returns:
        tuple: A tuple containing the priority queue and
        the next gate at each processing zone.
    """

    unique_sequence = OrderedDict()
    graph.next_gate_at_pz = {}
    for seq_elem in graph.sequence:
        # 1-qubit gate
        if len(seq_elem) == 1:
            elem = seq_elem[0]

            # add first gate of pz to next_gate_at_pz (if not already there)
            if graph.map_to_pz[elem] not in graph.next_gate_at_pz or graph.next_gate_at_pz[graph.map_to_pz[elem]] == []:
                graph.next_gate_at_pz[graph.map_to_pz[elem]] = seq_elem

            # add ion to unique_sequence
            if elem not in unique_sequence.keys():
                unique_sequence[elem] = graph.map_to_pz[elem]
                if len(unique_sequence) == max_length:
                    break

        # 2-qubit gate
        elif len(seq_elem) == 2:
            if seq_elem not in graph.locked_gates:
                # pick processing zone for 2-qubit gate
                pz_for_2_q_gate = pick_pz_for_2_q_gate(
                    graph, seq_elem[0], seq_elem[1]
                )
                graph.locked_gates[seq_elem] = pz_for_2_q_gate
            else:
                pz_for_2_q_gate = graph.locked_gates[seq_elem]

            # add first gate of pz to next_gate_at_pz (if not already there)
            if pz_for_2_q_gate not in graph.next_gate_at_pz or graph.next_gate_at_pz[pz_for_2_q_gate] == []:
                graph.next_gate_at_pz[pz_for_2_q_gate] = seq_elem
                # lock the processing zone for the 2-qubit gate for later iterations
                # otherwise maybe pz changes if both move in a way, that favors a new pz
                # -> could result in a bug, if the very next iterations
                # change state back to old pz
                graph.locked_gates[seq_elem] = pz_for_2_q_gate

            # add ions to unique_sequence
            for elem in seq_elem:
                if elem not in unique_sequence.keys():
                    unique_sequence[elem] = pz_for_2_q_gate
                    if len(unique_sequence) == max_length:
                        break
        else:
            raise ValueError("len gate 0 or > 2? - can only process 1 or 2-qubit gates")

        # at the end fill all empty pzs with []
        for pz in graph.pzs:
            if pz.name not in graph.next_gate_at_pz:
                graph.next_gate_at_pz[pz.name] = []

    # NEW: add ions in connections to priority queue as below for in move_list
    ions_edges = get_ions(graph)
    
    # order pz according to gates
    # First, extract the unique order from pz_gate_order
    ordered_keys = list(dict.fromkeys(pz_executing_gate_order))

    # Then, add any pzs from graph.pzs that are not already in ordered_keys
    ordered_pzs = ordered_keys + [pz.name for pz in graph.pzs if pz.name not in ordered_keys]

    for pz_name in ordered_pzs:
        pz = graph.pzs_name_map[pz_name]
        
        # EXIT
        # NEW: add ions in exit connections to priority queue as below for in move_list (need both?)
        ions_in_exit_connections = []
        for edge_idx in pz.path_to_pz_idxs:  # order follows path_to_pz_idxs
            for ion, ion_edge_idc in ions_edges.items():
                ion_edge_idx = get_idx_from_idc(graph.idc_dict, ion_edge_idc)
                if ion_edge_idx == edge_idx:
                    ions_in_exit_connections.append(ion)

        if len(ions_in_exit_connections) > 0:
            for ion in ions_in_exit_connections:
                with contextlib.suppress(Exception):
                    unique_sequence.remove(ion)
                unique_sequence[ion] = pz.name
                unique_sequence.move_to_end(ion, last=False)

        # ENTRY
        # get ions in all entry edges and place in front
        # ion in entry must move out TODO for multiple pzs could be the case that out of entry moves block each other -> can't move out of some pzs -> need blocks in these pzs?
        ions_in_entry_connections = []
        for edge_idx in pz.path_from_pz_idxs:  # order follows path_from_pz_idxs
            for ion, ion_edge_idc in ions_edges.items():
                ion_edge_idx = get_idx_from_idc(graph.idc_dict, ion_edge_idc)
                if ion_edge_idx == edge_idx:
                    ions_in_entry_connections.append(ion)

        if len(ions_in_entry_connections) > 0:
            for ion in ions_in_entry_connections:
                with contextlib.suppress(Exception):
                    unique_sequence.remove(ion)
                unique_sequence[ion] = pz.name
                unique_sequence.move_to_end(ion, last=False)
    
    return unique_sequence, graph.next_gate_at_pz


def get_partitioned_priority_queues(graph, priority_queue, partition):
    # partitioned_priority_queue is a dictionary with pzs as keys and ions as values
    # represents priority queue for each processing zone individually
    # is just priority reversed (keys and values exchanged)
    partitioned_priority_queues = defaultdict(list)
    for key, value in priority_queue.items():
        # Append the key to the list of the corresponding value
        partitioned_priority_queues[value].append(key)

    # partitioned_priority_queues = {}
    # for pz in graph.pzs:
    #     partitioned_priority_queues[pz.name] = [
    #         elem for elem in priority_queue if elem in partition[pz.name]
    #     ]
    return partitioned_priority_queues


def create_gate_info_list(graph):
    # create list of next gate at each processing zone
    gate_info_list = {pz.name: [] for pz in graph.pzs}
    for seq_elem in graph.sequence:
        if len(seq_elem) == 1:
            elem = seq_elem[0]
            pz = graph.map_to_pz[elem]
            if gate_info_list[pz] == []:
                gate_info_list[pz].append(elem)
        elif len(seq_elem) == 2:
            # only pick processing zone for 2-qubit gate if not already locked -> then lock it (add to locked_gates)
            if seq_elem not in graph.locked_gates:
                pz = pick_pz_for_2_q_gate(graph, seq_elem[0], seq_elem[1])
                graph.locked_gates[seq_elem] = pz
            else:
                pz = graph.locked_gates[seq_elem]
            if gate_info_list[pz] == []:
                gate_info_list[pz].append(seq_elem[0])
                gate_info_list[pz].append(seq_elem[1])
        # break if all pzs have a gate
        if all(gate_info_list.values()):
            break
    return gate_info_list


def create_move_list(graph, partitioned_priority_queue, pz):
    """
    Create a move list based on a given graph and partitioned priority queue.
    Move list is specific to a processing zone.

    That is also why the single processing zone logic can be implemented here (entry move has path_length 0, since it is the first ion in priority queue).
    """
    ion_chains = get_ions(graph)
    path_length_sequence = {}
    move_list = []
    for i, rotate_ion in enumerate(partitioned_priority_queue):
        edge_idc = ion_chains[rotate_ion]
        # shortest path is also 1 edge if already at pz -> set to 0
        # if edge_idc == pz.parking_edge:
        #     path_length_sequence[rotate_ion] = 0

        # parking, exit and entry moves are set to 0, since they should always move and not block other moves (and will be placed in front of move_list below)
        if graph.get_edge_data(edge_idc[0], edge_idc[1])['edge_type'] in ['exit', 'entry', 'first_entry_connection', 'parking_edge']:
            path_length_sequence[rotate_ion] = 0
        else:
            path_to_go = find_path_edge_to_edge(graph, edge_idc, pz.parking_edge, exclude_first_entry_connection=True)
            path_length_sequence[rotate_ion] = len(path_to_go)

        # if first ion or all paths are 0 (all ions to move are in pz already) or current path is longer than all other paths 
        # or ion in exit (always pushes through now, to simplify movement in pz, but may be not optimal - may move more important ion out of pz)
        if (
            i == 0
            or sum(path_length_sequence.values()) == 0
            or sum(
                np.array([path_length_sequence[rotate_ion]] * len(move_list))
                > np.array([path_length_sequence[ion] for ion in move_list])
            )
            == len(move_list)
            or graph.get_edge_data(edge_idc[0], edge_idc[1])['edge_type'] == 'exit'
        ):
            move_list.append(rotate_ion)
    
    # same thing done in priority queue (needed later while looping over priority queue in find_movable_cycles so now also added in priority queue)
    
    ions_edges = get_ions(graph)
    # exit logic from single shuttler (now also needed since ions may be pushed into exit, then wait and others move out of pz but are in front of it in move_list)
    # -> for rare case that dag changed -> chain in exit connection was placed in front
    # -> overwrote other chain in exit connection which was still in move list but now later than the one that was inserted in front
    chains_in_exit_connections = []
    for edge_idx in pz.path_to_pz_idxs:  # Ensure order follows path_to_pz_idxs
        for ion, ion_edge_idc in ions_edges.items():
            ion_edge_idx = get_idx_from_idc(graph.idc_dict, ion_edge_idc)
            if ion_edge_idx == edge_idx:
                chains_in_exit_connections.append(ion)

    if len(chains_in_exit_connections) > 0:
        for ion in chains_in_exit_connections:
            with contextlib.suppress(Exception):
                move_list.remove(ion)
            move_list = [ion, *move_list]
    
    # get chains in all entry edges and place in front
    # chain in entry must move out
    ions_in_entry_connections = []
    for edge_idx in pz.path_from_pz_idxs:  # Ensure order follows path_from_pz_idxs
        for ion, ion_edge_idc in ions_edges.items():
            ion_edge_idx = get_idx_from_idc(graph.idc_dict, ion_edge_idc)
            if ion_edge_idx == edge_idx:
                ions_in_entry_connections.append(ion)

    if len(ions_in_entry_connections) > 0:
        for ion in ions_in_entry_connections:
            with contextlib.suppress(Exception):
                move_list.remove(ion)
            move_list = [ion, *move_list]

    return move_list

def bfs_free_edge(graph, node, other_next_edges):
    state_idcs = get_ions(graph).values()
    state_idxs = [get_idx_from_idc(graph.idc_dict, state_idc) for state_idc in state_idcs]
    other_next_edges_idxs = [get_idx_from_idc(graph.idc_dict, next_edge) for next_edge in other_next_edges]
    for edge_idc in nx.edge_bfs(graph.mz_graph, node):
        if (
            get_idx_from_idc(graph.idc_dict, edge_idc) not in state_idxs
            and get_idx_from_idc(graph.idc_dict, edge_idc) not in other_next_edges_idxs
        ):
            return edge_idc
    return None

def create_cycles_for_moves(graph, move_list, prio_queue, cycle_or_paths, pz, other_next_edges=None):
    print('start with ', pz.name)
    pz.rotate_entry = False
    pz.ion_to_park = find_ion_in_edge(graph, pz.path_to_pz[-1])
    pz.ion_to_move_out_of_pz = None
    if len(find_ions_in_parking(graph, pz)) < pz.max_num_parking:# or gate_execution_finished: # TODO gate_execution_finished
        parking_open = True
    else:
        parking_open = False

    all_cycles = {}
    ion_chains = graph.state
    # TODO other_next_edges - would need to calculate other_next_edges for each pz before?

    in_and_into_exit_moves = {}
    for rotate_ion in move_list:
        edge_idc = ion_chains[rotate_ion]
        # change pz to pz of position (if ion is in a pz and moves to other, pz would otherwise be the target pz, need position pz for following logic)
        if graph.get_edge_data(edge_idc[0], edge_idc[1])['edge_type'] != 'trap':
            # save target pz
            target_pz = pz # for finding next edge
            edge_idx = get_idx_from_idc(graph.idc_dict, edge_idc)
            position_pz = graph.edge_to_pz_map[edge_idx]
        else:
            target_pz = pz
            position_pz = pz
            
        # move from entry to memory zone
        # following checks are done with ion position (edge_idc)
        if get_idx_from_idc(graph.idc_dict, edge_idc) in position_pz.path_from_pz_idxs:
            # if in entry (last edge of entry connections, connected to mz)
            if get_idx_from_idc(graph.idc_dict, edge_idc) == get_idx_from_idc(
                graph.idc_dict, position_pz.entry_edge
            ):  # in graph.pzgraph_creator.path_from_pz_idxs:
                other_next_edges = [] #TODO
                starting_search_node = position_pz.entry_node
                target_edge = bfs_free_edge(graph, starting_search_node, other_next_edges)
                # calc path to target edge (node path does not include target edge and in this case - also not the start node)
                start_node = [node for node in edge_idc if graph.nodes(data=True)[node]['node_type'] in ["entry_connection_node", "processing_zone_node"]][0]
                node_path = find_path_edge_to_edge(graph, edge_idc, target_edge, exclude_exit=True, exclude_first_entry_connection=False)
                end_node = [node for node in target_edge if node not in node_path][0]
                node_path = [start_node, *node_path, end_node]
                edge_path = []
                for edge in pairwise(node_path):
                    edge_path.append(edge)
                all_cycles[rotate_ion] = edge_path
            
            # if in entry connection -> move to next edge (is now done here instead of in find_next_edge)
            # this is basically an else here but need to loop over path_from_pz to find next edge
            for i, edge_idx in enumerate(position_pz.path_from_pz_idxs[:-1]):
                if get_idx_from_idc(graph.idc_dict, edge_idc) == edge_idx:
                    next_edge = get_idc_from_idx(graph.idc_dict, position_pz.path_from_pz_idxs[i + 1])
                    edge_idc, next_edge = find_ordered_edges(graph, edge_idc, next_edge)
                    all_cycles[rotate_ion] = [edge_idc, next_edge]
        
        # moves inside mz
        # and moves into pz (exit) and out of parking edge
        # following checks are done with next_edge
        else:
            next_edge = find_next_edge(graph, edge_idc, target_pz.parking_edge)
            edge_idc, next_edge = find_ordered_edges(graph, edge_idc, next_edge)
            print('next, ', rotate_ion, next_edge, target_pz.parking_edge)
            # moves in exit and into parking_edge (blocks if parking is full)
            if get_idx_from_idc(graph.idc_dict, next_edge) in [
                *position_pz.path_to_pz_idxs,
                get_idx_from_idc(graph.idc_dict, position_pz.parking_edge),
            ]:
                in_and_into_exit_moves[position_pz.name] = {rotate_ion: edge_idc}
                all_cycles[rotate_ion] = [edge_idc, next_edge]
                # block moves to pz if parking is full (now blocks if parking not open and ion moving in exit and its next edge is in state_idxs)
                if parking_open is False and (get_idx_from_idc(graph.idc_dict, next_edge) in get_state_idxs(graph).values()):
                    print(f'blocked ion {rotate_ion} since parking open {parking_open} and {(get_idx_from_idc(graph.idc_dict, next_edge) in get_state_idxs(graph).values())}')
                    all_cycles[rotate_ion] = [edge_idc, edge_idc]

            # covers all shuttling in memory zone (does not cover entry connections anymore, see above)
            elif not check_if_edge_is_filled(graph, next_edge) or edge_idc == next_edge or get_idx_from_idc(graph.idc_dict, edge_idc) == get_idx_from_idc(graph.idc_dict, position_pz.parking_edge):
                all_cycles[rotate_ion] = [edge_idc, next_edge]
            else:
                
                if cycle_or_paths == "Cycles":
                    all_cycles[rotate_ion] = create_cycle(graph, edge_idc, next_edge)
                else:
                    all_cycles[rotate_ion] = create_path_via_bfs_directional(graph, edge_idc, next_edge, other_next_edges=None) #TODO other next edges for pz

    return all_cycles, in_and_into_exit_moves

def update_entry_and_exit_cycles(graph, pz, all_cycles, in_and_into_exit_moves_pz, out_of_entry_moves_pz, prio_queue):
    # move ion out of parking edge if needed
    ions_in_parking = find_ions_in_parking(graph, pz)
    ion_in_entry = find_ion_in_edge(graph, pz.entry_edge)

    # always create an out of entry move (need to coordinate them since multiple now - can't just use next free edge)
    if ion_in_entry is not None:
        all_cycles[ion_in_entry] = out_of_entry_moves_pz

    # if pz full and no ion is moving out (not in state_idxs entry edge) but ion is moving in
    if len(ions_in_parking) >= pz.max_num_parking and pz.ion_to_park is not None:
        print('ion to park', pz.ion_to_park)
        # if gate finished -> new space could be in parking
        print('gate finished =', pz.gate_execution_finished)
        if pz.gate_execution_finished:
            # find least important chain in parking edge
            pz.ion_to_move_out_of_pz = find_least_import_ion_in_parking(
                prio_queue, [*ions_in_parking]#, pz.ion_to_park]    # ion to park now always moves to pz (simplifies the pz shuttling, may be not optimal though), for that also changed: ion in exit always in move_list
            )
            print('ion to move out', pz.ion_to_move_out_of_pz)
            # if ion moves out of parking edge -> check if ion in entry edge -> make sure it can move into MZ -> then move ion in parking to entry (if not, stop moves in exit) 
            # now only if the ion moving to pz is higher in prio queue or ion moving out not in prio queue (-> new logic of finding least important)
            # and new addition: if ion is moving in and can't move to pz because of the other clauses here -> if then no gate is executable with pz ions + ion moving in -> move ion in and least important out (worst case but to avoid complete blockages)
            if pz.ion_to_move_out_of_pz != pz.ion_to_park and (
                ion_in_entry is None or (ion_in_entry is not None and out_of_entry_moves_pz is not None)
                    ) and ((pz.ion_to_move_out_of_pz not in prio_queue or (pz.ion_to_park in prio_queue and prio_queue.index(pz.ion_to_park) < prio_queue.index(pz.ion_to_move_out_of_pz)))
                        ) or any(ion not in [pz.ion_to_park, *ions_in_parking] for ion in graph.next_gate_at_pz[pz.name]):
                # move it to entry (later through rotate_entry flag in rotate_free_cycles)
                pz.rotate_entry = True
                # change its path/circle to a stop move -> will be later placed into entry
                all_cycles[pz.ion_to_move_out_of_pz] = [
                    pz.path_from_pz[0],
                    pz.path_from_pz[0],
                ]
                if ion_in_entry is not None:
                    print('changed path of ion in entry: ', ion_in_entry)
                    print('from: ', all_cycles[ion_in_entry])
                    all_cycles[ion_in_entry] = out_of_entry_moves_pz
                    print('to: ', all_cycles[ion_in_entry])


            # TODO 3
            # # else -> pz.ion_to_park not needed right now
            # # if new gate can be executed -> bring everything to a stop in exit
            # elif pz.new_gate_starting:
            #     # if chain to park should move to entry -> all chains with next edge in exit -> stop move
            #     for chain, edge_idc in in_and_into_exit_moves.items():
            #         all_cycles[chain] = [edge_idc, edge_idc]
            #     # maybe already covered above
            #     all_cycles[pz.ion_to_move_out_of_pz] = (
            #         pz.path_to_pz[-1],
            #         pz.path_to_pz[-1],
            #     )

            # else -> no new gate possible with only parking chains
            # -> but also chain can't move to parking since it is least important
            # -> should be rare edge case -> chain moves from exit to entry (but what if gate is in process then for 3 time steps? Is ion still blocked or moves to entry?)
            else:
                # TODO 3* also changed this (could be problem if dag changes? now ion to park waits, does that work?)
                # pz.rotate_entry = True
                # # change its path/circle to a stop move -> will be later placed into entry
                # all_cycles[pz.ion_to_move_out_of_pz] = [
                #     pz.path_from_pz[0],
                #     pz.path_from_pz[0],
                # ]
                # new stop move in exit? -> instead of moving to entry and not to pz, ion_to_park now waits
                # if this works -> should be better solution than whole thing above, but could be rare bug cases
                pz.rotate_entry = False
                all_cycles[pz.ion_to_park] = [
                    pz.path_to_pz[-1],
                    pz.path_to_pz[-1],
                ]
        else:
            # else bring everything to a stop in exit
            # same as above
            for chain, edge_idc in in_and_into_exit_moves_pz.items():   # into exit move? dont move into exit if not all in prio queue (partioned?) in front of you are in exit, exitconn or parking?
                #if graph.get_edge_data(edge_idc[0], edge_idc[1])['edge_type'] == 'trap' and 
                all_cycles[chain] = [edge_idc, edge_idc]
                

            # maybe already covered above
            # all_circles[pz.ion_to_move_out_of_pz] = (
            #     memorygrid.pzgraph_creator.path_to_pz[-1],
            #     memorygrid.pzgraph_creator.path_to_pz[-1],
            # )
    
    # remove ions if they would move to exit but their next gate is a 2-qubit gate and their gate is not yet in locked_gates (maybe not needed, TODO check)
    for ion, edge_idc in in_and_into_exit_moves_pz.items():
        if graph.get_edge_data(edge_idc[0], edge_idc[1])['edge_type'] == 'trap':    # if edge_idc is in MZ -> moves into exit
            for gate in graph.sequence:
                if ion in gate:
                    if len(gate) == 1:
                        break
                    elif len(gate) == 2:
                        print(f'checked if ion {ion} moving to exit will execute 2 qubit gate')
                        if gate in graph.locked_gates.keys():
                            print(f'{ion} can move to exit since gate was locked')
                            break
                        else:
                            all_cycles.pop(ion)
                            print(f'stopped ion {ion} from entering exit since gate {gate} is not locked yet')
                            break

    return all_cycles


def find_movable_cycles(graph, all_cycles, priority_queue, cycle_or_paths):
    print('\nall_cycles')
    for key, value in all_cycles.items():
        if isinstance(value, (int, float)):
            print(f"{key}: {value:.2f}")
        elif isinstance(value, (list, tuple)):
            # Format each number inside list/tuple, leave non-number items as-is
            formatted = []
            for x in value:
                if isinstance(x, (int, float)):
                    formatted.append(f"{x:.2f}")
                else:
                    formatted.append(x)
            print(f"{key}: {formatted}")
        else:
            # Just print non-number values as-is
            print(f"{key}: {value}")
            
    if cycle_or_paths == "Cycles":
        nonfree_cycles = find_conflict_cycle_idxs(graph, all_cycles)
    else:
        nonfree_cycles = find_nonfree_paths(graph, all_cycles)
    # start with first ion in priority queue
    free_cycle_seq_idxs = [list(priority_queue.keys())[0]]

    # check if ion can move while first ion is moving and so on
    for seq_cyc in list(priority_queue.keys())[1:]:
        # skip ion of priority_queue if it is not in all_cycles
        # -> was removed before in individual move_list
        if seq_cyc not in all_cycles.keys():
            continue
        nonfree = False
        for mov_cyc in free_cycle_seq_idxs:
            if (seq_cyc, mov_cyc) in nonfree_cycles or (
                mov_cyc,
                seq_cyc,
            ) in nonfree_cycles:
                nonfree = True
                break
        if nonfree is False:
            free_cycle_seq_idxs.append(seq_cyc)
    print(f"free_cycle_seq_idxs: {free_cycle_seq_idxs}")

    return free_cycle_seq_idxs


def rotate(graph, ion, cycle_idcs):
    # ion {ion} along cycle {cycle_idcs}")
    state_dict = get_edge_state(graph)
    first = True
    last_ion = -1
    for current_edge, new_edge in pairwise(cycle_idcs):
        current_edge = tuple(sorted(current_edge, key=sum))
        new_edge = tuple(sorted(new_edge, key=sum))

        # reverse edge if necessary
        if state_dict.get((current_edge[1], current_edge[0])) is not None:
            current_edge = (current_edge[1], current_edge[0])

        # if ion already rotated via previous cycle
        # (now checks directly in state_dict, in case two ions on one edge)
        if first:
            # first ion is ion to move (resolved bug with taking first ion of multiple in edge, should not be necessary anymore, see below)
            current_ion = ion
            if ion not in state_dict[current_edge]:  # current_ion != ion:
                print(f"Ion {ion} already rotated via previous cycle")
                return

        else:
            current_ion = state_dict.get(current_edge)
            
            # take first ion in list to rotate (should be unnecessary since multiple ions only in pz, which is handled elsewhere)
            try:
                current_ion = current_ion[0]
            except IndexError:
                pass

        first = False
        if (
            current_ion != []
            and current_ion != last_ion
            and current_ion not in graph.in_process
        ):  # and not ion in pz and needed in 2-qubit gate
            graph.edges[current_edge]["ions"].remove(current_ion)
            graph.edges[new_edge]["ions"].append(current_ion)
            print(f"Rotated ion {current_ion} from {current_edge} to {new_edge}")

        # save last ion so each ion only rotates once
        last_ion = current_ion


def rotate_free_cycles(graph, all_cycles, free_cycles_idxs):
    rotate_cycles_idcs = {}
    for cycle_ion in free_cycles_idxs:
        try:
            rotate_cycles_idcs[cycle_ion] = all_cycles[cycle_ion]
        except KeyError:
            pass
    # skip stop moves
    for ion, indiv_cycle_idcs in rotate_cycles_idcs.items():
        if len(indiv_cycle_idcs) == 2:
            if indiv_cycle_idcs[0] == indiv_cycle_idcs[1]:
                # cycle {indiv_cycle_idcs}, since it is a stop move")
                continue
            else:   # TODO das hier checkt ob ein out of parking move sowieso durch geht. pz.out_of_parking_move würde aber eher gebraucht werden -> check ob man pz.out_of_parking_move auch noch später nutzen könnte? Jetzt Lösung bei rotate_entry implementiert
                for pz in graph.pzs:
                    if get_idx_from_idc(graph.idc_dict, indiv_cycle_idcs[0]) == get_idx_from_idc(graph.idc_dict, pz.parking_edge
                    ) and get_idx_from_idc(graph.idc_dict, indiv_cycle_idcs[1]) == get_idx_from_idc(graph.idc_dict, pz.first_entry_connection_from_pz):
                        pz.out_of_parking_move = ion
        rotate(graph, ion, indiv_cycle_idcs)

    # extra clause if an ion is moving out of pz anyway (on its way to another pz) -> do not need to move out another ion
    for pz in graph.pzs:
        if pz.rotate_entry:    
            if pz.out_of_parking_move is not None and pz.ion_to_move_out_of_pz != pz.out_of_parking_move:
                    print(f'ion {ion} does not move out of pz since {pz.out_of_parking_move} is moving out anyway')
            else:
                print(f'moving {pz.ion_to_move_out_of_pz} out of {pz.name}')
                graph.edges[pz.parking_edge]["ions"].remove(pz.ion_to_move_out_of_pz)
                graph.edges[pz.path_from_pz[0]]["ions"].append(pz.ion_to_move_out_of_pz)
                
def find_out_of_entry_moves(graph, other_next_edges):
    free_edges = {}
    for pz in graph.pzs:
        if check_if_edge_is_filled(graph, pz.entry_edge):
            # if ion in entry edge -> find space in memory zone
            free_edges[pz] = bfs_free_edge(graph, pz.entry_node, other_next_edges=other_next_edges)
            other_next_edges.append(free_edges[pz])
    # TODO implement priorization of out of entry moves

    out_of_entry_moves = {}
    mz_graph_copy = graph.mz_graph.copy()
    for pz, edge in free_edges.items():
        # start from entry node (=not entry_connection_node or processing_zone_node) -> needed this way since entry node can be exit node of other pz TODO check if exit/entry of two pzs cause other problems
        start_node = pz.entry_edge[1] if graph.nodes[pz.entry_edge[0]]['node_type'] in ["entry_connection_node", "processing_zone_node"] else pz.entry_edge[0]
        other_start_node = pz.entry_edge[1] if start_node == pz.entry_edge[0] else pz.entry_edge[0]

        node_path = find_path_node_to_edge(mz_graph_copy, start_node, edge)
        if node_path == None:
            out_of_entry_moves[pz] = ((other_start_node, start_node), (other_start_node, start_node))
        else:
            other_end_node = [node for node in edge if node not in node_path][0]
            node_path = [other_start_node, *node_path, other_end_node]
            edge_path = []
            for edge in pairwise(node_path):
                edge_path.append(edge)
            out_of_entry_moves[pz] = edge_path

            # # remove edges from mz_graph_copy -> other pzs can't use the same junctions to move ions out of entry
            # for edge in edge_path[1:]:
            #     mz_graph_copy.remove_edge(edge[0], edge[1])
            # remove nodes from mz_graph_copy -> other pzs can't use the same junctions to move ions out of entry
            for node in node_path[1:-1]:
                mz_graph_copy.remove_node(node)

    return out_of_entry_moves
    
        
            




