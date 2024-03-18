import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
from more_itertools import distinct_combinations, pairwise


# create dictionary to swap from idx to idc and vice versa
def create_idc_dictionary(nx_g):
    edge_dict = {}
    for edge_idx, edge_idc in enumerate(nx_g.edges()):
        edge_dict[edge_idx] = tuple(sorted(edge_idc, key=sum))
    return edge_dict


def get_idx_from_idc(edge_dictionary, idc):
    idc = tuple(sorted(idc, key=sum))
    return list(edge_dictionary.values()).index(idc)


def get_idc_from_idx(edge_dictionary, idx):
    return edge_dictionary[idx]


def get_path_to_node(nx_g, src, tar, exclude_exit=False, exclude_first_entry_connection=True):
    edge_path = []
    if exclude_first_entry_connection is True:
        # lambda function to give path over processing zone huge weight -> doesn't take that path if not necessary - now only encludes entry edge -> can use exit (in MemGrid was != trap before and then to exit node -> not PZ node)
        node_path = nx.shortest_path(
            nx_g,
            src,
            tar,
            lambda _, __, edge_attr_dict: (edge_attr_dict["edge_type"] == "first_entry_connection") * 1e8 + 1,
        )
        # also exclude exit edge if necessary
        if exclude_exit is True:
            node_path = nx.shortest_path(
                nx_g,
                src,
                tar,
                lambda _, __, edge_attr_dict: (edge_attr_dict["edge_type"] in ("first_entry_connection", "exit")) * 1e8
                + 1,
            )

    # only exclude exit edge
    elif exclude_exit is True:
        node_path = nx.shortest_path(
            nx_g, src, tar, lambda _, __, edge_attr_dict: (edge_attr_dict["edge_type"] == "exit") * 1e8 + 1
        )

    else:
        node_path = nx.shortest_path(nx_g, src, tar)
    # shortest path should always be the correct path in a grid -> care for changes

    for edge in pairwise(node_path):
        edge_path.append(edge)

    return edge_path


def calc_dist_to_pz(nx_g_creator, edge_idx):
    edge_idc = get_idc_from_idx(nx_g_creator.idc_dict, edge_idx)
    node1, node2 = edge_idc[0], edge_idc[1]

    path1 = get_path_to_node(
        nx_g_creator.networkx_graph, node1, nx_g_creator.processing_zone, exclude_first_entry_connection=True
    )
    path2 = get_path_to_node(
        nx_g_creator.networkx_graph, node2, nx_g_creator.processing_zone, exclude_first_entry_connection=True
    )
    if edge_idx == get_idx_from_idc(nx_g_creator.idc_dict, nx_g_creator.parking_edge):
        return 0
    if edge_idx == get_idx_from_idc(nx_g_creator.idc_dict, nx_g_creator.first_entry_connection_from_pz):
        return max(len(path1), len(path2)) + 1
    return min(len(path1), len(path2)) + 1


# def circle_is_contained_in_other_circle(subseq, seq):
#     subseq_len = len(subseq)
#     seq_len = len(seq)

#     if subseq_len > seq_len:
#         return False

#     # Duplicate the sequence to handle wrap-around cases (remove last element, circles constructed with first element added at the end)
#     seq_extended = seq[:-1] + seq
#     return any(seq_extended[i : i + subseq_len] == subseq for i in range(seq_len))


class MZGraphCreator:
    def __init__(self, m, n, ion_chain_size_vertical, ion_chain_size_horizontal):
        self.m = m
        self.n = n
        self.ion_chain_size_vertical = ion_chain_size_vertical
        self.ion_chain_size_horizontal = ion_chain_size_horizontal
        self.networkx_graph = self.create_graph()

        self.idc_dict = create_idc_dictionary(self.networkx_graph)

    def create_graph(self):
        self.m_extended = self.m + (self.ion_chain_size_vertical - 1) * (self.m - 1)
        self.n_extended = self.n + (self.ion_chain_size_horizontal - 1) * (self.n - 1)

        networkx_graph = nx.grid_2d_graph(self.m_extended, self.n_extended)
        self._set_trap_nodes(networkx_graph)
        self._remove_horizontal_edges(networkx_graph)
        self._remove_vertical_edges(networkx_graph)
        self._remove_horizontal_nodes(networkx_graph)
        self._set_junction_nodes(networkx_graph)
        nx.set_edge_attributes(networkx_graph, "trap", "edge_type")

        return networkx_graph

    def _set_trap_nodes(self, networkx_graph):
        for node in networkx_graph.nodes():
            networkx_graph.add_node(node, node_type="trap_node", color="b")

    def _remove_horizontal_edges(self, networkx_graph):
        for i in range(0, self.m_extended - self.ion_chain_size_vertical, self.ion_chain_size_vertical):
            for k in range(1, self.ion_chain_size_vertical):
                for j in range(self.n_extended - 1):
                    networkx_graph.remove_edge((i + k, j), (i + k, j + 1))

    def _remove_vertical_edges(self, networkx_graph):
        for i in range(0, self.n_extended - self.ion_chain_size_horizontal, self.ion_chain_size_horizontal):
            for k in range(1, self.ion_chain_size_horizontal):
                for j in range(self.m_extended - 1):
                    networkx_graph.remove_edge((j, i + k), (j + 1, i + k))

    def _remove_horizontal_nodes(self, networkx_graph):
        for i in range(0, self.m_extended - self.ion_chain_size_vertical, self.ion_chain_size_vertical):
            for k in range(1, self.ion_chain_size_vertical):
                for j in range(0, self.n_extended - self.ion_chain_size_horizontal, self.ion_chain_size_horizontal):
                    for s in range(1, self.ion_chain_size_horizontal):
                        networkx_graph.remove_node((i + k, j + s))

    def _set_junction_nodes(self, networkx_graph):
        for i in range(0, self.m_extended, self.ion_chain_size_vertical):
            for j in range(0, self.n_extended, self.ion_chain_size_horizontal):
                networkx_graph.add_node((i, j), node_type="junction_node", color="g")

    def get_graph(self):
        return self.networkx_graph


class GraphCreator:
    def __init__(self, m, n, ion_chain_size_vertical, ion_chain_size_horizontal):
        self.m = m
        self.n = n
        self.ion_chain_size_vertical = ion_chain_size_vertical
        self.ion_chain_size_horizontal = ion_chain_size_horizontal
        self.networkx_graph = self.create_graph()

        self.idc_dict = create_idc_dictionary(self.networkx_graph)
        self.path_to_pz_idxs = [get_idx_from_idc(self.idc_dict, edge) for edge in self.path_to_pz]
        self.path_from_pz_idxs = [get_idx_from_idc(self.idc_dict, edge) for edge in self.path_from_pz]

        # create lookup dictionaries for rest of path to and from processing zone
        self.rest_of_path_to_pz = {edge: self.path_to_pz[i + 1 :] for i, edge in enumerate(self.path_to_pz)}
        self.rest_of_path_from_pz = {edge: self.path_from_pz[i + 1 :] for i, edge in enumerate(self.path_from_pz)}

        self.pz_edges_idx = [
            get_idx_from_idc(self.idc_dict, edge)
            for edge in self.networkx_graph.edges()
            if nx.get_edge_attributes(self.networkx_graph, "edge_type")[edge] != "trap"
        ]

    def create_graph(self):
        self.m_extended = self.m + (self.ion_chain_size_vertical - 1) * (self.m - 1)
        self.n_extended = self.n + (self.ion_chain_size_horizontal - 1) * (self.n - 1)
        self.num_edges = self.n // 2

        networkx_graph = nx.grid_2d_graph(self.m_extended, self.n_extended)
        self._set_trap_nodes(networkx_graph)
        self._remove_horizontal_edges(networkx_graph)
        self._remove_vertical_edges(networkx_graph)
        self._remove_horizontal_nodes(networkx_graph)
        self._set_junction_nodes(networkx_graph)
        nx.set_edge_attributes(networkx_graph, "trap", "edge_type")
        self._set_processing_zone(networkx_graph)

        return networkx_graph

    def _set_trap_nodes(self, networkx_graph):
        for node in networkx_graph.nodes():
            networkx_graph.add_node(node, node_type="trap_node", color="b")

    def _remove_horizontal_edges(self, networkx_graph):
        for i in range(0, self.m_extended - self.ion_chain_size_vertical, self.ion_chain_size_vertical):
            for k in range(1, self.ion_chain_size_vertical):
                for j in range(self.n_extended - 1):
                    networkx_graph.remove_edge((i + k, j), (i + k, j + 1))

    def _remove_vertical_edges(self, networkx_graph):
        for i in range(0, self.n_extended - self.ion_chain_size_horizontal, self.ion_chain_size_horizontal):
            for k in range(1, self.ion_chain_size_horizontal):
                for j in range(self.m_extended - 1):
                    networkx_graph.remove_edge((j, i + k), (j + 1, i + k))

    def _remove_horizontal_nodes(self, networkx_graph):
        for i in range(0, self.m_extended - self.ion_chain_size_vertical, self.ion_chain_size_vertical):
            for k in range(1, self.ion_chain_size_vertical):
                for j in range(0, self.n_extended - self.ion_chain_size_horizontal, self.ion_chain_size_horizontal):
                    for s in range(1, self.ion_chain_size_horizontal):
                        networkx_graph.remove_node((i + k, j + s))

    def _set_junction_nodes(self, networkx_graph):
        for i in range(0, self.m_extended, self.ion_chain_size_vertical):
            for j in range(0, self.n_extended, self.ion_chain_size_horizontal):
                networkx_graph.add_node((i, j), node_type="junction_node", color="g")

    def _set_processing_zone(self, networkx_graph):
        # Define the key nodes
        self.exit = (self.m_extended - 1, self.n_extended - 1)
        self.processing_zone = (self.m_extended + self.num_edges - 1, self.n_extended + self.num_edges - 1)
        self.entry = (self.m_extended - 1, 0)
        self.parking_node = (self.processing_zone[0] + 1, self.processing_zone[1])
        self.parking_edge = (self.processing_zone, self.parking_node)

        # differences
        dy_exit = self.exit[1] - self.processing_zone[1]
        dy_entry = self.processing_zone[1] - self.entry[1]

        self.path_to_pz = []
        self.path_from_pz = []

        # Add exit edges
        for i in range(self.num_edges):
            exit_node = (self.exit[0] + (i + 1), self.exit[1] - (i + 1) * dy_exit / self.num_edges)

            if i == 0:
                networkx_graph.add_node(exit_node, node_type="exit_node", color="y")
                previous_exit_node = self.exit
                self.exit_edge = (previous_exit_node, exit_node)

            networkx_graph.add_node(exit_node, node_type="exit_connection_node", color="y")
            networkx_graph.add_edge(previous_exit_node, exit_node, edge_type="exit", color="k")
            self.path_to_pz.append((previous_exit_node, exit_node))
            previous_exit_node = exit_node

        # Add entry edges
        for i in range(self.num_edges):
            entry_node = (self.entry[0] + (i + 1), self.entry[1] + (i + 1) * dy_entry / self.num_edges)

            if i == 0:
                networkx_graph.add_node(entry_node, node_type="entry_node", color="orange")
                previous_entry_node = self.entry
                self.entry_edge = (previous_entry_node, entry_node)

            networkx_graph.add_node(entry_node, node_type="entry_connection_node", color="orange")
            # first entry connection is first edge after pz
            # entry is edge connected to memory grid, so last entry connection
            # if entry is one edge only -> first entry connection is the same as entry edge
            if entry_node == self.processing_zone:
                self.first_entry_connection_from_pz = (entry_node, previous_entry_node)
                networkx_graph.add_edge(previous_entry_node, entry_node, edge_type="first_entry_connection", color="k")
            else:
                networkx_graph.add_edge(previous_entry_node, entry_node, edge_type="entry", color="k")
            self.path_from_pz.insert(0, (entry_node, previous_entry_node))

            previous_entry_node = entry_node

        assert exit_node == entry_node, "Exit and entry do not end in same node"
        assert exit_node == self.processing_zone, "Exit and entry do not end in processing zone"

        # Add the processing zone node
        networkx_graph.add_node(self.processing_zone, node_type="processing_zone_node", color="r")

        # new parking edge
        networkx_graph.add_node(self.parking_node, node_type="parking_node", color="r")
        networkx_graph.add_edge(self.parking_edge[0], self.parking_edge[1], edge_type="parking_edge", color="g")

    def get_graph(self):
        return self.networkx_graph

    # plotting
    def plot_state(self, ion_moves, labels, plot_ions=True, show_plot=False, save_plot=False, filename=""):
        # idc_dict = create_idc_dicitonary(nx_G)
        pos = {(x, y): (y, -x) for i, (x, y) in enumerate(list(self.networkx_graph.nodes()))}
        if plot_ions is True:
            pass
            # edge_labels = nx.get_edge_attributes(self.networkx_graph,'ion_chain')
        else:
            edge_labels = {}
            for idc in self.networkx_graph.edges():
                # pass
                edge_labels[idc] = "$e_{%s}$" % get_idx_from_idc(self.idc_dict, idc)

        for edge_idc in self.networkx_graph.edges():
            # color all edges black
            self.networkx_graph.add_edge(edge_idc[0], edge_idc[1], color="k")

            ion_holder = {}
            colors = []
            np.random.seed(0)
            for _ in range(len(ion_moves)):
                r = np.round(np.random.rand(), 1)
                g = np.round(np.random.rand(), 1)
                b = np.round(np.random.rand(), 1)

                colors.append((r, g, b))
            np.random.seed()

            for i, ion_place in enumerate(ion_moves):
                ion_edge_idc = get_idc_from_idx(self.idc_dict, ion_place)
                try:
                    ion_holder[ion_place].append(i)
                except KeyError:
                    ion_holder[ion_place] = [i]
            for i, ion_place in enumerate(ion_moves):
                ion_edge_idc = get_idc_from_idx(self.idc_dict, ion_place)
                self.networkx_graph.add_edge(
                    ion_edge_idc[0], ion_edge_idc[1], ion_chain=ion_holder[ion_place], color=colors[i]
                )

        edge_color = nx.get_edge_attributes(self.networkx_graph, "color").values()
        node_color = list(nx.get_node_attributes(self.networkx_graph, "color").values())
        edge_labels = nx.get_edge_attributes(self.networkx_graph, "ion_chain")

        # plt.figure(figsize=(25, 15))
        plt.figure(figsize=(15, 8))
        nx.draw_networkx(
            self.networkx_graph,
            pos=pos,
            with_labels=True,
            node_size=300,
            node_color=node_color,
            width=8,
            edge_color=edge_color,
            font_size=6,
        )
        nx.draw_networkx_edge_labels(self.networkx_graph, pos, edge_labels)

        # reset edge labels
        for i, ion in enumerate(ion_moves):
            ion_edge_idc = get_idc_from_idx(self.idc_dict, ion)
            self.networkx_graph.add_edge(ion_edge_idc[0], ion_edge_idc[1], ion_chain="", color=colors[i])

        labels0, labels1 = labels
        plt.plot([], [], label=labels0)
        plt.plot([], [], label=labels1)
        plt.legend()

        if show_plot is True:
            plt.show()

        if save_plot is True:
            plt.savefig(filename)
        plt.close()


class MemoryZone:
    def __init__(
        self,
        m,
        n,
        v,
        h,
        starting_config,
        max_timestep,
        max_num_parking,
        time_2qubit_gate=2,
        time_1qubit_gate=1,
    ):
        # new graph MZ
        self.mz_Graph_creator = MZGraphCreator(m, n, v, h)
        self.mz_graph = self.mz_Graph_creator.get_graph()

        self.graph_creator = GraphCreator(m, n, v, h)
        self.graph = self.graph_creator.get_graph()
        self.starting_config = starting_config
        self.max_timestep = max_timestep
        self.max_num_parking = max_num_parking
        self.time_2qubit_gate = time_2qubit_gate
        self.time_1qubit_gate = time_1qubit_gate
        self.num_ion_chains = len(starting_config)
        self.idc_dict = self.graph_creator.idc_dict

        # create dictionary with all distances to entry
        self.dist_dict = {}
        for edge_idc in self.graph.edges():
            # keep node ordering consistent:
            edge_idx = get_idx_from_idc(self.idc_dict, edge_idc)
            self.dist_dict[get_idc_from_idx(self.idc_dict, edge_idx)] = calc_dist_to_pz(
                self.graph_creator, get_idx_from_idc(self.idc_dict, edge_idc)
            )

        # create dictionary with all distances to entry for all nodes
        self.dist_dict_nodes = {}
        for node in self.graph.nodes():
            self.dist_dict_nodes[node] = len(get_path_to_node(self.graph, node, self.graph_creator.processing_zone))

        # # create dictionary with all paths to entry (TODO artifact?)
        # self.path_dict = {}
        # for edge_idc in self.graph.edges():
        #     self.path_dict[edge_idc] = calc_dist_to_pz(self.graph_creator, get_idx_from_idc(self.idc_dict, edge_idc))

        self.ion_chains = self.starting_config.copy()

        self.junction_nodes = [
            node
            for node in self.graph.nodes()
            if (
                nx.get_node_attributes(self.graph, "node_type")[node]
                in ("junction_node", "exit_node", "exit_connection_node", "entry_node", "entry_connection_node")
            )
        ]

        self.path_entry_to_exit = get_path_to_node(
            self.graph, self.graph_creator.entry, self.graph_creator.exit, exclude_first_entry_connection=True
        )

        # precalulculate bfs for top left and exit
        self.bfs_top_left = nx.edge_bfs(self.mz_graph, (0, 0))
        self.bfs_exit = nx.edge_bfs(self.mz_graph, self.graph_creator.exit)

    # get edge idxs of ion chains
    def get_state_idxs(self):
        ion_chains_idx = []
        for chain in self.ion_chains.values():
            ion_chains_idx.append(get_idx_from_idc(self.idc_dict, chain))
        self.state_idxs = ion_chains_idx
        return ion_chains_idx

    # calc distance to parking edge for all ion chains
    def update_distance_map(self):
        self.distance_map = {}
        for ion_chain, edge_idx in enumerate(self.get_state_idxs()):
            self.distance_map[ion_chain] = self.dist_dict[get_idc_from_idx(self.idc_dict, edge_idx)]
        return self.distance_map

    def count_chains_in_pz(self):
        return len([chain_idx for chain_idx in self.get_state_idxs() if chain_idx in self.graph_creator.pz_edges_idx])

    def count_chains_in_exit(self):
        return len(
            [chain_idx for chain_idx in self.get_state_idxs() if chain_idx in self.graph_creator.path_to_pz_idxs]
        )

    def count_chains_in_parking(self):
        return len(
            [
                chain_idx
                for chain_idx in self.get_state_idxs()
                if chain_idx == get_idx_from_idc(self.idc_dict, self.graph_creator.parking_edge)
            ]
        )

    def find_chain_in_edge(self, edge_idc):
        chains = [
            ion
            for ion, chain_idx in enumerate(self.get_state_idxs())
            if chain_idx == get_idx_from_idc(self.idc_dict, edge_idc)
        ]
        assert (
            len(chains) <= 1
        ), f"more than one chain ({chains}) in edge {edge_idc}, if parking edge -> use find_chains_in_parking()"
        if len(chains) == 0:
            return None
        return chains[0]

    def find_chains_in_parking(self):
        return [
            ion
            for ion, chain_idx in enumerate(self.get_state_idxs())
            if chain_idx == get_idx_from_idc(self.idc_dict, self.graph_creator.parking_edge)
        ]

    def find_next_edge(self, edge_idc, towards=(0, 0)):
        ### find next edge given edge_idc of ion chain

        # if in entry connection -> move to next edge
        for i, edge_idx in enumerate(self.graph_creator.path_from_pz_idxs[:-1]):
            if get_idx_from_idc(self.idc_dict, edge_idc) == edge_idx:
                return get_idc_from_idx(self.idc_dict, self.graph_creator.path_from_pz_idxs[i + 1])

        if get_idx_from_idc(self.idc_dict, edge_idc) == get_idx_from_idc(self.idc_dict, self.graph_creator.entry_edge):
            if towards == (0, 0):
                next_edge = next(
                    edge
                    for edge in self.graph.edges(self.graph_creator.entry)
                    if edge not in (self.graph_creator.entry_edge, self.path_entry_to_exit[0])
                )
            elif towards == "exit":
                next_edge = self.path_entry_to_exit[0]
            else:
                msg = "towards must be (0,0) or 'exit'"
                raise ValueError(msg)

            # assert that next edge after entry is not entry or exit
            assert get_idx_from_idc(self.idc_dict, next_edge) != get_idx_from_idc(
                self.idc_dict, self.graph_creator.exit_edge
            )
            assert get_idx_from_idc(self.idc_dict, next_edge) != get_idx_from_idc(
                self.idc_dict, self.graph_creator.entry_edge
            )
            return next_edge

        # if in parking space
        if get_idx_from_idc(self.idc_dict, edge_idc) == get_idx_from_idc(
            self.idc_dict, self.graph_creator.parking_edge
        ):
            return self.graph_creator.parking_edge

        # find shortest path from both sides for all other edges
        path0 = nx.shortest_path(
            self.graph,
            edge_idc[0],
            self.graph_creator.parking_node,
            lambda _, __, edge_attr_dict: (edge_attr_dict["edge_type"] == "first_entry_connection") * 1e8 + 1,
        )
        path1 = nx.shortest_path(
            self.graph,
            edge_idc[1],
            self.graph_creator.parking_node,
            lambda _, __, edge_attr_dict: (edge_attr_dict["edge_type"] == "first_entry_connection") * 1e8 + 1,
        )

        # create chains in correct order -> chains are path from outer node to other outer node
        # start with outer node that is farther away from processing zone (longer path)
        # end with closer outer node (shorter path)
        if len(path1) < len(path0):
            next_edge = (path1[0], path1[1])
            # make next_edge_idc consistent
            next_edge_idx = get_idx_from_idc(self.idc_dict, next_edge)
            return get_idc_from_idx(self.idc_dict, next_edge_idx)

        next_edge = (path0[0], path0[1])
        # make next_edge_idc consistent
        next_edge_idx = get_idx_from_idc(self.idc_dict, next_edge)
        return get_idc_from_idx(self.idc_dict, next_edge_idx)

    def find_ordered_edges(self, edge1, edge2):
        # Find the common node shared between the two edges
        common_node = set(edge1).intersection(set(edge2))

        if len(common_node) != 1 and edge1 != edge2:
            msg = f"The input edges are not connected. Edges: {edge1}, {edge2}"
            raise ValueError(msg)

        common_node = common_node.pop()
        if edge1[0] == common_node:
            edge1_in_order = (edge1[1], common_node)
            edge2_in_order = (common_node, edge2[1]) if edge2[0] == common_node else (common_node, edge2[0])
        else:
            edge1_in_order = (edge1[0], common_node)
            edge2_in_order = (common_node, edge2[1]) if edge2[0] == common_node else (common_node, edge2[0])

        # new if same edge twice don't change order
        if get_idx_from_idc(self.idc_dict, edge1_in_order) == get_idx_from_idc(self.idc_dict, edge2_in_order):
            edge2_in_order = edge1_in_order

        return edge1_in_order, edge2_in_order

    def have_common_junction_node(self, edge1, edge2):
        # Extract nodes from the edges
        nodes_edge1 = set(edge1)
        nodes_edge2 = set(edge2)
        # TODO can change to self.junction_nodes?
        all_junctions = [
            *self.junction_nodes,
            self.graph_creator.processing_zone,
            self.graph_creator.entry,
            self.graph_creator.exit,
        ]

        # Check if the edges have any common junction nodes
        common_junction_nodes = nodes_edge1.intersection(nodes_edge2).intersection(all_junctions)

        return len(common_junction_nodes) == 1

    def create_outer_circle(self, edge_idc, next_edge, other_next_edges, towards=(0, 0)):
        # move from entry to memory zone
        if get_idx_from_idc(self.idc_dict, edge_idc) == get_idx_from_idc(
            self.idc_dict, self.graph_creator.entry_edge
        ):  # in self.graph_creator.path_from_pz_idxs:
            target_edge = self.bfs_free_edge(towards, other_next_edges)
            # calc path to target edge
            path0 = get_path_to_node(
                self.graph,
                self.graph_creator.processing_zone,
                target_edge[0],
                exclude_exit=True,
                exclude_first_entry_connection=False,
            )
            path1 = get_path_to_node(
                self.graph,
                self.graph_creator.processing_zone,
                target_edge[1],
                exclude_exit=True,
                exclude_first_entry_connection=False,
            )
            if len(path1) > len(path0):
                edge_path = [*path0, (target_edge[0], target_edge[1])]
            else:
                edge_path = [*path1, (target_edge[1], target_edge[0])]
            return edge_path

        # circles within memory zone
        node_path = nx.shortest_path(
            self.graph,
            next_edge[1],
            edge_idc[0],
            lambda node0, node1, _: [
                1e8
                if (
                    get_idx_from_idc(self.idc_dict, (node0, node1)) == get_idx_from_idc(self.idc_dict, edge_idc)
                    or get_idx_from_idc(self.idc_dict, (node0, node1)) == get_idx_from_idc(self.idc_dict, next_edge)
                    or get_idx_from_idc(self.idc_dict, (node0, node1))
                    == get_idx_from_idc(self.idc_dict, self.graph_creator.entry_edge)
                )
                else 1
            ][0],
        )
        edge_path = []
        for edge in pairwise(node_path):
            edge_path.append(edge)

        return [edge_idc, next_edge, *edge_path, edge_idc]

    def check_if_edge_is_filled(self, edge_idc):
        return get_idx_from_idc(self.idc_dict, edge_idc) in self.get_state_idxs()

    def find_nonfree_and_free_circle_idxs(self, circles_dict):
        # careful! If only two edges -> second edge must always be free! (could also be parking edge, because PE can store multiple)
        junction_nodes = [*self.junction_nodes, self.graph_creator.processing_zone]
        combinations_of_circles = list(distinct_combinations(circles_dict.keys(), 2))

        def get_circle_nodes(circle):
            # if next edge is free -> circle is just two edges -> can skip first and last node
            if len(circles_dict[circle]) == 2:
                if circles_dict[circle][0] != circles_dict[circle][1]:
                    circle_or_path = [(circles_dict[circle][0][1], circles_dict[circle][1][0])]
                    assert (
                        circles_dict[circle][0][1] == circles_dict[circle][1][0]
                    ), "circle is not two edges? Middle node should be the same"
                    # if middle node is exit or exit connection -> skip also middle node -> can always push through to parking edge
                    # TODO unskip
                    # if (
                    #     nx.get_node_attributes(self.graph, "node_type")[circles_dict[circle][0][1]] in ("exit_node", "exit_connection_node")
                    # ):
                    #     circle_or_path = []

                # new if same edge twice is parking edge -> skip completely
                elif get_idx_from_idc(self.idc_dict, circles_dict[circle][0]) == get_idx_from_idc(
                    self.idc_dict, self.graph_creator.parking_edge
                ):
                    # if self.count_chains_in_parking() >= self.max_num_parking:
                    #     circle_or_path = [(circles_dict[circle][0][0], circles_dict[circle][0][0])]
                    # else:
                    circle_or_path = []
                else:  # else if path is same edge twice skip (but of course keep first node -> no movement into this edge)
                    circle_or_path = [(circles_dict[circle][0][0], circles_dict[circle][0][0])]
            # if circle is real circle -> need to check all nodes
            elif circles_dict[circle][0] == circles_dict[circle][-1]:
                circle_or_path = circles_dict[circle]

            # if circle is only a path -> can skip first and last node
            # extra clause for when there is a stop at the end? -> if last two edges are same -> skip last two edges?
            elif circles_dict[circle][-1] == circles_dict[circle][-2]:
                circle_or_path = circles_dict[circle][1:-2]
            else:
                circle_or_path = circles_dict[circle][1:-1]
            nodes = set()
            for edge in circle_or_path:
                node1, node2 = edge
                if node1 == node2:
                    nodes.add(node1)
                else:
                    nodes.add(node1)
                    nodes.add(node2)
            return nodes

        junction_shared_pairs = []
        for circle1, circle2 in combinations_of_circles:
            nodes1 = get_circle_nodes(circle1)
            nodes2 = get_circle_nodes(circle2)

            # following clause is only to allow some cases that the assert would stop
            if (
                len(nodes1) != 0
                and len(nodes2) != 0
                and get_idx_from_idc(self.idc_dict, circles_dict[circle1][-1])
                != get_idx_from_idc(self.idc_dict, self.graph_creator.parking_edge)
                and not (
                    len(nodes1.intersection(nodes2).intersection(junction_nodes)) > 0
                    and self.graph_creator.processing_zone not in nodes1.intersection(nodes2)
                )
                # added new clause that it is allowed if they block each other (this is the if statement below, that adds the circles to junction_shared_pairs and thus the first blocks the second)
            ):
                # assert that circles don't end in same edge (just sanity check at this point - exceptions are added in if clause above after being checked)
                assert get_idx_from_idc(self.idc_dict, circles_dict[circle1][-1]) != (
                    get_idx_from_idc(self.idc_dict, circles_dict[circle2][-1])
                ), "circles end in same edge, was in if statement below as: or (get_idx_from_idc(self.idc_dict, circles_dict[circle1][-1]) == (get_idx_from_idc(self.idc_dict, circles_dict[circle2][-1]))), -> problem with circles: {}, {}".format(
                    circles_dict[circle1], circles_dict[circle2]
                )
            # new: exclude processing zone node -> if pz node in circles -> can both be executed (TODO check again for moves out of pz)
            if len(
                nodes1.intersection(nodes2).intersection(junction_nodes)
            ) > 0 and self.graph_creator.processing_zone not in nodes1.intersection(nodes2):
                junction_shared_pairs.append((circle1, circle2))

        # free_circle_combs = [
        #     circle_idx_pair
        #     for circle_idx_pair in combinations_of_circles
        #     if (circle_idx_pair not in junction_shared_pairs and (circle_idx_pair[1], circle_idx_pair[0]) not in junction_shared_pairs)
        # ]
        return junction_shared_pairs  # , free_circle_combs

    # change: if list in other list -> take longer list, delete other
    # if list can be connected to other list -> combine and delete both

    def rotate(self, full_circle_idxs, plot=False):
        # create dictionary of state
        # convert keys to values and vice versa, so one can iterate over the path (positions need to be keys here)
        edge_state_dict = {}
        for ion, edge in self.ion_chains.items():
            edge_state_dict[get_idx_from_idc(self.idc_dict, edge)] = ion

        new_edge_state_dict = {}
        for edge_bef, edge_aft in pairwise(full_circle_idxs):
            try:
                new_edge_state_dict[edge_aft] = edge_state_dict[edge_bef]
                del edge_state_dict[edge_bef]
            except KeyError:
                continue

        # change ion chains
        for idx, ion in new_edge_state_dict.items():
            self.ion_chains[ion] = get_idc_from_idx(self.idc_dict, idx)

        if plot is True:
            self.graph_creator.plot_state(
                [get_idx_from_idc(self.idc_dict, edge_idc) for edge_idc in self.ion_chains.values()], show_plot=True
            )

        return new_edge_state_dict

    def bfs_free_edge(self, node, other_next_edges):
        state_idxs = self.get_state_idxs()
        other_next_edges_idxs = [get_idx_from_idc(self.idc_dict, next_edge) for next_edge in other_next_edges]
        for edge_idc in nx.edge_bfs(self.mz_graph, node):
            if (
                get_idx_from_idc(self.idc_dict, edge_idc) not in state_idxs
                and get_idx_from_idc(self.idc_dict, edge_idc) not in other_next_edges_idxs
            ):
                return edge_idc
        return None

    def find_least_import_chain_in_parking(self, seq, ions_in_parking):
        for num in seq:
            if num in ions_in_parking:
                ions_in_parking.remove(num)
                if len(ions_in_parking) == 1:
                    return ions_in_parking[0]
        return ions_in_parking[-1]
