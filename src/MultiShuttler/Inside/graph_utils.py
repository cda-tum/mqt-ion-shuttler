import networkx as nx


# create dictionary to swap from idx to idc and vice versa
def create_idc_dictionary(graph):
    edge_dict = {}
    for edge_idx, edge_idc in enumerate(graph.edges()):
        edge_dict[edge_idx] = tuple(sorted(edge_idc, key=sum))
    return edge_dict


def get_idx_from_idc(edge_dictionary, idc):
    idc = tuple(sorted(idc, key=sum))
    return list(edge_dictionary.values()).index(idc)


def get_idc_from_idx(edge_dictionary, idx):
    return edge_dictionary[idx]


class GraphCreator:
    def __init__(self, m, n, ion_chain_size_vertical, ion_chain_size_horizontal):
        self.m = m
        self.n = n
        self.ion_chain_size_vertical = ion_chain_size_vertical
        self.ion_chain_size_horizontal = ion_chain_size_horizontal

        self.networkx_graph = self.create_graph()

    def create_graph(self):
        self.m_extended = self.m + (self.ion_chain_size_vertical - 1) * (self.m - 1)
        self.n_extended = self.n + (self.ion_chain_size_horizontal - 1) * (self.n - 1)

        networkx_graph = nx.grid_2d_graph(
            self.m_extended,
            self.n_extended,
        )
        self._set_trap_nodes(networkx_graph)
        self._remove_horizontal_edges(networkx_graph)
        self._remove_vertical_edges(networkx_graph)
        self._remove_horizontal_nodes(networkx_graph)
        networkx_graph.junction_nodes = []
        self._set_junction_nodes(networkx_graph)
        nx.set_edge_attributes(networkx_graph, "trap", "edge_type")
        nx.set_edge_attributes(networkx_graph, 1, "weight")

        return networkx_graph

    def _set_trap_nodes(self, networkx_graph):
        for node in networkx_graph.nodes():
            networkx_graph.add_node(
                node, node_type="trap_node", color="k", node_size=150
            )

    def _remove_horizontal_edges(self, networkx_graph):
        for i in range(
            0,
            self.m_extended - self.ion_chain_size_vertical,
            self.ion_chain_size_vertical,
        ):
            for k in range(1, self.ion_chain_size_vertical):
                for j in range(self.n_extended - 1):
                    networkx_graph.remove_edge((i + k, j), (i + k, j + 1))

    def _remove_vertical_edges(self, networkx_graph):
        for i in range(
            0,
            self.n_extended - self.ion_chain_size_horizontal,
            self.ion_chain_size_horizontal,
        ):
            for k in range(1, self.ion_chain_size_horizontal):
                for j in range(self.m_extended - 1):
                    networkx_graph.remove_edge((j, i + k), (j + 1, i + k))

    def _remove_horizontal_nodes(self, networkx_graph):
        for i in range(
            0,
            self.m_extended - self.ion_chain_size_vertical,
            self.ion_chain_size_vertical,
        ):
            for k in range(1, self.ion_chain_size_vertical):
                for j in range(
                    0,
                    self.n_extended - self.ion_chain_size_horizontal,
                    self.ion_chain_size_horizontal,
                ):
                    for s in range(1, self.ion_chain_size_horizontal):
                        networkx_graph.remove_node((i + k, j + s))

    def _remove_mid_part(self, networkx_graph):
        for i in range(self.ion_chain_size_vertical):
            networkx_graph.remove_node((self.m_extended // 2, self.n_extended // 2 + i))
        for i in range(1, self.ion_chain_size_vertical):
            networkx_graph.remove_node((self.m_extended // 2, self.n_extended // 2 - i))
        for i in range(1, self.ion_chain_size_horizontal):
            networkx_graph.remove_node((self.m_extended // 2 + i, self.n_extended // 2))
        for i in range(1, self.ion_chain_size_horizontal):
            networkx_graph.remove_node((self.m_extended // 2 - i, self.n_extended // 2))

    def _set_junction_nodes(self, networkx_graph):
        for i in range(0, self.m_extended, self.ion_chain_size_vertical):
            for j in range(0, self.n_extended, self.ion_chain_size_horizontal):
                networkx_graph.add_node(
                    (i, j), node_type="junction_node", color="g", node_size=400
                )
                networkx_graph.junction_nodes.append((i, j))

    def get_graph(self):
        return self.networkx_graph
