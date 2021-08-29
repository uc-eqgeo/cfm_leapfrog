from typing import List, Set

import networkx
from networkx.algorithms.components.connected import connected_components


def to_graph(node_list: List[List[str]]):
    """
    Answer from https://stackoverflow.com/a/4843408
    :param node_list:
    :return:
    """
    graph = networkx.Graph()
    for part in node_list:
        # each sublist is a bunch of nodes
        graph.add_nodes_from(part)
        # it also implies a number of edges:
        graph.add_edges_from(zip(part[:-1], part[1:]))
    return graph


def connected_nodes(node_list: List[List[str]]):
    graph = to_graph(node_list)
    return list(connected_components(graph))

def suggest_combined_name(connected_nodes: set):
    pass


