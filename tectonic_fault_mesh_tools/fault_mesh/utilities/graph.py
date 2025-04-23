"""
Submodule that uses networkx to find connections between segments in a fault map

Functions:
    - to_graph: Converts a list of node lists into a graph
    - connected_nodes: Finds all connected components in a graph created from node lists
    - suggest_combined_name: Suggests a combined name for a set of connected nodes based on their names
"""
from collections import Counter
from typing import List

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
    """Find all connected components in a graph created from node lists.

    :param node_list: Lists of nodes where each list represents a path or chain of connected nodes
    :type node_list: List[List[str]]
    :return: List of sets, where each set contains nodes that are connected to each other
    :rtype: list[set]
    """
    graph = to_graph(node_list)
    return list(connected_components(graph))


def suggest_combined_name(connected_node_set: list):
    """Suggest a combined name for a set of connected nodes.
    
    This function analyzes a set of connected node names and attempts to find
    a meaningful combined name based on common elements or patterns in the names.
    
    :param connected_node_set: A set of node names that are connected in a graph
    :type connected_node_set: list
    :return: A suggested name for the combined nodes
    :rtype: str
    """
    connected_node_list = list(connected_node_set)
    split_list = [set(name.split()) for name in connected_node_list]
    intersections = split_list[0].intersection(*split_list)

    if len(intersections):
        if len(intersections) > 1:
            ordered_example = None
            for item in connected_node_list:
                if all([x in item for x in intersections]):
                    ordered_example = item

            if ordered_example is not None:
                ordered_name = " ".join([x for x in ordered_example.split() if x in intersections])
                out_name = ordered_name + " combined"

            else:
                out_name = " ".join(list(intersections)) + " combined"

        else:
            out_name = " ".join(list(intersections)) + " combined"

    else:
        colons = [name.split(":")[0] for name in connected_node_list if ":" in name]
        if len(colons):
            counted = Counter(colons)
            out_name = counted.most_common()[0][0] + " combined"

        elif len(connected_node_list) == 2:

            if connected_node_list[0] in connected_node_list[1]:
                out_name = connected_node_list[1] + " combined"
            elif connected_node_list[1] in connected_node_list[0]:
                out_name = connected_node_list[0] + " combined"
            else:
                out_name = " - ".join(connected_node_list) + " combined"

        else:
            out_name = connected_node_list[0] + " combined"

    return out_name
