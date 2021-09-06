from typing import List, Set
from collections import Counter

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


def suggest_combined_name(connected_node_set: list):

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

    # colons = [name.split(":")[0] for name in connected_node_list if ":" in name]
    # print(colons)
    # if len(colons):
    #     counted = Counter(colons)
    #     out_name = counted.most_common()[0][0] + " combined"
    # elif len(connected_node_list) == 2:
    #
    #     if connected_node_list[0] in connected_node_list[1]:
    #         out_name = connected_node_list[1] + " combined"
    #     elif connected_node_list[1] in connected_node_list[0]:
    #         out_name = connected_node_list[0] + " combined"
    #     else:
    #         out_name = " ".join(connected_node_list) + " combined"
    # else:
    #     out_name = connected_node_list[0] + " combined"
    #
    return out_name




