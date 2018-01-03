#!/usr/bin/env python

if __name__ == '__main__':
    # argument parser
    import argparse
    parser = argparse.ArgumentParser(description='Acquire phylogenetic data for PhyCLIP.')
    parser.add_argument('-t', '--tree', type=str,  default='toy_tree_small.nwk', help='Rooted input phylogenetic tree (NEWICK format).')
    parser.add_argument('-c', '--min_clade_size', type=int, default=3, help='Minimum number of leaves to be considered a clade (default: %(default)s)')
    parser.add_argument('-m', '--mad', type=int, default=3, help='Admissible multiple of positive-skew median absolute deviation of median patristic distance of a clade (default: %(default)s)')
    parser.add_argument('-o', '--outfname', type=str, help='Output filename (optional).')
    params = parser.parse_args()

    import ete3, itertools
    # parse tree using ete3
    tree = ete3.Tree(params.tree)
    tree.ladderize() # ladderize tree

    # acuqire data
    nindex_to_node, node_to_nindex = {}, {}
    nindex_to_list_of_leaf_nodes = {}
    node_to_ancestral_nodes = {}
    for n, node in enumerate(tree.traverse()):
        # proceed only if node is an internal node (i.e. not a leaf) and has at least specified number of leaves to be considered a clade
        if node.is_leaf() == False and len(node.get_leaves()) >= params.min_clade_size:
            # get list of leaf nodes for internal node
            nindex_to_list_of_leaf_nodes[n] = node.get_leaves()
            # enumerate internal nodes - better for integer linear (ILP) programming solver than using ete3 node
            nindex_to_node[n] = node
            node_to_nindex[node] = n

    # list of internal nodes (indices)
    list_of_internal_node = nindex_to_node.keys()
    # dictionary of ancestral lineage to each internal node - required for determining inter-cluster divergence
    node_to_ancestral_nodes = {n:[node_to_nindex[anc_node] for anc_node in nindex_to_node[n].iter_ancestors()] for n in list_of_internal_node}

    # pairwise distances between all leaves in tree
    leafpair_to_patristic_distance = {}
    for x, y in itertools.combinations(tree.get_leaves(), 2):
        leafpair_to_patristic_distance[(x,y)] = leafpair_to_patristic_distance[(y,x)] = x.get_distance(y)

    # dictionary of pairwise distance distribution for each internal node
    node_to_pairwise_leaf_distance_distribution = {n:[leafpair_to_patristic_distance[(x,y)] for x,y in itertools.combinations(nindex_to_list_of_leaf_nodes[n], 2)] for n in list_of_internal_node}
    print nindex_to_node[2].get_leaf_names()
    print  node_to_pairwise_leaf_distance_distribution[2]
    print nindex_to_node[12].get_leaf_names()
    print  node_to_pairwise_leaf_distance_distribution[12]
    from scipy import stats
    print stats.ks_2samp(node_to_pairwise_leaf_distance_distribution[12], node_to_pairwise_leaf_distance_distribution[2])
 
