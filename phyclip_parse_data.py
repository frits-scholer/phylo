#!/usr/bin/env python

if __name__ == '__main__':
    # argument parser
    import argparse
    parser = argparse.ArgumentParser(description='Acquire phylogenetic data for PhyCLIP.')
    parser.add_argument('-t', '--tree', type=str, required=True, help='Rooted input phylogenetic tree (NEWICK format).')
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

    import numpy as np
    ### Determine list of nodes which pairwise distance distribution falls under limits as partially inferred from tree ###
    # mean pairwise leaf distance distribution for each node
    node_to_mean_pairwise_dist = {n: np.mean(node_to_pairwise_leaf_distance_distribution[n]) for n in list_of_internal_node}
    # based on pariwise leaf distant distritbution of all nodes, calculate upper-limit of within-clade pairwise distance = med_x + user-defined-multiple*mad_x
    med_x = np.median(node_to_mean_pairwise_dist.values())
    mad_x = np.median([x - med_x for x in node_to_mean_pairwise_dist.values() if x >= med_x])  # because global distribution is likely not normal
    # update list of internal nodes to be considered for clade delineation to those whose mean pairwise distance < upper-limit
    list_of_internal_node = [n for n in list_of_internal_node if node_to_mean_pairwise_dist[n] <= (med_x + (params.mad*mad_x))]

    from scipy import stats
    ### Calculate p-values to test for inter-clade divergence ###
    # Apply two-sample Kolmogorov-Smirnov (KS) test only on remaining nodes
    nodepair_to_KS_pvalue = {}
    for i, j in itertools.combinations(list_of_internal_node, 2):
        # if either i or j is an ancestral node of the other
        if (i in node_to_ancestral_nodes[j]) or (j in node_to_ancestral_nodes[i]):
            nodepair_to_KS_pvalue[(i,j)] = stats.ks_2samp(node_to_pairwise_leaf_distance_distribution[i], node_to_pairwise_leaf_distance_distribution[j]).pvalue
        else:
            # find common ancestor node (ca) linking i and j
            pairwise_dist_distribution_common_ancestor = [leafpair_to_patristic_distance[(x, y)] for x, y in itertools.combinations(list(set(nindex_to_list_of_leaf_nodes[i])|set(nindex_to_list_of_leaf_nodes[j])), 2)]
            # perform KS-tests for (i vs ca) and (j vs ca), then take the conservative (higher) p-value
            nodepair_to_KS_pvalue[(i,j)] = max([stats.ks_2samp(node_to_pairwise_leaf_distance_distribution[i], pairwise_dist_distribution_common_ancestor).pvalue, stats.ks_2samp(node_to_pairwise_leaf_distance_distribution[j], pairwise_dist_distribution_common_ancestor).pvalue])

    import statsmodels.api as sm
    # multiple-testing correction using BH procedure
    nodepair_list = nodepair_to_KS_pvalue.keys()
    qval_list = sm.stats.multipletests([nodepair_to_KS_pvalue[nodepair] for nodepair in nodepair_list], method='fdr_bh')[1].tolist()
    nodepair_to_KS_pvalue = {nodepair: qval_list[i] for i, nodepair in enumerate(nodepair_list)}
    for i, j in nodepair_to_KS_pvalue.keys():
        nodepair_to_KS_pvalue[(j,i)] = nodepair_to_KS_pvalue[(i,j)]

    # write output
    print ('\nPrinting output...')
    if params.outfname:
        outfname = params.outfname
    else:
        outfname = 'phyclip_output.txt'
    with open(outfname, 'w') as output:
        # print index of eligible nodes and the leaves subtended to be considered for clades
        for n in list_of_internal_node:
            output.write('{}\t{}\n'.format(n, ','.join(nindex_to_node[n].get_leaf_names())))
        # print KS p-values for node pairs
        for (i,j) in nodepair_to_KS_pvalue.keys():
            output.write('{}\t{}\t{}\n'.format(i, j, nodepair_to_KS_pvalue[(i,j)]))

    print ('\n...done.\n')
    exit(0)
