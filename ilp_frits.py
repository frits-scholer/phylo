from __future__ import division
from os import listdir
from gurobipy import *
from scipy import stats
from argparse import ArgumentParser, FileType
from os.path import expanduser
import statsmodels.api as sm
import sys, ete3, subprocess
import numpy as np

# parse arguments
def parse_args():
    parser = ArgumentParser(description='')
    # required inputs
    parser.add_argument('-t', '--tree', type=str, required=True, help='Rooted NEWICK tree file.')
    parser.add_argument('-cs', '--cs', type=int, default=3, help='Minimum cluster size (Default: 3)')
    parser.add_argument('-fdr', '--fdr', type=float, default=0.2, help='FDR (Default: 0.2)')
    parser.add_argument('-mad', '--mad', type=float, default=3.0, help='MAD (Default: 3.0)')
    return parser.parse_args()

# global
if __name__ == '__main__':
    # get annotation info
    tag_to_annotation = {}
    tag_to_unswitched = {}
    tag_to_patient = {}
    for file in listdir(expanduser('~/Dropbox/antibody_literature/malaria_annotations')):
        fhandle = filter(None, open(expanduser('~/Dropbox/antibody_literature/malaria_annotations/{}'.format(file)), 'rU'))
        fhandle.pop(0)
        patient = re.search('(1017|1019|2207)', file).group()
        for line in fhandle:
            unknown, tag, label, vgene, jgene, mutV, mutJ, cdr3len, sequence = line.strip().split(',')
            tag = re.search('MS6_\d+', tag).group()
            tag_to_annotation[tag] = '_'.join([tag, label.replace('"',''), 'PT{}'.format(patient), 'mutV{}'.format(mutV), 'mutJ{}'.format(mutJ)])
            tag_to_patient[tag] = patient
            if re.search('(IGHM|IGHD)', label):
                tag_to_unswitched[tag] = 0

    # parse arguments
    args = parse_args()

    # read newick tree
    tree = ete3.Tree(args.tree, format=0)
    # tree info
    taxon_list = tree.get_leaf_names()
    # ladderize
    tree.ladderize()

    # data
    # enumerate nodes
    nindex_to_node, node_to_nindex = {}, {}
    ancestral_node_to_list_of_leaf_nodes = {}  # list of leaves to each ancestral node
    for n, node in enumerate(tree.traverse()):
        if node.is_leaf() == False:
            ancestral_node_to_list_of_leaf_nodes[n] = node.get_leaves()
            nindex_to_node[n] = node
            node_to_nindex[node] = n

    list_of_ancestral_node = nindex_to_node.keys()  # list of ancestral nodes (index)
    node_to_ancestral_nodes = {n:[node_to_nindex[anc_node] for anc_node in nindex_to_node[n].iter_ancestors()] for n in list_of_ancestral_node} # ancestral lineage of each node

    # pairwise patristic distances between all leaf nodes in tree
    leafpair_to_patristic_distance = {}
    for x, y in itertools.combinations(tree.get_leaves(), 2):
        leafpair_to_patristic_distance[(x,y)] = leafpair_to_patristic_distance[(y,x)] = x.get_distance(y)
    # record pairwise leaf distance distribution for each node
    node_to_pairwise_leaf_distance_distribution = {n:[leafpair_to_patristic_distance[(x,y)] for x,y in itertools.combinations(ancestral_node_to_list_of_leaf_nodes[n], 2)] for n in list_of_ancestral_node}
    # mean pairwise leaf distance distribution for each node
    node_to_mean_pairwise_dist = {n:np.mean(node_to_pairwise_leaf_distance_distribution[n]) for n in list_of_ancestral_node}
    # setting max limit of mean within-clade pairwise distance
    med_x = np.median(node_to_mean_pairwise_dist.values())
    mad_x = np.median([abs(x-med_x) for x in node_to_mean_pairwise_dist.values() if x >= med_x]) # because global distribution is likely not normal

    print ('...cs{}, fdr{}, mad{}...'.format(args.cs, args.fdr, args.mad))

    # update current list of ancestral nodes
    curr_list_of_ancestral_node = ''
    curr_list_of_ancestral_node = [n for n in list_of_ancestral_node if len(ancestral_node_to_list_of_leaf_nodes[n]) >= args.cs and node_to_mean_pairwise_dist[n] <= (med_x + (args.mad * mad_x))]

    # KS-test for remaining nodes
    nodepair_to_KS_pvalue = ''
    nodepair_to_KS_pvalue = {}
    for i,j in itertools.combinations(curr_list_of_ancestral_node, 2):
        if (i in node_to_ancestral_nodes[j]) or (j in node_to_ancestral_nodes[i]):
            nodepair_to_KS_pvalue[(i,j)] = stats.ks_2samp(node_to_pairwise_leaf_distance_distribution[i], node_to_pairwise_leaf_distance_distribution[j]).pvalue
        else:
            pairwise_dist_distribution_common_ancestor = [leafpair_to_patristic_distance[(x,y)] for x,y in itertools.combinations(list(set(ancestral_node_to_list_of_leaf_nodes[i])|set(ancestral_node_to_list_of_leaf_nodes[j])), 2)]
            # take the conservative (max) p-value comparing node i/j individually to i+j
            nodepair_to_KS_pvalue[(i,j)] = max([stats.ks_2samp(node_to_pairwise_leaf_distance_distribution[i], pairwise_dist_distribution_common_ancestor).pvalue, stats.ks_2samp(node_to_pairwise_leaf_distance_distribution[j], pairwise_dist_distribution_common_ancestor).pvalue])

    # multiple-testing correction using BH procedure
    nodepair_list = nodepair_to_KS_pvalue.keys()
    qval_list = sm.stats.multipletests([nodepair_to_KS_pvalue[nodepair] for nodepair in nodepair_list], method='fdr_bh')[1].tolist()
    nodepair_to_KS_pvalue = {nodepair:qval_list[i] for i, nodepair in enumerate(nodepair_list)}
    for i,j in nodepair_to_KS_pvalue.keys():
        nodepair_to_KS_pvalue[(j,i)] = nodepair_to_KS_pvalue[(i,j)]

    # GUROBI SOLVER
    # model
    model = Model()

    # variable
    node_decision = ''
    node_decision = model.addVars(curr_list_of_ancestral_node, vtype=GRB.BINARY)

    # objective function - maximize number of strains clustered
    model.ModelSense = GRB.MAXIMIZE
    model.setObjective(quicksum(node_decision[n]*len(ancestral_node_to_list_of_leaf_nodes[n]) for n in curr_list_of_ancestral_node))
    model.update()

    # constraints
    # must select at least one node
    model.addConstr(quicksum(node_decision[n] for n in curr_list_of_ancestral_node) >= 1)

    # ancestral node of selected descendant cannot be chosen
    model.addConstrs(node_decision[n] + quicksum(node_decision[m] for m in node_to_ancestral_nodes[n] if m in curr_list_of_ancestral_node) <= 1 for n in curr_list_of_ancestral_node)

    # inter-clade pairwise leaf distance distribution must be distinct from that should they be combined as a single clade
    model.addConstrs((2 - node_decision[n] - node_decision[m])*100 >= nodepair_to_KS_pvalue[(n,m)]-args.fdr for n,m in tuplelist([(n,m) for n,m in itertools.combinations(curr_list_of_ancestral_node, 2)]))

    # update and optimize
    model.update()
    model.optimize()

    if model.status == GRB.Status.OPTIMAL:
        # get solution
        solution = model.getAttr('x', node_decision)
        selected_nodes = [n for n in curr_list_of_ancestral_node if solution[n] == 1]
    else:
        sys.exit('\nERROR: NO SOLUTION FOUND.\n')

    print ('\nWrite output...')
    with open('ilp_output.txt', 'w') as output:
        for n in selected_nodes:
            node = nindex_to_node[n]
            leaves = node.get_leaf_names()
            output.write('>Cluster_Node_{}\n{}\n\n'.format(n, '\n'.join(leaves)))

    print ('\n...done.\n')
    exit(0)
