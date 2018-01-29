import networkx as nx


def distribute_cross_validation(network_annotated_gen):
    """
    @deprecated DOT NOT USE :distribute all annotated gene to different cross validation group
    """
    group_genes = {}
    index = 0
    for gene in network_annotated_gen:
        group = index % config.CV
        if not group in group_genes:
            group_genes[group] = [gene]
        else:
            group_genes[group].append(gene)
        index += 1
    return group_genes


def cross_validation(network, network_annotated_gene, gene_annotation, go_num):
    """ @deprecated: DO NOT USE: Put annotated gene into different cross validation groups"""
    group_genes = distribute_cross_validation(network_annotated_gen)

    recall_avg = 0.0
    precision_avg = 0.0

    for i in range(0, config.CV):
        print('Cross Validation %d.... ' % i)
        annotated_gene_cv = {}
        removed_genes = group_genes[i]
        for gene in network_annotated_gene:
            if not gene in removed_genes:
                annotated_gene_cv[gene] = network_annotated_gene[gene]
        #...


def remove_unannotated_genes(network, gene_annotation):
    """
    @deprecated: DO NOT USE:Remove unannotated gene from network
    """
    for node in network.nodes():
        if not node in gene_annotation:
            network.remove_nodes(node)
    print("Number of nodes in network after removing unannonated genes %d" %
          network.number_of_nodes())


def write_test2_prediction(proteins, predictions):
    """
    Just a method to take progeins, and prdictions and save those into a file
    """

    f = open('results/PredictionResultsTest2.txt', 'w')
    for p, pred in zip(proteins, predictions):
        string = ','.join(["Prot_" + str(p), pred]) + '\n'
        f.write(string)
    f.close()


def predict_cps(proteins, hgraph, fgraph):
    """
    Using the weight functions to predict 2 kind of graphs and make predictions.
    Returns list of predictions, i.e. ['cancer','nonCancer', etc.]

    Parameters
    ----------
    proteins: list of Protein instances.
    hgraph: networkx.Graph() instance
    fgraph: networkx.Graph() instance
    """

    prediction = []
    for p in proteins:
        if p.cancerweight(hgraph, fgraph) > 15000:
            prediction.append('cancer')
        else:
            prediction.append('nonCancer')
    return prediction


def get_shortest_path(hgraph):
    """
    Returns a dictionary where keys represent all the proteins in the graph
    and values represent the number of occurences of non-cancer proteins
    along the shortest paths between any two cancer proteins.
    """
    cps = [p for p in hgraph.nodes() if p.is_cancer_protein()]
    aps = [p for p in hgraph.nodes()]

    paths = []

    for i in range(len(cps)):
        p1 = cps[i]
        for j in range(i, len(cps)):
            p2 = cps[j]
            paths.append(nx.shortest_path(hgraph, p1, p2))
    flagged = {}
    for p in aps:
        flagged[p] = 0
    for path in paths:
        for p in path:
            if not p.is_cancer_protein():
                flagged[p] += 1
    return flagged
