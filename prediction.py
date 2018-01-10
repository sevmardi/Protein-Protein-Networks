import networkx as nx
from operator import *
from mix_bag import *

def evaluate(predicted_genes, removed_genes, network_annotated_gen):
    pass


def distribute_cross_validation(network_annotated_gen):
    pass


def cross_validation(network, network_annotated_gene, gene_annotation, go_num):
    pass


def remove_unannotated_genes(network, gene_annotation):
    pass


def write_test2_prediction(proteins, predictions):
    f = open('PredictionResultsTest2.txt', 'w')
    for p, pred in zip(proteins, predictions):
        string = ','.join(["Prot_" + str(p), pred]) + '\n'
        f.write(string)
    f.close()


def predict_cps(proteins, hgraph, fgraph):
    """
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



	