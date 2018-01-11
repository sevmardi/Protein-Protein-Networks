import h5py
import argparse
import numpy as np
import networkx as nx

from main import load_data


def main():
    data = load_data()
    humanppi_G = nx.Graph(data.humanppi, name="HumanPPI")
    write_to_hdf5(humanppi_G)
    return 0


def write_to_hdf5(graph):
    pagerankings = np.array(nx.pagerank(graph, max_iter=1000).items())
    closeness = np.array(nx.closeness_centrality(graph).items())
    degrees = np.array(graph.degree().items())
    
    hdf5file = h5py.File(args.outputfilepath, 'w')
    hdf5file['protein_nr'] = graph.nodes()
    hdf5file['degree'] = degrees[:, 1]
    hdf5file['pagerank'] = pagerankings[:, 1]
    hdf5file['closeness'] = closeness[:, 1]

    hdf5file.close()
    del hdf5file


def get_options():
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', '--outputfilepath', type=str, default='nodeinfo.hdf5',
                        help='output filepath for hdf5 file.')
    args = parser.parse_args()

    return args

if __name__ == '__main__':
    args = get_options()
    print(args)
    main()
