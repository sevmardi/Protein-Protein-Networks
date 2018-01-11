import numpy as np
import networkx as nx
import matplotlib as mpl
import matplotlib.pyplot as plt
from misc import zeros
from misc import load_data


def main():
    data = load_data()
    hgraph = nx.Graph(data.humanppi)
    fgraph = nx.Graph(data.functions)

    plot_cp_degree(hgraph)
    # plot_fn_degree(hgraph, fgraph)
    # plot_fn_cp_weight(hgraph, fgraph)

def plot_cp_degree(hgraph):
    cps = [p for p in hgraph.nodes() if p.is_cancer_protein()]
    ps = [p for p in hgraph.nodes() if not p.is_cancer_protein()]

    cp_degrees_cp = [p.cp_degree(hgraph) for p in cps]
    avg_cp = round(np.mean(cp_degrees_cp), 2)
    std_cp = round(np.std(cp_degrees_cp), 2)

    cp_degrees = [p.cp_degree(hgraph) for p in ps]
    avg = round(np.mean(cp_degrees), 2)
    std = round(np.std(cp_degrees), 2)

    fig = plt.figure(figsize=(10, 10))
    ax1 = fig.add_subplot(211)
    ax2 = fig.add_subplot(212)

    kwargs = {'bins': 50,
              'normed': True,
              'histtype': 'step',
              'color': 'C0'}
    ax1.hist(cp_degrees_cp, **kwargs)
    ax1.set_title('Cancer Proteins Only')
    ax1.set_xlim(0, 15)
    ax1.set_ylabel('Normalized frequency')
    ax1.text(5, 0.3, 'Avg: %s   Std: %s' % (str(avg_cp), str(std_cp)))

    ax2.hist(cp_degrees, **kwargs)
    ax2.set_title('Non-Cancer Proteins Only')
    ax2.set_xlim(0, 15)
    ax2.set_ylabel('Normalized frequency')
    ax2.set_xlabel('Number of Cancer Protein Neighbors')
    ax2.text(5, 0.3, 'Avg: %s   Std: %s' % (str(avg), str(std)))

    filepathname = 'report/plots/cp_degrees.png'
    plt.savefig(filepathname,  bbox_inches='tight')
    print('Done!' + ' Check this folder => ' + filepathname)
    fig.clf()
    plt.close()

    del fig


def plot_fn_degree(hgraph, fgraph):
    proteins_in_fgraph = [p for p in fgraph.nodes()
                          if not p.is_function()]
    cps = [p for p in proteins_in_fgraph if p.is_cancer_protein()]
    ps = [p for p in proteins_in_fgraph if not p.is_cancer_protein()]

    cps_fn_degrees = [p.degree(fgraph) for p in cps]
    avg_cp = round(np.mean(cps_fn_degrees), 2)
    std_cp = round(np.std(cps_fn_degrees), 2)

    ps_fn_degrees = [p.degree(fgraph) for p in ps]
    avg = round(np.mean(ps_fn_degrees), 2)
    std = round(np.std(ps_fn_degrees), 2)

    fig = plt.figure(figsize=(10, 10))
    ax1 = fig.add_subplot(211)
    ax2 = fig.add_subplot(212)

    kwargs = {'bins': 50,
              'normed': True,
              'histtype': 'step',
              'color': 'C0'}

    ax1.hist(cps_fn_degrees, **kwargs)
    ax1.set_title('Cancer Proteins Only')
    ax1.set_xlim(0, 80)
    ax1.set_ylabel('Normalized frequency')
    ax1.text(50, 0.03, 'Avg: %s   Std: %s' % (str(avg_cp), str(std_cp)))

    ax2.hist(ps_fn_degrees, **kwargs)
    ax2.set_title('Non-Cancer Proteins Only')
    ax2.set_xlim(0, 80)
    ax2.set_ylabel('Normalized frequency')
    ax2.set_xlabel('Number of Functions')
    ax2.text(50, 0.03, 'Avg: %s   Std: %s' % (str(avg), str(std)))

    filepathname = 'report/plots/fn_degrees.png'
    plt.savefig(filepathname,  bbox_inches='tight')
    print('Done!' + ' Check this folder => ' + filepathname)
    fig.clf()
    plt.close()

    del fig


def plot_fn_cp_weight(hgraph, fgraph):
    proteins_in_fgraph = [p for p in fgraph.nodes()
                          if not p.is_function()]

    cps = [p for p in proteins_in_fgraph if p.is_cancer_protein()]
    ps = [p for p in proteins_in_fgraph if not p.is_cancer_protein()]

    cps_fn_cp_weight = [p.fn_cp_weight(hgraph, fgraph) for p in cps]
    avg_cp = round(np.mean(cps_fn_cp_weight), 2)
    std_cp = round(np.std(cps_fn_cp_weight), 2)

    ps_fn_cp_weight = [p.fn_cp_weight(hgraph, fgraph) for p in ps]
    avg = round(np.mean(ps_fn_cp_weight), 2)
    std = round(np.std(ps_fn_cp_weight), 2)

    fig = plt.figure(figsize=(10, 10))
    ax1 = fig.add_subplot(211)
    ax2 = fig.add_subplot(212)

    kwargs = {'bins': 50,
              'normed': True,
              'histtype': 'step',
              'color': 'k'}
    ax1.hist(cps_fn_cp_weight, **kwargs)

    ax1.set_title('Cancer Proteins Only')
    ax1.set_ylabel('Normalized frequency')
    ax2.hist(ps_fn_cp_weight, **kwargs)
    ax2.set_title('Non-Cancer Proteins Only')
    ax2.set_ylabel('Normalized frequency')
    ax2.set_xlabel('Number of Functions')


    filepathname = 'report/plots/fn_cp_weight.png'
    plt.savefig(filepathname,  bbox_inches='tight')
    print('Done!' + ' Check this folder => ' + filepathname)
    fig.clf()
    plt.close()

    del fig





def make_function_graphs(data):
    pass


if __name__ == '__main__':
    main()
