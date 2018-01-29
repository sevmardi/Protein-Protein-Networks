import numpy as np
import networkx as nx
import matplotlib as mpl
import matplotlib.pyplot as plt
from misc import zeros
from misc import load_data

# Helper class: To do some plotting


def main():
    data = load_data()
    hgraph = nx.Graph(data.humanppi)
    fgraph = nx.Graph(data.functions)

    # plot_degree(hgraph)
    # plot_cp_degree(hgraph)
    #plot_fn_degree(hgraph, fgraph)
    plot_fn_cp_weight(hgraph, fgraph)


def plot_degree(hgraph):

    cps = [p for p in hgraph.nodes() if p.is_cancer_protein()]
    ps = [p for p in hgraph.nodes() if not p.is_cancer_protein()]

    cp_degrees_cp = [p.degree(hgraph) for p in cps]
    avg_cp = round(np.mean(cp_degrees_cp), 2)
    std_cp = round(np.std(cp_degrees_cp), 2)

    cp_degrees = [p.degree(hgraph) for p in ps]
    avg = round(np.mean(cp_degrees), 2)
    std = round(np.std(cp_degrees), 2)

    lims = 50
    b = 20
    bins = np.linspace(0, lims, b)
    kwargs = {'bins': bins,
              'normed': True,
              'histtype': 'bar',
              'color': ['g', 'r'],
              'label': ['Normal Proteins', 'Cancer Proteins'],
              'alpha': 0.7
              }

    plt.hist([cp_degrees, cp_degrees_cp], **kwargs)
    plt.legend(loc='upper right')
    #plt.title('Distribution of the number of edges')
    plt.xlim(0, lims)
    plt.xticks(np.linspace(0, lims, lims / 5 + 1))
    plt.ylabel('Normalized frequency')
    plt.xlabel('Number of edges')
    plt.text(15, 0.08, 'Normal proteins: Avg: %s   Std: %s' %
             (str(avg), str(std)))
    plt.text(15, 0.07, 'Cancer proteins: Avg: %s   Std: %s' %
             (str(avg_cp), str(std_cp)))

    filepathname = 'report/degrees_combined.png'
    plt.savefig(filepathname,  bbox_inches='tight')
    print('Done!' + ' Check this folder => ' + filepathname)
    plt.show()


def plot_cp_degree(hgraph):
    cps = [p for p in hgraph.nodes() if p.is_cancer_protein()]
    ps = [p for p in hgraph.nodes() if not p.is_cancer_protein()]

    cp_degrees_cp = [p.cp_degree(hgraph) for p in cps]
    avg_cp = round(np.mean(cp_degrees_cp), 2)
    std_cp = round(np.std(cp_degrees_cp), 2)

    cp_degrees = [p.cp_degree(hgraph) for p in ps]
    avg = round(np.mean(cp_degrees), 2)
    std = round(np.std(cp_degrees), 2)

    lims = 15
    b = 10
    bins = np.linspace(0, lims, b)
    kwargs = {'bins': bins,
              'normed': True,
              'histtype': 'bar',
              'color': ['g', 'r'],
              'label': ['Normal Proteins', 'Cancer Proteins'],
              'alpha': 0.7
              }

    plt.hist([cp_degrees, cp_degrees_cp], **kwargs)
    plt.legend(loc='upper right')
    #plt.title('Distribution of the number of cancer protein neighbors')
    plt.xlim(0, lims)
    plt.xticks(np.linspace(0, lims, lims / 5 + 1))
    plt.ylabel('Normalized frequency')
    plt.xlabel('Number of cancer protein neighbors')
    plt.text(5, 0.25, 'Normal proteins: Avg: %s   Std: %s' %
             (str(avg), str(std)))
    plt.text(5, 0.22, 'Cancer proteins: Avg: %s   Std: %s' %
             (str(avg_cp), str(std_cp)))

    filepathname = 'report/cp_degrees_combined.png'
    plt.savefig(filepathname,  bbox_inches='tight')
    print('Done!' + ' Check this folder => ' + filepathname)
    plt.show()


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

    lims = 60
    b = 30
    bins = np.linspace(0, lims, b)
    kwargs = {'bins': bins,
              'normed': True,
              'histtype': 'bar',
              'color': ['g', 'r'],
              'label': ['Normal Proteins', 'Cancer Proteins'],
              'alpha': 0.7
              }

    plt.hist([ps_fn_degrees, cps_fn_degrees], **kwargs)
    plt.legend(loc='upper right')
    plt.xlim(0, lims)
    plt.xticks(np.linspace(0, lims, lims / 5 + 1))
    plt.ylabel('Normalized frequency')
    plt.xlabel('Number of functions')
    plt.text(20, 0.05, 'Normal proteins: Avg: %s   Std: %s' %
             (str(avg), str(std)))
    plt.text(20, 0.045, 'Cancer proteins: Avg: %s   Std: %s' %
             (str(avg_cp), str(std_cp)))

    filepathname = 'report/fn_degrees_combined.png'
    plt.savefig(filepathname,  bbox_inches='tight')
    print('Done!' + ' Check this folder => ' + filepathname)
    plt.show()


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

    lims = 80000
    b = 15
    bins = np.linspace(0, lims, b)
    kwargs = {'bins': bins,
              'normed': True,
              'histtype': 'bar',
              'color': ['g', 'r'],
              'label': ['Normal Proteins', 'Cancer Proteins'],
              'alpha': 0.7
              }

    plt.hist([ps_fn_cp_weight, cps_fn_cp_weight], **kwargs)
    plt.legend(loc='upper right')
    plt.xlim(0, lims)
    plt.xticks(np.linspace(0, lims, 11))
    plt.ylabel('Normalized frequency')
    plt.xlabel('Assigned weight')
    plt.text(15000, 0.00007, 'Normal proteins: Avg: %s   Std: %s' %
             (str(avg), str(std)))
    plt.text(15000, 0.00006, 'Cancer proteins: Avg: %s   Std: %s' %
             (str(avg_cp), str(std_cp)))

    filepathname = 'report/fn_w_combined.png'
    plt.savefig(filepathname,  bbox_inches='tight')
    print('Done!' + ' Check this folder => ' + filepathname)
    plt.show()


if __name__ == '__main__':
    main()
