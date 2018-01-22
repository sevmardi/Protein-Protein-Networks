import time
import networkx as nx
import matplotlib.pyplot as plt
import argparse
import math
from misc import zeros
from misc import load_data


def main():
    data = load_data()

    if args.edges == 0:
        ppiGraph = nx.Graph(data.humanppi, name="HumanPPI")
    else:
        ppiGraph = nx.Graph(data.humanppi[0:args.edges], name="HumanPPI")
    draw_kwargs = dict(node_size=5, font_size=6,
                       with_labels=False, linewidths=0.1, width=0.2)
    layout_kw = dict(iterations=args.iterations)
    img_title = 'Nodes: %i  Edges: %i  Iterations: %i  Runtime: %s'

    if args.iterations == None:
        # Number of iters will be set to 10, check the get_iterations method.
        iters = get_iterations()

        for iteration in iters:
            start_time = time.time()

            fig = plt.figure(figsize=(10, 10))
            ax = fig.add_subplot(111)
            nx.draw(ppiGraph, pos=nx.spring_layout(
                ppiGraph, iterations=iteration), **draw_kwargs)

            run_time = time.time() - start_time
            runtime_str = str(round(run_time, 2)) + ' seconds'

            print("Iterations: %i  Runt-time: %s" % (iteration, runtime_str))

            subs = (len(ppiGraph), args.edges, iteration, runtime_str)
            ax.set_title(img_title % subs)
            filepath = 'report/graphs' + '/' + \
                zeros(iteration, padlength=4) + '.png'
            plt.savefig(filepath,  bbox_inches='tight')
            print("Done! Saved at : " + filepath)

            fig.clf()
            plt.close()
    else:
        start_time = time.time()
        fig = plt.figure(figsize=(10, 10))
        ax = fig.add_subplot(111)
        nx.draw(ppiGraph, pos=nx.spring_layout(
            ppiGraph, **layout_kw), **draw_kwargs)

        runtime = time.time() - starttime
        runtime_str = str(round(runtime, 2)) + ' seconds'

        subs = (len(ppiGraph), args.edges, args.iterations, runtime_str)
        ax.set_title(img_title % subs)

        filepath = 'report/graphs' + '/' + \
            zeros(args.iterations, padlength=4) + '.png'
        plt.savefig(filepath,  bbox_inches='tight')
        print("Done! Saved at " + filepath)
        fig.clf()
        plt.close()


def get_iterations():
    """ Define the list of iterations here for now."""
    x1 = range(1, 10)
    x2 = [int(str(i) + '0') for i in x1]
    x3 = [int(str(i) + '0') for i in x2]
    x4 = [int(str(i) + '0') for i in x3]
    iters = x1 + x2 + x3 + x4
    return iters


def get_options():
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', '--outputdir', type=str, default='graphs',
                        help='output directory for images.')
    parser.add_argument('-i', '--iterations', type=int, default=None,
                        help='Do a single graph with i iterations')
    parser.add_argument('-e', '--edges', type=int, default=0,
                        help='cutoff graph at e edges. ')
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    args = get_options()
    main()
