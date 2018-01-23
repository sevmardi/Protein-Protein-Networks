from misc import *
from prediction import *
from mix_bag import *
import argparse
import networkx as nx


def main():
    # Load the data first
    data = load_data()
    # Create two graph using networkx
    hgraph = nx.Graph(data.humanppi)
    fgraph = nx.Graph(data.functions)
    # If chosen to test, run the below.
    if args.test2:
        print('run predictions on test2 set')
        proteins_test2 = data.test2
        # Run the prediction on test2
        predictions = predict_cps(proteins_test2, hgraph, fgraph)
        print(zip(proteins_test2, predictions))
        # Write results to file.
        write_test2_prediction(proteins_test2, predictions)
    else:
        # take sample if no arguments were provided and test this on test1
        # and predict the cancer-relatedness of each of the proteins using
        # methods/algorithms such as Precision, Recall and F-measures
        sampler = Sampler(data.test1, args.samplesize)
        results = []

        for i in range(args.kfoldcross):
            sample = zip(*sampler.remainder)[0]
            predictions = predict_cps(sample, hgraph, fgraph)
            result = Results(sampler.remainder, predictions)
            results.append(result)
            sampler.resample()
            if args.verbose:
                result.print_results()
            result.print_confusion_matrix()
        f_measures = [i.f_measure for i in results]
        avg_f_measure = sum(f_measures) / float(len(f_measures))
        print("\n\nAverage F-Measure:", avg_f_measure)

        pause = raw_input("Presss any key to continue..")

        return 0


def get_options():
    parser = argparse.ArgumentParser()

    parser.add_argument('-v', '--verbose',
                        action='store_true',
                        help='Increases verbosity. ')
    parser.add_argument('-k', '--kfoldcross', type=int, default=10,
                        help='Cross validation in K times')

    parser.add_argument('-s', '--samplesize', type=float, default=0.9,
                        help='Percentage of test1 to train on. Example:\
                   -s 0.9 will use 90 percent of test1 to train on.')

    parser.add_argument('-t', '--test2', default=False,
                        action='store_true',
                        help='Run predictions on trainingset Test2.txt.')

    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = get_options()
    assert 0 < args.samplesize < 1
    print(args)
    main()
