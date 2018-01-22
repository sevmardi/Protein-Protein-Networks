from colorama import Fore, Style
import random


class Bunch(object):

    def __init__(self, **kwds):
        self.__dict__.update(kwds)


class Sampler(object):
    """Sampler class - Takes a sample from a dataset.

    Parameters
    ----------
    dataset: a list of datapoints
    percentage: relative size of sample w.r.t. size of the dataset.

    """

    def __init__(self, set_to_sample_from, percentage):
        self.dataset = set(set_to_sample_from)
        self.samplesize = int(len(set_to_sample_from) * percentage)

        self.take_sample_from_set()

        self.resample = self.take_sample_from_set

    def take_sample_from_set(self):
        """
        Return a sample from the dataset.

        Parameters
        ----------
        testset: list of tuples.
                 i.e.: [(0, 'nonCancer'), (1, 'nonCancer'),
                        (1208, 'cancer'), (2431, 'cancer'),
                        (2, 'nonCancer')]

        percentage: sample size in terms of a percentage of the testset.

        Return
        ------
        random subset of testset.

        """
        self.sample = set(random.sample(self.dataset, self.samplesize))
        self.remainder = set(self.dataset) - self.sample
        assert len(self.sample.intersection(self.remainder)) == 0
        assert len(self.sample) + len(self.remainder) == len(self.dataset)


class Node(object):

    def __init__(self, name):
        self.name = name

    def __str__(self):
        return str(self.name)

    def __hash__(self):
        return hash(self.name)

    def __eq__(self, other):
        return self.__dict__ == other.__dict__

    def is_function(self):
        """
        Returns
        -------
        False (this method is overridden in the subclass Function)
        """
        return False

    def is_cancer_protein(self):
        """
        Returns
        -------
        False (this method is overridden in the subclass CancerProtein)
        """
        return False

    def degree(self, graph):
        """
        Parameters
        ----------
        graph: A <class 'networkx.classes.graph.Graph'> object

        Returns
        -------
        Degree of the node in the graph

        """
        return graph.degree(self)

    def neighbors(self, graph):
        """
        Parameters
        ----------
        graph: A <class 'networkx.classes.graph.Graph'> object

        Returns
        -------
        Neighbors of the node.

        """
        return graph.neighbors(self)

    def neighbors_iter(self, graph):
        """
        Parameters
        ----------
        graph: A <class 'networkx.classes.graph.Graph'> object

        Returns
        -------
        Neighbors of the node.

        """
        return graph.neighbors_iter(self)

    def cp_neighbors(self, graph):
        """
        Returns
        -------
        cps: a list of cancerous protein neighbors

        """
        nbrs = self.neighbors(graph)
        cps = [p for p in nbrs if p.is_cancer_protein()]
        return cps

    def cp_degree(self, hgraph):
        """
        Returns
        -------
        number of cancerous protein neighbors

        """
        return len(self.cp_neighbors(hgraph))


class Function(Node):
    """
    Subclass of Node.
    """

    def __init__(self, name):
        self.name = name

    def __repr__(self):
        string = "%s"
        subst = (str(self.name),)
        return Fore.GREEN + string % subst + Style.BRIGHT + Fore.WHITE

    def color(self):
        return 'g'

    def is_function(self):
        """
        Returns
        -------
        True. (overrrides the Node.is_function() method.)
        """
        return True

    def has_cps(self, fgraph):
        """
        Returns
        -------
        True if the function contains cancerous proteins.

        """
        if self.cp_neighbors(fgraph) == 0:
            return False
        else:
            return True

    def cp_p_ratio(self, fgraph):
        """
        Returns
        -------
        The ratio of (cancerous proteins)/(proteins) of the function.
        """
        number_of_proteins = fgraph.degree(self)
        return float(self.cp_degree(fgraph)) / number_of_proteins

    def cancerweight(self, hgraph, fgraph):
        """
        Return
        ------
        sum of cp_degrees of all cp's in this function
        """
        weight = 0
        for cp in self.cp_neighbors(fgraph):
            weight = weight + cp.cp_degree(hgraph)
        return weight


class Protein(Node):

    def __repr__(self):
        string = "P%s"
        subst = (str(self.name),)
        return string % subst

    def color(self):
        return 'w'

    def functions(self, fgraph):

        if self in fgraph.nodes():
            return self.neighbors(fgraph)
        else:
            return []

    def fn_degree(self, fgraph):
        return len(self.functions(fgraph))

    def cpp_ratio(self, hgraph):
        """
        Return The ratio of (neighbor cancerous proteins)/(neighbor proteins)
        of the protein.
        """
        return float(self.cp_degree(hgraph)) / number_of_nbrs

    def shared_cp_functions(self, hgraph, fgraph):
        """
        Returns a list of functions that shares with neighbors
        cancerous proteins.
        """
        shared_fns = set()
        for cp_nbr in self.cp_neighbors(hgraph):
            if cp_nbr in fgraph.nodes():
                cp_nbr_fns = set(cp_nbr.neighbors(fgraph))
                shared_fns = shared_fns.union(cp_nbr_fns)
            else:
                pass

        return list(shared_fns)

    def fn_cp_degree(self, hgraph, fgraph):
        """Returns the number of functions that protein shares with neighbors cancerous proteins."""
        return len(self.shared_cp_functions(hgraph, fgraph))

    def fn_cp_weight(self, hgraph, fgraph):
        """
        Return a weight based on the cancerweight of the functions that
        the protein shares with cancerous neighbors

        """
        weight = 0
        for fn in self.shared_cp_functions(hgraph, fgraph):
            weight = weight + fn.cancerweight(hgraph, fgraph)
        return round(weight, 2)

    def cp_degree_of_neighbors(self, hgraph):
        """Returns the sum of cp_degrees of each cancerous protein neighbors """
        nbrs = self.neigbors(hgraph)
        total = 0
        for i in nbrs:
            total += n.cp_degree(hgraph)
        return total

    def cancerweight(self, hgraph, fgraph):
        """
        Returns a weight to determine the likelihood a protein is a cancerous
        protein """

        return self.fn_cp_weight(hgraph, fgraph)


class CancerProtein(Protein):
    """
    Subclass of Protein"""

    def __repr__(self):
        string = "P%s"
        subst = (str(self.name),)
        return Fore.RED + string % subst + Fore.WHITE

    def color(self):
        return 'r'

    def is_cancer_protein(self):
        """
        True (overrides protein.is_cancer_protein)
        """

        return True


class Results(object):
    """ Results class - Container for a prediction set. Analyses
    predictions and calculates Precision, Recall and F-Measure.

    Parameters
    ----------
    protein_answer = list of tuples = [ (protein, answer),
                                       (123, 'cancer'),
                                       ... , ]
    predictions: list of predictions = ['cancer','nonCancer', ... ]
    """

    def __init__(self, protein_answer, predictions):
        self.true_positives = 0
        self.true_negatives = 0
        self.false_positives = 0
        self.false_negatives = 0

        self.proteins, self.answers = zip(*list(protein_answer))
        self.predictions = predictions

        for i, prediction in enumerate(predictions):
            if self.answers[i] == prediction and \
                    prediction == 'nonCancer':
                outcome = "True  Negative"
                self.true_negatives += 1
            elif self.answers[i] == prediction and \
                    prediction == 'cancer':
                outcome = "True Positive"
                self.true_positives += 1
            elif self.answers[i] != prediction and \
                    prediction == 'nonCancer':
                outcome = "False Negative"
                self.false_negatives += 1
            else:
                outcome = "False Positive"
                self.false_positives += 1

    def print_results(self):
        print("Protein    Answer     Prediction     Outcome \n" +
              "--------------------------------------------")
        line = "%s %s %s  %s"
        for i, prediction in enumerate(self.predictions):
            if self.answers[i] == prediction and \
                    prediction == 'nonCancer':
                outcome = "True  Negative"
            elif self.answers[i] == prediction and \
                    prediction == 'cancer':
                outcome = "True  Positive"
            elif self.answers[i] != prediction and \
                    prediction == 'nonCancer':
                outcome = "False Negative"
            else:
                outcome = "False Positive"

            substitutions = (str(self.proteins[i]).ljust(10), str(
                self.answers[i]).ljust(10), prediction.ljust(13), outcome)

            print(line % substitutions)

    def print_confusion_matrix(self):
        string = "\n" + \
                 "Confusion Matrix \n\n" + \
                 "  Tp  |  Fn        %s  |  %s  \n" + \
                 "-------------    ------------------\n" + \
                 "  Fp  |  Tn        %s  |  %s  \n\n" + \
                 "Precision: %f  Recall: %f  F-Measure: %f"

        substitutions = (str(self.true_positives).ljust(3),
                         str(self.false_negatives),
                         str(self.false_positives).ljust(3),
                         str(self.true_negatives),
                         self.precision, self.recall, self.f_measure)
        print(string % substitutions)

    def write_to_file(self):
        """Write the results into a file"""
        pass

    @property
    def precision(self):
        """
        Parameters
        ----------
        true_positives: number of true positives (int)
        false_positives: number of false positives (int)


        Returns
        -------
        Precision, which is defined as:

                                   # True_pos
               Precision =  --------------------------
                             (#True_pos + #False_pos)

        """
        Tp = self.true_positives
        Fp = self.false_positives
        if Tp == 0:
            return 0
        else:
            return float(Tp) / (Tp + Fp)

    @property
    def recall(self):
        """
        Parameters
        ----------
        true_positives: number of true positives (int)
        false_negatives: number of false negatives (int)


        Returns
        -------
        Recall, which is defined as:

                                 # True_pos
               Recall  = --------------------------
                          (#True_pos + #False_neg)

        """
        Tp = self.true_positives
        Fn = self.false_negatives
        if Tp == 0:
            return 0
        else:
            return float(Tp) / (Tp + Fn)

    @property
    def f_measure(self):
        """
        Parameters
        ----------
        p: precision (float)
        r: recall (float)


        Returns
        -------
        F-Measure, which is defined as:

                         2 * precision * recall
          F-Measure  =  ------------------------
                         ( precision + recall )

        """
        p = self.precision
        r = self.recall
        if (p or r) == 0:
            return 0
        else:
            return float(2 * p * r) / (p + r)


class Table(object):
    """Class to facilitate table formatting"""

    def __init__(self, name="Table_name"):
        self.name = name
        self.operator = ', '
        self.colwidth = 15
        self.arrays = {}

    def add(self, array_name, elem):
        """ Add element to array """
        if array_name in self.arrays:
            self.arrays[array_name].append(elem)
        else:
            self.arrays[array_name] = [elem]

    def get_column_order(self):
        """ Order the columns and return it. """
        try:
            columns = self.order
        except AttributeError:
            columns = self.arrays.keys()
        return columns

    def columns_names(self):
        """ Helper function: Join the ordered columns and join it toegether """
        columns = self.get_column_order()
        sep = self.separator
        string = sep.join(c.ljust(self.colwidth) for c in columns) + '\n'
        return string

    def lines_iter(self):
        """ 
        Iterator for lines to be written. Assumes all arrays are
        of equal length.
        Return Iterator
        """
        width = self.colwidth
        sep = self.separator
        arrays = self.arrays
        columns = self.get_column_order()
        arraylength = len(self.arrays.values()[0])

        i = 0
        if i < arraylength:
            line = sep.join([str(arrays[c][i]).ljust(width)
                for c in columns] + "\n"
        i += 1
        yield line
