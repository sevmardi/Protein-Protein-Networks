import os
from mix_bag import Sampler
from mix_bag import Results
from mix_bag import Bunch
from mix_bag import Protein
from mix_bag import CancerProtein
from mix_bag import Function


def zeros(number, pad_length):
    """Pad with zeros"""
    return str(number).zfill(pad_length)


def setup_dirs(*dirs):
    """ Creates directories if they don't already exist."""
    for dir in dirs:
        if os.path.exists(dir):
            pass
        else:
            os.makedires(dir)


def read_file(file_name):
    file_handler = open(file_name, 'r')
    raw = file_handler.readlines()
    file_handler.close()
    return raw


def load_data():
    """ Loads data from the data folder"""

    cancer_text = read_file('data/Cancer.txt')
    cancer = [line.strip() for line in cancer_text]
    cancer = [int(line[5:]) for line in cancer]

    humanppi_txt = read_file('data/humanPPI.txt')
    humanppi = [line.strip().split(",") for line in humanppi_txt]
    humanppi = [tuple([int(elem[5:]) for elem in line]) for line in humanppi]

    temp = []
    for p1, p2 in humanppi:
        if p1 in cancer:
            a = CancerProtein(p1)
        else:
            a = Protein(p1)
        if p2 in cancer:
            b = CancerProtein(p2)
        else:
            b = Protein(p2)
        temp.append((a, b))
    del humanppi
    humanppi = temp

    functions_txt = read_file('data/Functions.txt')
    functions = [line.strip().split(',') for line in functions_txt]
    functions = [tuple([int(line[0][5:]), "F" + line[1][5:]]) for line
                 in functions]

    f = []
    for p, fn in functions:
        if p in cancer:
            a = CancerProtein(p)
        else:
            a = Protein(p)
        b = Function(fn)
        f.append((a, b))
    del functions
    functions = f

    test1_txt = read_file('data/Test1.txt')
    test1 = [line.strip().split(',') for line in test1_txt]
    test1 = [tuple([int(line[0][5:]), line[1]]) for line in test1]
    test1 = [(Protein(p), answer) for p, answer in test1]

    test2_txt = read_file('data/Test2.txt')
    test2 = [line.strip().split(',') for line in test2_txt]
    test2 = [int(line[0][5:]) for line in test2]
    test2 = [Protein(p) for p in test2]

    data = Bunch(humanppi=humanppi,
                 functions=functions,
                 cancer=cancer,
                 test1=test1,
                 test2=test2)

    return data



