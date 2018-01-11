# Predict Cancer-Related Proteins in Protein-Protein Interaction Networks
# Intro
Computational protein function prediction is a very important and actively studied research area. Current
state of the art methods use techniques such as homology modeling (for related proteins), sequence analysis
and weighted and annotated Protein-Protein Interaction (PPI) networks. 

# Dataset
PPI-network is modeled as a graph where each node is a protein and each edge is a physical interaction between two proteins. There are two types of annotation information for each protein p in the PPI network. A ’functional annotation’ shows the biological functions of p in the PPI network and a ’cancer-relatedness annotation’ shows if p is involved in cancer or not. The task of predicting cancer related proteins in PPI networks is trying to predict which of the new proteins are involved in cancer.

Note: Dataset is is not allowed to be published. 

# Task
You are asked to design, implement, and evaluate a strategy for predicting cancer-related proteins in PPI networks. A PPI network as described above is given in a dataset. You will be asked to evaluate your method with one test-set and give your predictions for a second test-set.


# Dataset 
- HumanPPI.txt: In this file each row < Prot i , Prot j > represents one edge (modeling a physical
interaction) between two proteins Prot i and Prot j in the PPI network. Please note, the edges are
undirected, this means every edge in the PPI network will be represented in the file twice, by
< Prot i , Prot j > and < Prot j , Prot i >.

- Functions.txt: Here each row < Prot i , Func j > adds the function Func j to the function set of protein
Prot i . Please be careful that, there exist some proteins in the dataset which are not functionally
annotated yet.

- Cancers.txt: Each row < Prot i > in this file indicates that protein Prot i is involved in cancer.


# Evaluation

Evalution should of the algorithms based on the following files:

Test1.txt: This file contains proteins that are known to be involved in cancer (as given in the file
Cancers.txt) or not. The file should be used to evaluate the performance of your proposed prediction
algorithm. You should predict the cancer-relatedness of each of the proteins in the file Test1.txt using your methods/algorithms, and report the evaluation results according to the Precision, Recall and F-measures:

Precision = tp/(tp+fp)
Recall = tp/(tp+fn)
F-measure = 2 ∙ Precision ∙ Recall / Precision + Recall

where proteins involved in cancer are considered as the positive class, and tp, fp and fn denote the
number of true positives, false positives, and false negatives, respectively.

Test2.txt: Determine for each protein p in Test2.txt file, if p is involved in cancer or not and report
your predictions in a file PredicionResultsTest2.txt. For example if there are three proteins p1 , p2
and p3 in Test2.txt and your algorithm predicts just p1 as a cancer-related protein, then
PredictionResultsTest2.txt should be EXACTLY as follows:

p1 , Cancer
p2 , nonCancer
p3 , nonCancer


# References
- Xuebing Wu and Shao Li, Cancer Gene Prediction Using a Network Approach. Chapter 11
Mathematical and Computational Biology. Cancer Systems Biology (Ed. Edwin Wang). Series:
Chapman and Hall/CRC, 2010. 

- Renzhi Cao, and Jianlin Cheng, Integrated protein function prediction by mining function associations,
sequences, and protein–protein and gene–gene interaction networks, Methods, Vol. 93, January 2016,
pp. 84 – 91. (Available here: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4894840/ )

- Predrag Radivojac et al., A large-scale evaluation of computational protein function prediction, Nat
Methods, Vol. 10(3), March 2013, pp 221 - 227.
