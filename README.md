# satute
Satute is a Python module to test for branch saturation in a phylogeny. In most of cases, a branch is saturated if it is too long, although this is not always true.
Importantly, Satute is highly dependent on IQ-TREE (http://www.iqtree.org), which we assumed the user has installed. Let us explain the workflow of Satute.

[workflow.pdf](https://github.com/cassiusma/satute/files/11662009/workflow.pdf)

The minimum input required to run Satute is an alignment file in phylip or fasta format. IQ-TREE reconstructs an evolutionary tree and then Satute will tell us, for every branch of this tree, whether this branch is saturated or not. If there are different evolutionary speeds (e.g. due to the Gamma model), then we need to separate the alginment sites into the category that likely generated them.

The final output is a results file listing, for every branch, its saturation status. If there are different categories, then one results file is output for each category.


The statistic we use to test for branch saturation is the coherence of the branch, which has a value typically between 0 (saturated) and 1 (unsaturated).


