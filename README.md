# satute

---- What is Satute?

Satute is a Python module to test for branch saturation in a phylogeny. In most of cases, a branch is saturated if it is too long, although sometimes the cause is more technical than just branch length.

Importantly, Satute is highly dependent on IQ-TREE (http://www.iqtree.org), which we assume the user has installed.

---- Workflow of Satute (see the image attached)

[workflow.pdf](https://github.com/cassiusma/satute/files/11662009/workflow.pdf)

The minimum input required to run Satute is an alignment file in phylip or fasta format. IQ-TREE reconstructs an evolutionary tree and then Satute will tell us, for every branch of this tree, whether this branch is saturated or not. If there are different evolutionary speeds (e.g. due to the Gamma model), then Satute will separate the alignment sites into the categories that likely generated them.

The final output is a results file listing, for every branch, its saturation status. If there are different categories, then one results file is output for each category.

The statistic we use to test for branch saturation is the coherence of the branch, which has a value typically between 0 (saturated) and 1 (unsaturated).

---- Ingredients

You will need:

--An alignment file alignment.phy for satute, with path e.g. /path/alignment.phy

--IQ-TREE installed, whose executable has path e.g. /path/iqtree2

--The main of this project, with path e.g. /path/main.py.  Installing the modules of the requirements.txt into your python3.

Having run IQ-TREE before running Satute is optional. BUT IMPORTANTLY: I assume you added no prefix nor suffix to the IQ-TREE output files.


---- Examples using the command line

python3 /path/main.py -s /path/alignment.phy -iq /path/iqtree2 







