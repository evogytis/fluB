## Reassortment between influenza B lineages and the emergence of a co-adapted PB1-PB2-HA gene complex
### Scripts

The following python scripts were used in the manuscript to extract information from trees:
- `Hrothvitnir_eLife.py`
- `Hydra_eLife.py`
- `LD_calculator_ChiSqdf_eLife.py`

Most figures were made using an IPython notebook:
- `fluB_eLife.ipynb`

## Data
To replicate findings you will need two things:
* Sequences used in the manuscript. The [sequences we used](https://github.com/evogytis/fluB/blob/master/acknowledgement%20tables/fluB_gisaid_acknowledge_table.tsv) can be downloaded from [GISAID](http://platform.gisaid.org). Sequences should be codon-aligned and contain just the coding sequences.
* [BEAUti and BEAST](https://code.google.com/p/beast-mcmc/). The XML files we used are [available here](https://github.com/evogytis/fluB/tree/master/data/BEAST%20XML%20files) and only require the BEAUti-parsed sequences.

## Tree processing
Both `Hrothvitnir_eLife.py` (Hróðvitnir means "fame-wolf" in Old Norse) and `Hydra_eLife.py` are distant descedants of a [short linked list script](http://stackoverflow.com/questions/280243/python-linked-list/280286#280286) from StackOverflow. Both were made to analyze and summarize the properties of trees drawn from a posterior distribution of trees and have been written to be as general as possible. It should take minimal effort to change the code to suit almost any BEAST-related need.

### Hrothvitnir
`Hrothvitnir_eLife.py` is used to extract information about diversity of phylogenies [sampled from a posterior distribution of trees](https://code.google.com/p/beast-mcmc/). You will require **python** and **numpy** to run this script.

To run type:
``python Hrothvitnir_eLife.py -i [path to input] -m [mode] 1> [path to output]``
ignoring the square brackets.

Input should be posterior distribution of trees made by BEAST.

Mode takes one of 4 arguments:
* __diversityOT__ - returns the date of most recent common ancestor of all lineages existing at different time points.
* __FstOT__ - returns mean pairwise date of most recent common ancestor between branches under different trait values at time slices. Only works if trait values are "V" and "Y".
* __stateTime__ - returns the nucleotide, synonymous and non-synonymous rates of evolution and the total amount of time spent under "VVV", "YYY" or other combinations of trait values of PB1, PB2 and HA traits.
* __X__ - returns the number of objects parsed by the script. Completely unnecessary, since the script has inbuilt checks to see whether the tree has been parsed correctly.

### Hydra
`Hydra_eLife.py` is more complicated and takes in 4 trees at a time (hence the name) and focuses on quantifying tree-to-tree distances. It requires **python**, **numpy** and, for SPR distances, [**RSPR**](http://kiwi.cs.dal.ca/Software/RSPR).

To run you will have to edit the script itself under the bit where it says "OPTIONS". There you can choose whether to give it a directory and files to work on automatically or whether you would like to pick the files for analysis manually. There are additional options and explanations within the script itself.
After editing simply type in:
``python Hydra_eLife.py``

Input is 4 trees - a set of trees sampled from the posterior distribution of alignment A, another set for alignment B and independent analyses of alignment A and B, which we refer to as trees A' and B'.

Briefly, it normalizes comparisons between two different trees by comparisons between replicate trees.


## Linkage disequilibrium
`LD_calculator_ChiSqdf_eLife.py` calculates the [Chi-squared df](http://www.genetics.org/content/112/1/135) statistic of linkage disequilibrium. It requires **python**, **biopython** and **numpy** to run, as well as nucleotide or amino acid alignments. To run edit the script to point to the directory and the file name format and type in:
`python LD_calculator_ChiSqdf_eLife.py`

## Figures
`fluB_eLife.ipynb` makes most of the data figures from the manuscript. In order to make the figures you will require **IPython**, **matplotlib**, **numpy** and **scipy** (the latter is optional). After everything has been installed simply put `fluB_eLife.ipynb` where the IPython notebook can find it. You will need to point it to the output files which are [available here](https://github.com/evogytis/fluB/tree/master/data/) if you can't be bothered running all other scripts.

The notebook in its intended form can be viewed [here](http://nbviewer.ipython.org/github/evogytis/fluB/blob/master/scripts/fluB_eLife.ipynb?create=1) (give it some time, it's quite big).
