## MAAMOUL: A method for detecting microbiome-metabolome alterations in disease using metabolic networks

**Table of contents:**
 - [Method overview](#ch1)
 - [Installation](#ch2)
 - [Instructions - Running MAAMOUL on your own data](#ch3)
 - [Usage example](#ch4)
 - [FAQs](#ch5)

<a id="ch1"></a>
## Method overview

MAAMOUL is a knowledge-based computational method that integrates metagenomic and metabolomic data to identify custom data-driven microbial metabolic modules associated with disease states. Unlike traditional statistical approaches, MAAMOUL leverages prior biological knowledge about bacterial metabolism to link genes to metabolites through a global, microbiome-wide metabolic network, and then projects genes' and metabolites' disease-association scores onto this network. The identified 'modules' are sub-networks in this graph that are significantly enriched with disease-associated features, both metagenomic and metabolomic.

For further details see: Muller E, Baum S, and Borenstein E. __"Detecting Microbiome-Metabolome Alterations in Disease Using Metabolic Networks."__ _In preparation_.

<img src="man/Figure 1.png" width="700">

***

<a id="ch2"></a>
## Installation

MAAMOUL can be installed directly from GitHub, by running the following:

```
install.packages("devtools")  
library(devtools)   
install_github("efratmuller/MAAMOUL")   
library(MAAMOUL)
```

Note: The MAAMOUL package is dependant on the installation of the 'BioNet' package [1]. See installation instructions [here](https://www.bioconductor.org/packages/release/bioc/html/BioNet.html).

***
   
<a id="ch3"></a>
## Instructions - Running MAAMOUL on your own data

_Coming soon..._

***

<a id="ch4"></a>
## Usage example

```
library(MAAMOUL)
write_test_files()
maamoul(global_network_edges = 'test_input/enzyme_compound_edges_kegg.csv',
  ec_pvals = 'test_input/ec_pvals.tsv',
  metabolite_pvals = 'test_input/mtb_pvals.tsv',
  out_dir = 'test_outputs',
  N_REPEATS = 100,
  N_VAL_PERM = 9,
  N_THREADS = 4
)
```

***

<a id="ch5"></a>
## FAQs

#### What does it mean if MAAMOUL identified a module but it did not come out significant?
MAAMOUL first identifies significant nodes (i.e., ECs or metabolites that are significantly associated with a phenotype) that are grouped closely in a global metabolic network. These are termed 'modules'. A p-value is then computed per module to describe the likelihood of identifying such a module given the topology of the network and the specific dataset un hand (e.g., the number of significant ECs, the number of significant metabolites, etc.). Intuitively, a non-significant p-value means that the identified module could have resulted from data characteristics and random chance rather than actual biological signal.

#### MAAMOUL has identified a microbiome-metabolome module in my data, however I see there are additional significant nodes (ECs/metabolite) in proximity to the identified module. Why weren't these nodes added to the identified module?
Most likely, the significant nodes that weren't included in the identified module were initially assigned to another module, which eventually did not meet the module significance threshold. Future method versions will include an option for post-processing of modules, potentially merging them with nearby nodes as long as modules' significance isn't harmed.

*** 

For questions about the pipeline, please open an issue (https://github.com/efratmuller/MAAMOUL/issues) or contact Prof. Elhanan Borenstein at elbo@tauex.tau.ac.il.

***

__References__

1. Beisser D, Klau GW, Dandekar T, Mueller T, Dittrich M (2010). “BioNet: an R-package for the Functional Analysis of Biological Networks.” Bioinformatics, 26(8), 1129-1130.
