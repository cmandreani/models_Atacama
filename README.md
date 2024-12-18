# Modeling the emergent metabolic potential of soil microbiomes in Atacama landscapes

This script runs the analyses and plots included in the article by CM Andreani-Gerard et al (2025).
Briefly, our work simulates the metabolic potential of prokaryotic communities from six constrasting sites of the Talabre-Lejía transect in the Atacama Desert.
For this, we used sequence data from whole communities and MAGs that enabled us to build community-wide and genome-resolved metabolic models.
In our systems biology framework, five nutritional conditions were defined to dinamically simulate organic sources of carbon, nitrogen, and sulfur. 
We refer to the results of the simulations (combining six sites, two datasets, and five seeds) as "scopes" which consist of lists of metabolites predicted to be producible.
Finally, we associate metabolic and environmental data through a regression-based approach.
Note that Fig. 5 is generated directly by Metage2Metabo (Belcour et al. 2020. eLife. 10.7554/eLife.61968).

If you find this methodology useful, please cite our work: XXXX.

Zenodo: DOIXXXXX

## 1) BACKGROUND OF THE TLT
####    1.1 Environmental metadata
####    1.2 Alpha diversity
####    1.3 Taxonomic and functional profiling
## 2) METABOLISM OF THE TLT
####    2.1 Extraction of scopes for metagenomic data (MetaG-GEMs)
####    2.2 Extraction of scopes for genomic data (MAG-GEMs)
####    2.3 Conversion of scopes into MXs
####    2.4 Draw: boxplot
####    2.5 Draw: flowchart
####    2.6 Draw: PCoA
####    2.7 Draw: heatmap
## 3) Selection of key metabolites with an elastic net regression 
####    3.1 Creation of elastic net
####    3.2 Draw: heatmap
## 
/!\ Don't forget to adapt the codes hereby provided to your own data structure ;)
