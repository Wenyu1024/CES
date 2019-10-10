# Combined Essentiality Score (CES)

Interrogating the genetic dependencies of cancer cells provides important first evidence for targetbased drug discovery. RNAi-based shRNA and CRISPR-Cas9-based sgRNA have been commonly utilized in functional genetic screens to derive cancer dependence map. However, previous studies suggested limited overlap of the essentiality profile based on two technologies. Existing computational methods mainly focused on the estimation of true gene essentiality from genetic screens using single technologies, while there is a lack of integrative methods to combine the gene essentiality profiles from both CRISPR and shRNA screens.

We developed a computational approach called combined gene essentiality score (CES) to integrate CRISPR and shRNA gene essentiality profiles as well as molecular features of the cancer cells. CES improved significantly the performance of the gene essentiality prediction, not only for shared genetic dependencies across multiple cell lines, but also for therapeutic targets that are selective for a specific cancer cell line.

# Code 
R markdown for replicating the analysis in the manuscript is provided (analysis_report_revised2.Rmd).

# Data
Preapred gene-level data, gene name annotation as well as reference gene sets used for the code are avaliable at https://drive.google.com/file/d/1wDxQLiLja6r-OdiSQX6HcfVg9CIGrmwm/view?usp=sharing



