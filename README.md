# Combined Essentiality Score (CES)

Probing the genetic dependencies of cancer cells helps understand the tumor biology and identify potential drug targets. RNAi-based shRNA and CRISPR/Cas9-based sgRNA have been commonly utilized in functional genetic screens to identify essential genes affecting growth rates in cancer cell lines. However, questions remain whether the gene essentiality profiles determined using these two technologies are comparable. 
We collected 42 cell lines representing a variety of 10 tissue types, which had been screened both by shRNA and CRISPR techniques. We observed poor consistency of the essentiality scores between the two screens for the majority of the cell lines. The consistency did not improve after correcting the off-target effects in the shRNA screening, suggesting a minimal impact of off-target effects. 
We considered a linear regression model where the shRNA essentiality score is the predictor and the CRISPR essentiality score is the response variable. We showed that by including molecular features such as mutation, gene expression and copy number variation as covariates, the predictability of the regression model greatly improved, suggesting that molecular features may provide critical information in explaining the discrepancy between the shRNA and CRISPR-based essentiality scores.
We provided a Combined Essentiality Score (CES) based on the model prediction and showed that the CES greatly improved the consensus of common essential genes. Furthermore, the CES also identified novel essential genes that are specific to individual cell types. 


# Code 
R code for condcuting the analysis is provided.

# Data
Preapred gene-level data, gene name annotation as well as a housekeeping gene list used for the code are avaliable at https://drive.google.com/open?id=1QyrAxUzNEvbvYwHg4_xXh-dCgWB7a0A2



