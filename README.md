
## scRNAzirf: extending random forests to handle zero-inflated count data

<!-- badges: start -->
### Daniel Conn and S&#252;nd&#252;z Kele&#351;
<!-- badges: end -->

## scRNAzirf

scRNAzirf is an implementation of the zero-inflated random forests ZIRFs
algorithm. ZIRFs is an extension of random forests designed to explicitly
account for the high rates of zero-inflation found in certain types of
single-cell RNA-seq (scRNA-seq) data. This package serves as

## Installation
To install and use scRNAzirf you'll need download SCENIC. Directions for downloading
SCENIC can be found at this link:
# https://htmlpreview.github.io/?https://github.com/aertslab/SCENIC/blob/master/inst/doc/SCENIC_Setup.html#installation.

## Example

We first need to load the required packages to run SCENIC with ZIRFs variable importance measures
(VIMs). For guidance on downloading the requisite packages to use SCENIC see SCENIC's
documentation on installation and set-up:
https://github.com/aertslab/SCENIC/blob/master/vignettes/SCENIC_Setup.Rmd.
```{r setup}
suppressPackageStartupMessages({
  library(scRNAzirf)
  library(SCENIC)
  library(foreach)
  library(BiocParallel)
  library(doMC)
})
knitr::opts_knit$set(root.dir="zirfsSCENIC_vignette")
```

We will go through a toy example using diaphragm cell data from the Tabula Muris Consortium
(https://tabula-muris.ds.czbiohub.org/). During this vignette we will be downloading databases,
creating intermediate files, and output files. We recommend creating new directories to keep
these files organized.

## Directories
```{r dirsetup, eval = FALSE}
#directory which we will be primarily working in
dir.create("zirfsSCENIC_vignette")
setwd("zirfsSCENIC_vignette")
#directory where databases for motif scoring in SCENIC will be held
dir.create("mousedb")
```

To run this example you'll need to download the following data sets which allow motifs to be 
scored in SCENIC: mm9-500bp-upstream-7species.mc9nr.feather and mm9-tss-centered-10kb-7species.mc9nr.feather,
which are located at https://resources.aertslab.org/cistarget/databases/mus_musculus/mm9/refseq_r45/mc9nr/gene_based.

We then initialize our SCENIC run. 
```{r}
scenicOptions <- initializeScenic(org = "mgi", dbDir = "mousedb",  nCores=1)
```

We will then load toy data sets derived from diaphragm cells from the Tabula Muris Consortium.
The data sets consist of the scRNA-seq count data,`r countData`, and normalized
scRNA-seq data, `r normData`.
They also include the vector of size factors, `r sizeFactors`, which will serve
as covariates in the zero-inflation part of the model.
The names of the transcription factors are found in `r TF_names`. These 
transcription factors serve as covariates in Poisson part of the model.
The names of the genes used in ZIRFs are given in the vector `r target_names`.
```{r load_toy_data, eval = TRUE}
data("countDat")
data("normDat")
data("sizeFactors")
data("TF_names")
data("target_names")
```
We then fit ZIRFs using TFs as covariates and gene expression as the outcome.
We do this for all genes in `r target_names`.
```{r run_zirf_genie3, eval = TRUE}
zirf_genie3_fit <- zirf_genie3(normDat, countDat, regulators = TF_names,
                               targets = union(target_names, TF_names),
                               z = as.matrix(sizeFactors),
                               nCores = 2)
```

We use the output of ZIRFs to run the SCENIC pipeline.
```{r run_SCENIC}
SCENIC_results <- scRNAzirf::run_SCENIC(scenicOptions = scenicOptions,
                                        weightMatrix_list = zirf_genie3_fit,
                                        log_scale = F,
                                        exprMatrix = normDat)
```

We can then look at the results of SCENIC by accessing the "output" folder.
The regulons and AUCell data can also be read in as follows. 
```{r extract_AUCell}
regulons <- readRDS("int/3.1_regulons_forAUCell.Rds")
AUCell_dat <- readRDS("int/3.4_regulonAUC.Rds")
```
