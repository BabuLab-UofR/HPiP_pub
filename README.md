<!-- badges: start -->
[![](https://img.shields.io/badge/lifecycle-stable-green.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
[![](https://img.shields.io/github/last-commit/BabuLab-UofR/HPiP_pub.svg)](https://github.com/BabuLab-UofR/HPiP_pub/commits/main)
[![License: MIT (&gt;=
3)](https://img.shields.io/badge/license-MIT-blue.svg)](https://cran.r-project.org/web/licenses/MIT)
<!-- badges: end -->

# Scripts for "HPiP: an R/Bioconductor package for predicting hostpathogen protein-protein interactions from protein sequences using an ensemble machine learning"
Despite arduous and time-consuming experimental efforts, protein-protein interactions (PPIs) for many pathogenic microbes with their human host are still unknown, limiting our ability to understand the intricate interactions during infection and the identification of therapeutic targets. Since computational tools offer a promising alternative, we developed an R/Bioconductor package, HPiP (Host-Pathogen Interaction Prediction) software toolkit with a series of amino acid sequence property descriptors and an ensemble machine-learning (ML) classifiers to define the yet unmapped interactions between pathogen and host proteins.

# Installation
This package required R version 4.1 or higher. If you are using an older version of R you will be prompted to upgrade when installing the package. 

The official release of HPiP is available on Bioconductor. You can install the `HPiP` from bioconductor using:

```r
if(!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager") 
}
BiocManager::install("HPiP")
```
To install the development version in `R`, run:
  
```r
if(!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools") 
}
devtools::install_github("mrbakhsh/HPiP")
```
# HPiP Help Page
- The main page: https://github.com/BabuLab-UofR/HPiP

- Vignette 
```r
browseVignettes("HPiP")
```

# This Repository Contains:
1. [Scripts](https://github.com/mrbakhsh/HPiP_pub/blob/main/R/Scripts_Model1.R) to run case study 1 
2. [Scripts](https://github.com/mrbakhsh/HPiP_pub/blob/main/R/Scripts_Model2.R) to run case study 2 
3. [Scripts](https://github.com/mrbakhsh/HPiP_pub/blob/main/R/Scripts_Model3.R) to run case study 3 

# Sample Data Description:

## Case study 1:
- `Trainingset_priotFC.csv` includes Sars-CoV-1-human PPIs. To construct this set, positive set was retrieved from  [Gordon, Hiatt, *et al.,* 2020](https://www.science.org/doi/10.1126/science.abe9403?url_ver=Z39.88-2003&rfr_id=ori%3Arid%3Acrossref.org&rfr_dat=cr_pub++0pubmed&). Negative sampling was used to construct negative PPIs from the positive ones using the `get_negativePPI` function.Both sets were mixed to construct labelled training set.

- `Testset_priotFC.csv` includes Sars-CoV-2-human PPIs. This set contains high-confidence PPI interactions between SARS-CoV-2 and human proteins were extracted from [Gordon, Hiatt, *et al.,* 2020](https://www.science.org/doi/10.1126/science.abe9403?url_ver=Z39.88-2003&rfr_id=ori%3Arid%3Acrossref.org&rfr_dat=cr_pub++0pubmed&). Negative sampling was then used to construct negative PPIs from the positive ones using the `get_negativePPI` function. Both sets were mixed to construct test set.

## Case study 2:
- `Trainingset_priotFC.csv` includes Sars-CoV-2-human PPIs. To construct this set, all the postive interactions including CoV-2-human PPIs deposited in [BioGRID](https://thebiogrid.org/) was retrived using `get_postivePPI` function provided in the HPiP package. Negative sampling was used to construct negative PPIs from the positive ones. Both sets were then mixed to construct labelled training set. Only (20%) of the data was used for model construction. This set was further split into a training set (70%) for model training and a test set (30%) for performance assessment of the classifiers.

## Case study 3:
- `Trainingset_priotFC.csv` includes Mtb-human PPIs. To construct this set, all the postive interactions including Mtb-human PPIs was retrived from
[Penn, Bennett H, *et al.,* 2018](https://pubmed.ncbi.nlm.nih.gov/30118682/), while negative instanses were constructed usig negative sampling. Both sets were mixed to construct labelled training set. This set was further split into a training set (70%) for model training and a test set (30%) for performance assessment of the classifiers.

# Citation 
Rahmatbakhsh,M. et al. (2022) HPiP: an R/Bioconductor package for predicting host–pathogen protein–protein interactions from protein sequences using ensemble machine learning approach. Bioinforma. Adv., 2, vbac038.



