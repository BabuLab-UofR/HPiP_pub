<!-- badges: start -->
[![](https://img.shields.io/badge/lifecycle-stable-green.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
[![](https://img.shields.io/github/last-commit/BabuLab-UofR/HPiP_pub.svg)](https://github.com/BabuLab-UofR/HPiP_pub/commits/main)
[![License: MIT (&gt;=
3)](https://img.shields.io/badge/license-MIT-blue.svg)](https://cran.r-project.org/web/licenses/MIT)
<!-- badges: end -->

# Scripts for "HPiP: an R/Bioconductor package for predicting hostpathogen protein-protein interactions from protein sequences using an ensemble machine learning"
Despite arduous and time-consuming experimental efforts, protein-protein interactions (PPIs) for many pathogenic microbes with their human host are still unknown, limiting our ability to understand the intricate interactions during infection and the identification of therapeutic targets. Since computational tools offer a promising alternative, we developed an R/Bioconductor package, HPiP (Host-Pathogen Interaction Prediction) software toolkit with a series of amino acid sequence property descriptors and an ensemble machine-learning (ML) classifiers to define the yet unmapped interactions between pathogen and host proteins.

# This Repository Contains:
1. [Scripts](https://github.com/mrbakhsh/HPiP_pub/blob/main/R/Scripts_Model1.R) to run case study 1 (using SARS-CoV-1-human host PPIs as a training set)
2. [Scripts](https://github.com/mrbakhsh/HPiP_pub/blob/main/R/Scripts_Model2.R) to run case study 2 (using SARS-CoV-2-human host PPIs as a training set)
3. [Scripts](https://github.com/mrbakhsh/HPiP_pub/blob/main/R/Scripts_Model3.R) to run case study 3 (using Mtb-human host PPIs as a training set)

# Sample Data Description:
## Case study 1:
- `Trainingset_priotFC.csv` includes Sars-CoV-1-human PPIs. To construct this set, positive set was retrieved from  [Gordon, Hiatt, *et al.,* 2020](https://www.science.org/doi/10.1126/science.abe9403?url_ver=Z39.88-2003&rfr_id=ori%3Arid%3Acrossref.org&rfr_dat=cr_pub++0pubmed&). Negative sampling was used to construct negative PPIs from the positive ones using the `get_negativePPI` function.Both sets were mixed to construct labelled training set.

- `Testset_priotFC.csv` includes Sars-CoV-2-human PPIs. This set contains high-confidence PPI interactions between SARS-CoV-2 and human proteins were extracted from [Gordon, Hiatt, *et al.,* 2020](https://www.science.org/doi/10.1126/science.abe9403?url_ver=Z39.88-2003&rfr_id=ori%3Arid%3Acrossref.org&rfr_dat=cr_pub++0pubmed&). Negative sampling was then used to construct negative PPIs from the positive ones using the `get_negativePPI` function. Both sets were mixed to construct test set.

## Case study 2:
- `Trainingset_priotFC.csv` includes Sars-CoV-2-human PPIs. To construct this set, all the postive interactions including CoV-2-human PPIs deposited in BioGRID was retrived using `get_postivePPI` function provided in the HPiP package. Negative sampling was used to construct negative PPIs from the positive ones. Both sets were then mixed to construct labelled training set. Only (20%) of the data was used for model construction. This set was further split into a training set (70%) for model training and a test set (30%) for performance assessment of the classifiers.

## Case study 3:
- `Trainingset_priotFC.csv` includes Mtb-human PPIs. To construct this set, all the postive interactions including Mtb-human PPIs was retrived from [Penn, Bennett H, *et al.,* 2018](chrome-extension://efaidnbmnnnibpcajpcglclefindmkaj/viewer.html?pdfurl=https%3A%2F%2Fwww.cell.com%2Fmolecular-cell%2FpdfExtended%2FS1097-2765(18)30557-4&clen=4948760&pdffilename=mmc11.pdf), while negative instanses were constructed usig negative sampling.Both sets were mixed to construct labelled training set. This set was further split into a training set (70%) for model training and a test set (30%) for performance assessment of the classifiers.

