# Scripts for "HPiP: an R/Bioconductor package for predicting hostpathogen protein-protein interactions from protein sequences using an ensemble machine learning"
Despite arduous and time-consuming experimental efforts, protein-protein interactions (PPIs) for many pathogenic microbes with their human host are still unknown, limiting our ability to understand the intricate interactions during infection and the identification of therapeutic targets. Since computational tools offer a promising alternative, we developed an R/Bioconductor package, HPiP (Host-Pathogen Interaction Prediction) software toolkit with a series of amino acid sequence property descriptors and an ensemble machine-learning (ML) classifiers to define the yet unmapped interactions between pathogen and host proteins.

# This Repository Contains:
1. [Scripts](https://github.com/mrbakhsh/HPiP_pub/blob/main/R/Scripts_Model1.R) to run case study 1 (using SARS-CoV-1-human host PPIs as a training set)
2. [Scripts](https://github.com/mrbakhsh/HPiP_pub/blob/main/R/Scripts_Model2.R) to run case study 2 (using SARS-CoV-2-human host PPIs as a training set)
3. [Scripts](https://github.com/mrbakhsh/HPiP_pub/blob/main/R/Scripts_Model3.R) to run case study 3 (using Mtb-human host PPIs as a training set)

# Sample Data Description:
## Case study 1:
- `Trainingset_priotFC.csv` includes Sars-CoV-1-human PPIs
- `Testset_priotFC.csv` includes Sars-CoV-2-human PPIs
## Case study 2:
- `Trainingset_priotFC.csv` includes Sars-CoV-2-human PPIs
## Case study 3:
- `Trainingset_priotFC.csv` includes Mtb-human PPIs

