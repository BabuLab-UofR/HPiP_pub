install.packages("seqinr")
library(readr)
library(HPiP)
library(dplyr)
library(seqinr)
library(data.table)

setwd("~/Library/CloudStorage/OneDrive-Personal/Desktop/cran_packages/HPiP_newAnalysis/githubFile/data/Model_1")
tdata <- fread("Trainingset_priotFC.csv")
pdata <- fread("predict_cov2_hPPI.csv")
uniq_p <- setdiff(pdata$H_unip, tdata$host_prot)
pdata <- 
   filter(pdata, H_unip %in% uniq_p) %>%
   .$H_unip




if (!require("BiocManager", quietly = TRUE))
   install.packages("BiocManager")
BiocManager::install("Biostrings")
library(Biostrings)
library()
setwd("~/Library/CloudStorage/OneDrive-Personal/Desktop/cran_packages/HPiP_newAnalysis/githubFile/data/Model_1/FASTA")
host_train <- read.fasta("host_TrainingSet.fasta")
host_train <- host_train[1:200]
host_pred <- read.fasta("host_PredictedSet.fasta")
nms <- setdiff(names(host_pred), names(host_train))
host_pred <- host_pred[1:300]
host_pred <- 
   host_pred[names(host_pred) %in% nms]



seq1string <- lapply(host_train, function(x) toupper(c2s(x)))
names(seq1string) <- paste(names(seq1string), "~")
seq2string <- lapply(host_pred, function(x) toupper(c2s(x)))


library(Biostrings)
data("BLOSUM62")

globalAlign<-
  unlist(lapply(seq1string, function(X) {
    lapply(seq2string, function(Y) {
      pairwiseAlignment(X, Y, substitutionMatrix = "BLOSUM62",
                        gapOpening = 10.0,
                        gapExtension=0.5)
    })
  }), recursive = FALSE)


dd <- lapply(globalAlign, function(x) pid(x, type = "PID3"))

dd <-
  do.call(rbind, dd) %>%
  as.data.frame(.) %>%
  tibble::rownames_to_column("id")
filter(V1 >= 90)
dd1 <- separate(dd, id, c("cov","mars"), sep = "~", remove = F)

