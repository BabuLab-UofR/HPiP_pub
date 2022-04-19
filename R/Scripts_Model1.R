library(readr)
library(HPiP)
library(dplyr)

#### 1. open training data , extract fasta sequences ####
t_data <- 
    read_csv("Trainingset_priotFC.csv")


# open pre-constructed FASTA seqeunces for viral
v_seq <- 
    read_csv("SARS_COV_1_FASTA.csv")

# retrieve fasta sequences for host proteins 
hostid <- unique(t_data$host_prot)
hostseq <-getFASTA(hostid)
hostseq <- do.call(rbind, hostseq) %>%
    as.data.frame(.) %>%
    tibble::rownames_to_column("UniprotKB")

#### 2. compute the features ####

#viral side
feat_v <- calculateAAC(v_seq)
feat_v <- 
    feat_v[order(match(feat_v$identifier,unique(t_data$viral_prot))),]
x_viral1 <- as.matrix(feat_v[, -1])
row.names(x_viral1) <- feat_v$identifier


#host side
feat_h <- calculateAAC(hostseq)
feat_h <- 
    feat_h[order(match(feat_h$identifier,hostid)),]
x_host1 <- as.matrix(feat_h[, -1])
row.names(x_host1) <- feat_h$identifier

#### 3. map the training to those protiens with available fasta sequences ####

t_data <- 
    filter(t_data, host_prot %in% feat_h$identifier)
t_data <- 
    filter(t_data, viral_prot %in% feat_v$identifier)

#### 5. map features to training data ####

x1_viral <- 
    matrix(NA, nrow = nrow(t_data), ncol = ncol(x_viral1))
for (i in 1:nrow(t_data)) x1_viral[i, ] <- x_viral1[which(t_data$viral_prot[i] == feat_v$identifier), ]

x1_host <- matrix(NA, nrow = nrow(t_data), ncol = ncol(x_host1))
for (i in 1:nrow(t_data)) x1_host[i, ] <- x_host1[which(t_data$host_prot[i] == feat_h$identifier), ]

#### 6. Generate the HPI descriptors ####
x <- getHPI(x1_viral,x1_host, type = "combine")
x_train <- as.data.frame(x)
x_train <- cbind(t_data$PPI, t_data$class, x_train)
colnames(x_train)[1:2] <- c("PPI", "class")


#### 7. open test set, , extract fasta sequences ####
e_data <- 
    read_csv("Testset_priotFC.csv")

# open pre-constructed FASTA seqeunces for viral
v_seq <- 
    read_csv("SARS_COV_2_FASTA.csv")

# retreive fasta sequences for host proteins 
hostid <- unique(e_data$host_prot)
hostseq <-getFASTA(hostid)
hostseq <- do.call(rbind, hostseq) %>%
    as.data.frame(.) %>%
    tibble::rownames_to_column("UniprotKB")

#### 8. compute the features (test set) ####

#viral side
feat_v <- calculateAAC(v_seq)
feat_v <- 
    feat_v[order(match(feat_v$identifier,unique(e_data$viral_prot))),]
x_viral1 <- as.matrix(feat_v[, -1])
row.names(x_viral1) <- feat_v$identifier


#host side
feat_h <- calculateAAC(hostseq)
feat_h <- 
    feat_h[order(match(feat_h$identifier,hostid)),]
x_host1 <- as.matrix(feat_h[, -1])
row.names(x_host1) <- feat_h$identifier

#### 9. map the testset to those protiens with available fasta sequences (test_set) ####

e_data <- 
    filter(e_data, host_prot %in% feat_h$identifier)
e_data <- 
    filter(e_data, viral_prot %in% feat_v$identifier)

#### 10. map features to test data ####

x1_viral <- 
    matrix(NA, nrow = nrow(e_data), ncol = ncol(x_viral1))
for (i in 1:nrow(e_data)) x1_viral[i, ] <- x_viral1[which(e_data$viral_prot[i] == feat_v$identifier), ]

x1_host <- matrix(NA, nrow = nrow(e_data), ncol = ncol(x_host1))
for (i in 1:nrow(e_data)) x1_host[i, ] <- x_host1[which(e_data$host_prot[i] == feat_h$identifier), ]

#### 11. Generate the HPI descriptors (test_set) ####
x <- getHPI(x1_viral,x1_host, type = "combine")
x_test <- as.data.frame(x)
x_test <- cbind(e_data$PPI, e_data$class, x_test)
colnames(x_test)[1:2] <- c("PPI", "class")

#### 12. ML prediction ####

#create data with features
features <- 
    rbind(x_test[, -2], x_train[,-2])
#create label data
gd <- 
    x_train[,c(1,2)]


set.seed(101)
ml_output <- pred_ensembel(features,
                           gd,
                           classifier = c("glm","svmRadial", "ranger"),
                           resampling.method = "cv",
                           ncross = 5,
                           verboseIter = TRUE,
                           plots = TRUE)



