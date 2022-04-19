library(readr)
library(HPiP)
library(dplyr)

#### 1. open training data , extract fasta sequences ####
setwd("C:/Users/tom/OneDrive/Desktop/OneDrive/Desktop/cran_packages/HPiP_newAnalysis/githubFile/data/Model_2")
t_data <- 
    read_csv("Trainingset_priotFC.csv")


# open pre-constructed FASTA seqeunces for viral
v_seq <- 
    read_csv("SARS_COV_2_FASTA.csv")

# retreive fasta sequences for host proteins 
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


#### 7. Split data into training & test set ####

## 75% of the sample size
set.seed(101)
smp_size <- floor(0.70 * nrow(x))
## set the seed to make your partition reproducible
x_ind <- sample(seq_len(nrow(x_train)), size = smp_size)

gd <- x_train[x_ind, ]
gd <- gd[, c(1,2)]
#gd <- gd[c(1:50),]


features <- 
    x_train[,-2]

set.seed(101)
ml_output <- pred_ensembel(features,
                           gd,
                           classifier = c("glm","svmRadial", "ranger"),
                           resampling.method = "cv",
                           ncross = 10,
                           verboseIter = TRUE,
                           plots = TRUE)


predfd <- ml_output[["predicted_interactions"]]


df.bench <- left_join(predfd,x_train)
df.bench$pred_label <- 
    ifelse(df.bench$ensemble >= 0.5, "Pos", "Neg")
table(df.bench$class,df.bench$pred_label)

df.bench <- 
    filter(df.bench, class == "Positive" & pred == "Pos" )


pos <- 
    filter(df.bench, class == "Positive") %>%
    dplyr::select(4,5)
table(pos$class)




negat <- 
    filter(df.bench, class == "Negative") %>%
    dplyr::select(4,5)
table(negat$class)



dff <- rbind(pos,negat)





ggplot() +
    geom_density(dff, aes(x=ensemble, color=class),size = 1) 
    theme_classic() 
    scale_color_manual(values=c("#999999", "#E69F00")) +
    geom_line(dff2, aes(y=ratio))  # Divide by 10 to get the same range than the temperature
    




size.ranges <- pos %>% 
    mutate(cuts = cut(ensemble, seq(0,0.9, by = 0.1))) %>% 
    dplyr::group_by(cuts) %>%
    dplyr::summarise(n=n()) 


size.ranges2 <- negat %>% 
    mutate(cuts = cut(ensemble, seq(0,0.9, by = 0.1))) %>% 
    dplyr::group_by(cuts) %>%
    dplyr::summarise(n=n()) 


dff2<- 
    left_join(size.ranges,size.ranges2, by = c("cuts"))

dff2$ratio <- 
    dff2$n.x/dff2$n.y
ro2 <- c(0,0.2,0.4,0.6,0.8,1.0)

d <-cbind(ro2,dff2)



library(ggplot2)
ggplot(dff, aes(x=ensemble, color=class)) +
    geom_density(size = 1) +
    scale_color_manual(values=c("#999999", "#E69F00")) +
    theme_classic() +
    theme(text = element_text(size=14, color = "black")) +
    theme(axis.text.x = element_text(angle = 1, hjust = 1, color = "black")) +
    theme(axis.ticks.length=unit(0.25, "cm")) +
    geom_point(data=d, aes(x = as.numeric(ro2), y=ratio/4),  shape = 16, col = "#006600", size = 4) +
    geom_line(data=d, aes(x = as.numeric(ro2), y=ratio/4),  shape = 16, col = "#006600", size = 1) +
    scale_x_continuous(expand = c(0,0), limits = c(0,1))+
    scale_y_continuous(name="data1", sec.axis = sec_axis(~ 2*., name="data2"),
                       expand = c(0,0), limits = c(0,4)) +
    theme(axis.title = element_text()) +
    theme(aspect.ratio = 1.8/2.2) +
    theme(text = element_text(size=15)) 
    



    








geom_path() +
    geom_point() +
    scale_y_continuous(name="data1", sec.axis = sec_axis(~ 2*., name="data2")) +
    scale_color_manual(name="z", values = mycolors) +
    theme(
        axis.title.y = element_text(color = mycolors["data1"]),
        axis.text.y = element_text(color = mycolors["data1"]),
        axis.title.y.right = element_text(color = mycolors["data2"]),
        axis.text.y.right = element_text(color = mycolors["data2"])
    )


