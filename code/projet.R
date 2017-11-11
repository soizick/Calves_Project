# *********************************************************************************************
# ***************************************** Projet Veaux **************************************
# *********************************************************************************************

library(mixOmics)
library(metagenomeSeq)

####################################################

df_abundances <- read.delim("../data/abondances.csv", sep = ",", 
                            stringsAsFactors = FALSE)
summary(df_abundances[ ,1:15])

abundances <- df_abundances[ ,grep("A_", colnames(df_abundances))] +
  df_abundances[ ,grep("B_", colnames(df_abundances))]
dim(abundances)

condition <- rep("LBA", ncol(abundances))
condition[grep("EN", colnames(abundances))] <- "EN"
condition<-as.factor(condition)
table(condition)

id_abundances <- as.character(colnames(abundances))
id_abundances <- sapply(id_abundances, function(ac) 
  substr(ac, nchar(ac) - 19, nchar(ac) - 18))
id_abundances <- gsub("A", "0", id_abundances)
id_abundances <- gsub("N", "0", id_abundances)
table(id_abundances)

####################################################
# Normalisations :

df <- abundances
names(df) <- paste0("sample", 1:ncol(df))
ggplot(df, aes(x = sample1)) + geom_histogram(bins = 50) + scale_y_log10() +
  theme_bw() + xlab("counts (sample 1)") + ggtitle("Histogram of sample1 column with log scale")

boxplot(df,main="Boxplot of sample1 column with no transfo")

# log sur donnees brutes

abundances_log<-apply(abundances, 2,function(x){log(x+1)})
df <- as.data.frame(abundances_log)
names(df) <- paste0("sample", 1:ncol(df))
ggplot(df, aes(x = sample1)) + geom_histogram(bins = 50) +
  theme_bw() + xlab("counts (sample 1)") + ggtitle("Histogram of sample1 column with log transfo")

boxplot(df,main="Boxplot of sample1 column with log transfo")

# TSS

abundances_TSS <- apply(abundances, 2, function(asample)
  asample / sum(asample))
df <- as.data.frame(abundances_TSS)
names(df) <- paste0("sample", 1:ncol(df))
ggplot(df, aes(x = sample1)) + geom_histogram(bins = 50) + scale_y_log10() +
  theme_bw() + xlab("counts (sample 1)")+ ggtitle("Histogram of sample1 column in TSS with log scale")

boxplot(df,main="Boxplot of sample1 column in TSS")

# TSS + CLR

abundances_CLR <- logratio.transfo(abundances_TSS, logratio = "CLR", 
                                   offset = 1)
class(abundances_CLR) <- "matrix"
df <- data.frame(abundances_CLR)
names(df) <- paste0("sample", 1:ncol(df))
ggplot(df, aes(x = sample1)) + geom_histogram(bins = 50) + theme_bw() + 
  xlab("counts (sample 1)")+ ggtitle("Histogram of sample1 column in TSS+CLR with log scale")

boxplot(df,main="Boxplot of sample1 column in TSS+CLR")

# TSS + ILR

abundances_ILR <- logratio.transfo(abundances_TSS, logratio = "ILR", 
                                   offset = 1)
class(abundances_ILR) <- "matrix"
df <- data.frame(abundances_ILR)
names(df) <- paste0("sample", 1:ncol(df))
ggplot(df, aes(x = sample1)) + geom_histogram(bins = 50) + theme_bw() + 
  xlab("counts (sample 1)")+ ggtitle("Histogram of sample1 column in TSS+ILR with log scale")

boxplot(df,main="Boxplot of sample1 column in TSS+ILR")

# CSS

abundances_CSS <- newMRexperiment(abundances)
abundances_CSS <- cumNorm(abundances_CSS)

df <- data.frame(MRcounts(abundances_CSS))
names(df) <- paste0("sample", 1:ncol(df))
ggplot(df, aes(x = sample1)) + geom_histogram(bins = 50) + theme_bw() + 
  xlab("counts (sample 1)")+ ggtitle("Histogram of sample1 column in CSS with log scale")

boxplot(df,main="Boxplot of sample1 column in CSS")

####################################################
# ACP

pca_raw <- pca(t(abundances_log), ncomp = ncol(abundances), 
               logratio = 'none')
plot(pca_raw)

plotIndiv(pca_raw, 
          comp = c(1,2),
          pch = 16, 
          ind.names = FALSE, 
          group = condition, 
          col.per.group = color.mixo(1:2),
          legend = TRUE,
          title="ACP on raw data with log transfo")

pca_CLR <- pca(t(abundances_TSS)+1, ncomp = ncol(abundances_TSS), 
               logratio = 'CLR')
plot(pca_CLR)

plotIndiv(pca_CLR, 
          comp = c(1,2),
          pch = 16, 
          ind.names = FALSE, 
          group = condition, 
          col.per.group = color.mixo(1:2),
          legend = TRUE,
          title="ACP with CLR transfo")

pca_ILR <- pca(t(abundances_TSS) + 1, ncomp = ncol(abundances_TSS) - 1, # pourquoi '-1' ?
               logratio = 'ILR')
plot(pca_ILR)

plotIndiv(pca_ILR, 
          comp = c(1,2),
          pch = 16, 
          ind.names = FALSE, 
          group = condition, 
          col.per.group = color.mixo(1:2),
          legend = TRUE,
          title="ACP with ILR transfo")

counts_CSS <- MRcounts(abundances_CSS)
pca_CSS <- pca(t(counts_CSS), ncomp = ncol(counts_CSS))
plot(pca_CSS)

plotIndiv(pca_CSS, 
          comp = c(1,2),
          pch = 16, 
          ind.names = FALSE, 
          group = condition, 
          col.per.group = color.mixo(1:2),
          legend = TRUE,
          title="ACP with CSS transfo")

####################################################
# PLS-DA :
# sur TSS+CLR
diverse.plsda = plsda(X = t(abundances_TSS)+1, condition, ncomp = nlevels(condition), logratio = 'CLR')
diverse.perf.plsda = perf(diverse.plsda, validation = 'Mfold', folds = 5,
                          progressBar = FALSE, nrepeat = 60)
plot(diverse.perf.plsda, overlay = 'measure', sd = TRUE)

plotIndiv(diverse.plsda , comp = c(1,2), ind.names = FALSE, 
          ellipse = TRUE, legend = TRUE, title = 'HMP Most Diverse, PLSDA comp 1 - 2')

# this chunk takes ~2 min to run with 6.3.0
set.seed(33)  # for reproducible results for this code
diverse.tune.splsda = tune.splsda(t(abundances_TSS)+1, condition, ncomp = nlevels(condition), 
                                  logratio = 'CLR',
                                  test.keepX = c(seq(5,150, 5)), validation = 'Mfold', 
                                  folds = 5, dist = 'mahalanobis.dist', nrepeat = 10,
                                  progressBar = TRUE)

plot(diverse.tune.splsda)

#log

diverse.plsda = plsda(X = t(abundances_log)+1, condition, ncomp = nlevels(condition))
diverse.perf.plsda = perf(diverse.plsda, validation = 'Mfold', folds = 5,
                          progressBar = FALSE, nrepeat = 60)
plot(diverse.perf.plsda, overlay = 'measure', sd = TRUE)

plotIndiv(diverse.plsda , comp = c(1,2), ind.names = FALSE, 
          ellipse = TRUE, legend = TRUE, title = 'HMP Most Diverse, PLSDA comp 1 - 2')

set.seed(33)  
diverse.tune.splsda = tune.splsda(t(abundances_log)+1, condition, ncomp = nlevels(condition),
                                  test.keepX = c(seq(5,150, 5)), validation = 'Mfold', 
                                  folds = 5, dist = 'mahalanobis.dist', nrepeat = 10,
                                  progressBar = TRUE)

plot(diverse.tune.splsda)

select.keepX = diverse.tune.splsda$choice.keepX[1:2]
select.keepX

# on doit enlever le veau 29 : -which(id_abundances=="29")
diverse.splsda = splsda(t(abundances_log[,-which(id_abundances=="29")])+1,
                        condition[-which(id_abundances=="29")], 
                        ncomp = nlevels(condition), 
                        multilevel = id_abundances[-which(id_abundances=="29")], keepX = select.keepX) 

plotIndiv(diverse.splsda, comp = c(1,2),
          ind.names = FALSE, 
          ellipse = TRUE, legend = TRUE,
          title = 'HMP Most Diverse, sPLSDA comp 1 - 2')

#Interpretation: on regarde quelles bacteries influencent le plus l'axe 1:

head(selectVar(diverse.splsda, comp = 1)$value)
list_var<-head(selectVar(diverse.splsda, comp = 1)$name)
indices<-substr(list_var,2,nchar(list_var))
names<-df_abundances[indices,1]
cat("* * * Nom des variables qui expliquent le mieux la différences entre les 2 clusters : * * *",names,fill=T)

plotLoadings(diverse.splsda, comp = 1, method = 'mean', contrib = 'max',size.title = 1)
plotLoadings(diverse.splsda, comp = 2, method = 'mean', contrib = 'max',size.title = 1)

# heatmap

diverse.splsda = splsda(t(abundances_log[,-which(id_abundances=="29")])+1,
                        condition[-which(id_abundances=="29")], 
                        ncomp = 1, 
                        multilevel = id_abundances[-which(id_abundances=="29")], keepX = select.keepX) 

cim(diverse.splsda, row.sideColors = color.mixo(condition[-which(id_abundances=="29")]),row.names=id_abundances[-which(id_abundances=="29")])

par(mfrow=c(1,1))
plotArrow(diverse.splsda, legend=T)


#****************************************************************************************************#
#*********************************** REGRESSION LOGISTIQUE ******************************************#
#****************************************************************************************************#

df_pathogenes<-read.delim("../data/pathogenes.csv",sep = ",")

ID_veaux<-as.character(df_pathogenes[1:paste(dim(df_pathogenes)[1]/2),1])
ID_veaux_totalp <- as.character(df_pathogenes[,1])
ID_veaux_totalp <- gsub(" ", "", ID_veaux_totalp)

# on prend que les LBA:

df_pathogenes_LBA<-df_pathogenes[which(df_pathogenes[,9]=="LBA"),]
rownames(df_pathogenes_LBA)<-df_pathogenes_LBA[,1]
df_pathogenes_LBA<-df_pathogenes_LBA[,-c(1,9)]
df_pathogenes_LBA<-ifelse(df_pathogenes_LBA=='p',1,0)

df_abundances_LBA<-df_abundances[,which(condition=="LBA")]

reg <- glm(df_pathogenes_LBA[,1]~ ., 
           data = as.data.frame(t(df_abundances_LBA)), family = binomial(logit))
summary(reg)

length(unique(df_abundances[,1]))
