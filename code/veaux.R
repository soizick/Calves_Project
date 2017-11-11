# *********************************************************************************************
# ****************************** Test de MixOmics sur les donnees *****************************
# *********************************************************************************************

setwd("C:/Users/Soizick/Documents/5AGMM/Projet")
#setwd("C:/Users/Admin/Desktop/Projet_Innovation")

df_abundances <- read.delim("abondances.csv", sep = ",", stringsAsFactors = FALSE)
summary(df_abundances[ ,1:15])

df_pathogenes<-read.delim("pathogenes.csv",sep = ",")
ID_veaux<-as.character(df_pathogenes[1:paste(dim(df_pathogenes)[1]/2),1])

abundances <- df_abundances[ ,grep("A_", colnames(df_abundances))] + df_abundances[ ,grep("B_", colnames(df_abundances))]
abundances <- abundances[,-grep("LBA29",colnames(abundances))]
dim(abundances)

condition <- rep("LBA", ncol(abundances))
condition[grep("EN", colnames(abundances))] <- "EN"
condition <- as.factor(condition)
table(condition)

names<-colnames(abundances)

# ********************** ACP on individuals ************************

library(mixOmics)

pca_raw <- pca(t(abundances) + 1, ncomp = ncol(abundances), logratio = 'CLR')
plot(pca_raw)

plotIndiv(pca_raw, 
          comp = c(1,2), # the components to plot
          pch = 16, 
          ind.names = F, 
          group = condition, 
          col.per.group = color.mixo(1:2),
          legend = TRUE,
          title="PCA on individuals")

# ************************** Histograms ****************************
# work on the first column

h<-hist(as.numeric(abundances[,1])+1, breaks = 100,plot=T, main="Without normalisation")
plot(h$count, log="y", type='h', lwd=10, lend=2, main="Normalisation with log scale",ylab="log(frequences)")  # for having logaritmic scale
# the histogram is skewed, a normalisation is necessary

# ************************* Normalisation **************************

# Cumulative Sum Scaling normalisation (CSS)

library(metagenomeSeq)
data<-newMRexperiment(t(abundances))
p<- cumNormStat(data) 
data.cumnorm <- cumNorm(data, p=p)
h<-hist(as.numeric(MRcounts(data.cumnorm, norm=TRUE, log=TRUE)), breaks = 200,plot=F)
plot(h$count, log="y", type='h', lwd=10, lend=2, main="Normalisation with CSS (log scale)",ylab="log(frequences)")

# CLR normalisation

h<-hist(logratio.transfo(t(abundances[,1])+1, logratio = "CLR", offset = 0), breaks = 100, main="Normalisation with CLR",plot=F)
plot(h$count, log="y", type='h', lwd=10, lend=2, main="Normalisation with CLR (log scale)",ylab="log(frequences)")


# ************************* PLSDA **************************


diverse.plsda = plsda(X = t(abundances)+1, condition, ncomp = nlevels(condition), logratio = 'CLR')
diverse.perf.plsda = perf(diverse.plsda, validation = 'Mfold', folds = 5,
                          progressBar = FALSE, nrepeat = 60)
plot(diverse.perf.plsda, overlay = 'measure', sd = TRUE)


plotIndiv(diverse.plsda , comp = c(1,2), ind.names = FALSE, 
          ellipse = TRUE, legend = TRUE, title = 'HMP Most Diverse, PLSDA comp 1 - 2')


# this chunk takes ~2 min to run with 6.3.0
set.seed(33)  # for reproducible results for this code
diverse.tune.splsda = tune.splsda(t(abundances)+1, condition, ncomp = 2, 
                                  logratio = 'CLR',
                                  test.keepX = c(seq(5,150, 5)), validation = 'Mfold', 
                                  folds = 5, dist = 'mahalanobis.dist', nrepeat = 10,
                                  progressBar = TRUE)
# may show some convergence issues for some of the CV-runs, it is ok for tuning



# ID_veaux[11]<-substr(ID_veaux[11],1,7)
# 
# id<-match(substr(ID_veaux,6,nchar(ID_veaux)),substr(names,12,nchar(test)-18))
# 
# names[id]
