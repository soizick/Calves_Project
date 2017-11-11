library(mixOmics)

# function to perform pre-filtering
low.count.removal = function(
  data, # OTU count data frame of size n (sample) x p (OTU)
  percent=0.01 # cutoff chosen
){
  keep.otu = which(colSums(data)*100/(sum(colSums(data))) > percent)
  data.filter = data[,keep.otu]
  return(list(data.filter = data.filter, keep.otu = keep.otu))
}

data(diverse.16S)
data.raw = diverse.16S$data.raw #the diverse.16s raw data include an offset of 1
head(data.raw)

result.filter = low.count.removal(data.raw, percent=0.01)
data.filter = result.filter$data.filter
length(result.filter$keep.otu) # check the number of variables kept after filtering


data("diverse.16S")

# the 16S normalised data
data.mixMC = diverse.16S$data.TSS

# the outcome  
Y = diverse.16S$bodysite

# unique ID of each individual for multilevel analysis
sample = diverse.16S$sample

pca.res = pca(data.mixMC, ncomp = 10, logratio = 'CLR', multilevel = sample)
#pca.res
plot(pca.res)


plotIndiv(pca.res, 
          comp = c(1,2), # the components to plot
          pch = 16, 
          ind.names = F, 
          group = Y, 
          col.per.group = color.mixo(1:3),
          legend = TRUE)


plsda <- plsda(X = data.mixMC, Y, ncomp = nlevels(Y), logratio = 'CLR', multilevel = sample)

perf.plsda <- perf(plsda, validation = 'Mfold', folds = 5,
                   progressBar = FALSE, nrepeat = 10)

plot(perf.plsda, overlay = 'measure', sd = TRUE)

plotIndiv(plsda , comp = c(1,2), ind.names = FALSE, 
          ellipse = TRUE, legend = TRUE, title = 'HMP Most Diverse, PLSDA comp 1 - 2')

plotIndiv(plsda , comp = c(1,3), ind.names = FALSE, 
          ellipse = TRUE, legend = TRUE, title = 'HMP Most Diverse, PLSDA comp 1 - 3')


# this chunk takes ~ 12 min to run
splsda.tune = tune.splsda(data.mixMC, Y, ncomp = 3, 
                          logratio = 'CLR',
                          multilevel = sample,
                          test.keepX = c(seq(5,150, 5)), validation = 'Mfold', 
                          folds = 5, dist = 'max.dist', nrepeat = 10)
# may show some convergence issues for some of the cases, it is ok for tuning