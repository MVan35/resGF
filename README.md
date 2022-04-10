# resGF

By design, resGF has the ability to handle multiple environmental predicators to generate resistance surface.

# install
```R
library(devtools)
install_github("MVan35/resGF")
```
# Description
Used in community ecology, Gradient Forest is an extension of random forest (RF), and has been implemented in genomic studies to model species genetic offset under future climatic scenarios. 

Gradient forest (GF), an extension of random forest that uses regression trees to fit a model of associations between individual responses variables (e.g. single-nucleotide polymorphisms (SNP)) to predicator variables (e.g. climate variables). 

Moving along environmental gradients, few changes in allelic frequencies might be observed while in some places large changes might occur. In the GF approach, partition of the data is distributed on either side of a split produced in the random forest. Split values are obtained from the ensemble of regression models and explain changes in allelic frequencies across each environmental predicators. The amount of variations, known as the ‘raw split importance’ (Fitzpatrick & Keller, 2015), is cumulatively added in a step-wise manner and a cumulative frequency curve,  is computed. Similarly as building a staircase, the importance values are cumulatively summed and each step strength is proportional to the importance split at that location. 

The predictive performance for each SNP is quantified using the out-of-bag proportion (R2) which provides a cross-validated estimate of the generalization error (Ellis et al., 2012). These goodness-of-fit R2 estimates allow to assess the relative contribution of each environmental variables in explaining changes in allelic frequencies. Lastly, the accuracy importance of predicators is determined as the decrease in performance when a predictor is randomly permuted.

In the resGF approach, these large steps along the gradient observed on the cumulative frequency curve  are considered as evidence of resistance impeding gene flow whereas flatter regions represent regions where gene flow is facilitated. Each environmental predicator can be weighted by its R2-weighted importance to build a multilayer resistance surface. 

# Usage

Generate resistance surface using allellic frequencies and landscape variables.


```R
## Gradient Forest
library(resGF)
# extract points from
clim.points <- raster::extract(clipped_stack, sample.coord.sp)
#generates the PCNMs
pcnm <- pcnm(dist(clim.points[,c('x', 'y')]))  
keep <- round(length(which(pcnm$value > 0))/2)
pcnm.keep <- scores(pcnm)[,1:keep]  #keep half of positive ones as suggested by some authors
clim.points <- cbind(p2, clim.points)

nbvar <- length(colnames(clim.points))
env.gf <- cbind(clim.points[ , c(seq(3, nbvar, by=1)) ], pcnm.keep)
library(gradientForest)
maxLevel <- log2(0.368*nrow(env.gf)/2)
env.gf <- as.data.frame(env.gf)
snp <- as.data.frame(snp)

gf <- gradientForest(cbind(env.gf, snp), predictor.vars=colnames(env.gf), response.vars=colnames(snp), ntree=500, maxLevel=maxLevel, trace=T, corr.threshold=0.50)

single_r <- resGF(gf_final, clipped_stack, save.image = T, GF_env_variables)
```
