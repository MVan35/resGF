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
# get data - (https://github.com/MVan35/resGF_method_paper)
"""
Using data obtain from the following study:
Bekkevold, D., Piper, A., Campbell, R., Rippon, P., Wright, R. M., Crundwell, C., ... & Maltby, A. (2021). Genetic stock identification of sea trout (Salmo trutta L.) along the British North Sea Coast shows prevalent long-distance migration. ICES Journal of Marine Science, 78(3), 952-966.
"""

library(raster)
## create a genind object
geno <- read.table("genotypes.txt", header = T, sep = '\t', row.names=1) # no tamar
x <- pegas::as.loci(geno)
x <-pegas::loci2genind(x, ploidy = 2, na.alleles = c("9"), unphase = TRUE)

## obtain spatial data
sample.coord  <-read.table("lat_lon.txt", header=T, stringsAsFactors=F, sep="\t", row.names=1)
colnames(sample.coord) <- c("x", "y")
sample.coord.sp <- SpatialPointsDataFrame(sample.coord[,c('y','x')], proj4string=CRS(crs.wgs), data=sample.coord)

## Obtain some environmental data
library(sdmpredictors)
e <- c(-3, 16, 50, 60) # study area

# load layers using sdmpredictors
#environment.bottom <- load_layers( layercodes = c("BO_bathymean",
#                                                  "BO_calcite",
#                                                  "BO_ph",
#                                                  "BO_dissox",
#                                                  "BO_chlorange",
#                                                  "BO2_ironltmin_bdmax",
#                                                  "BO2_carbonphytomean_bdmean",
#                                                  "BO2_chlomax_bdmax",
#                                                  "BO2_curvelmean_ss",
#                                                  "BO2_curvelmean_bdmin"
#                                                  "MS_sss09_5m"), equalarea=FALSE, rasterstack=TRUE)
# or - load tif filesdirectly
BO_bathymean <- raster("BO_bathymean_lonlat.tif")
BO_calcite <- raster("BO_calcite_lonlat.tif")
BO_chlorange <- raster("BO_chlorange.tif")
BO_dissox <- raster("BO_dissox_lonlat.tif")
BO_ph <- raster("BO_ph_lonlat.tif")
environment.bottom <- stack(BO_bathymean, BO_calcite, BO_chlorange, BO_dissox, BO_ph)

environment.bottom.crop <- crop(environment.bottom, e, snap="out")
# obtain a raster stack
clipped_stack <- crop(environment.bottom.crop, e, snap="out")
## create a resistance surface using resistantGF
crs.wgs <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

# Gradient Forest
library(gradientForest)
# extract points from
clim.points <- raster::extract(clipped_stack, sample.coord.sp)
#generates the PCNMs
library(vegan)
pcnm <- pcnm(dist(sample.coord.sp@data))
keep <- round(length(which(pcnm$value > 0))/2)
pcnm.keep <- scores(pcnm)[,1:keep]  #keep half of positive ones as suggested by some authors

nbvar <- length(colnames(clim.points))
env.gf <- cbind(clim.points[ , c(seq(1, nbvar, by=1)) ], pcnm.keep)
maxLevel <- log2(0.368*nrow(env.gf)/2)
env.gf <- as.data.frame(env.gf)
snp <- as.data.frame(x$tab)

# run gradient forest
gf <- gradientForest(cbind(env.gf, snp), predictor.vars=colnames(env.gf), response.vars=colnames(snp), ntree=500, maxLevel=maxLevel, trace=T, corr.threshold=0.50)
# generate resistance surface
single_r <- resGF(gf, clipped_stack, save.image = T)

```
