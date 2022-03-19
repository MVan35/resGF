# In general, try to group together related functions into the same .R file
# https://tinyheero.github.io/jekyll/update/2015/07/26/making-your-first-R-package.html

# tag above your function to indicate this function to be “exposed” to users to use.


#' Generate resistance surface from GF analysis
#'
#' This function generate a resistance surface based on gradient forest analysis.
#'
#' @param obj A GF object
#' @param raster_stack A raster stack oject containing the same variable used in the GF analysis
#' @param save.image Save individual conversion of landscape variables into resistance surfaces
#' @param results_dir Directory to save the transformation images
#' @return A resistance surface
#' @export
resGF <- function(obj,
                  raster_stack,
                  save.image = TRUE,
                  results_dir) {
  #GF_R_Total <- GF_total_importance(obj) # simplify this?
  mylength <- raster::ncell(raster_stack[[1]])
  myvector <- vector(mode="numeric", length=mylength)
  Env_pc <- 0
  s <- raster::stack()
  for (i in names(raster_stack)) {
    importance.df <- obj$res[obj$res$var==i,] # me
    mycumu <- getCU(importance.df, extendedForest::importance(obj, "Weighted")[i], i, obj)
    CU <- gradientForest::cumimp(obj, i)
    ymax <- max(CU$y)
    imp <- extendedForest::importance(obj)[i]
    resA <- obj$res[obj$res$var == i, ]
    splits <- resA$split
    w <- pmax(resA$improve.norm, 0)
    X <- na.omit(obj$X[,i])
    rX <- range(X)
    dX <- diff(rX)

    # importance density I(x)
    dImp <- density(splits, weight = w/sum(w), from = rX[1], to = rX[2])
    if((dX/dImp$bw) > 50)
      dImp <- density(splits, weight = w/sum(w), from = rX[1], to = rX[2], bw=dX/50)
    dImpNorm <- tryCatch( { normalize.density(dImp,imp,integrate=T)}, error = function(e){dImpNorm <- normalize.histogram(dImp,imp)} )

    # data density d(x)
    dObs <- density(X, from = rX[1], to = rX[2])
    if((dX/dObs$bw) > 50)
      dObs <- density(X, from = rX[1], to = rX[2], bw=dX/50)
    dObs <- whiten(dObs,lambda=0.9)
    # dObsNorm <- normalize.density(dObs,imp,integrate=T)
    dObsNorm <-tryCatch( { normalize.density(dObs,imp,integrate=T)}, error = function(e){dObsNorm <- normalize.histogram(dObs,imp)} )

    # standardized density f(x) = I(x)/d(x)
    dStd <- dImp
    dStd$y <- dImp$y/dObs$y
    # dStdNorm <- tryCatch( { normalize.density(dStd,imp,integrate=T)}, error = function(e){dStdNorm <- dStd} )
    dStdNorm <- try(normalize.density(dStd,imp,integrate=T), silent=TRUE) # this sometimes does not converge, added the silence
    if (class(dStdNorm) == "try-error")
      dStdNorm <- normalize.histogram(dStd,imp)

    # getting the y axis normalised to get the individual contribution of each variables
    dStd2 <- dStd
    dStd2$y <- dStd2$y * extendedForest::importance(obj, "Weighted")[i] / max(dStd2$y)
    dImpNorm2 <- dImpNorm
    dImpNorm2$y <- dImpNorm$y * extendedForest::importance(obj, "Weighted")[i] / max(dImpNorm$y)
    dObsNorm2 <- dObsNorm
    dObsNorm2$y <- dObsNorm$y * extendedForest::importance(obj, "Weighted")[i] / max(dObsNorm$y)
    dStdNorm2 <- dStdNorm
    dStdNorm2$y <- dStdNorm$y * extendedForest::importance(obj, "Weighted")[i] / max(dStdNorm$y)

    # dinv <- dStdNorm
    dinv <- dStdNorm2
    dinv$call <- "fp(x)"
    r <- raster_stack[[i]]
    #print(r)
    v0 <- raster::values(r)
    v1 <- as.matrix(v0)
    df1 <- as.data.frame(v0)
    df1$ID <- seq.int(nrow(df1))
    v1b <- v1[!rowSums(!is.finite(v1)),]
    uniq <- unique(v1b)

    # Transofomation
    v2 <- approx(dinv$x, dinv$y, uniq, rule=2)$y
    df2 <- data.frame(uniq,v2)
    names(df2)[names(df2) == "uniq"] <- "v0"
    myderivDF5 <- merge(df1, df2, by='v0', all = TRUE)
    myderivDF5 <- myderivDF5[order(myderivDF5$ID),]
    # Getting individual raster
    r_GF_tempmin  <- raster::raster()
    names(r_GF_tempmin)=i
    raster::extent(r_GF_tempmin) <- raster::extent(r)
    dim(r_GF_tempmin) <- dim(r)
    raster::crs(r_GF_tempmin) <- raster_stack
    # myderivDF5$v2[myderivDF5$v2<0] <- 0 # change $y into $v2
    # stdvalues <- myderivDF5$v2 * gf$overall.imp[v] / GF_R_Total # method 1
    stdvalues <- myderivDF5$v2  # method 2
    # stdvalues <- myderivDF5$v2 * imp / sum(extendedForest::importance(obj)) # method 2
    # stdvalues <- myderivDF5$v2 * extendedForest::importance(obj, "Weighted") / extendedForest::importance(obj, "Weighted")[i]

    myvector <- myvector + stdvalues # change $y into $v2
    raster::values(r_GF_tempmin) <- myderivDF5$v2 # change $y into $v2 - non-standardized values

    variable_name <- paste0(i, ".pdf")
    #Plotting
    bin=F
    nbin=101
    leg.panel=1
    barwidth=1
    leg.posn="topright"
    cex.legend=0.8
    line.ylab=1.5
    ci <- gradientForest::cumimp(obj,i,standardize=T, standardize_after=T)
    ci$y <- diff(c(0,ci$y))
    nice.names <- i
    # ci <- tryCatch( {normalize.histogram(ci, imp, bin=bin | !is.binned(obj), nbin=nbin)}, error = function(e){ci <- ci} )
    ci <-normalize.histogram(ci, imp, bin=bin | !is.binned(obj), nbin=nbin)

    ci$y <- ci$y * extendedForest::importance(obj, "Weighted")[i] / max(ci$y)


    if (save.image) {
      pdf(paste0(results_dir, variable_name)) # to get it save to file
      par(mfrow=c(2,3), oma=c(0,0,2,0))
      # plot(mycumu$x, mycumu$y, main= "Cummulative values")
      # lines(mycumu, col=2)
      plot(gradientForest::cumimp(obj,i,standardize=FALSE),type='n',main="",ylab="Cum. importance",xlab=i)
      # lines(gradientForest::cumimp(gf,i),type='l',col="black")
      # lines(gradientforest::cumimp(gf,i,standardize_after=TRUE),type='l',col="blue")
      lines(gradientForest::cumimp(obj,i,standardize=FALSE),type='l',col="black")
      # legend(par("usr")[1],par("usr")[4],legend=c("Standardize then Normalize (default)",
      #                                             "Normalize then Standardize","Normalize only"),col=c("black","blue","red"),lty=1,cex=0.8)
      raster::plot(r, main= "original raster")
      hist(v1, breaks=50, col = 'grey', xlab = i, main = "intial values" ) # value of the initial raster ; main = paste("intial values for" , v)
      plot(ci, type='h', col="grey60", xlim=range(splits), lwd=barwidth,
           ylim = c(0, max(dStdNorm2$y)*1.1), lend=2, xlab = i, ylab = "density")
      lines(dStdNorm2, col = "blue", lwd = 2)
      abline(h = mean(dStdNorm2$y)/mean(dStd2$y), lty = 2, col = "blue")
      # legend(leg.posn, legend = c("Density of splits", "Density of data", "Ratio of densities", "Ratio=1"), lty = c(1, 1, 1, 2), col = c("black", "red", "blue", "blue"), cex = cex.legend, bty = "n", lwd = 1)
      raster::plot(r_GF_tempmin, main= "resistance surface")
      hist(myderivDF5$v2, breaks=50, col = 'grey', xlab = i,  main = "transformed values") # value of the initial raster
      mtext(i, line=0, side=3, outer=TRUE, cex=1.2)
      dev.off() # to get it save to file

    }
    # percent = round(gf$overall.imp[v] / GF_R_Total * 100, 2)
    percent = round(imp / sum(extendedForest::importance(obj)) * 100, 2)
    Env_pc <- Env_pc + percent
    par(mfrow=c(2,3), oma=c(0,0,2,0))
    # plot(mycumu$x, mycumu$y, main= "Cummulative values")
    # lines(mycumu, col=2)
    plot(gradientForest::cumimp(obj,i,standardize=FALSE),type='n',main="",ylab="Cum. importance",xlab=i)
    # lines(gradientForest::cumimp(gf,i),type='l',col="black")
    # lines(gradientForest::cumimp(gf,i,standardize_after=TRUE),type='l',col="blue")
    lines(gradientForest::cumimp(obj,i,standardize=FALSE),type='l',col="black")
    raster::plot(r, main= "original raster")
    hist(v1, breaks=50, col = 'grey', xlab = i, main = "intial values" ) # value of the initial raster ; main = paste("intial values for" , v)
    plot(ci, type='h', col="grey60", xlim=range(splits), lwd=barwidth,
         ylim = c(0, max(dStdNorm2$y)*1.1), lend=2, xlab = i, ylab = "density")
    lines(dStdNorm2, col = "blue", lwd = 2)
    abline(h = mean(dStdNorm2$y)/mean(dStd2$y), lty = 2, col = "blue")
    raster::plot(r_GF_tempmin, main= "resistance surface")
    hist(myderivDF5$v2, breaks=50, col = 'grey', xlab = i,  main = "transformed values") # value of the initial raster
    mtext(paste0(i,': ',percent,'%'), line=0, side=3, outer=TRUE, cex=1.2)
    s <- raster::stack(s, r_GF_tempmin)
    print(paste0(i,': ',percent,'%'))
  }
  # creating a single resistant raster
  par(mfrow=c(1,1))
  single_r  <- raster::raster()
  names(single_r)='final_resistance'
  raster::extent(single_r) <- raster::extent(r)
  dim(single_r) <- dim(r)
  raster::values(single_r) <- myvector
  raster::crs(single_r) <- raster::crs(raster_stack)
  single_r <- climateStability::rescale0to1(single_r)
  # raster::plot(single_r)
  print(paste0('overall predicators: ',Env_pc,'%'))
  print(paste0('positive snp: ', obj$species.pos.rsq))
  return(single_r)
}



#' Return gradient forest importance
#'
#' This function generate a resistance surface based on gradient forest analysis.
#'
#' @param obj A GF object
#' @param raster_stack A raster stack oject containing the same variable used in the GF analysis
#' @return importance of the different variables
#' @export
get_importance <- function(obj, raster_stack) {
  Env_pc <- 0
  myimp <- NULL;
  for (i in names(raster_stack)) {
    importance.df <- obj$res[obj$res$var==i,] # me
    # https://r-forge.r-project.org/scm/viewvc.php/pkg/gradientForest/R/whiten.R?view=markup&revision=2&root=gradientforest&pathrev=16
    imp <- extendedForest::importance(obj)[i]
    percent = round(imp / sum(extendedForest::importance(obj)) * 100, 2)
    print(paste0(i,': ',percent,'%'))
    myobj <- c(i, percent)
    myimp <- rbind(myimp, myobj)
    Env_pc <- Env_pc + percent
  }
  print(paste0('overall predicators: ',Env_pc,'%'))
  myobj <- c('overall', Env_pc)
  myimp <- rbind(myimp, myobj)
  rownames(myimp) <- myimp[,1]
  myimp <- as.data.frame(myimp)
  myimp[,1] <- NULL
  return(myimp)
}

######################
### Other function
######################
lower <- function(matrix) {
  if (is.vector(matrix) == TRUE ||
      dim(matrix)[1] != dim(matrix)[2]) {
    warning("Must provide square distance matrix with no column or row names")
  }
  lm <- matrix[lower.tri(matrix)]
  return(lm)
}

getCU <- function(importance.df, Rsq, i, gf) {  # ME - made using both standardization
  agg <- with(importance.df, agg.sum(improve.norm, list(split), sort.it=TRUE))
  cum.split <- agg[,1]
  height <- agg[,2]
  dinv <- normalize(inverse(gf$dens[[i]])) # crucial to normalize this case
  dinv.vals <- approx(dinv$x, dinv$y, cum.split, rule=2)$y
  #par(mfrow=c(1,2), oma=c(0,0,2,0))
  #plot(dinv)
  #plot(dinv.vals)
  par(mfrow = c(1, 1))
  if (any(bad <- is.na(height))) {
    cum.split <- cum.split[!bad]
    height <- height[!bad]
    dinv.vals <- dinv.vals[!bad]
  }
  # height <- height/sum(height)*Rsq # if (standardize & !standardize_after
  height <- height * dinv.vals
  res <- list(x=cum.split, y=cumsum(height))
}

# aggregate
agg.sum <- function(x, by, sort.it = F)
{
  if(!is.data.frame(by))
    by <- data.frame(by)
  if(sort.it) {
    ord <- do.call("order", unname(by))
    x <- x[ord]
    by <- by[ord,  , drop=F]
  }
  logical.diff <- function(group)
    group[-1] != group[ - length(group)]
  change <- logical.diff(by[[1]])
  for(i in seq(along = by)[-1])
    change <- change | logical.diff(by[[i]])
  by <- by[c(T, change),  , drop = F]
  by$x <- diff(c(0, cumsum(x)[c(change, T)]))
  by
}

normalize <- function(f) {
  integral <- integrate(approxfun(f,rule=2),lower=min(f$x),upper=max(f$x),  stop.on.error = FALSE)$value
  f$y <- f$y/integral*diff(range(f$x));
  f
}

inverse <- function(dens) {dens$y <- 1/dens$y; dens}

normalize.histogram <- function(ci,integral=1,bin=F,nbin=101) {
  # scale y values so that histogram integrates to integral
  # optionally aggregate the y's into binned x ranges
  if (bin) {
    brks <- seq(min(ci$x),max(ci$x),len=nbin)
    xx <- cut(ci$x,breaks=brks,inc=T)
    yy <- tapply(ci$y,xx,sum)
    yy[is.na(yy)] <- 0
    ci <- list(x=0.5*(brks[-1]+brks[-nbin]),y=yy)
  }
  dx <- min(diff(ci$x));
  Id <- sum(ci$y*dx);
  ci$y <- ci$y/Id*integral;
  ci
}

normalize.density <- function(d,integral=1,integrate=T) {
  # scale y values so that density integrates to integral
  Id <- if(integrate) integrate.density(d) else 1;
  d$y <- d$y/Id*integral;
  d
}

integrate.density <- function(d) {
  integrate(approxfun(d,rule=2),lower=min(d$x),upper=max(d$x))$value
}




scale.density <- function(d,scale=1/mean(d$y)) {
  d$y <- d$y*scale
  d
}

whiten <-
  function(dens, lambda=0.9)
  {
    # add a small uniform value to the density to avoid zeroes when taking inverse
    dens$y <- lambda*dens$y + (1-lambda)/diff(range(dens$x))
    dens
  }

is.binned <- function(obj) {
  compact <- obj$call$compact
  if (is.null(compact))
    FALSE
  else
    eval(compact)
}


gf_importance <- function(obj) {
  Env_pc <- 0
  myimp <- NULL;
  variable <- NULL;
  for (i in names(obj$overall.imp) ) {
    imp <- extendedForest::importance(obj)[i]
    percent = round(imp / sum(extendedForest::importance(obj)) * 100, 2)
    #print(paste0(i,': ',percent,'%'))
    myimp <- rbind(myimp, percent)
    variable <- rbind(variable, i)
    Env_pc = Env_pc + imp
  }
  #print(paste0('overall predicators: ',Env_pc,'%'))
  myimp <- rbind(myimp, Env_pc)
  variable <- rbind(variable, 'overall predicators')
  rownames(myimp) <- variable
  colnames(myimp) <- "imp"
  return(myimp)
}

