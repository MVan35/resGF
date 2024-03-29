% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/resGF.R
\name{get_importance}
\alias{get_importance}
\title{Return gradient forest importance}
\usage{
get_importance(obj, raster_stack)
}
\arguments{
\item{obj}{A GF object}

\item{raster_stack}{A raster stack oject containing the same variable used in the GF analysis}
}
\value{
importance of the different variables
}
\description{
This function generate a resistance surface based on gradient forest analysis.
}
