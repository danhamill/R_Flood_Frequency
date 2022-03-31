##From blog http://www.headwateranalytics.com/blog/flood-frequency-analysis-in-r 
#' FrequencyAnalysis Fits a given extreme value distribution to an extreme value series
#' @param series A vector representing an extreme value series (e.g., annual maximum flood)
#' @param distribution A three-character name of the extreme value distribution (see ?dist.list())
#' @param nep A vector of non-exceedance probabilities
#' @return A list object containing: (1) distribution information and (2) output
#'  (quantile estimates at various non-exceedance probabilities)
#' @export
#' @import lmomco
FrequencyAnalysis <- function( series, distribution, nep = nonexceeds() ) {
  
  library(lmomco)
  
  distribution <- tolower(distribution)
  transformed <- FALSE
  
  # add log Pearson Type 3 to list of distributions supported
  # by lmomco package
  base.dist <- c('lp3', dist.list())
  
  if( any(distribution %in% base.dist) ) {
    
    # log transform series 
    if( distribution == 'lp3' ) {
      series <- log10(series)
      transformed <- TRUE
      distribution <- 'pe3'
    }
    
    # compute L-moments
    samLmom <- lmom.ub(series)
    
    # estimate distribution parameters
    distPar <- lmom2par(samLmom, type = distribution)
    
    # compute quantiles for nonexceedances
    quant <- par2qua(f = nep, para = distPar)
    
    if( distribution == 'pe3' & transformed ) {
      distribution <- 'lp3'
      quant <- 10^quant
    }
    
    # return result as list object
    return(
      list(
        distribution = list(
          name = distribution,
          logTransformed = transformed,
          parameters = distPar),
        output = data.frame(nep = nep, rp = prob2T(nep), estimate = quant) 
      ) )
    
  } else {
    stop(
      sprintf('Distribution \'%s\' not recognized!', distribution))
  }
}
