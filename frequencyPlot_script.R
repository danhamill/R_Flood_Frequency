#' FrequencyPlot Generates a nice-looking (hydrologist-centric) frequency plot
#' @param series A vector representing an extreme value series (e.g., annual maximum flood)
#' @param ci A data frame containing confidence intervals (lower, true, upper) derived from 
#'  calling BootstrapCI()
#' @export ggplot
#' @import ggplot2
#' @import scales
frequencyPlot <- function(series, ci) {
  
  library(ggplot2)
  library(scales)
  
  # determine plotting positions
  bwpeaks <- data.frame(PROB = pp(series, sort = FALSE), FLOW = series)
  xbreaks <- c(0.002,0.01,0.10,0.25,0.5,0.8,0.9,0.95,0.975,0.99,0.995, 0.998)
  log.range <- log10(range(series, ci[,-1], na.rm = TRUE))
  lower <- 10^floor(log.range[1])
  upper <- 10^ceiling(log.range[2])
  cap <- lower
  ybreaks <- NULL
  while(cap < upper) {
    ybreaks <- c(ybreaks, seq(cap, cap*5, by = cap))
    cap <- cap * 10
  }
  
  # now plot
  ggplot(bwpeaks) + 
    geom_point(aes(x=PROB, y=FLOW)) + 
    theme_bw() + 
    scale_y_continuous(trans="log10", breaks=ybreaks, name="Discharge (cfs)", labels=comma) +
    scale_x_continuous(trans=probability_trans(distribution="norm"),
                       breaks=xbreaks, labels=signif(prob2T(xbreaks), digits=3),
                       name="Return period (yrs)") +
    geom_line(data=ci, aes(x=nonexceed_prob, y=true), color="blue", size=1.0) +
    geom_ribbon(data=ci, aes(x=nonexceed_prob, ymin=lower, ymax=upper), fill="blue",alpha=0.2)+
    theme(axis.text.x=element_text(size=10, vjust=0.5,colour="black", face="bold"), 
          axis.text.y=element_text(size=10, hjust=0.5, colour="black",face="bold"),
          axis.title.x=element_text(color="black", vjust=-0.35, face="bold",size=12),
          axis.title.y=element_text(color="black", vjust=0.35, face="bold", size=12))
    
  
}




