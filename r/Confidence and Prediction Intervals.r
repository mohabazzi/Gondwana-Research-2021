#########################################
### Confidence & Prediction Intervals ###
#########################################

##  Author: Nicolas E. Campione

##  DESCRIPTION
#   This function calculates the confidence and prediction intervals of the mean, based on a numeric vector.

##  ARGUMENTS
#   xvar  - a numeric vector representing the sample for which intervals will be calculated.
#   alpha - a numeric value specifying the significance level for the confidence intervals. 95% confidence is the default
#           (i.e., alpha=0.975)

##  VALUES
#   Returns a list of length=8, including:
#   $mean     - mean of 'xvar'
#   $st.dev   - standard deviation of 'xvar'
#   $N        - sample size
#   $t.stat   - t-statistic, give 'alpha' and n-1 degrees of freedom
#   $Lower.CI - lower confidence interval
#   $Upper.CI - uuper confidence interval
#   $Lower.PI - lower prediction interval
#   $Upper.PI - upper prediction interval

intervals <- function(xvar,alpha = 0.975) {
  mean <- mean(xvar,na.rm = TRUE)
  std <- sd(xvar,na.rm = TRUE)
  n <- length(na.omit(xvar))
  ci.s.x <- sqrt((std^2)/n)
  pi.s.x <- std*sqrt(1 + 1/n)
  t <- qt(alpha,df = n - 1)
  ci.uci <- mean + t*ci.s.x
  ci.lci <- mean - t*ci.s.x
  pi.uci <- mean + t*pi.s.x
  pi.lci <- mean - t*pi.s.x
  res <- as.list(c(round(mean,4),round(std,4),n,round(t,4),round(ci.lci,4),round(ci.uci,4),round(pi.lci,4),round(pi.uci,4)))
  names(res) <- c('mean','st.dev','N','t.stat','Lower.CI','Upper.CI','Lower.PI','Upper.PI')
  return(res)
}