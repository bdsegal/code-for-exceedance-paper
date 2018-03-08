# functions -------------------------------------------------------------------

exceedProb <- function(cutoff.vec, theta.hat, sd.hat, m) {
  # This function returns a point estimate for the exceedance probability

  p <- pnorm(q = sqrt(m) * (cutoff.vec - theta.hat) / sd.hat, lower.tail = FALSE)
  return(p)
}

deltaFun <- function(delta, k.stand, df, q) {
  # This function returns the difference between the upper tail
  # probability of a non-central t-distribution and a quantile q, 
  # where the t-distribution has n-1 degrees of freedom and 
  # non-centrality parameter delta.
  #
  # Arguments:
  #   delta (num): non-centrality parameter
  #   k.stand (num): quantile at which to evaluate the t distribution
  #   df (num): degrees of freedom
  #   q (num): confidence level (usually alpha/2 or 1-alpha/2)
  # Returns:
  #   num: t(k)_{n-1, delta} - q

  pt(q = k.stand, df = df, ncp = delta, lower.tail = FALSE) -  q
}

getDeltaCI <- function(cutoff, 
                       theta.hat, 
                       sd.hat, 
                       n, 
                       d,
                       alpha, 
                       interval = c(-100, 100)) {
  # This function obtains confidence intervals for the non-centrality
  # parameter of a t-distribution.
  #
  # Arguments:
  #   cutoff (num scalar or vector): cutoff value
  #   theta.hat (num): point estimate of mean parameter
  #   sd.hat (num): point estimate of standard deviation
  #   n (num): number of observations used to estimate theta.hat and sd.hat
  #   d (num): number of mean parameters
  #   alpha: confidence level of interval
  #   interval (num vector length 2): interval within which to search
  #            for roots; passed to the uniroot function
  # Returns:
  #   CI (vector or matrix): Ff scalar cutoff passed to getDeltaCI,
  #      returns a vector of length two giving lower and upper values
  #      If vector cutoff passed to getDeltaCI, returns a matrix with
  #      2 rows and number of columns equal to length of cutoff vector.

  k.stand <- sqrt(n) * (cutoff - theta.hat) / sd.hat

  delta.lower <- uniroot(f = deltaFun, 
                         interval = interval, 
                         k.stand = k.stand, 
                         df = n - d,
                         q = alpha / 2)$root

  delta.upper <- uniroot(f = deltaFun, 
                         interval = interval, 
                         k.stand = k.stand, 
                         df = n - d,
                         q = 1 - alpha / 2)$root

  out <- c(delta.lower, delta.upper)
  names(out) <- c("lower", "upper")
  return(out)
}
getDeltaCI <- Vectorize(getDeltaCI, vectorize.args = "cutoff")
