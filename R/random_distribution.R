
#------ random.distribution ----------

#' Create and plot a random distribution
#' 
#' This function computes a random truncated distribution which is a mixture of
#' one up to 4 different distributions that can be uniform, exponential,
#' normal or lognormal distributed. It can be used to test the evectiveness
#' of function \code{\link{is.unimodal}}.
#' @param lower lower bound of the support of the distribution. Defaults to 0.
#' @param upper upper bound of the support of the distribution. Defaults to 10.
#' @param plot boolean variable. Should the distribution be plotted?
#' @return A list with the distributions that are components of the mixture
#'   distribution. It has the same structure as parameter \code{distrib} in
#'   function \code{\link{dmix}}.
#' @examples
#' par(mfrow = c(1, 2))
#' d <- random.distribution(plot = TRUE)
#' curve(dmix(0, 10, distrib = d)(x), -1, 10, type = "l", col = "red")
#' @importFrom stats runif
#' @export
random.distribution <- function(lower = 0, upper = 10, plot = TRUE){
  x <- NULL # to avoid warning 'no visible binding...' in R CMD check
  n <- sample(1:4, 1)                         # random number of components
  specs <- c("norm", "unif", "exp", "lnorm")  # possible shape of components
  weights <- sample(1:3, n, replace=TRUE)     # random weights
  
  distributions <- as.list(NULL) # random values for distribution parameters
  for(i in 1:n){
    unifMin <- sample(lower:((lower+upper)/2),1)
    t <- switch(sample(specs,1),
                "norm"  = list(spec="norm", mean = sample(i:(2*i),1)),
                "unif"  = list(spec="unif", min = unifMin, max = unifMin+1),
                "exp"   = list(spec="exp", rate = runif(1, min = 0.1, max = 2)),
                "lnorm" = list(spec="lnorm", meanlog = sample(i:(2*i),1)))
    distributions[[i]] <- t
  }
  if(plot) graphics::curve(dmix(distrib=distributions, lower = lower, upper = upper)(x), 
                           from = lower-1, to = upper+1, ylab = "distribution")
  return(distributions)
}