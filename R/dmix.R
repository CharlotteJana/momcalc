#' Mixtures of truncated distributions
#' 
#' Compute the distribution and raw moments of a mixture of different, 
#' one dimenisonal, eventually truncated distributions.
#' 
#' @param lower numeric giving the lower bound of the support. Defaults to -Inf.
#' @param upper numeric giving the upper bound of the support. Defaults to Inf.
#' @param distrib a list. Every element is itself a list representing a
#'   distribution. This list should have one named element \code{spec} which
#'   describes the distribution, i. e. "exp" for the exponential distribution or
#'   "norm" for the normal distribution. The other elements are optional
#'   additional parameters for the specfied distribution.
#' @param weights numeric vector with the same length as \code{distrib}.
#' Provides weights for every distribution given in \code{distrib}.
#' @examples 
#' distributions <- list(list(spec="exp", rate = 2),
#'                      list(spec="norm", mean = 0, sd = 0.5),
#'                      list(spec="unif", min = 2, max = 3))
#' d <- dmix(-1, 3, weights = c(.1,.3,.1), distrib = distributions)
#' curve(d(x), -2, 5)
#' @name dmix
#' @aliases mmix
#' @export
dmix <- function(lower = -Inf, upper = Inf, distrib, weights){
  
  n <- length(distrib)
  if(missing(weights)) weights <- rep(1/n, n)
  
  weights <- weights/sum(weights)
  function(x){
    h <- vapply(1:n, 
                function(i){
                  weights[i]*do.call("dtrunc", 
                                     c(x = list(x), distrib[[i]], 
                                     a = lower, b = upper))
                },
                numeric(length(x))
    )
    rowSums(h)
  }
}

#' @param order order of the moment
#' @rdname dmix
#' @export
mmix <- function(order = 1:4, lower = -Inf, upper = Inf, distrib, weights){
  
  n <- length(distrib)
  if(missing(weights)) weights <- rep(1/n, n)
  
  weights <- weights/sum(weights)
  h <- vapply(1:n, 
              function(i){ 
                weights[i]*do.call("mtrunc", c(list(order = order), 
                                               a = lower, 
                                               b = upper, 
                                               distrib[[i]]))
              },
              numeric(length(order))
  )
  if(length(order) > 1)
    return(rowSums(h))
  else
    return(sum(h))
}
