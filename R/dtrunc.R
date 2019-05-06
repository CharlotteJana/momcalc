#===============================================================================
#v1 Zur Dokumentation von dtrunc: Zitat aus http://r-pkgs.had.co.nz/check.html
# "If the licenses are compatible you can copy and paste the exported function 
#  into your own package. If you do this, remember to update Authors@R."
#t1 Die Argumente upper, a, lower, b einheitlich verwenden!
#t1 Brauche ich dtrunc wirklich? (mtrunc gibt es in truncdist nicht!)
#t1 vielleicht fehler testen und dann durch Nullen ersetzen?

#------ dtrunc & mtrunc ----------

#' Probability density function of truncated random variables
#' 
#' This function computes values for the probability density function of a
#' truncated random variable. It was originally implemented in package
#' \pkg{truncdist} and slightly modified to return zeros in case that the
#' trunction interval [a, b] is not inside the support of the density function.
#' @inherit truncdist::dtrunc
#' @param ... other arguments are passed to the corresponding quantile function
#' @examples 
#' x <- seq(0, 3, 0.1)
#' dtrunc(x, spec = "norm", a = 1, b = 2)
#' curve(dtrunc(x, spec = "norm", a = -Inf, b = 1), -10, 2)
#' 
#' \dontrun{
#' # different results for intervals outside the support of the density function:#' 
#' truncdist::dtrunc(x, spec = "norm", a = 20, b = 30) # gives error
#' momcalc::dtrunc(x, spec = "norm", a = 20, b = 30) # gives only a warning}
#' @importFrom utils str
#' @export
dtrunc <- function (x, spec, a = -Inf, b = Inf, ...) {
  if (a >= b) 
    stop("argument a is greater than or equal to b")
  tt <- rep(0, length(x))
  g <- get(paste("d", spec, sep = ""), mode = "function")
  G <- get(paste("p", spec, sep = ""), mode = "function")
  if(G(b, ...)-G(a, ...) != 0){ # additional to truncdist::dtrunc
    tt[x >= a & x <= b] <- g(x[x >= a & x <= b], ...)/(G(b, ...) - G(a, ...))
  }
  else  
   warning("Truncation interval is outside the support of the density function")
  
  return(tt)
}

#' Moments of truncated random variables
#' 
#' This function computes the raw moments of a truncated random variable.
#' @inheritParams truncdist::dtrunc
#' @param order numeric vector giving the order of the moments
#' @param ... other arguments are passed to the corresponding moment function
#' @examples
#' mtrunc(1:6, spec = "norm")
#' mtrunc(1:6, spec = "norm", b = 0)
#' @seealso \code{\link{dtrunc}} for the probability distribution of a truncated variable.
#' @importFrom actuar mnorm mlnorm mexp munif
#' @importFrom stats integrate
#' @export
mtrunc <- function (order, spec, a = -Inf, b = Inf , ...){
  argList <- list(...)
  
  if(spec == "norm" & a == -Inf & b == Inf) 
    return(actuar::mnorm(order, ...))
  else if(spec == "lnorm" & a <= 0 & b == Inf)   
    return(actuar::mlnorm(order, ...))
  else if(spec == "exp" & a <= 0 & b == Inf)     
    return(actuar::mexp(order, ...))
  else if(spec == "unif") {
    if(a <= argList$min & b >= argList$max) 
      return(actuar::munif(order, ...))
    if(a > argList$max | b < argList$min) 
      stop("this is not a distribution")
  }
  else{
  #message("Direct integration is used.")
  sapply( order, 
    function(i) integrate( 
      Vectorize(function(x) x^i * dtrunc(x, spec = spec, a = a, b = b, ...)), 
      lower=a, upper=b)$value)
  }
}