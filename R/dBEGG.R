#======== todo =================================================================
#t2 Note weg oder wirklich beispiele für is.unimodal mit dBEGG

#' Bimodal extension of the generalized Gamma-Distribution
#' 
#' The bimodal extension of the generalized Gamma-Distribution (BEGG) was first
#' introduced by Bulut et. al. in 2015. It is a scale mixture of the
#' generalized gamma distribution that is almost always bimodal.
#' The two modes can have different shapes, depending on the parameters 
#' \eqn{\alpha, \beta, \delta_0, \delta_1, \eta, \epsilon, \mu}{α, β, δ0, δ1, η, ε, μ} 
#' and \eqn{\sigma}{σ}.
#' @param x vector of quantiles.
#' @param alpha a positive number. Controls the kurtosis of the distribution.
#'  The distribution is leptokurtic for \eqn{\alpha \in (0, 2)}{α ϵ (0,2)} and
#'   \eqn{\beta = 1}{β = 1}. It is platikurtic for \eqn{\alpha > 2}{α > 2} and
#'   \eqn{\beta = 1}{β = 1}.
#' @param beta a positive number. Controls the kurtosis of the distribution.
#' @param delta0 a positive number. Controls the bimodality. If
#'   \eqn{\delta_0 = \delta_1}{δ0 = δ1}, the density function will have two modes
#'   with the same height. If \eqn{\delta_0 = \delta_1 = 1}{δ0 = δ1 = 0}, the
#'   distribution will be unimodal.
#' @param delta1 a positive number. Controls the bimodality.
#' @param eta a positive number. Controls the tail thickness of the
#'   distribution.
#' @param eps numeric. Controls the skewness of the distribution. When
#'   \eqn{\epsilon = 0}{ε = 0}, the distribution will be symmetric with two
#'   modes with different height.
#' @param mu numeric. The location parameter of the distribution. Defaults to 0.
#' @param sigma a positive number. The scaling parameter of the distribution.
#'   Defaults to 1.
#' @param order integer vector. Specifies all orders for which the raw moments
#' shall be calculated.
#' @note This distribution is included in package \pkg{momcalc} because it is a 
#' good test case for function \code{\link{is.unimodal}} and the raw moments are
#' known.
#' @aliases mBEGG dBEGG begg dbegg mbegg
#' @examples  # The first 3 examples are the same as in the paper:
#' par(mfrow=c(2, 2))
#' x <- seq(-2, 2, .01)
#' y <- dBEGG(x, alpha = 2, beta = 2, delta0 = 1, delta1 = 4, eta = 1, eps = 0)
#' plot(x, y, type = "l", ylab = "")
#' x <- seq(-2, 3, .01)
#' y <- dBEGG(x, alpha = 2, beta = 1, delta0 = 0, delta1 = 2, eta = 1, eps = -0.5)
#' plot(x, y, type = "l", ylab = "")
#' x <- seq(-2.5, 1.5, .01)
#' y <- dBEGG(x, alpha = 3, beta = 2, delta0 = 4, delta1 = 2, eta = 2, eps = 0.3)
#' plot(x, y, type = "l", ylab = "")
#' x <- seq(-2, 2, .01)
#' y <- dBEGG(x, alpha = 2, beta = 1, delta0 = 0, delta1 = 0, eta = 1, eps = 0.7)
#' plot(x, y, type = "l", ylab = "")
#'
#' @references Çankaya, M. N.; Bulut, Y. M.;  Doğru, F. Z. & Arslan, O.(2015). A
#' bimodal extension of the generalized gamma distribution. \emph{Revista
#' Colombiana de Estadística}, 38(2), 371-384.
#' @name BEGG
#' @aliases dBEGG mBEGG
NULL

#' @describeIn BEGG density function
#' @export
dBEGG <- function(x, alpha, beta, delta0, delta1, eta, eps, mu = 0, sigma = 1){
  a0 <- (delta0+1)/alpha
  a1 <- (delta1+1)/alpha
  f <- function(x){
    if (x >= mu) {
      return ((alpha * beta) / 
        (2 * sigma * eta ^ a0 * (1 - eps) ^ delta0 * gamma(a0 / beta)) * 
        ((x - mu) / sigma) ^ delta0 * exp(-((x - mu) ^ (alpha * beta)) / 
        (eta ^ beta * ((1 - eps) * sigma) ^ (alpha * beta))))
    } else {
      return ((alpha * beta) / 
        (2 * sigma * eta ^ a1 * (1 + eps) ^ delta1 * gamma(a1 / beta)) * 
        ((mu - x) / sigma) ^ delta1 * exp(-((mu - x) ^ (alpha * beta)) / 
        (eta ^ beta * ((1 + eps) * sigma) ^ (alpha * beta))))
    } 
  }
  vapply(x, f, numeric(1))
}

#' @describeIn BEGG raw moments
#' @export
mBEGG <- function(alpha, beta, delta0, delta1, eta, 
                  eps, mu = 0, sigma = 1, order = 1:4){
  
  ((-1) ^ order * eta ^ (order / alpha) * (1 + eps) ^ (order + 1))/2 * 
  gamma((delta1 + order + 1) / (alpha * beta)) / 
  gamma((delta1 + 1) / (alpha * beta)) +  
  (eta ^ (order / alpha) * (1 - eps) ^ (order + 1))/2 * 
  gamma((delta0 + order + 1) / (alpha * beta)) / 
  gamma((delta0 + 1) / (alpha * beta))  
  
}

