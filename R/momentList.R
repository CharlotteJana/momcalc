
#' Class momentList
#' 
#' Create and validate objects of class \code{momentList}. Class \code{momentList}
#' is needed as argument for function \code{\link{transformMoment}}. The returned value
#' of \code{transformMoment} is of class \code{momentList}, too. The class provides a
#' convenient structure to store moments of a multidimensional distribution. \cr \cr
#' Class \code{momentList} consists of four elements: Element
#' \code{centralMoments} contains all known central moments of the distribution,
#' where as \code{rawMoments} contains all raw moments of the distribution. Both
#' are stored as list. The elements \code{centralMomentOrders} and
#' \code{rawMomentOrders} contain the corresponding orders of the moments. They
#' are stored as matrix or data.frame, each row represents one order of the
#' moment. The number of columns of these matrices should be equal to the
#' dimension of the distribution.
#' 
#' There are five different functions available for this class: \cr
#' Function \code{new_momentList} is a low level constructor that does not check for
#' correct types of the arguments or if they fit together. It should be only used
#' with care. \cr
#' Function \code{momentList} creates an object of class \code{momentList} and performs
#' various checks on the elements. Most of them are executed within function 
#' \code{validate_momentList}. This function particularly checks for moments of order
#' one and zero. If they contain false values, a warning will be thrown. If moments of these 
#' orders are not contained in \code{momentList}, they are added automatically. \cr
#' The functions \code{mean} and \code{cov} don't do calculations. They pick the
#' appropriate values out of momentList in case they exist. Function \code{cov} may
#' transform raw moments of order two into covariances, if the corresponding central
#' moments are not given.
#' 
#' @param centralMomentOrders matrix or data.frame. Every row gives the order of a central moment that is already known.
#' @param centralMoments list. The i-th entry is the central Moment of order \code{centralMomentOrders[i, ]}.
#' @param rawMomentOrders matrix or data.frame. Every row gives the order of a raw moment that is already known.
#' @param rawMoments list. The i-th entry is the raw Moment of order \code{rawMomentOrders[i, ]}.
#' @param x object of class \code{momentList}
#' @param ... additional arguments to function mean. They are currently not used.
#' @name momentList
#' @aliases validate_momentList new_momentList momentlist cov.momentList mean.momentList
NULL

####### constructors ########

#' @export
#' @rdname momentList
new_momentList <- function(rawMomentOrders = NULL,
                           rawMoments = list(),
                           centralMomentOrders = NULL,
                           centralMoments = list()){
  
  structure(list(rawMomentOrders = rawMomentOrders,
                 rawMoments = rawMoments,
                 centralMomentOrders = centralMomentOrders,
                 centralMoments = centralMoments,
                 class = "momentList"))
}


#' @export
#' @rdname momentList
momentList <- function(rawMomentOrders = NULL,
                       rawMoments = list(),
                       centralMomentOrders = NULL,
                       centralMoments = list()){
  
  if(is.null(rawMomentOrders) & is.null(centralMomentOrders))
    stop("Please provide either values for 'rawMomentOrders' and 'rawMoments'
          or values for 'centralMomentOrders' and 'centralMoments'")
  
  if(is.null(rawMomentOrders)){
    rawMomentOrders <- t(rep(0, ncol(centralMomentOrders)))
    rawMoments <- list(1)
  }
  if(is.null(centralMomentOrders)){
    centralMomentOrders <- rbind(rep(0, ncol(rawMomentOrders)), 
                                 diag(ncol(rawMomentOrders)))
    centralMoments <- append(1, as.list(rep(0, ncol(rawMomentOrders))))
  }
  
  mList <- new_momentList(rawMomentOrders = rawMomentOrders,
                          rawMoments = rawMoments,
                          centralMomentOrders = centralMomentOrders,
                          centralMoments = centralMoments)
  
  return(validate_momentList(mList))
}

####### validate_momentList ######

#' @rdname momentList 
#' @export
validate_momentList <- function(x){
  
  stopifnot(is.list(x$rawMoments))
  stopifnot(is.list(x$centralMoments))
  stopifnot(is.matrix(x$rawMomentOrders) | is.data.frame(x$rawMomentOrders))
  stopifnot(is.matrix(x$centralMomentOrders) | is.data.frame(x$centralMomentOrders))
  
  if(length(unique(x$rawMomentOrders)) != length(x$rawMomentOrders)){
    stop("Some rows in 'rawMomentOrders' appear several times. 
         They should only appear once.")
  }
  if(length(unique(x$centralMomentOrders)) != length(x$centralMomentOrders)){
    stop("Some rows in 'centralMomentOrders' appear several times. 
         They should only appear once.")
  }
  if(nrow(x$rawMomentOrders) != length(x$rawMoments)){
    stop("The number of elements in 'rawMoments' should be equal to
         the number of rows in 'rawMomentOrders'.")
  }
  if(nrow(x$centralMomentOrders) != length(x$centralMoments)){
    stop("The number of elements in 'centralMoments' should be equal to
         the number of rows in 'centralMomentOrders'.")
  }
  if(nrow(x$rawMomentOrders) != 0 &
     nrow(x$centralMomentOrders) != 0 & 
     ncol(x$rawMomentOrders) != ncol(x$centralMomentOrders)){
    stop("The number of columns in 'rawMomentOrders' and 'centralMomentOrders'
         should be identical.")
  }
  
  # ------- check moments of order 1 -------
  
  n <- ncol(x$rawMomentOrders)
  for(i in seq_len(n)){
    unitVector <- rep(0, n)
    unitVector[n-i+1] <- 1
    rowIndex <- prodlim::row.match(unitVector, x$centralMomentOrders)
    if(!is.na(rowIndex)){
      if(x$centralMoments[[rowIndex]] != 0) 
        warning("Central moments of order 1 should be 0.")
    }
    else{
      x$centralMomentOrders <- rbind(unitVector, x$centralMomentOrders)
      x$centralMoments <- append(0, x$centralMoments)
    }
  }
  rownames(x$centralMomentOrders) <- NULL
  
  #------- check moments of order 0 -------
  
  indexRaw <- prodlim::row.match(rep(0, n), x$rawMomentOrders)
  indexCentr <- prodlim::row.match(rep(0, n), x$centralMomentOrders)
  
  if(!is.na(indexRaw)){
    if(x$rawMoments[[indexRaw]] != 1)
      warning("Moments of order 0 should have value 1.")
  }
  if(!is.na(indexCentr)){
    if(x$centralMoments[[indexCentr]] != 1)
      warning("Moments of order 0 should have value 1.")
  }
  
  if(is.na(indexRaw)){
    x$rawMomentOrders <- rbind(rep(0, n), x$rawMomentOrders)
    x$rawMoments <- append(1, x$rawMoments)
  }
  if(is.na(indexCentr)){
    x$centralMomentOrders <- rbind(rep(0, n), x$centralMomentOrders)
    x$centralMoments <- append(1, x$centralMoments)
  }
  
  return(x)
}

####### mean #######


#' @rdname momentList
#' @export
mean.momentList <- function(x, ...){
  n <- ncol(x$rawMomentOrders)
  mean <- list()
  
  for(i in 1:n){
    row <- rep(0, n)
    row[i] <- 1
    mean[[i]] <- x$rawMoments[[prodlim::row.match(row, x$rawMomentOrders)]]
  }
  
  return(mean)
}

###### cov ######

cov <- function(x, y = NULL, use = "everything", 
                method = c("pearson", "kendall", "spearman"), ...){
  UseMethod("cov")
}

cov.default <- function(x, y = NULL, use = "everything", 
                        method = c("pearson", "kendall", "spearman"), ...){
   stats::cov(x, y = NULL, use = "everything", 
              method = c("pearson", "kendall", "spearman"), ...)
}

#' @rdname momentList
#' @export
cov.momentList <- function(x){
  
  n <- ifelse(length(x$centralMoments) > 0, 
              ncol(x$centralMomentOrders), 
              ncol(x$rawMomentOrders))

  covOrder <- expand.grid(lapply(1:n, function(i) 0:2))
  covOrder <- covOrder[which(rowSums(covOrder) == 2),]
  covOrderMissing <- is.na(prodlim::row.match(covOrder, x$centralMomentOrders))
  
  for(i in seq_along(covOrderMissing)){
    if(covOrderMissing[i]){
      x <- transformMoment(order = as.numeric(covOrder[i, ]),
                           type = "central",
                           momentList = x)  
    }
  }
  
  cov <- list()
  for(i in 1:n){
    cov[[i]] <- list()
    for(j in 1:i){
      
      row <- rep(0, n)
      row[i] <- row[i] + 1
      row[j] <- row[j] + 1
      
      cov[[i]][[j]] <- x$centralMoments[[prodlim::row.match(row, x$centralMomentOrders)]]
      cov[[j]][[i]] <- cov[[i]][[j]]
    }
  }
  
  return(cov)  
}