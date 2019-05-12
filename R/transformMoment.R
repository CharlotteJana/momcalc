##### transformMoment ####

#' Transform raw moments into central moments and vice versa
#'
#' Let X be a random variable and Y the corresponding centered variable, i.e. \eqn{Y = X - \mu}. 
#' Let p be the \code{order} of the desired moment. \cr \cr
#' If \bold{type = 'raw'}, this function returns
#' \deqn{
#' E(X^p) = \sum_{k_1=0}^{p_1}...\sum_{k_n=0}^{p_n} {p \choose k} \mu^{p-k} E(Y^k).}{
#' E(X^p) = \sum choose(p,k)*\mu^(p-k)*E(Y^k),
#' }
#' \ifelse{latex}{}{where the sum is taken over any row \code{k} in \code{expand.grid(lapply(1:n, function(i) 0:p[i]))}. \cr \cr}
#' The values of \eqn{E(Y^k)} are replaced by
#' \itemize{
#'   \item \code{centralMoment[[i]]}, if \code{k} = \code{centralMomentOrders[i, ]}
#'   \item \code{\link{transformMoment}(type = 'central', ...)}, if \code{k} is a row in \code{rawMomentOrders}
#'   \item \code{\link{symbolicMoments}(closure, k)}, if \code{k} is neither a row in \code{centralMomentOrders} nor in \code{rawMomentOrders}
#' }
#' If \bold{type = 'central'}, this function returns
#' \deqn{
#' E(Y^p) = \sum_{k_1=0}^{p_1}...\sum_{k_n=0}^{p_n} (-1)^{p-k} {p \choose k} \mu^{p-k} E(X^k).}{
#' E(Y^p) = \sum (-1)^(p-k)*choose(p,k)*\mu^(p-k)*E(X^k),
#' }
#' \ifelse{latex}{}{where the sum is taken over any row \code{k} in \code{expand.grid(lapply(1:n, function(i) 0:p[i]))}. \cr \cr} 
#' The values of \eqn{E(X^k)} are replaced by
#' \itemize{
#'   \item \code{rawMoment[[i]]}, if \code{k} = \code{rawMomentOrders[i, ]}
#'   \item \code{\link{transformMoment}(type = 'raw', ...)}, if \code{k} is not a row in \code{rawMomentOrders}
#' }
#' 
#' @return This function returns an object of class 'momentList'.
#' @aliases transformmoment
#' @param order numeric vector giving the order of the desired moment.
#' @param type string, either 'central' or 'raw'.
#' @param momentList object of class \code{\link{momentList}}.
#' @param closure string giving the closure method to use if a central moment is unknown. 
#' Possible values are the same as for argument \code{distribution} of \code{\link{symbolicMoments}}.
#' @param simplify bool indiciating if the resulting expressions should be simplified.
#' Function \code{\link[Deriv]{Simplify}} from package \pkg{Deriv} is used for simplification.
#' @importFrom Deriv Simplify
#' @export
transformMoment <- function(order, type, momentList, 
                            closure = "zero", simplify = TRUE){
  
  momentList <- validate_momentList(momentList, warnings = FALSE)
  
  p <- order
  n <- length(p)
  otherType <- setdiff(c("central", "raw"), type) # if type = 'raw' then otherType = 'central' and vice versa
  k_indexes <- as.matrix(expand.grid(lapply(1:n, function(i) 0:p[i])))
  
  # to avoid the R CMD Check NOTE 'no visible binding for global variable ...'
  typeOrders <- typeMoments <- otherOrders <- otherMoments <- NULL
  k_in_typeOrders <- k_in_otherOrders <- NULL

  # a help function to assign some values:
  readMoments <- function(mList, type){ 
    
    if(type == "raw"){
      assign("typeOrders", mList$rawMomentOrders, envir = parent.frame())
      assign("typeMoments", mList$rawMoments, envir = parent.frame())
      assign("otherOrders", mList$centralMomentOrders, envir = parent.frame())
      assign("otherMoments", mList$centralMoments, envir = parent.frame())
    }
    if(type == "central"){
      assign("typeOrders", mList$centralMomentOrders, envir = parent.frame())
      assign("typeMoments", mList$centralMoments, envir = parent.frame())
      assign("otherOrders", mList$rawMomentOrders, envir = parent.frame())
      assign("otherMoments", mList$rawMoments, envir = parent.frame())
    }
    
    k_in_typeOrders <- apply(k_indexes, 1, function(i){
      if(is.null(typeOrders)) NA
      else prodlim::row.match(i, typeOrders)
    })
    assign("k_in_typeOrders", k_in_typeOrders, envir = parent.frame())
    
    k_in_otherOrders <- apply(k_indexes, 1, function(i){
      if(is.null(otherOrders)) NA
      else prodlim::row.match(i, otherOrders)
    })
    assign("k_in_otherOrders", k_in_otherOrders, envir = parent.frame())
  }
  readMoments(momentList, type)
  
  # before calculation: check if calculation is needed
  if(!is.na(prodlim::row.match(p, typeOrders))) 
    return(momentList) # if there is already an entry in momentList
  if(type == "central" & is.na(prodlim::row.match(p, otherOrders))){
    moment <- symbolicMoments(distribution = closure, 
                              missingOrders = t(p),
                              mean = rep(0, n),
                              cov = extractCov(momentList))[[1]]
    
    momentList$centralMomentOrders <- rbind(momentList$centralMomentOrders, p)
    momentList$centralMoments <- append(momentList$centralMoments, moment)
    return(momentList)
  }
  
  # mu = vector with expected values
  muRowIndex <- lapply(1:n, function(i){
    muRow <- rep(0, n)
    muRow[i] <- 1
    prodlim::row.match(muRow, momentList$rawMomentOrders)
  })
  mu <- lapply(muRowIndex, function(i) momentList$rawMoments[[i]])
  if(sum(vapply(mu, is.null, logical(1))) > 0){
    stop("momentList$rawMoments should contain all expected values.")
  }
  
  #### calculation #####
    
  sum <- list()
  for(k_index in seq_len(nrow(k_indexes))){
    
    k <- k_indexes[k_index, ]
    
    # if k is neither a row in typeMomentOrders nor in otherMomentOrders
    if(is.na(k_in_typeOrders[k_index]) & is.na(k_in_otherOrders[k_index])){
      momentCentral <- symbolicMoments(distribution = closure, 
                                       missingOrders = t(k),
                                       mean = extractMean(momentList),
                                       cov = extractCov(momentList))[[1]]
      momentList$centralMomentOrders <- rbind(momentList$centralMomentOrders, k)
      momentList$centralMoments <- append(momentList$centralMoments, 
                                          momentCentral)
      readMoments(momentList, type)
    }
    # if k is a row in typeMomentOrders but not in otherMomentOrders
    if(!is.na(k_in_typeOrders[k_index]) & is.na(k_in_otherOrders[k_index])){
      momentList <- transformMoment(order = k, 
                                    type = otherType, 
                                    momentList = momentList, 
                                    closure = closure,
                                    simplify = simplify)
      readMoments(momentList, type)
    }
    # if k is a row in otherMomentOrders
    if(!is.na(k_in_otherOrders[k_index])){
      momentK <- otherMoments[[k_in_otherOrders[k_index]]]
    }
    
    # compute summand
    pChooseK <- as.numeric(Reduce("*", lapply(1:n, function(i) choose(p[i], k[i]))))
    muPower <- lapply(1:n, function(i) bquote(.(mu[[i]])^.(as.numeric(p[1:n]-k)[i])))
    muPower <- Reduce(function(a,b) bquote(.(a)*.(b)), muPower)
    if(type == "raw"){
      summand <- bquote(.(pChooseK)*.(muPower)*(.(momentK)))
    }
    if(type == "central"){
      sign <- (-1)^(sum(p)-sum(k))
      summand <- bquote(.(sign*pChooseK)*.(muPower)*(.(momentK)))
    }
    sum[[k_index]] <- summand
  }
  sum <- Reduce(function(a,b) bquote(.(a)+.(b)), sum)
  if(simplify){
    sum <- stringr::str_remove_all(pattern = "\"", string = format(sum))
    sum <- str2lang(Deriv::Simplify(sum))
  }
  
  if(type == 'raw'){
    momentList$rawMomentOrders <- rbind(typeOrders, p)
    momentList$rawMoments <- append(typeMoments, sum)
  }
  if(type == 'central'){
    momentList$centralMomentOrders <- rbind(typeOrders, p)
    momentList$centralMoments <- append(typeMoments, sum)
  }
  return(momentList)
}  
