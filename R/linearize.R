linearize <- function (y, weights = NULL, offset = NULL, family = myfamily, mustart=NULL) 
  {
    nobs <- NROW(y)
    if (is.null(weights)) 
      weights <- rep.int(1, nobs)
    if (is.null(offset)) 
      offset <- rep.int(0, nobs)
    variance <- family$variance
    linkinv <- family$linkinv
    if (!is.function(variance) || !is.function(linkinv)) 
      stop("'family' argument seems not to be a valid family object", 
           call. = FALSE)
    dev.resids <- family$dev.resids
    aic <- family$aic
    mu.eta <- family$mu.eta
    
    unless.null <- function(x, if.null) if (is.null(x)) if.null else x
    valideta <- unless.null(family$valideta, function(eta) TRUE)
    validmu <- unless.null(family$validmu, function(mu) TRUE)
    
     if (is.null(mustart)) {
        eval(family$initialize)
     }
     else {
         mukeep <- mustart
         eval(family$initialize)
         mustart <- mukeep
     }
    
    eta <-family$linkfun(mustart)
    mu <- linkinv(eta)
    if (!(validmu(mu) && valideta(eta))) 
      stop("cannot find valid starting values: please specify some", 
           call. = FALSE)
    dev<- sum(dev.resids(y, mu, weights))

    good <- weights > 0
    varmu <- variance(mu)[good]
    if (any(is.na(varmu))) 
      stop("NAs in V(mu)")
    if (any(varmu == 0)) 
      stop("0s in V(mu)")
    mu.eta.val <- mu.eta(eta)
    if (any(is.na(mu.eta.val[good]))) 
      stop("NAs in d(mu)/d(eta)")
    good <- (weights > 0) & (mu.eta.val != 0)
    if (all(!good)) {
      conv <- FALSE
      stop("no informative  observations found")
    }
    z <- (eta - offset)[good] + (y - mu)[good]/mu.eta.val[good]
    w <- sqrt((weights[good] * mu.eta.val[good]^2)/variance(mu)[good])
    #aic.model <- aic(y[good], length(y[good]), mu, weights, dev) + 2 * rank

  return(list(z=z, w=w , good=good, dev=dev))     
    
  }   




