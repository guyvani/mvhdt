glmnetX=function(formula, family=NULL, varpower=NULL, data, na.action, weights=NULL, offset=NULL,alpha=1.0,nlambda=100,
                lambda.min.ratio=ifelse(nobs<nvars,1e-2,1e-4),lambda=NULL,standardize=TRUE,intercept=TRUE,thresh=1e-7,dfmax=nvars+1,
                pmax=min(dfmax*2+20,nvars),exclude,penalty.factor=rep(1,nvars),lower.limits=-Inf,upper.limits=Inf,maxit=100000,
                type.gaussian=ifelse(nvars<500,"covariance","naive"),type.logistic=c("Newton","modified.Newton"),standardize.response=FALSE,
                type.multinomial=c("ungrouped","grouped"),  method = "glmnet.fit",  x = FALSE, y = TRUE, 
                contrasts = NULL, maxit=100, epsilon = 0.00001,  CV=FALSE,  parallel=FALSE 
                ... ){
  # for monotasklearning , formula:  Y ~ X1 + X2 + X3 + X4 ...,
  # for mutlitasklearning , formula:  cbind(Y1,Y2, Y3, ...) ~ X1 + X2 + X3 + X4 ...,
  # for mutlitasklearning , character vector or matrix, same for offset,
  
  
  call <- match.call()
   if (is.null(family)) {
       stop(" family must be specified")
   } else  {
      if(!(all(family %in% c("gaussian","binomial","poisson","multinomial","cox","mgaussian", "gamma", "tweedie"))))
       stop("'family' not recognized")
   }
  
  if (missing(data)) 
    data <- environment(formula)
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "weights", "na.action", "offset"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- as.name("model.frame")
   mf <- eval(mf, parent.frame())
   if (identical(method, "model.frame")) 
    return(mf)
  if (!is.character(method) && !is.function(method)) 
  if (method != "glmnet.fit")   
    stop("invalid 'method' argument")
  #if (identical(method, "glmnet.fit")) 
  #  control <- do.call("glm.control", control)
  mt <- attr(mf, "terms")
  Y <- model.response(mf, "any")
  
  #weights and offsets
  weights <- model.weights(mf)
  if (!is.null(weights) && !is.numeric(weights)) 
    stop("'weights' must numeric values")
  if (!is.null(weights) && any(weights < 0)) 
    stop("negative weights not allowed")
  offset <- model.offset(mf))

 #matrix of predictors
  X <- if (!is.empty.model(mt)) 
  model.matrix(mt, mf, contrasts)
  else matrix(, NROW(Y), 0L)
  
  
#monotasklearning  
  if (length(dim(Y)) == 1L) {
    nm <- rownames(Y)
    dim(Y) <- NULL
    if (!is.null(nm)) 
      names(Y) <- nm
  

  
  weights <- as.vector(model.weights(mf))
  offset <- as.vector(model.offset(mf))
  
  if (!is.null(weights)) {
    if (length(weights) != NROW(Y)) 
      stop(gettextf("number of weights is %d should equal %d (number of observations)", 
                    length(weights), NROW(Y)), domain = NA)
  }
  
  if (!is.null(offset)) {
    if (length(offset) != NROW(Y)) 
      stop(gettextf("number of offsets is %d should equal %d (number of observations)", 
                    length(offset), NROW(Y)), domain = NA)
  }
  
  fit<- glmnet.fit(x=X,y=Y,family=family, weights=weights,offset=offset,alpha=alpha,nlambda=nlambda,
                        lambda.min.ratio=lambda.min.ratio,lambda=lambda,standardize=standardize,intercept=intercept,
                        thresh=thresh,dfmax=dfmax, pmax=pmax,exclude=exclude,penalty.factor=penalty.factor,lower.limits=lower.limits,
                        upper.limits=upper.limits, maxit=maxit, type.gaussian=type.gaussian, type.logistic=type.logistic,
                        standardize.response=standardize.response,type.multinomial=type.multinomial)
  
  #multitasklearning 
  } else if (length(dim(Y)) > 1L)  {
    
    if ( dim(Y)[2] ! = length(family)) stop("Specify a family for each response")
    
    if (!is.null(offset)) {
      if (is.character(offset)) {
        if ( dim(Y)[2] ! = length(offset)) stop("length of offset should be equal to number of responses")
      } else {
        if ( dim(Y)[2] ! = dim(offset)[2]) stop("Specify an offset for each response")
        if ( dim(Y)[1] ! = dim(offset)[1]) stop("length of offsets should be equal to the number of observations")
      }
    }
    
    if (!is.null(weights)) {
      if (is.character(weights)) {
        if ( dim(Y)[2] ! = length(weights)) stop("length of weights should be equal to number of responses")
        
      }  else {
        if ( dim(Y)[2] ! = dim(weights)[2]) stop("Specify weights for each response")
        if ( dim(Y)[1] ! = dim(weights)[1]) stop("length of weights should be equal to the number of observations")
      }
    }
    
    #get tweedie variance power
    tweedie.response = family[which(family %in% "tweedie")]
    if (length(tweedie.response) >= 1) {
       if (is.null(varpower)) { 
           stop("must specify variance power for tweedie response")
       }  else {
          if ( length(family) != length(varpower) ) stop("specify a variance power for each tweedie response and NULL for non-tweedie")
          if (any(varpower < 1) | any(varpower > 3)) stop("invalid tweedie variance power")
       }
    }
    
    #c("gaussian","binomial","poisson","multinomial","cox","mgaussian", "gamma", "tweedie")
    family.obj<-list()
    z0<-as.matrix(Y)
    w0<-as.matrix(Y)
    good0<-as.matrix(Y)
    #dev0<-rep(NA, dim(Y)[2])    

    for (i in 1:length(family)) {
      
      if (family[i] %in%  c("gaussian","binomial","poisson", "gamma")) {
           family.obj[i]=switch(family[i],
                  "gaussian"=gaussian(link = "identity") ,
                  "poisson"=poisson(link = "log") ,
                  "binomial"=binomial(link = "logit") ,
                  "gamma"=Gamma(link = "log") 
                )
      } else if (family[i] %in%  c("tweedie") ) {
            require(statmod) # tweedie family
            if(is.null(varpower[i])) {
             stop("variance power must be specified for tweedie response")
            } else {
              if (varpower[i] < 1 | varpower[i] > 3) stop("invalid tweedie variance power")
            }
          family.obj[i]=tweedie(var.power=varpower[i],link.power=0) #tweedie with log link function
        
      } else {
        stop("family not allowed in the multitasklearning framework")
      }
       
     
      linearize <-(y=Y[,i], weights=weights[,i], offset=offset[,i], family=family.obj[i], mustart=NULL) 
      z0[,i]<-linearize$z
      w0[,i]<-linearize$w
      good0[,i]<-linearize$good
      #dev0[i]<-linearize$dev  
      
    }
    
    iter <- 1
    dev0 <- 0.00001
    epsilon <- 0.01
    diff<-rep(1000000000000, length(Y))
    
    while(iter <= maxit) {
      
      z=z0
      w=w0
      good<-good0
      #dev<-dev0
      
      #X=X[good, , drop = FALSE]*w 
      #X=X[good, ]*w 
      #Z=z[good,]*w
      
      fit<-  glmnet.fit(x=X[good, ]*w ,y=z[good,]*w,family="mgaussian", weights=NULL,offset=NULL,alpha=alpha,nlambda=nlambda,
                        lambda.min.ratio=lambda.min.ratio,lambda=lambda,standardize=standardize,intercept=intercept,
                        thresh=thresh,dfmax=dfmax, pmax=pmax,exclude=exclude,penalty.factor=penalty.factor,lower.limits=lower.limits,
                        upper.limits=upper.limits, maxit=maxit, type.gaussian=type.gaussian, type.logistic=type.logistic,
                        standardize.response=standardize.response,type.multinomial=type.multinomial)
      
      idx<-which.max(fit$dev.ratio)
      dev<-fit$dev.ratio[idx]
      
      if ( (abs(dev - dev0)/(0.0000001 + abs(dev0)) <= epsilon) & (iter > 1) ) break 
       dev0 <- dev
      
      lbd<-fit$lambda[idx]
      pred<-predict(fit, newx = X[good, ]*w , s = lbd, type="response")
      nrcds <- nrow(z[good,])
      nrsp <- ncol(z[good,])
      Zhat<-matrix( rep(NA), nrcds, nrsp )
       for (j in 1:nrsp){
         Zhat[,j] <- family.obj$linkinv(pred[1:nrcds, j, ]) 
       }
      
      for (i in 1:length(family)) {
        beta<- fit$coef[,i] # call properly
        linearize <-(y=Y[,i], weights=weights[,i], offset=offset[,i], family=family.obj[i], mustart=Zhat[,j]) 
        z0[,i]<-linearize$z
        w0[,i]<-linearize$w
        good0[,i]<-linearize$good
        }
      
    }

  }
  
###Output
  if(is.null(lambda))fit$lambda=fix.lam(fit$lambda)##first lambda is infinity; changed to entry point
  fit$call <- call 
  if (length(dim(Y)) > 1L) {     #get family information for predictions in mutlitask learning
    fit$family<- family.obj
    }
  fit$nobs <- dim(Y)[1]
  class(fit)=c(class(fit),"glmnet")
  
  return(fit) 
}
