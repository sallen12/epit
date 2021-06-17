#' Determine optimal number of weights by CN criterion
#'
#' @param tData input data (transformed to fit into the given interval).
#' @param L diagonal weight matrix.
#'
#' @return
#' The optimal number of weights (single integer).
#'
#' @author
#' Code by Bradley C. Turnbull, NC State University, Dr. Sujit K. Ghosh, NC
#' State University, August 25, 2013.
#'
#' @keywords internal
mOpt.CN <- function(tData, L){

  #set starting point for m
  n <- length(tData)
  m <- floor(  n^(1/3) ) - 1

  #set start value for while loop
  logratio <- 1

  while( logratio < sqrt(n) ){
    m <- m+1

    #construct D matrix using m value
    B <- NULL
    for(k in 1:m){
      B <- cbind(B, stats::pbeta(tData, shape1=k, shape2=m-k+1))
    }
    Dmat <- t(B) %*% L %*% B

    #take spectral decomposition to find eigenvalues
    spec <- eigen(Dmat, symmetric=TRUE)
    d <- spec$values
    min.eigenValue <- max( min(d), 0 )
    max.eigenValue <- max(d)

    logratio <- log10(max.eigenValue) - log10(min.eigenValue)
  }
  m-1 #return number of weights
}

#' Generate constraint matrix
#'
#' Generates constraint matrix for estimating monotone density
#'
#' @param m maximal number of degrees for Bernstein polynomial.
#' @param maxplace location of the mode.
#'
#' @return
#' The constraint matrix.
#'
#' @author
#' Code by Bradley C. Turnbull, NC State University, Dr. Sujit K. Ghosh, NC
#' State University, August 25, 2013.
#'
#' @keywords internal
constraintMat <- function( m, maxplace){
  A <- suppressWarnings(
    rbind( rep(1,m), diag(rep(1,m)),
           matrix( rep( c(-1,1, rep(0,m - 1)) , maxplace-1), maxplace-1, m, byrow=TRUE),
           matrix( rep( c( rep (0,maxplace-1),1,-1,rep(0, m-maxplace)),m-maxplace),
                   m-maxplace,m,byrow=TRUE))
  )
  Amat <- t(A)
}

#' Compute weight vector
#'
#' @param m maximum degree in Bernstein polynomial.
#' @param Fn ecdf of transformed data.
#' @param lower,upper lower and upper bound of support.
#' @param Dmat,dvec matrix and vector for QP problem.
#' @param monotone direction of monotonicity.
#' @param settings settings for \code{\link[osqp]{solve_osqp}}.
#'
#' @author
#' Original code by Bradley C. Turnbull, NC State University, Dr. Sujit K.
#' Ghosh, NC State University, August 25, 2013, suitably adapted by Sebastian
#' Arnold and Alexander Henzi.
#'
#' @return
#'
#'
#' @importFrom osqp solve_osqp
#'
#' @keywords internal
solveWeights <- function(m, Fn, lower, upper, Dmat, dvec,monotone, settings){
  #find the location of the maximum weight
  if(monotone=="Iso"){max.place <- m}
  if(monotone=="Decr"){max.place <- 1}
  #make the constraint matrix
  Amat <- constraintMat(m, max.place)
  #make bvec vector of constraints
  bvec=c(1,rep(0,2*m-1))
  #find weights using solve.QP function
  w.hat = osqp::solve_osqp(
    P=Dmat,
    q=-dvec,
    A=t(Amat), l= bvec, u=c(bvec[1],rep(Inf,length(bvec)-1)),pars = settings)$x
  #function to find max of an element and 0
  max0 <- function(x){max(x,0)}
  #make sure no weights are < 0
  w.hat <- sapply( w.hat, max0)
  #make sure the weights sum to 1
  wsum <- sum(w.hat)
  w.hat / wsum
}

#' Monotone density estimation with Bernstein polynomials
#'
#' Estimates a smooth monotone mixture of Beta densities.
#'
#' @param data observations (in [0,1]).
#' @param monotone chose \code{"Iso"} for monotone increasing or \code{"Decr"}
#'     for monotone decreasing.
#' @param m maximum degree of the Bernstein polynomials in the mixture. If
#'     \code{NULL}, an optimal \code{m} is estimated with the AIC, BIC or CN
#'     criterion specified.
#' @param crit the type of criterion to use for selecting the number of weights.
#'     Defaults to \code{"CN"}.
#' @param settings options for \code{\link[osqp]{osqp}}.
#'
#' @author
#' Original code by Bradley C. Turnbull, NC State University, Dr. Sujit K.
#' Ghosh, NC State University, August 25, 2013, suitably adapted by Sebastian
#' Arnold and Alexander Henzi.
#'
#' @return
#' A list containing the estimated weights, \code{m}, the density, distribution,
#' random number generation and quantile function for the estimate density,
#' and the bounds of the support.
#'
#' @importFrom stats dbeta
#' @importFrom stats pbeta
#' @importFrom stats rbeta
#' @importFrom stats runif
#' @importFrom stats sd
#' @importFrom stats uniroot
#' @importFrom stats ecdf
#'
#' @export
umd <- function(data,monotone="Iso", crit="CN",
                m=NULL, settings = list()){

  #delta definition
  n <- length(data)
  delta <- stats::sd(data)/ sqrt(n)
  sort.data <- sort(data)

  lower <- 0
  upper <- 1

  #Make the transformed data set
  tdata <- (data-lower)/(upper-lower)	# redudant line since we work with [0,1]-valued data

  #Construct the L matrix and make vector of ecdf values
  Fn <- stats::ecdf(tdata)
  ep <- 3/( 8*length(data) )
  ecdf.vec <- Fn(tdata)

  #construct L
  L <- diag( n / ( ( ecdf.vec + ep)*(1 + ep - ecdf.vec ) ) )


  #Find the optimal number of weights, depends on the given crit
  #then return the desired functions
  if( crit == "CN"){

    if( !is.null(m)){

      #make B, Dmat, and dvec
      B <- NULL
      for(k in 1:m){
        B <- cbind(B, stats::pbeta(tdata, shape1=k, shape2=m-k+1))
      }
      Dmat <- t(B) %*% L %*% B
      dvec <- t(B)%*% L %*% ecdf.vec

      #Make sure there are no zero eigen values
      spec <- eigen(Dmat, symmetric=TRUE)
      d <- spec$values
      Q <- spec$vectors

      #find which values are < 10e-6, and set to smallest
      #eigen value >= 10e-6
      if( min(d) < 10e-6){
        tooSmall <- which( d < 10e-6)
        d[tooSmall] <- d[ min(tooSmall) - 1]
        #Recreate pos. def. Dmat matrix
        Dmat <- Q %*% diag(d) %*% t(Q)
      }

    } #if no m provided then use the optimal m
    else{ m <- mOpt.CN(tdata, L)
    #make B, Dmat, and dvec
    B <- NULL
    for(k in 1:m){
      B <- cbind(B, stats::pbeta(tdata, shape1=k, shape2=m-k+1))
    }
    Dmat <- t(B) %*% L %*% B
    dvec <- t(B)%*% L %*% ecdf.vec

    #Make sure there are no zero eigen values
    spec <- eigen(Dmat, symmetric=TRUE)
    d <- spec$values
    Q <- spec$vectors

    #find which values are < 10e-6, and set to smallest
    #eigen value >= 10e-6
    if( min(d) < 10e-6){
      tooSmall <- which( d < 10e-6)
      d[tooSmall] <- d[ min(tooSmall) - 1]
      #Recreate pos. def. Dmat matrix
      Dmat <- Q %*% diag(d) %*% t(Q)
    }
    }

    #Solve for the weights
    weights <- solveWeights(m,Fn,lower,upper,Dmat,dvec,monotone, settings)

    #use the weights to create the 4 distribution functions
    #pdf
    dumd = function(x){
      mix.pdf = function(x){ifelse( (x-lower)*(upper-x) >= 0,
                                    sum(weights*stats::dbeta((x-lower)/(upper-lower),1:m,m:1))/(upper-lower),0)}
      sapply(x, mix.pdf)
    }
    #cdf
    pumd = function(x){
      mix.cdf= function(j){sum(weights*stats::pbeta( (j-lower)/(upper-lower), 1:m, m:1)) }
      sapply(x, mix.cdf)
    }
    #rand generator
    rumd = function(n=1){
      rsample = function(){
        k <- sample(1:m,size=1,prob=weights)
        return(  lower+( stats::rbeta(1, k, m-k+1)*(upper-lower) )   )
      }
      replicate(n, rsample() )
    }
    #quantile function
    qumd = function(q){
      mix.quantile = function(q){ g = function(x){ pumd(x)-q }
      return( uniroot(g, interval=c(lower,upper) )$root )
      }
      sapply(q, mix.quantile)
    }

    #Return list (ends function)
    return(list(weights=weights,m.hat=m,dumd=dumd,pumd=pumd,rumd=rumd,qumd=qumd,
                lower=lower, upper=upper))
  }
  else if(crit %in% c("AIC","BIC")) {

    #Find Optimal number of weights using AIC and BIC criterion
    m.test <- seq( 2, ceiling(n/log(n)) )

    #lists and vectors to hold the values in the loop
    AICs <- rep(0, length(m.test) )
    BICs <- rep(0, length(m.test) )
    weights.list <- list()

    #loop through m values and calculate the AIC
    for( i in 1:length(m.test)){

      #make B, Dmat, and dvec for this m
      B <- NULL
      for(k in 1:m.test[i]){
        B <- cbind(B, stats::pbeta(tdata, shape1=k, shape2=m.test[i]-k+1))
      }
      Dmat <- t(B) %*% L %*% B
      dvec <- t(B)%*% L %*% ecdf.vec

      #Make sure there are no zero eigen values
      spec <- eigen(Dmat, symmetric=TRUE)
      d <- spec$values
      Q <- spec$vectors
      if( min(d) < 10e-6){
        tooSmall <- which( d < 10e-6)
        d[tooSmall] <- d[ min(tooSmall) - 1]
        #Recreate pos. def. Dmat matrix
        Dmat <- Q %*% diag(d) %*% t(Q)
      }

      #Calculate weights
      weights.list[[i]] <- solveWeights(m.test[i],Fn,lower,upper,Dmat,dvec,monotone)

      #Make the pdf
      dumd.now = function(x){
        mix.pdf = function(x){ifelse( (x-lower)*(upper-x) >= 0,
                                      sum(weights.list[[i]]*stats::dbeta((x-lower)/(upper-lower),
                                                                  1:m.test[i],m.test[i]:1))/(upper-lower),0)}
        sapply(x, mix.pdf)
      }

      #Calculate AIC and BIC
      AICs[i]<- -2*sum( log(dumd.now(data))) + 2*(m.test[i]-2)
      BICs[i]<- -2*sum( log(dumd.now(data))) + log(n)*(m.test[i]-2)
    }

    #Which AIC or BIC is the smallest
    if(crit == "AIC"){
      best.m <- which.min(AICs)
    }else if(crit == "BIC"){
      best.m <- which.min(BICs)
    }

    #Set weights and mOpt
    weights <- weights.list[[best.m]]
    mOpt <- m.test[best.m]

    #see if user supplied their own m value
    if( is.na(m) == FALSE){

      if( mOpt < m && warning==TRUE){
        cat("WARNING: given number of weights is larger than optimal number,\n"
            , "\t \t optimal number of weights =", mOpt, "\n")
      }
      if( mOpt > m && warning==TRUE){
        cat("WARNING: given number of weights is less than optimal number,\n"
            ,"\t \t optimal number of weights =", mOpt, "\n")
      }

      #make B, Dmat, and dvec
      B <- NULL
      for(k in 1:m){
        B <- cbind(B, stats::pbeta(tdata, shape1=k, shape2=m-k+1))
      }
      Dmat <- t(B) %*% L %*% B
      dvec <- t(B)%*% L %*% ecdf.vec

      #Make sure there are no zero eigen values
      spec <- eigen(Dmat, symmetric=TRUE)
      d <- spec$values
      Q <- spec$vectors

      #find which values are < 10e-6, and set to smallest
      #eigen value >= 10e-6
      if( min(d) < 10e-6){
        tooSmall <- which( d < 10e-6)
        d[tooSmall] <- d[ min(tooSmall) - 1]
        #Recreate pos. def. Dmat matrix
        Dmat <- Q %*% diag(d) %*% t(Q)
      }

      #Calculate weights
      weights <- solveWeights(m,Fn,lower,upper,Dmat,dvec)
    }
    else{
      #if no m is provided then just set m to mOpt
      #weights are already set
      m <- mOpt
    }

    #Calculate items to return to user:
    #pdf
    dumd = function(x){
      mix.pdf = function(x){ifelse( (x-lower)*(upper-x) >= 0,
                                    sum(weights*stats::dbeta((x-lower)/(upper-lower),1:m,m:1))/(upper-lower),0)}
      sapply(x, mix.pdf)
    }
    #cdf
    pumd = function(x){
      mix.cdf= function(j){sum(weights*stats::pbeta( (j-lower)/(upper-lower), 1:m, m:1)) }
      sapply(x, mix.cdf)
    }
    #rand generator
    rumd = function(n=1){
      rsample = function(){
        k <- sample(1:m,size=1,prob=weights)
        return(  lower+( stats::rbeta(1, k, m-k+1)*(upper-lower) )   )
      }
      replicate(n, rsample() )
    }
    #quantile function
    qumd = function(q){
      mix.quantile = function(q){ g = function(x){ pumd(x)-q }
      return( stats::uniroot(g, interval=c(lower,upper) )$root )
      }
      sapply(q, mix.quantile)
    }

    #Return list (ends function)
    return(list(weights=weights,m.hat=m,dumd=dumd,pumd=pumd,rumd=rumd,qumd=qumd,
                lower=lower, upper=upper))
  }
  #end the function
}
