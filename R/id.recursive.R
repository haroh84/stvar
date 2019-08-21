#--------------------------------------------------#
## Recursive identification                       ##
#--------------------------------------------------#

id.chol <- function(x){

  # if(is.null(residuals(x))){
  #   stop("No residuals retrieved from model")
  # }
  u <- Tob <- p <- k <- residY <- coef_x <- yOut <- type <- y <-  NULL
  get_var_objects(x)
  ########### starting the computations ------------------------------------------------------------------------

  sigg <- crossprod(u) / (Tob- 1 - k * p)

  # Choleski decomposition
  P_chol <- t(chol(sigg))

  # minimize dCov with 'steadyICA'
  u_chol <- t(solve(P_chol)%*%t(u))


  # obtaining VAR parameter
  if(inherits(x, "var.boot")){
    A_hat <- coef_x
  }else{
    A <- matrix(0, nrow = k, ncol = k * p)
    for(i in 1:k){
      A[i,] <- coef_x[[i]][1:(k * p),1]
    }

    A_hat <- A

    if(type == "const"){
      v <- rep(1, k)

      for(i in 1:k){
        v[i] <- coef_x[[i]][(k*p+1), 1]
      }

      A_hat <- cbind(v, A)
    }else if (type == "trend"){
      trend <- rep(1, k)

      for(i in 1:k){
        trend[i] <- coef_x[[i]][(k*p+1), 1]
      }

      A_hat <- cbind(trend, A)
    }else if(type == "both"){
      v <- rep(1, k)

      for(i in 1:k){
        v[i] <- coef_x[[i]][(k*p+1), 1]
      }

      trend <- rep(1, k)

      for(i in 1:k){
        trend[i] <- coef_x[[i]][(k*p+2), 1]
      }

      A_hat <- cbind(v, trend, A)
    }
  }

  result <- list(B = P_chol,       # estimated B matrix (unique decomposition of the covariance matrix)
                 A_hat = A_hat,  # estimated VAR parameter
                 method =        "Recursive",
                 n = Tob,        # number of observations
                 type = type,    # type of the VAR model e.g 'const'
                 y = yOut,       # Data
                 p = unname(p),  # number of lags
                 K = k          # number of time series

  )
  class(result) <- "svars"
  return(result)
}
