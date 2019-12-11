##########################################################
### modify the ident dCov function from svars package: ###
###____________________________________________________###


id.dc <- function(x, PIT = FALSE){
  require("steadyICA")
  # if(is.null(residuals(x))){
  #   stop("No residuals retrieved from model")
  # }
  R <- u <- Tob <- p <- k <- residY <- coef_x <- yOut <- type <- y <-  NULL
  get_var_objects(x)
  ########### starting the computations ------------------------------------------------------------------------

  sigg <- crossprod(u) / (Tob- 1 - R*k * p)

  # Choleski decomposition
  P_chol <- t(chol(sigg))

  # minimize dCov with 'steadyICA'
  u_chol <- t(solve(P_chol)%*%t(u))
  ICA <- suppressMessages(steadyICA(u_chol, symmetric=TRUE, PIT=PIT))

  # structural matrix Sigma_u = BB'
  P <- P_chol%*%ICA$W

  # obtaining VAR parameter
  if(inherits(x, "var.boot") | inherits(x, "stvar")){ ## modified for stvar
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

  result <- list(B = P,       # estimated B matrix (unique decomposition of the covariance matrix)
                 A_hat = A_hat,  # estimated VAR parameter
                 method =        "Distance covariances",
                 n = Tob,        # number of observations
                 type = type,    # type of the VAR model e.g 'const'
                 y = yOut,       # Data
                 p = unname(p),  # number of lags
                 K = k,          # number of time series
                 PIT=PIT,        #
                 R = R           # number of regimes, only for ST-VARs (otherwise always 1)
  )
  if (inherits(x, "stvar")){
    class(result) <- c("svars", "stvar")
  } else {
    class(result) <- "svars"
  }
  return(result)
}



