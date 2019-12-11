####################################################################################################
##### function to estimate the ST-VAR model with regime-dependent AR paras and det. terms only  ####
##________________________________________________________________________________________________##


stvar_fun <- function(x, tv, nlags, MLE = TRUE, NLS = TRUE ,op){

  require(zoomgrid)

  # Transition variable
  z_t <- scale(tv)          # standardize transition variable
  ## ts properties
  tprop <- tsp(op$daten$year)

  op$meanZ <- mean(tv)
  op$stdZ  <- sd(tv)

  op$nlags <- nlags



  # Create regressor matrix, LHS data, and decide the det. terms:
  temp <- lag.create(inp = t(x), p = op$nlags, det_term = op$det_term )
  X    <- temp$X
  Y    <- temp$y
  Tobs <- dim(Y)[1] # # of observations

  if(ncol(x)>=nrow(x)){
    x <- t(x)
  }


  op$Y     <- Y
  op$Tobs  <- Tobs

  #### count parameters:


  number.paras <- para_count(op, R = 2, rd_covar = FALSE)
  print(paste("The ratio of parameters to observations per equation is",number.paras$percent_of_obs, sep = ""))
  if(number.paras$percent_of_obs > 0.25 ){
    stop("Too many parameters to be estimated")
  }

  # Construct/specify transition variable and use lagged values:
  if (op$zlagged==TRUE){
    ends    <- length(z_t)
    diff_zt <- length(z_t) - dim(X)[1]
    Z       <- z_t[diff_zt: (ends-1)]
    Zts     <- ts(Z, start = c(tprop[1] + nlags/tprop[3]), frequency = tprop[3])

  } else {
    begin   <- nlags + 1
    Z       <- z_t[begin:length(z_t)]
    Zts     <- ts(Z, start = c(tprop[1] + nlags/tprop[2]), frequency = tprop[3])
  }


  op$Z    <- Z
  op$Zts  <- Zts
  #### specify intercept switching or constant intercept model in the first order expansion:
  if(op$swint == TRUE){
    swint_flag  = 0
    swint_flag2 = 2  ## for the number of constants to set up the storage matrices later on
    # Xm <- cbind(X, X*Z, sign(Z)*(X * (Z^2)))
    X2 <- X[,-1] ## here we need the X matrix for estimation to contain the intercept column, in order to get the switchung intercept estimators
    Xm <- cbind(rep(1, Tobs), X2, X2*Z, sign(Z)*(X2 * (Z^2)))
  } else {
    swint_flag  = 1
    swint_flag2 = 1
    X <- X[,-1]
    # full set of regressors
    Xm <- cbind(rep(1, Tobs), X, X*Z, sign(Z)*(X * (Z^2)))
  }


  # Build grid: -------------------------------------
  op$X            <- X
  op$Xm           <- Xm
  op$swint_flag   = swint_flag
  op$swint_flag2  = swint_flag2

  # Support for theta and c -------------------------------------------------

  #### use now a transformation of gamma:
  ## gamma = exp(theta) ==> log(gamma) = theta
  op$mintheta    <- log(2) - log(max(abs(Z)))   ## 1/10 of Max theta
  op$maxtheta    <- log(2) - log(max(abs(Z))/op$gbound)
  if(op$step_length_theta == 0){
    op$step_length_theta   <- (op$maxtheta - op$mintheta)/(op$Ntg)
  }

  op$gamma_range <- c(from = op$mintheta, to = op$maxtheta, length = op$Ntg)
  ## support for c:
  op$cmin       <- unname(quantile(Z, probs = op$c_low))
  op$cmax       <- unname(quantile(Z, probs = op$c_high) )
  op$c_range    <- c(from= op$cmin, to= op$cmax, length = op$Ncg)
  ##-------------------------------

  netz <- build_grid(op$gamma_range, op$c_range)


  ###################################################
  ### Grid search for starting values gamma_0 and C_0 for the NLS estimation:

  op$ar_calc <- TRUE

  gs_res <- grid_search(FUN = grid_SSR, netz, MoreArgs = list(paras2 = op), silent = FALSE)


  init_val_nls <- gs_res$par #list(gs_test$par[1], gs_test$par[2])


  ar_init <- get_ar_paras(gs_res$par, paras2 = op, smallsample = op$smallsample, vcv = diag(op$K))
  op$A_hat <- ar_init$all_ar


  #init_val_nls <-  c(1, 0)
  ###************************************************************************###
  ##  Minimize SSR w.r.t. Gamma and c; get starting values for ML estimation  ##

  #### test optim:

  if(NLS==TRUE){

    op$ar_calc == FALSE
    eps = 10^-3
    dSSR = 1000
    iterNLS = 1

    SSRv <- vector("numeric")
    SSRv[[1]] <- ar_init$SSR

    convvec <- vector("numeric")


    while(dSSR > eps & iterNLS <= op$maxiter){

      para_nls <- optim(init_val_nls, grid_SSR, paras2 = op, method = "L-BFGS-B",
                        lower = c(op$mintheta, op$cmin), upper = c(op$maxtheta, op$cmax),
                        control = list(factr = 1e-7, maxit = 50*length(init_val_nls)))

      init_val_nls       <- para_nls$par
      convvec[[iterNLS]] <- para_nls$convergence

      ar_temp <- get_ar_paras(init_val_nls, paras2 = op, smallsample = op$smallsample, vcv = diag(op$K))
      SSRv[[iterNLS+1]] <- ar_temp$SSR

      ##### check for improving SSR compared to Grid search:

      ssr_comp <- SSRv[[1]] < SSRv[[iterNLS +1]]
      print(paste("Overshooting?:", ssr_comp, sep = ""))


      dSSR <- abs(SSRv[[iterNLS + 1]] - SSRv[[iterNLS]])
      print(dSSR)
      iterNLS = iterNLS + 1


    } # while
    print(SSRv)

    op$NLSactive        <- TRUE
    op$convNLS          <- para_nls$convergence
    op$parasNLS         <- para_nls$par
    op$SSRNLS           <- ar_temp$SSR

    stparas <- para_nls

  } # if


  if(MLE==TRUE){

    ###*************************###
    ##  ML estimation step:      ##
    #_____________________________#

    ### numerical optimization loop:
    # first LL value
    llv <- vector("numeric")
    llv[[1]] <- ll_optim(para_nls$par, paras2 = op, pi_mat = ar_temp$all_ar, covar = ar_temp$covar)

    eps = 10^-3
    dLL = 1000
    iter = 1
    while(dLL > eps & iter <= op$maxiter){

      stparas <- optim(stparas$par, ll_optim, paras2 = op, covar = ar_temp$covar, pi_mat = ar_temp$all_ar , method = "L-BFGS-B",
                       lower = c(op$mintheta, op$cmin), upper = c(op$maxtheta, op$cmax),
                       control = list(maxit = 50*2))
      #  print(stparas$value)
      #  print(stparas$convergence)

      ar_temp    <- get_ar_paras(stparas$par, paras2 = op, smallsample = op$smallsample, vcv = ar_temp$covar)

      llv[[iter + 1]] <- ll_optim(stparas$par, paras2 = op, pi_mat = ar_temp$all_ar, covar = ar_temp$covar)

      dLL <- abs(llv[[iter + 1]] - llv[[iter]])
      iter = iter + 1
    }

    op$MLEactive          <- TRUE
    op$convMLE            <- stparas$convergence
    op$bestLL             <- -llv[iter]
    op$LLik               <- llv
    op$iterML             <- iter             ### later: include messages for convergence.
    print(llv)
  }


  ### save output of stvar est:
  op$ratio_paras_obs  <- number.paras$percent_of_obs
  op$A_hat            <- ar_temp$all_ar
  op$A_low            <- ar_temp$arparas_low
  op$A_high           <- ar_temp$arparas_high
  op$sigma_u          <- ar_temp$covar
  op$u                <- ar_temp$kack
  op$SSRMLE           <- ar_temp$SSR
  op$Xm               <- ar_temp$Xm  ## (pseudo) regressor mat.
  op$c                <- stparas$par[2]
  op$theta            <- stparas$par[1]

  class(op) <- "stvar"

  return(op)

} ## for main function



###########################################
### Likelihood function to be optimized ###
#_________________________________________#


## minimize only the neg. LL given the AR-paras and the Covariance mat.

ll_optim <- function(param, paras2, pi_mat, covar){

  # Input: theta and c
  # Output: SSR
  # Calls: Function unvec
  # Outside object: X (extended Regressor mat), Y, K, Tobs, swint, swint_flag, Z
  theta        <- param[[1]] # 1
  c            <- param[[2]] # 0.5

  ###
  Z            <- paras2$Z
  X            <- paras2$X
  Y            <- paras2$Y
  K            <- paras2$K
  Tobs         <- paras2$Tobs
  swint        <- paras2$swint
  swint_flag   <- paras2$swint_flag
  #print(param[[1]])

  # transition equation
  F_Z_t    <- exp(-exp(theta)*(Z-c))/(1 + exp(-exp(theta)*(Z-c)))

  high     <- (1 - F_Z_t) ## wurde vertauscht: wenn z_t -> infinity --> F_zt -> 0
  low      <- F_Z_t ## daher ist 1-F_z_t HIGH
  if(swint == TRUE){
    XM     <- cbind(X*low, X*high)
  } else {
    XM     <- cbind(rep(1, nrow(X)), X*low, X*high)
  }

  # calculate residuals (here gamma and c play a role)
  residsM <- Y - XM%*%pi_mat

  # loglikelihoodvalue
  log_value <- 0
  wgtM <- covar
  for (i in 1: Tobs){
    log_value <- log_value-0.5*((log(det(wgtM))) +
                                  (residsM[i,] %*%solve(wgtM) %*% residsM[i,]))  ## positive LLik
  }

  return(-log_value)
}


# Likelihood and parameter function ---------------------------------------

# 1. Likelihood function which also calculates the AR parameters is needed!

### Function calculating the SSR vbased in Gamma/theta and c.
grid_SSR <- function(param, paras2){

  # Input: theta and c
  # Output: SSR
  # Calls: Function unvec
  # Outside object: X (extended Regressor mat), Y, K, Tobs, swint, swint_flag, Z
  #if(class(param)=="list"){
  #print(param)
  # print(str(paras2, max.level = 1))
  theta        <- param[[1]] # 1
  c            <- param[[2]] # 0.5

  #print(c(theta, c))

  Z            <- paras2$Z
  X            <- paras2$X
  Y            <- paras2$Y
  K            <- paras2$K
  Tobs         <- paras2$Tobs
  swint        <- paras2$swint
  swint_flag   <- paras2$swint_flag
  #print(param[[1]])


  # transition equation
  F_Z_t    <- exp(-exp(theta)*(Z-c))/(1 + exp(-exp(theta)*(Z-c)))

  high     <- (1 - F_Z_t) ## wurde vertauscht: wenn z_t -> infinity --> F_zt -> 0
  low      <- F_Z_t ## daher ist 1-F_z_t HIGH
  if(swint == TRUE){
    XM     <- cbind(X*low, X*high)
  } else {
    XM     <- cbind(rep(1, nrow(X)), X*low, X*high)
  }


  ### new flag: decide if the AR mat should be computed or be pre-specified:

  if(paras2$ar_calc){

    kronsum  <- 0
    vecsum   <- 0

    wgtM <- diag(rep(K))
    for(i in 1:Tobs){
      krons <- kronecker(wgtM, (XM[i,])%*% t(XM[i,]))
      kronsum <- kronsum + krons
      vecs <- as.vector(XM[i,] %*% t(Y[i,]) %*% wgtM)
      vecsum <- vecsum + vecs
    }

    #print(kronsum)

    pis <- solve(kronsum)%*%vecsum
    pi_mat <- unvec(pis, dim(X)[2]*2+swint_flag, K) ## returns a matrix with each column containing all parameters belonging to a variable + possibly the exogenous variables
  } else {

    pi_mat <- paras2$A_hat
  }

  residsM <- Y - XM%*%pi_mat #U(T x K)%*%U' (K x T)
  vecU    <- matrix(as.vector(t(residsM)), ncol = 1)
  SSR     <- t(vecU)%*%vecU
  # print(SSR)
  return(SSR)
}


### Function calculating AR parameters and SSR:

get_ar_paras <- function(param, paras2, smallsample = FALSE, vcv){

  # Input: theta and c
  # Output: SSR
  # Calls: Function unvec
  # Outside object: X (extended Regressor mat), Y, K, Tobs, swint, swint_flag, Z
  if(class(param)=="list"){
    theta        <- param[[1]] # 1
    c            <- param[[2]] # 0.5
  } else {
    theta        <- param[1] # 1
    c            <- param[2] # 0.5
  }

  p            <- paras2$nlags
  Z            <- paras2$Z
  X            <- paras2$X
  Y            <- paras2$Y
  K            <- paras2$K
  Tobs         <- paras2$Tobs
  swint        <- paras2$swint
  swint_flag   <- paras2$swint_flag
  #print(param[[1]])

  # transition equation
  F_Z_t    <- exp(-exp(theta)*(Z-c))/(1 + exp(-exp(theta)*(Z-c)))

  high     <- (1 - F_Z_t) ## wurde vertauscht: wenn z_t -> infinity --> F_zt -> 0
  low      <- F_Z_t ## daher ist 1-F_z_t HIGH
  if(swint == TRUE){
    XM     <- cbind(X*low, X*high)
  } else {
    XM     <- cbind(rep(1, nrow(X)), X*low, X*high)
  }

  # wgtM     <- list()
  dims     <- c(K*dim(XM)[2], K*dim(XM)[2], Tobs)
  krons    <- matrix(0, dims[1], dims[2])
  kronsum  <- krons
  vecs     <- rep(0, dims[1])
  vecsum   <- vecs

  wgtM <- solve(vcv)
  for(i in 1:Tobs){
    krons <- kronecker(wgtM, (XM[i,])%*% t(XM[i,]))
    kronsum <- kronsum + krons
    vecs <- as.vector(XM[i,] %*% t(Y[i,]) %*% wgtM)
    vecsum <- vecsum + vecs
  }

  pis <- solve(kronsum)%*%vecsum
  pi_mat <- unvec(pis, dim(X)[2]*2+swint_flag, K) ## returns a matrix with each column containing all parameters belonging to a variable + possibly the exogenous variables
  Kp    <- K*nlags
  if(swint==TRUE){
    pi_l <- t(pi_mat[c(2:(Kp+1)),])
    pi_h <- t(pi_mat[c((Kp + 3):(2*Kp+2)),])
    const <- t(pi_mat[c(1, Kp+2),])
  } else {
    pi_l <- t(pi_mat[c(2:(Kp+swint_flag)),])
    pi_h <- t(pi_mat[(Kp + 2):(2*Kp+1),])
    const <- t(pi_mat[1,])
  }

  residsM <- Y - XM%*%pi_mat # T x K
  UU      <- t(residsM)%*%residsM
  vecU    <- matrix(as.vector(t(residsM)), ncol = 1) #as.vector() is = vec(), need to get resids as K x T mat!
  SSR     <- t(vecU)%*%vecU
  if(smallsample==TRUE){
    covar <- UU/(Tobs - 1 - 2*K*p)
  } else {
    covar   <- UU/Tobs
  }

  ### positive log L value:
  log_value <- 0
  #  wgtM <- covar
  #  for (i in 1: Tobs){
  #    log_value <- log_value-0.5*((log(det(wgtM))) +
  #                                  (residsM[i,] %*%solve(wgtM) %*% residsM[i,]))  ## positive LLik
  #  }
  return(list(SSR = SSR,
              arparas_low     = pi_l, # only AR paras
              arparas_high    = pi_h, # only AR paras
              all_ar          = pi_mat, # all parameters (AR + det.) in rowform/long form, i.e., (dt + 2*K*p x K) where dt = number of det. terms!
              covar           = covar,
              negLL           = -log_value, ## not available (yet)
              kack            = residsM,
              detterms        = const,     # rowform matrix of det. terms, i.e., (#of det terms x K)
              Xm              = XM))
}
