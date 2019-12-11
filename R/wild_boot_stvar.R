######################################################################
### For ST-VAR adjusted wild bootstrap function from svars package ###
####______________________________________________________________####



wild.boot <- function(x, design = "fixed", distr = "rademacher", n.ahead = 20,
                      nboot = 500, nc = 1, dd = NULL, signrest = NULL, itermax = 300,
                      steptol = 200, iter2 = 50, rademacher = "deprecated", redfstvar = NULL){

  #require(pbapply)
  #require(steadyICA)
  #require(expm)

  # x: vars object
  # B: estimated covariance matrix from true data set
  # distr: whether the bootstraop works with gaussian, rademacher or mammen distribution
  # design: choice of recursive or fixed design for bootstrap
  # n.ahead: Time n.ahead for Irf
  # nboot: number of bootstrap replications
  if(x$method == "Cramer-von Mises distance" & is.null(dd)){
    dd <- copula::indepTestSim(x$n, x$K, verbose=F)
  }

  sqrt.f <- function(Pstar, Sigma_u_star){
    yy <- suppressMessages(sqrtm(Sigma_u_hat_old))%*%solve(suppressMessages(sqrtm(Sigma_u_star)))%*%Pstar
    return(yy)
  }


  # gathering informations from vars object
  y <- x$y  ## (TxK)
  p <- x$p
  obs <- x$n
  k <- x$K
  B <- x$B
  R <- ifelse(is.null(x$R), 1, x$R)  ## NEW: number of Regimes (stvar)
  restriction_matrix = x$restriction_matrix
  restriction_matrix <- get_restriction_matrix(restriction_matrix, k)
  restrictions <- length(restriction_matrix[!is.na(restriction_matrix)])

  if(length(signrest) > k){
    stop('too many sign restrictions')
  }

  # calculating covariance from actual VAR

  ##__________ change here___________________##
  if(inherits(x, "stvar")){
    A <-  x$A_hat   ## A_hat is (2*Kp x K)  --> A (K x 2*Kp) not anymore!
    #print(dim(A))
    Z <-  t(redfstvar$Xm)  ## Xm is (T x 2*Kp)
    #print(dim(Z))
    u <-  t(redfstvar$u)
    Sigma_u_hat_old <- tcrossprod(u)/(obs - 1 - R*k * p)
    #print(Sigma_u_hat_old)
    swint_flag <- redfstvar$swint_flag
  } else {
    #####___________________end change _______________###
    A <- x$A_hat
    Z <- t(YLagCr(y, p))   ### YlagCr returns (T x K) reg. mat., while Z is (#no of reg. x T) reg. mat

    if(x$type == 'const'){
      Z <- rbind(rep(1, ncol(Z)), Z)
    }else if(x$type == 'trend'){
      Z <- rbind(seq(p + 1, ncol(Z)+ p), Z)
    }else if(x$type == 'both'){
      Z <- rbind(rep(1, ncol(Z)), seq(p + 1, ncol(Z) + p), Z)
    }else{
      Z <- Z
    }

    u <- t(y[-c(1:p),]) - A %*% Z
    Sigma_u_hat_old <- tcrossprod(u)/(obs - 1 - k * p)
  }


  ub <- u

  # creating new error terms
  errors <- list()

  if(rademacher != "deprecated"){
    if(rademacher == "TRUE"){
      warning("The argument 'rademacher' is deprecated and may not be supported in the future. Please use the argument 'distr' to decide upon a distribution.",
              call. = TRUE, immediate. = FALSE, noBreaks. = FALSE, domain = NULL)
    } else if(rademacher == "FALSE"){
      distr <- "gaussian"
      warning("The argument 'rademacher' is deprecated and may not be supported in the future. Please use the argument 'distr' to decide upon a distribution.",
              call. = TRUE, immediate. = FALSE, noBreaks. = FALSE, domain = NULL)
    } else{
      warning("Invalid use of deprecated argument 'rademacher'. Please use the argument 'distr' to decide upon a distribution!",
              call. = TRUE, immediate. = FALSE, noBreaks. = FALSE, domain = NULL)
    }
  }

  for(i in 1:nboot){
    ub <- u
    #my <- rnorm(1)
    if (distr == "rademacher") {
      my <- rnorm(n = ncol(u))
      my <- (my > 0) - (my < 0)
    } else if (distr == "mammen") {
      cu <- (sqrt(5)+1)/(2*sqrt(5))
      my <- rep(1,ncol(u))*(-(sqrt(5)-1)/2)
      uni <- runif(n = ncol(u), min = 0, max = 1)
      my[uni > cu] <- (sqrt(5)+1)/2
    } else if (distr == "gaussian") {
      my <- rnorm(n = ncol(u))
    }
    errors[[i]] <- ub* my
  }

  # print(dim(errors[[1]]))

  # Ustar1 <- errors[[1]]

  # Bootstrapfunction
  bootf <- function(Ustar1){

    if(design == "fixed"){

      ###___________ change here ______________________###  should be fine now!

      Ystar <- t(A %*% Z + Ustar1)  ## (TxK)

      if(inherits(x, "stvar")){
        XM <- t(Z)
        ####
        kronsum  <- 0
        vecsum   <- 0
        if(is.null(x$MLEactive)){
          wgtM <- diag(k)
        } else {
          wgtM  <- solve(Sigma_u_hat_old)
        }
        for(i in 1:x$n){
          krons <- kronecker(wgtM, XM[i,]%*% t(XM[i,]))
          kronsum <- kronsum + krons
          vecs <- as.vector(XM[i,] %*% t(Ystar[i,]) %*% wgtM)
          vecsum <- vecsum + vecs
        }
        pis <- solve(kronsum)%*%vecsum
        Bstar <- t(unvec(pis, dim(redfstvar$X)[2]*2+redfstvar$swint_flag, k)) ## A_hat = [K x Kp] here, in stvar it is [Kp x K]
        ###
      } else {
        Bstar <- t(Ystar) %*% t(Z) %*% solve(Z %*% t(Z))
      }

      Ustar <- Ystar - t(Bstar %*% Z)
      print(dim(Ustar))
      Sigma_u_star <- crossprod(Ustar)/(ncol(Ustar1) - 1 - R*k*p)
      #print(Sigma_u_star)

      ###____________end _______________________________###

      varb <- list(y = Ystar,
                   coef_x = Bstar,
                   residuals = Ustar,
                   p = p,
                   type = x$type)

      if(inherits(x, "stvar")){
        class(varb)   <- c("var.boot", "stvar")  ## NEW
      } else {
        class(varb)   <- 'var.boot'
      }
      if(x$method == "Non-Gaussian maximum likelihood"){
        temp <- id.ngml_boot(varb, stage3 = x$stage3, Z = Z, restriction_matrix = x$restriction_matrix)
      }else if(x$method == "Changes in Volatility"){
        temp <- tryCatch(id.cv_boot(varb, SB = x$SB, Z = Z, restriction_matrix = x$restriction_matrix),
                         error = function(e) NULL)
      }else if(x$method == "Cramer-von Mises distance"){
        temp <- id.cvm(varb, itermax = itermax, steptol = steptol, iter2 = iter2, dd)
      }else if(x$method == "Distance covariances"){
        temp <- id.dc(varb, PIT=x$PIT)
      }else if(x$method == "GARCH"){
        temp <- tryCatch(id.garch(varb, restriction_matrix = x$restriction_matrix, max.iter = x$max.iter,
                                  crit = x$crit),
                         error = function(e) NULL)
      }else if(x$method == "Recursive"){
        temp <- id.chol(varb)
      } else {
        temp <- tryCatch(id.st_boot(varb, c_fix = x$est_c, transition_variable = x$transition_variable, restriction_matrix = x$restriction_matrix,
                                    gamma_fix = x$est_g, max.iter = x$iteration, crit = 0.01, Z = Z),
                         error = function(e) NULL)
      }
      temp$Ustar <- Ustar
    } else if (design == "recursive") {
      Ystar <- matrix(0, nrow(y), k)

      # adding pre sample values
      Ystar[1:p,] <- y[1:p,]

      if (x$type == 'const' | x$type == 'trend') {
        for (i in (p + 1):nrow(y)) {
          for (j in 1:k) {
            Ystar[i, j] <- A[j, 1] + A[j, -1] %*% c(t(Ystar[(i - 1):(i - p), ])) + Ustar1[j, (i - p)]
          }
        }
      } else if (x$type == 'both') {
        for (i in (p + 1):nrow(y)) {
          for (j in 1:k) {
            Ystar[i, j] <- A[j, 1] + A[j, 2] + A[j, -c(1, 2)] %*% c(t(Ystar[(i - 1):(i - p),])) + Ustar1[j, (i - p)]
          }
        }
      }else if (x$type == 'none') {
        for (i in (p + 1):nrow(y)) {
          for (j in 1:k) {
            Ystar[i, j] <- A[j, ] %*% c(t(Ystar[(i - 1):(i - p), ])) + Ustar1[j, (i - p)]
          }
        }
      }

      # Delete pre sample values
      Ystar <- Ystar[-c(1:p), ]

      ##____________________ change here:  first omit this.. only fixed design

      varb <- suppressWarnings(VAR(Ystar, p = x$p, type = x$type))
      Ustar <- residuals(varb)
      print(dim(Ustar))
      ##_________end change ___________##

      Sigma_u_star <- crossprod(Ustar)/(obs - 1 - k * p)

      if(x$method == "Non-Gaussian maximum likelihood"){
        temp <- id.ngml_boot(varb, stage3 = x$stage3, restriction_matrix = x$restriction_matrix)
      }else if(x$method == "Changes in Volatility"){
        if (length(x$SB) > 3) {
          SB <- x$SB[-c(1:p)]
        } else {
          SB <- x$SB
        }
        temp <- tryCatch(id.cv_boot(varb, SB = SB, restriction_matrix = x$restriction_matrix),
                         error = function(e) NULL)
      }else if(x$method == "Cramer-von Mises distance"){
        temp <- id.cvm(varb, itermax = itermax, steptol = steptol, iter2 = iter2, dd)
      }else if(x$method == "Distance covariances"){
        temp <- id.dc(varb, PIT=x$PIT)
      }else if(x$method == "Smooth transition"){
        temp <- id.st(varb, c_fix = x$est_c, transition_variable = x$transition_variable, restriction_matrix = x$restriction_matrix,
                      gamma_fix = x$est_g, max.iter = x$iteration, crit = 0.01)
      }else if(x$method == "GARCH"){
        temp <- tryCatch(id.garch(varb, restriction_matrix = x$restriction_matrix, max.iter = x$max.iter,
                                  crit = x$crit),
                         error = function(e) NULL)
      }else {
        temp <- id.chol(varb)
        # print(temp$B)
      }
      temp$Ustar <- Ustar
    } ## for design == recursive
    if(!is.null(temp)){
      if(x$method != "Recursive"){

        Pstar <- temp$B

        if(restrictions > 0){
          Pstar1 <- Pstar
          frobP <- frobICA_mod(Pstar1, B, standardize=TRUE)
        }else {
          Pstar1 <- sqrt.f(Pstar, Sigma_u_star)
          diag_sigma_root <- diag(diag(suppressMessages(sqrtm(Sigma_u_hat_old))))

          frobP <- frobICA_mod(t(solve(diag_sigma_root)%*%Pstar1), t(solve(diag_sigma_root)%*%B), standardize=TRUE)
        }
        Pstar <- Pstar1%*%frobP$perm
        temp$B <- Pstar
      } else {

        Pstar <- temp$B

      }
      ###______-Change here______: Use the irf_quant function??

      if(inherits(x, "stvar")){
        #class(temp) <- c("var.boot", "svars")
        ip <- irf_quant(x = temp, redf = redfstvar, quantiles = list(low = 0, high = 1), h = n.ahead) ## quantiles could be changed to wild.boot argument
        #class(temp) <- c("var.boot", "svars")
      } else {

        ip <- irf(temp, n.ahead = n.ahead)
      }
      ####____ end change ___________###

      return(list(ip, Pstar, temp$A_hat, temp))  ## NEW: temp is also returned (i.e., the svars BS object)
    }else{
      return(NA)
    }
  }

  bootstraps <- pblapply(errors, bootf, cl = nc)

  delnull  <-  function(x){
    x[unlist(lapply(x, length) != 0)]
  }

  bootstraps <- lapply(bootstraps, function (x)x[any(!is.na(x))])
  bootstraps <- delnull(bootstraps)

  Bs       <- array(0, c(k,k,length(bootstraps)))
  ipb      <- list()
  svarsobj <- list() ## NEW


  ## Obtaining Bootstrap estimates of VAR parameter
  Aboot <- array(0, c(nrow(A), ncol(A),length(bootstraps)))

  for(i in 1:length(bootstraps)){
    Bs[,,i] <- bootstraps[[i]][[2]]
    ipb[[i]] <- bootstraps[[i]][[1]]
    Aboot[, , i] <- bootstraps[[i]][[3]]
    svarsobj[[i]] <- bootstraps[[i]][[4]] ## NEW
  }

  A_hat_boot <- matrix(Aboot, ncol = nrow(A)*ncol(A), byrow = TRUE)
  A_hat_boot_mean <- matrix(colMeans(A_hat_boot), nrow(A), ncol(A))

  # calculating covariance matrix of vectorized bootstrap matrices
  v.b      <-  matrix(Bs, ncol = k^2, byrow = TRUE)
  cov.bs   <- cov(v.b)

  # Calculating Standard errors for LDI methods
  #if(x$method == "Cramer-von Mises distance" | x$method == "Distance covariances"){
  SE <- matrix(sqrt(diag(cov.bs)),k,k)
  rownames(SE) <- rownames(x$B)
  #}else{
  #  SE <- NULL
  # }

  # Calculating Bootstrap means
  boot.mean <- matrix(colMeans(v.b),k,k)
  rownames(boot.mean) <- rownames(x$B)

  # Checking for signs
  if(restrictions > 0 | x$method == "Recursive"){
    if(!is.null(signrest)){
      cat('Testing signs only possible for unrestricted model \n')
    }
    sign.part <- NULL
    sign.complete <- NULL
  }else{
    if(is.null(signrest)){
      sign.mat <- matrix(FALSE, nrow = k, ncol = k)
      sign.complete <- 0
      sign.part <- rep(0, times = k)

      for(i in 1:length(bootstraps)){

        pBs <- permutation(Bs[,,i])
        sign.mat <-lapply(pBs, function(z){sapply(1:k, function(ii){all(z[,ii]/abs(z[,ii])  == x$B[,ii]/abs(x$B[,ii])) | all(z[,ii]/abs(z[,ii])  == x$B[,ii]/abs(x$B[,ii])*(-1))})})
        #print(pBs[1:length(pBs)])

        if(any(unlist(lapply(sign.mat, function(sign.mat)all(sign.mat == TRUE))))){
          sign.complete <- sign.complete + 1
        }

        for(j in 1:k){
          check <- rep(FALSE, k)
          for(l in 1:k){
            check[l] <- any(all(pBs[[1]][,l]/abs(pBs[[1]][,l]) == x$B[,j]/abs(x$B)[,j]) | all(pBs[[1]][,l]/abs(pBs[[1]][,l]) == x$B[,j]/abs(x$B)[,j]*(-1)))
          }
          if(sum(check) == 1){
            sign.part[[j]] <- sign.part[[j]] + 1
          }
        }
      }
    }else{
      nrest <- length(signrest)
      sign.part <- rep(list(0), nrest )
      sign.complete <- 0
      for(j in 1:length(bootstraps)){
        check.full <- 0
        for(i in 1:nrest){
          check <- rep(FALSE, length(signrest[[i]][!is.na(signrest[[i]])]))
          for(l in 1:k){
            check[l] <- any(all(Bs[!is.na(signrest[[i]]),l,j]/abs(Bs[!is.na(signrest[[i]]),l,j]) == signrest[[i]][!is.na(signrest[[i]])]) |
                              all(Bs[!is.na(signrest[[i]]),l,j]/abs(Bs[!is.na(signrest[[i]]),l,j]) == signrest[[i]][!is.na(signrest[[i]])]*(-1)))
          }
          if(sum(check) == 1){
            sign.part[[i]] <- sign.part[[i]] + 1
            check.full <- check.full + 1
          }
        }
        if(check.full == nrest){
          sign.complete <- sign.complete + 1
        }
      }
      names(sign.part) <- names(signrest)
    }
  }

  ## Impulse response of actual model

  ####__________ change here: IRF_quant

  if(inherits(x, "stvar")){

    ip <- irf_quant(x = x, redf = redfstvar, quantiles = list(low = 0, high = 1), h = n.ahead) ## quantiles could be changed to wild.boot argument

  } else {
    ip <- irf(x, n.ahead = n.ahead)
  }

  ####_______________________________

  result <- list(true = ip,
                 bootstrap = ipb,
                 SE = SE,
                 nboot = nboot,
                 distr = distr,
                 point_estimate = x$B,
                 boot_mean = boot.mean,
                 signrest = signrest,
                 sign_complete = sign.complete,
                 sign_part = sign.part,
                 cov_bs = cov.bs,
                 A_hat = x$A_hat,
                 design = design,
                 A_hat_boot_mean = A_hat_boot_mean,
                 Omodel   = x,
                 boot_B   = Bs,
                 Allbs = svarsobj,  ## NEW: List of id.dc, svars objects, for later use with IRF_quant
                 method   = 'Wild bootstrap')
  if(inherits(x, "stvar")){
    class(result) <- 'sbootst'
  } else {
    class(result) <- 'sboot'
  }
  return(result)
}


