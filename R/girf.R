############# computation of GIRFs


girf_bs <- function(inp, redf, bsrep = 100, design = c("automatic", "semi-fixed"),nc = 1,
                    quant = NULL, period = NULL, shock, transvar_pos, size, n.ahead = 10, iter1, iter2){
  ## inp = sbootst object
  ## bsrep = number of bootstrap replications
  ## nc = number of cores used for parallel computing

  ## choose desing:
  design <- match.arg(design)
  if(design=="automatic"){

    design <- inp$design
  }
  ## set up modified stvar object:
  redfbs <- list()
  rftmp <- redf
  for(i in 1:bsrep){
    if(design%in%c("recursive", "semi-fixed")){
      #rftmp$Z      <-     ### semi-fixed ergibt wohl keinen sinn...
      #rftmp$Zts
    }
    if(design=="recursive"){
      #rftmp$X    <- ###
    }
    rftmp$A_hat  <- t(inp$Allbs[[i]]$A_hat) ## muss transponiert werden, da die GIRF fun. [Kp x K] braucht!
    rftmp$B      <-   inp$Allbs[[i]]$B    ## lieber B aus der rftmp beziehen, dann ist gesichert, dass später die zuordnung korrekt ist!
    rftmp$ystar  <- inp$Allbs[[i]]$y
    rftmp$u      <- inp$Allbs[[i]]$Ustar  #rftmp$ystar - t(rftmp$A_hat %*% rftmp$X)

    redfbs[[i]] <- rftmp
  }
  redfbs[[1]]$B
  ####______ Compute GIRF for the original model:

  girf_true <- girf(x = redf, B = inp$Omodel$B, quant = quant, period = period, shock = shock, transvar_pos = transvar_pos, size = size,
                    n.ahead = n.ahead, iter1 = iter1, iter2 = iter2 )

  ####______ BS parallel computing loop here:

  # set up girf fun para list:
  parasgirf <- list( B = NULL,
                     quant = quant,
                     period = period,
                     shock = shock,
                     transvar_pos = transvar_pos,
                     size = size,
                     n.ahead = n.ahead,
                     iter1 = iter1,
                     iter2 = iter2 )


  tmp <- pblapply(redfbs, function(x) girf(x = x, B = NULL, quant = quant, period = period,
                                           shock = shock, transvar_pos = transvar_pos, size, n.ahead = n.ahead, iter1 = iter1, iter2 = iter2 ), cl = nc)


  if(!is.null(quant)){
    girfbs  <- lapply(tmp, function(x) x[c("irflow", "irfhigh", "irfdiff", "Fresp_low", "Fresp_high")])
    method  <- "regimes"
  } else if(!is.null(period)){
    girfbs <- lapply(tmp, function(x) x[c("irflow", "Fresp")])
    method <- "irfsingle"
  } else {
    girfbs <- lapply(tmp, function(x) x[c("irf", "Fresp")])
    method <- "irf3D"
  }

  ####------ to use for output


  out <- list(true      = girf_true[c("irflow", "irfhigh", "irfdiff", "Fresp_low", "Fresp_high")],
              bootstrap      = girfbs,       ###  get from bootstrap => only keep irf$low, irf$high, Fresp_low etc...
              nboot          = bsrep,         ### from bsrep
              distr          = inp$distr,         ### from inp
              point_estimate = inp$Omodel$B,  ### from 0model
              design         = design,    ## take from design
              Omodel         = girf_true,      ### girf_true saves all GIRF related base info!
              method         = method)  ### change here the method depending on the input: history, quant or period!

  class(out) <- c("sbootst", "girf")

}

#### outerfun, for GIRFs without

girf <- function(x, B = NULL, quant = NULL, period = NULL, shock, transvar_pos, size, n.ahead = 10, iter1, iter2 ){
  # x is redf. stvar object, trying to use as much info from x as possible
  theta       <- x$theta
  c           <- x$c
  Zts         <- x$Zts
  Z           <- x$Z
  u           <- x$u
  A_hat       <- x$A_hat
  Xm          <- x$Xm
  X           <- x$X
  Y           <- x$Y
  Zmean       <- x$meanZ
  Zstd        <- x$stdZ

  print(B)
  if(is.null(B)){
    B  <- x$B
  } else {
    B <- B
  }

  ## needed object:
  F_Z_t     <- exp(-exp(theta)*(Z-c))/(1 + exp(-exp(theta)*(Z-c)))
  F_Z_t     <- ts(F_Z_t, start = start(Zts), frequency = tsp(Zts)[3])
  #high     <- (1 - F_Z_t) ## wurde vertauscht: wenn z_t -> infinity --> F_zt -> 0
  #low      <- F_Z_t ## daher ist 1-F_z_t HIGH
  histories <- list()
  #histories <- asplit(x$X, c(1))  ## get all histories
  for(i in 1:x$Tobs){
    temp <- list()
    temp$Z               <- Z[i]
    temp$Y               <- Y[i,]
    temp$X               <- X[i,]
    temp$F_z_t           <- F_Z_t[i]
    temp$Xm              <- Xm[i,]
    histories[[i]]       <- temp
  }
  epsilon   <- t(solve(B)%*%t(u))

  ##### gather info for output:
  out <- list(B               = B,
              Regime_quantile = quant,
              Period          = period,
              shock           = shock,
              TV_pos_in_Y     = transvar_pos,
              Shock_size      = size,
              n.ahead         = n.ahead,
              inner_iter      = iter1,
              outer_iter      = iter2)

  # ===> !!! Umdenken: X is needed not the extended regressor mat !!! <===


  ##________________inner function to compute the single GIRFs for a given history by MC integration________________##

  hgirf <- function(hist1, iter1){

    #hist1 <- histories[[1]]
    ### Draw shocks
    ### Simulate Y* for 1,...,H
    ### Calc GIRF
    ### repeat Iter1 times
    F_response   <- matrix(NA, nrow = (n.ahead+1), ncol = iter1)
    y_girf       <- array(NA, dim = c((n.ahead+1), ncol(epsilon), iter1))  ## store single IRFs
    eps_draws_wo <- array(NA, dim = c((n.ahead+1), ncol(epsilon), iter1))
    eps_draws_ws <- array(NA, dim = c((n.ahead+1), ncol(epsilon), iter1))
    ##____ draw Iter1 matrices of shocks with dimension (H x K) and store in list

    for(i in 1:iter1){
      eps_draws_wo[,,i]                            <- t(B%*%t(epsilon[sample(nrow(epsilon), size=(n.ahead+1), replace=T),]))
      eps_draws_ws[1,,i]                           <- epsilon[sample(nrow(epsilon), size=1, replace=T),]
      eps_draws_ws[1,shock,i]                      <- eps_draws_ws[1,shock,i] + size
      eps_draws_ws[1,,i]                           <- B%*%eps_draws_ws[1,,i]   ## get reduced form error for pertubated shock!
      eps_draws_ws[2:dim(eps_draws_ws)[1],,i]      <- t(B%*%t(epsilon[sample(nrow(epsilon), size=(n.ahead), replace=T),]))
    }

    ##____ for loop to iterate the system

    ## h = 0: baseline for history
    for(j in 1:iter1){
      #j = 1
      y_wo <- matrix(NA, ncol = dim(eps_draws_ws)[2], nrow = (n.ahead+1))
      y_ws <- matrix(NA, ncol = dim(eps_draws_ws)[2], nrow = (n.ahead+1))
      x_wo <- matrix(NA, ncol = ncol(X), nrow = (n.ahead+2))
      x_wo[1,] <- hist1$X
      x_ws <- matrix(NA, ncol = ncol(X), nrow = (n.ahead+2))
      x_ws[1,] <- hist1$X
      #F_wo  <- vector("numeric", (n.ahead+1))
      # F_ws  <- vector("numeric", (n.ahead+1))
      F_wo <- F_ws <- hist1$F_z_t



      for (i in 1:(n.ahead+1)){
        # i = 1
        y_wo[i,]        <- c(x_wo[i,]*F_wo[i], x_wo[i,]*(1-F_wo[i]))%*%A_hat + eps_draws_wo[i,,j]
        y_ws[i,]        <- c(x_ws[i,]*F_ws[i], x_ws[i,]*(1-F_ws[i]))%*%A_hat + eps_draws_ws[i,,j]
        Z1_wo           <- (y_wo[i,transvar_pos]- Zmean)/Zstd                 ### first only option to normalize!
        Z1_ws           <- (y_ws[i,transvar_pos]- Zmean)/Zstd
        F_wo[i+1]          <- exp(-exp(theta)*(Z1_wo-c))/(1 + exp(-exp(theta)*(Z1_wo-c)))
        F_ws[i+1]            <- exp(-exp(theta)*(Z1_ws-c))/(1 + exp(-exp(theta)*(Z1_ws-c)))

        #high     <- (1 - F_Z_t)
        #low      <- F_Z_t
        x_wo[i+1,]  <- c(1, y_wo[i,], x_wo[i,2:(1+ncol(y_wo))])                  ### account later for det. terms other than switching constant and exo vars..
        x_ws[i+1,]  <- c(1, y_ws[i,], x_ws[i,2:(1+ncol(y_ws))])
      }
      #### compute difference between shocked and non-shocked IRFs for each horizon:
      ## as in Caggiano et al.: Y* - Y*_shocked
      F_response[,j]  <- F_wo[-length(F_wo)] - F_ws[-length(F_ws)]
      y_girf[,,j]     <- y_wo - y_ws

    }

    ## calculate the average response of all variables
    GIRF_h     <- apply(y_girf, c(1,2), mean)
    Fresp_mean <- apply(F_response, 1, mean)
    #comb       <- cbind(Fresp_mean, GIRF_h)
    out        <- list(girf   = GIRF_h,
                       Fresp  = Fresp_mean) # list(girf   = comb)

    ###  possbility to add, e.g., TF response as list object, or explosiveness flag || Fresp  = Fresp_mean
    return(out)
  }


  ##______________________specification part and estimation part______________________________##

  if(is.null(quant)==FALSE){
    qsplit    <- unname(quantile(Z, probs = quant))
    F_q       <- exp(-exp(theta)*(qsplit - c))/(1 + exp(-exp(theta)*(qsplit - c)))
    F_ind     <- which(F_Z_t>= F_q) ## low regime
    hist_low  <- histories[F_ind]
    hist_high <- histories[which(F_Z_t < F_q)]

    ### sample iter2 histories from each "regime":
    hist_low_sample  <- hist_low[sample(length(hist_low), size=iter2, replace=T)]
    hist_high_sample <- hist_high[sample(length(hist_high), size = iter2, replace = TRUE)]
    ### calculate GIRFs for both "regimes" iter2 times:
    #temp             <- sapply(hist_low_sample, function(x) hgirf(x, iter1 = iter1), simplify = 'array')

    temp         <-  lapply(hist_low_sample, function(x) hgirf(x, iter1 = iter1))
    girfs_low    <-  sapply(temp, function(x) x$girf , simplify = "array")
    Fresp_low    <-  sapply(temp, function(x) x$Fresp , simplify = "array")

    ## use this, if hirf(.) output is a list:
    #girf_low
    ### store first in list and then use sapply "[[" to get arrays !??


    #girfs_high                <-  sapply(hist_high_sample, function(x) hgirf(x, iter1 = iter1), simplify = 'array')
    temp         <-  lapply(hist_high_sample, function(x) hgirf(x, iter1 = iter1))
    girfs_high    <-  sapply(temp, function(x) x$girf , simplify = "array")
    Fresp_high   <-  sapply(temp, function(x) x$Fresp , simplify = "array")

    ##____________Regime Difference, high - low

    girfs_diff <- girfs_high - girfs_low

    ##________ ADD EXPLOSIVE RESPONSE CONTROL HERE
    # War jetzt so schön vektorisiert... vielleicht einfach nachträglich noch die fehlenden mehr berechnen: mit flag: explosive
    # exp = TRUE/FALSE, dann zählen --> nachberechnen mit while? oder einfach 10% mehr berechnen und Warning falls mehr als 10% explosive waren?

    ####___________
    GIRF_low                  <- cbind(0:(n.ahead), apply(girfs_low, c(1,2), mean))
    Fresp_low_m               <- cbind(0:(n.ahead), apply(Fresp_low, 1, mean))
    dimnames(GIRF_low)        <- list(c(paste("h", 0:n.ahead, sep = "")),   c("V1",colnames(x$daten)[-1]))
    GIRF_high                 <-  cbind(0:(n.ahead), apply(girfs_high, c(1,2), mean))
    Fresp_high_m              <- cbind(0:(n.ahead), apply(Fresp_high, 1, mean))
    dimnames(GIRF_high)       <- list(c(paste("h", 0:n.ahead, sep = "")),   c("V1",colnames(x$daten)[-1]))
    GIRF_diff                 <- cbind(0:(n.ahead), apply(girfs_diff, c(1,2), mean))
    dimnames(GIRF_diff)       <- list(c(paste("h", 0:n.ahead, sep = "")),   c("V1",colnames(x$daten)[-1]))
    out$irflow                <- list(irf = as.data.frame(GIRF_low))
    out$irfhigh               <- list(irf = as.data.frame(GIRF_high))
    out$irfdiff               <- list(irf = as.data.frame(GIRF_diff))
    out$Fresp_low             <- as.data.frame(Fresp_low_m)
    out$Fresp_high            <- as.data.frame(Fresp_high_m)
    class(out)                <- "stvarirf"
    ##_____________
  } else if(is.null(period)==FALSE){
    pind                     <- which(F_Z_t %in% window(F_Z_t, start = period[[1]], end = period[[2]]))
    hist_p                   <- histories[pind]
    ##________ ADD EXPLOSIVE RESPONSE CONTROL HERE

    hist_p_sample            <- hist_p[sample(length(hist_p), size=iter2, replace=T)]
    ### inner function call here
    temp                     <- lapply(hist_p_sample, function(x) hgirf(x, iter1 = iter1))
    girfs_p                  <- sapply(temp, function(x) x$girf, simplify = "array")
    Fresp                    <- sapply(temp, function(x) x$Fresp, simplify = "array")
    GIRF_p                   <- cbind(0:(n.ahead), apply(girfs_p, c(1,2), mean))
    Fresp_m                  <- cbind(0:(n.ahead), apply(Fresp, 1, mean))
    dimnames(GIRF_p)         <- list(c(paste("h", 0:n.ahead, sep = "")),  c("V1", colnames(x$daten)[-1]))
    out$irflow               <- list(irf = as.data.frame(GIRF_p))
    out$Fresp                <- list(Fresp = as.data.frame(Fresp_m))
    out$method               <- "girf_single"
    class(out)               <- "girfsingle"

    ##_____________

  } else {
    ### calc GIRF for all histories
    temp                     <- lapply(histories, function(x) hgirf(x, iter1 = 2*iter1))
    girf_time                <- sapply(temp, function(x) x$girf, simplify = "array")   ### add cbind(0:(n.ahead), ...
    dimnames(girf_time)      <- list(c(paste("h", 0:n.ahead, sep = "")),   colnames(x$daten)[-1], time(Zts))
    Fresp                    <- sapply(temp, function(x) x$Fresp, simplify = "array")
    out$irf                  <- girf_time
    out$Fresp                <- as.data.frame(Fresp)
    out$method               <- "girf_time"
    class(out)               <- "irf3D"

  }

  #### Output addition?

  return(out)

}
#***** Notes on implementation **********#
## ->  The stboot object can be used to get Ystar (the fixed design bs y, the only thing that changed), u the bs redf. errors, B = bs struc. impact mat
##     and also A_hat the extended regressor mat! Problematic are Z, but could be recovered from Ystar if endogenous...
##
## the TVar starts either at t=0 or at t=-1, hence the data info can be used also for the regressors...
## use which to get the rows that belong to spec. regime/quantile or to specific time period.


