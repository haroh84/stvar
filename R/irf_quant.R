#############################################################################################
#### wrapper function that produced two AR matrices according to the specified quantiles ####
##-----------------------------------------------------------------------------------------##


irf_quant <- function(x, redf, quantiles, h = 20){

  ### wird dies benötigt? eigentlich nicht...
  u <- Tob <- p <- k <- residY <- coef_x <- yOut <- type <- y <-  NULL
  get_var_objects(redf)


  ## calculate quantile AR matrices:
  qlow      <- quantile(redf$Z, prob = quantiles$low)
  tf_low    <- exp(-exp(redf$theta)*(qlow - redf$c))/(1 + exp(-exp(redf$theta)*(qlow - redf$c)))
  qhigh     <- quantile(redf$Z, prob = quantiles$high )
  tf_high   <- exp(-exp(redf$theta)*(qhigh - redf$c))/(1 + exp(-exp(redf$theta)*(qhigh - redf$c)))

  ### use redform obj. from stvar_fun to get needed info on det. term etc.
  if(!any((class(x) %in% c("stvar", "sbootst" )))){
    stop("\nPlease provide an object of class 'stvar' or 'stbootst'.\n")
  }

  ir.fun <- function(y){
    # if(inherits(y, "stvar")){

    A_low  <- y$A_hat[,2:(1+k*p)]     ### adjust later for other det. term specifications (e.g., non-switching intercept)
    A_high <- y$A_hat[,c((k*p + 3):(2*k*p+2))]

    #}
    #else {
    # A_low  <- redf$A_low     ### adjust later for other det. term specifications (e.g., non-switching intercept)
    #  A_high <- redf$A_high
    #}
    Amat_low  <- (1-tf_low)*A_high + tf_low*A_low
    Amat_high  <- (1-tf_high)*A_high + tf_high*A_low

    ### create two objects of class "svars"

    ###  IMPORTANT: use A_l and A_h for the claculation of IRFs and change type to = none
    obj_low             <- y
    obj_low$A_hat       <- Amat_low
    obj_low$type        <- "none"
    obj_high            <- y
    obj_high$A_hat      <- Amat_high
    obj_high$type       <- "none"

    ### calculate IRFs:
    ir_low             <- irf(obj_low,  n.ahead = h)
    ir_low$tf          <- tf_low
    ir_high            <- irf(obj_high, n.ahead = h)
    ir_high$tf         <- tf_high
    class(ir_low)[2]   <- "stvarirf"   ### sollte hier die klasse svarirf entfernt werden also nur "class(x) <- "stvarirf" "?
    class(ir_high)[2]  <- "stvarirf"

    out <- list(
      irflow    =  ir_low,
      irfhigh   =  ir_high
    )

    return(out)

  }


  # (1 - F_Z_t) ## wurde vertauscht: wenn z_t -> infinity --> F_zt -> 0, daher ist (1-F_z_t) HIGH

  if(inherits(x, "sbootst")){  ### instead of confibands = easier
    newx <- x
    newx$bootstrap <- lapply(x$Allbs, ir.fun)
    newx$true      <- ir.fun(x$Omodel)
    newx           <- newx[names(newx) != "Allbs"]
    out            <- newx
    class(out)     <- "sbootst"
  } else {
    out <- ir.fun(x)
    class(out) <- "stvarirf"
  }


  return(out)
}
