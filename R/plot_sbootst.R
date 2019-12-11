####################################
###  stvar IRF plot functions:  ####
#
### plot st-var irfs with confidence bands
### supports plotting of multiple model's GIRFs




plot.sbootst <- function(x, regime = c("first","second" , "both"), modname = NULL, scales = "free_y", lowerq = 0.16, upperq = 0.84, percentile = 'standard', ...){

  strtonum         <- c(1, 2, 3)
  names(strtonum)  <-  c("first","second" , "both")
  regime <- match.arg(regime)
  regime <- strtonum[regime]

  irplot <- function(r){

    impulse <- reshape2::melt(x$true[[r]]$irf, id = 'V1')
    confidence <- x$bootstrap
    horizon <- nrow(confidence[[1]][[1]]$irf)
    kk <- ncol(confidence[[1]][[1]]$irf)
    nboot <- length(confidence)

    intervals <- array(0, c(horizon, kk, nboot))
    for(i in 1:nboot){
      intervals[,,i] <- as.matrix(confidence[[i]][[r]]$irf)
    }

    lower <- matrix(0, horizon, kk)
    for(i in 1:horizon){
      for(j in 1:kk){
        lower[i,j] <- quantile(intervals[i,j,], probs = lowerq)
      }
    }

    upper <- matrix(0, horizon, kk)
    for(i in 1:horizon){
      for(j in 1:kk){
        upper[i,j] <- quantile(intervals[i,j,], probs = upperq)
      }
    }
    lower <- as.data.frame(lower)
    upper <- as.data.frame(upper)
    if(percentile == 'hall'){
      lower <- 2*x$true[[r]]$irf - lower
      upper <- 2*x$true[[r]]$irf - upper
    }

    lower <- reshape2::melt(lower, id = 'V1')
    upper <- reshape2::melt(upper, id = 'V1')

    out <- cbind(impulse, lower = lower$value, upper = upper$value)

    return(out)

  }

  ##### call to irplot fun depending on the chosen regime:
  reg1   <- irplot(1)
  reg2   <- irplot(2)
  reg    <- list(reg1, reg2)
  combplo <- rbind(reg1, reg2)
  ####### change to one plot matrix to allow for legends:
  if(is.null(modname) == FALSE){
    Modind    <- rep(modname, each = nrow(reg1))
  } else {
    Modind      <- paste("Regime",rep(c("low", "high"), each = nrow(reg1)))
  }
  combplo       <- cbind(combplo, Modind)
  colnames(combplo)[ncol(combplo)] <- "Model"

  #####
  if(regime == 1 | regime == 2){

    gg <-  ggplot(reg[[regime]], aes_(x = ~V1, y = ~value)) +  geom_ribbon(aes(x = V1, ymin= lower, ymax= upper), alpha=.6, fill = "#74C2E1") +
      geom_line() + geom_hline(yintercept = 0, color = 'blue') +
      facet_wrap(~variable, scales = scales, labeller = label_parsed) +
      xlab("Observation Time") + ylab("Response") +
      theme_bw()
  } else{

    gg  <-  ggplot(combplo, aes_(x = ~V1, y = ~value, group = ~Model)) + geom_line(aes(linetype = Model, color = Model)) + geom_hline(yintercept = 0, color = 'blue') + geom_ribbon(aes(x = V1, ymin= lower, ymax= upper, fill = Model), alpha=0.6) +
      facet_wrap(~variable, scales = scales, labeller = label_parsed) +
      xlab("Observation Time") + ylab("Response") + scale_linetype_manual(values=c(1,2))+
      scale_color_manual(values=c("black", "red") ) + scale_fill_manual(values= c("grey", "#74C2E1")) +
      theme_bw()
  }
  return(gg)
  print(gg)

}
