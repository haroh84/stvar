####################################
###  stvar IRF plot functions:  ####
#
### plot st-var generalized irfs for specific threshold values of the transition variable
### supports plotting of multiple model's GIRFs


plot.girfsingle <- function(x, modname = NULL, scales = "free_y",...){

  cbp1 <- c("#E69F00", "#56B4E9", "#009E73",
            "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  lt   <- 2:6

  multplot <- list(...)


  Nplot <- length(multplot)

  if(Nplot!=0){
    impadd <- x$irflow$irf
    for(i in 1:Nplot){
      impadd      <- rbind(impadd, multplot[[i]]$irf)
    }
    if(is.null(modname) == FALSE){
      Modind    <- rep(modname, each = nrow(impadd)/(Nplot+1))
    } else {
      Modind      <- paste("M",rep(0:Nplot, each = nrow(impadd)/(Nplot+1)))
    }
    impadd                         <- cbind(impadd, Modind)
    colnames(impadd)[ncol(impadd)] <- "Model"
    impulse                        <- reshape2::melt(impadd, id = c("V1", "Model"))
    gg <-  ggplot(impulse, aes_(x = ~V1, y = ~value, group = ~Model)) + geom_line(aes(linetype = Model, color = Model)) + geom_hline(yintercept = 0, color = 'blue') +
      facet_wrap(~variable, scales = scales, labeller = label_parsed) +
      xlab("Observation Time") + ylab("Response") + scale_linetype_manual(values=c(1, lt[1:i]))+
      scale_color_manual(values=c("black", cbp1[1:i]) ) +
      theme_bw()
  } else {

    impulse      <- reshape2::melt(x$irflow$irf, id = 'V1')

    gg <-  ggplot(impulse, aes_(x = ~V1, y = ~value)) + geom_line() + geom_hline(yintercept = 0, color = 'blue') +
      facet_wrap(~variable, scales = scales, labeller = label_parsed) +
      xlab("Observation Time") + ylab("Response") +
      theme_bw()

  }

  print(gg)
  return(gg)
}





### function to extract the irf difference btw. generalized IRFs from trans. varibable threshold (above and below) or similar,
### such that it can be plotted by plot.girfsingle
get_irf_diff <- function(x){
  out <- list(irflow = x$irfdiff)
  class(out) <- "girfsingle"
  return(out)
}
