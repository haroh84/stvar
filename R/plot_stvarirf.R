####################################
###  stvar IRF plot functions:  ####

plot.stvarirf <- function(x, regime = c("first","second" , "both"),  scales = "free_y",...){

  strtonum         <- c(1, 2, 3)
  names(strtonum)  <-  c("first","second" , "both")
  regime <- match.arg(regime)
  regime <- strtonum[regime]

  multplot <- list(...)

  impulse <- list()

  impulse[[1]]   <- reshape2::melt(x$irflow$irf, id = 'V1')
  impulse[[2]]   <- reshape2::melt(x$irfhigh$irf, id = 'V1')

  if(regime %in% c(1,2)){
    gg <-  ggplot(impulse[[regime]], aes_(x = ~V1, y = ~value)) + geom_line() + geom_hline(yintercept = 0, color = 'blue') +
      facet_wrap(~variable, scales = scales, labeller = label_parsed) +
      xlab("Observation Time") + ylab("Response") +
      theme_bw()
  } else if(regime == 3){
    gg <-  ggplot(impulse[[1]], aes_(x = ~V1, y = ~value)) + geom_line() + geom_hline(yintercept = 0, color = 'blue') +
      facet_wrap(~variable, scales = scales, labeller = label_parsed) +
      xlab("Observation Time") + ylab("Response") +
      theme_bw()
    gg <- gg + geom_line(data = impulse[[2]], aes(x = V1, y = value), col = "red")
  } else {

    stop("Specify a regime to be ploted! Choose first, second or both.")
  }

  Nplot <- length(multplot)

  if(Nplot!=0){
    impadd <- list()
    if (regime == 3){
      for(i in 1:Nplot){
        impadd[[i]][[1]] <- melt(multplot[[i]][[1]]$irf)
        impadd[[i]][[2]] <- melt(multplot[[i]][[2]]$irf)

        gg <- gg +  geom_line(impadd[[i]][[1]], aes_(x = ~V1, y = ~value))
        gg <- gg +  geom_line(impadd[[i]][[2]], aes_(x = ~V1, y = ~value))
      }
    } else{
      for(i in 1:Nplot){
        impadd[[i]] <- reshape2::melt(multplot[[i]][[regime]]$irf)
        gg <- gg + geom_line(impadd[[i]], aes_(x = ~V1, y = ~value))
      }
    }
  }

  print(gg)

}






