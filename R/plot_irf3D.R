####################################
###  stvar IRF plot functions:  ####
#
### 3D plot of generalized IRFs over time
### .....



plot.irf3D <- function(y, vars = "all",phi = 25, theta = -48 ,...){
  ## vars is either all, which means that all responses to a specific shock are ploted. or a vector with entries denoting the variables to
  # be ploted!
  require("plot3D")

  impulse <- y$irf
  mat3d   <- list()
  M       <- list()

  old.par  <- par()
  if(vars == "all"){
    nvars    <- dim(impulse)[2]
    varsplo  <- c(1:nvars)
  } else {
    nvars   <- length(vars)
    varsplo <- vars
  }
  windmat  <- c(ceiling(nvars/2), ceiling(ifelse(nvars<= 2, 2, nvars/2)))
  varnames <- colnames(impulse)

  par(mar = c(0.5, 0.5, 0.5, 0.5), mfrow = windmat)

  for( i in varsplo){
    mat3d[[i]] <- as.data.frame(cbind(0:(dim(impulse)[1]-1), impulse[,i,] ))
    names(mat3d[[i]])[1] <- "V1"
    M[[i]] <- mesh(mat3d[[i]]$V1, as.numeric(colnames(mat3d[[i]])[-1]) )
    x <-  M[[i]]$x
    y <-  M[[i]]$y
    z <- as.matrix(mat3d[[i]][,-1])

    surf3D(x = x, y = y, z = z, colvar = z, col = alpha.col(gg.col(), alpha = 1), colkey = FALSE, facets = TRUE ,
           bty="b2", box = TRUE, lwd.grid = 0, phi = phi, theta = theta, xlab = "horizon", ylab = "time", zlab = "response", main = varnames[i],
           r = 10, d = 1, ticktype = "detailed", border = "black", nticks = 10, cex.axis = 0.7)

  }
  par(old.par)

}
