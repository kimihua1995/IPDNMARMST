ad.rmst <- function(data,tau){
  nt <- length(table(data$trial))
  ntrt <- length(table(data$trt))

  rmst_m <- rmst_sd <- matrix(NA, nrow = nt, ncol = ntrt)
  for (j in 1:nt){
    datj <- data[data$trial == j,]
    for (i in 1:ntrt){
      datij <- datj[datj$trt == i,]
      if (nrow(datij) > 0){
        res <- my.rmst(datij$time, datij$status, tau)
        rmst_m[j,i] <- res$rmst
        rmst_sd[j,i] <- res$se
      }
    }
  }

  return(list(rmst = rmst_m, se = rmst_sd))
}
