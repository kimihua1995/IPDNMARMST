cen.weight <- function(data,ntrt,nt,t){
  weight <- id <- NULL
  for (j in 1:nt){
    data.j <- data[data$trial == j,]
    y <- pmin(data.j$time, t)
    data.j$status[y == t] <- 1
    for (k in 1:ntrt){
      yk <- y[data.j$trt == k]
      dk <- data.j$status[data.j$trt == k]
      idk <- data.j$id[data.j$trt == k]
      if (length(yk) > 0){
        fitck <- func_surv(yk, 1-dk, idk)
        weightk <- dk[order(yk)]/rep(pmax(fitck$surv,0.001), table(yk[order(yk)]))
        idk_w <- fitck$id
        weight <- c(weight, weightk)
        id <- c(id, idk_w)
      }
    }
  }
  return(data.frame(weight = weight,
                    id = id))
}
