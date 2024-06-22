fit.NPF <- function(data, tau){
  res <- ad.rmst(data,tau)
  rmst_m <- res$rmst
  rmst_sd <- res$se
  nt <- length(table(data$trial))
  ntrt <- length(table(data$trt))

  theta <- as.vector(t(rmst_m))
  S <- matrix(NA, nrow = nt, ncol = ntrt*(ntrt+1)/2)
  for (j in 1:nt){
    s=1
    for (i in 1:ntrt){
      for (l in i:ntrt){
        S[j,s] <- ifelse(i==l,rmst_sd[j,i]*rmst_sd[j,l],0)
        s <- s + 1
      }
    }
  }

  fit <- mvmeta(rmst_m, S = S,
                method = "reml", bscov = "unstr", control = list(maxiter=1000))

  return(fit)
}
