fit.adj.two <- function(data,tau,covX,family){
  nt <- length(table(data$trial))
  ntrt <- length(table(data$trt))
  weights <- cen.weight(data,ntrt=ntrt,nt=nt,t=tau)
  data_w <- merge(data, weights, by = "id")
  data_w$Y <- pmin(data_w$time,tau)
  data_w <- data_w[data_w$weight != 0,]
  nt <- length(table(data_w$trial))
  rmst_m <- S <- NULL
  for (j in 1:nt){
    dataj <- data_w[data_w$trial==j,]
    fit <- glm(formula = paste0("Y ~ -1 + ", covX),
               #Y ~ -1 + trt1 + trt2 + trt3 + trt1:X + trt2:X + trt3:X,
               family = family, weights = weight, data = dataj)
    rmst_m <- rbind(rmst_m, coef(fit))
    vcov_fit <- vcov(fit)
    S <- rbind(S, as.vector(vcov_fit[lower.tri(vcov_fit,diag = T)]))
  }


  fit <- mvmeta(rmst_m, S = S,
                method = "reml", bscov = "unstr", control = list(maxiter=1000))
  return(fit)
}
