fit.adj.one <- function(data,tau,covX){
  nt <- length(table(data$trial))
  ntrt <- length(table(data$trt))
  weights <- cen.weight(data,ntrt=ntrt,nt=nt,t=tau)
  data_w <- merge(data, weights, by = "id")
  data_w$Y <- pmin(data_w$time,tau)
  data_fit <- data_w[data_w$weight != 0,]
  data_fit$trt <- factor(data_fit$trt)
  fit <- MASS::glmmPQL(fixed = as.formula(paste0("Y ~ -1 + ", covX)),
                       #fixed = Y ~ -1 + trt1 + trt2 + trt3 + trt1:X + trt2:X + trt3:X,
                       random = as.formula(paste0("~ -1 + ", covX, " | trial")),
                       #random = ~ -1 + trt1 + trt2 + trt3 + trt1:X + trt2:X + trt3:X | trial,
                       family = quasipoisson(),
                       data = data_fit,
                       weights = data_fit$weight,
                       control = list(opt = "optim"), niter = 100, verbose = F)


  return(fit)
}
