fit.NPB <- function(data,tau,Omega,df,jag_model){
  rmst <- ad.rmst(data,tau)
  rmst_m <- c(t(rmst$rmst))
  rmst_sd <- c(t(rmst$se))
  nt <- nrow(rmst$rmst)
  ntrt <- ncol(rmst$se)

  trt <- rep(1:ntrt,nt); trt <- trt[!is.na(rmst_m)]
  trial <- rep(1:nt, each=ntrt); trial <- trial[!is.na(rmst_m)]
  rmst_m <- as.vector(na.omit(rmst_m))
  rmst_sd <- as.vector(na.omit(rmst_sd))
  dat_model <- list("y" = rmst_m, "s" = trial, "t" = trt,
                    "sd" = rmst_sd, "NT" = ntrt, "NS" = nt,
                    "N" = length(rmst_m), "zero.AB" = rep(0,ntrt),
                    "Omega" = Omega, "df" = df)
  par_model <- c("mu","tau")
  fit <- jags(data = dat_model, inits = NULL,
              parameters.to.save = par_model,
              n.iter=5000, n.burnin=1000, n.chains=2, n.thin=1,
              DIC = TRUE, jags.seed=212,
              model.file=jag_model,
              quiet = TRUE, progress.bar = "none")
  return(fit)
}
