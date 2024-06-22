library(Matrix)
library(zoo)
library(MASS)
library(lme4)
library(nlme)
library(survival)
library(survminer)
library(R2jags)
library(MCMCvis)
library(mvmeta)
library(devtools)
#install_github("kimihua1995/IPDNMARMST")  # uncomment to install the R pacakge "IPDNMARMST"
library(IPDNMARMST)


# load simulated data
data(sim_data)
head(data)


covX <- "trt1 + trt2 + trt3 + trt1:X + trt2:X + trt3:X"
tau <- 4


# 1) Proposed Two-Stage IPD-NMA RMST Model
fit_two <- fit.adj.two(data, tau, covX, quasi())
summary(fit_two)


# 2) Proposed One-Stage IPD-NMA RMST Model
fit_one <- fit.adj.one(data, tau, covX, quasi())
summary(fit_one)


# 3) Nonparametric Frequentist (NPF) Two-Stage Model
## Note: this approach cannot incorporate covariate, instead we fit the model on the subgroup X=0 and X=1.
fit_NPF_X0 <- fit.NPF(data[data$X==0,], tau)
fit_NPF_X1 <- fit.NPF(data[data$X==1,], tau)
summary(fit_NPF_X0)
summary(fit_NPF_X1)

# 4) Nonparametric Bayesian (NPB) Two-Stage Model
## Note: this approach cannot incorporate covariate, instead we fit the model on the subgroup X=0 and X=1.
df=5
Omega=diag(1,3)
rmst_ad_model <- function(){
  for (i in 1:N) {
    y[i] ~ dnorm(mean[s[i],t[i]], prec[i])
    prec[i] <- 1/pow(sd[i],2)
    mean[s[i],t[i]] <- mu[t[i]] + v[s[i],t[i]]
  }


  for (j in 1:NS) { v[j, 1:NT] ~ dmnorm(zero.AB[1:NT], invR[1:NT, 1:NT]) }

  invR[1:NT, 1:NT] ~ dwish(Omega[ , ],df)
  R[1:NT, 1:NT] <- inverse(invR[1:NT, 1:NT])


  for (k in 1:NT) {
    mu[k] ~ dnorm(0, 0.001)
    tau[k] <- sqrt(R[k,k])
  }
}

fit_NPB_X0 <- fit.NPB(data[data$X==0,], tau, Omega, df, rmst_ad_model)
fit_NPB_X1 <- fit.NPB(data[data$X==1,], tau, Omega, df, rmst_ad_model)

MCMCsummary(fit_NPB_X0, params = c("mu", "tau"), round = 3)
MCMCsummary(fit_NPB_X1, params = c("mu", "tau"), round = 3)
MCMCtrace(fit_NPB_X0, params = c("mu", "tau"),
         pdf = FALSE)
MCMCtrace(fit_NPB_X1, params = c("mu", "tau"),
          pdf = FALSE)









