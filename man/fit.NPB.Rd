\name{fit.NPB}
\alias{fit.NPB}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Non-Parametric Bayesian Two-Stage Model
}
\description{
Fit the individual participant level network meta-analysis with time-to-event endpoint using the non-parametric Bayesian two-stage model. This model does not allow incorporating individual covariate data.
}
\usage{
fit.NPB(data,tau,Omega,df,jag_model)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{data}{a data frame including "time": time-to-event, "status": censoring indicator; "trial": study id, "trt": factor treatment groups, and other covariates of interests.}
\item{tau}{truncation time for RMST.}
\item{Omega}{scale matrix in the Wishart prior of inverse variance-covariance matrix for the random effects.}
\item{df}{degrees of freedom in the Wishart prior of inverse variance-covariance matrix for the random effects.}
\item{jag_model}{jag model}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A object of JAGS model. See more details in R2jags::jags().
}
\references{
Hua, K., Wang, X., & Hong, H. (2023). Network Meta-Analysis of Time-to-Event Endpoints with Individual Participant Data using Restricted Mean Survival Time Regression. arXiv preprint arXiv:2310.13162.

Tang, X., & Trinquart, L. (2022). Bayesian multivariate network meta‐analysis model for the difference in restricted mean survival times. Statistics in medicine, 41(3), 595-611.
}


\examples{
df=5
Omega=diag(1,3)
jag_model <- function(){
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
}
