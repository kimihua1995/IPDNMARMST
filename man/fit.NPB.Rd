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
}
