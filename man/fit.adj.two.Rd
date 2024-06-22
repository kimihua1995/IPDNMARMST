\name{fit.adj.two}
\alias{fit.adj.two}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Two-Stage RMST IPD-NMA Model
}
\description{
Fit the individual participant level network meta-analysis with time-to-event endpoint using the two-stage RMST IPD-NMA model. This model allows incorporating individual covariate data.
}
\usage{
fit.adj.two(data,tau,covX,family)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{data}{a data frame including "time": time-to-event, "status": censoring indicator; "trial": study id, "trt": factor treatment groups, and other covariates of interests.}
\item{tau}{truncation time for RMST.}
\item{covX}{a vector of variables in the formula, e.g., covX = "trt1 + trt2 + trt3 + trt1:X + trt2:X + trt3:X".}
\item{family}{a description of the link function to be used in the model. Use `quasi()` for identity link and `quasipoisson()` for log link.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A list object of class "mvmeta" representing the meta-analytical model fit. See more details in mvmeta::mvmeta().
}
\references{
Hua, K., Wang, X., & Hong, H. (2023). Network Meta-Analysis of Time-to-Event Endpoints with Individual Participant Data using Restricted Mean Survival Time Regression. arXiv preprint arXiv:2310.13162.
}


\examples{
}

