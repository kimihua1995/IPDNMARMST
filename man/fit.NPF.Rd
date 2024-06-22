\name{fit.NPF}
\alias{fit.NPF}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Non-Parametric Frequentist Two-Stage Model
}
\description{
Fit the individual participant level network meta-analysis with time-to-event endpoint using the non-parametric frequentist twp-stage model. This model does not allow incorporating individual covariate data.
}
\usage{
fit.NPF(data,tau)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{data}{a data frame including "time": time-to-event, "status": censoring indicator; "trial": study id, "trt": factor treatment groups, and other covariates of interests.}
\item{tau}{truncation time for RMST.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A list object of class "mvmeta" representing the meta-analytical model fit. See more details in mvmeta::mvmeta().
}
\references{
Hua, K., Wang, X., & Hong, H. (2023). Network Meta-Analysis of Time-to-Event Endpoints with Individual Participant Data using Restricted Mean Survival Time Regression. arXiv preprint arXiv:2310.13162.

Weir, I. R., Tian, L., & Trinquart, L. (2021). Multivariate meta-analysis model for the difference in restricted mean survival times. Biostatistics, 22(1), 82-96.
}


\examples{
}
