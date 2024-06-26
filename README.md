<h2>Individual participant level data (IPD) network meta-analysis (NMA) with time-to-event endpoint using restricted mean survival time (RMST) regression</h2>

1) Proposed methods (allow incorporating individual covariate data)
   - Two-stage IPD-NMA RMST model `fit.adj.two(data, tau, covX, family)`
   - One-stage IPD-NMA RMST model `fit.adj.one(data, tau, covX, family)`
2) Existing approaches (does not allow incorporating individual covariate)
   - Nonparametric frequentist (NPF) two-stage model `fit.NPF(data, tau)`
   - Nonparametric Bayesian (NPB) two-stage model `fit.NPB(data, tau, Omega, df, jag_model)`

<h3>Install</h3>

In your R console
<pre>
library(devtools)
install_github("kimihua1995/IPDNMARMST")
</pre>

<h3>Examples</h3>

see `Example.R`

<h3>Reference</h3>

Hua, K., Wang, X., & Hong, H. (2023). Network Meta-Analysis of Time-to-Event Endpoints with Individual Participant Data using Restricted Mean Survival Time Regression. arXiv preprint arXiv:2310.13162.

Weir, I. R., Tian, L., & Trinquart, L. (2021). Multivariate meta-analysis model for the difference in restricted mean survival times. Biostatistics, 22(1), 82-96.

Tang, X., & Trinquart, L. (2022). Bayesian multivariate network meta‐analysis model for the difference in restricted mean survival times. Statistics in medicine, 41(3), 595-611.
