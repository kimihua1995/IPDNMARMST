my.rmst=function(time, status, tau){
  ft <- survfit(Surv(time, status)~1)
  idx <- ft$time<=tau

  wk.time <- sort(c(ft$time[idx],min(tau,ft$time)))
  wk.surv <- ft$surv[idx]
  wk.n.risk  <- ft$n.risk[idx]
  wk.n.event <- ft$n.event[idx]

  time.diff <- diff(c(0, wk.time))
  areas <- time.diff * c(1, wk.surv)
  rmst <-  sum(areas)

  wk.var <- ifelse((wk.n.risk-wk.n.event)==0, 0,
                   wk.n.event /(wk.n.risk *(wk.n.risk - wk.n.event)))
  wk.var <- c(wk.var,0)
  rmst.var <- sum( cumsum(rev(areas[-1]))^2 * rev(wk.var)[-1])
  rmst.se  <- sqrt(rmst.var)

  return(list(rmst = rmst, se = rmst.se))
}
