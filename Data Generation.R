dat.sim.aft <- function(nt,n1,n2,alphas,a,betas,px,sigma,t_trt){
  ntrt <- length(alphas)
  alphaj <- matrix(NA, nrow = nt, ncol = ntrt)
  for (k in 1:ntrt){
    alphaj[,k] <- rnorm(nt, mean = alphas[k], sd = a)
  }
  n <- runif(nt,n1,n2)
  dat_comb <- list()
  for (j in 1:nt){
    prob <- rep(0,ntrt)
    prob[t_trt[[j]]] <- 1
    trt0 <- sample(x = 1:ntrt, n[j], replace = T, prob = prob)
    trt <- matrix(0,nrow = n[j],ncol = ntrt)
    for (i in 1:n[j]){
      trt[i,trt0[i]] <- 1
    }
    
    Xj <- rbinom(n[j],1,px)
    D <- matrix(NA,nrow = n[j],ncol = ntrt)
    #W <- revd(n[j])
    W <- rnorm(n[j])
    for (k in 1:ntrt){
      D[,k] <- exp(alphaj[j,k] + betas[k]*Xj + sigma[k]*W)
    }
    timee <- apply(D*trt, 1, sum)
    timec <- rexp(n[j],0.15)
    time <- pmin(timee, timec)
    status <- as.numeric(timee <= timec)
    
    dat <- as.data.frame(cbind(time, status, timee, timec, trt0, trt, Xj))
    colnames(dat) <- c("time","status","timee","timec",
                       "trt","trt1","trt2","trt3","X")
    dat$trial <- j
    dat_comb[[j]] <- dat
  }
  data <- do.call("rbind", dat_comb)
  data$trial <- factor(data$trial)
  data$id <- 1:nrow(data)
  return(data)
}


alphas = c(0.5,1.5,1)
betas = c(0.3,0.5,0.7)
sigma = c(1,1.5,2)
a.list = c(.1,.3,.5)
tau = 4
px = 0.5
nt = 10
n1 = 300
n2 = 700
t_trt = c(rep(list(1:3),4), rep(list(c(1,2)),2), rep(list(c(1,3)),2), rep(list(c(2,3)),2))
a = 0.2
set.seed(212)


data <- dat.sim.aft(nt,n1,n2,alphas,a,betas,px,sigma,t_trt)


save(data, file = "sim_data.RData")

