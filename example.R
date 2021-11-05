# -----------------------------------------------------------------------------------
#  More complicated example using the data generation model in Sec 7.1
#  where dropout leads to an upward-biased estimate of the incremental effects
# -----------------------------------------------------------------------------------
library(ranger)
library(reshape2)
library(latex2exp)
library(SuperLearner)
library(KernelKnn)
library(kernlab)
library(npcausal)
suppressWarnings(suppressMessages(library("npcausal")))

source("inc_dropout_utils.R")

expit <- function(x){ exp(x)/(1+exp(x))} 
logit <- function(x){ log(x/(1-x)) }

# data generation model
set.seed(321) 
n = 1000
T.max = 2
c.xi.1 <- 1; c.xi.2 <- 16; y.c <- 2; r.c <- 2; a.c <- 8; 
B = 25
delta.seq <- c(seq(0.25, 1-0.05, by=0.05), seq(1, 3, by=0.1))
nsplits = 2 # > 1
ipsis <- ours <- NULL
r.d = 0
sig.3 = 0.25

x1.tm1 <- rnorm(n)
x2.tm1 <- rnorm(n)
B <- rbinom(n,1,0.5)
x3.tm1 <- 3*(2*B - 1)
a.tm1 <- rbinom(n, 1, .5)
R.t <- rep(1,n)
Y.t <- rnorm(n, 1)
dat <- cbind(1:n, rep(1,n), x1.tm1, x2.tm1, x3.tm1, a.tm1, Y.t, R.t)
for (t in 2:T.max) {
  xi.t <- ifelse(abs(x1.tm1+x2.tm1)<c.xi.1, (2*a.tm1 - 1)*abs(x1.tm1+x2.tm1)/c.xi.2, (2*a.tm1 - 1)/c.xi.2/c.xi.1)
  x1.t <- rnorm(x1.tm1,x1.tm1)  
  x2.t <- rnorm(x2.tm1*a.tm1,x2.tm1*a.tm1)
  x3.t <- x3.tm1 + rnorm(n,sig.3)
  a.t <- rbinom(n,1,1/2 + xi.t/a.c)
  y.mu.t <- (x1.tm1)/y.c  + a.t*x3.t + 1
  Y.t <- rnorm(y.mu.t, y.mu.t)
  x3.ind <- ifelse(x3.t>0,1,-1)
  p <- expit(r.d + 5*x3.ind + a.t/a.c + (x1.t + x2.t)/c.xi.2)
  Rt <- rbinom(n = n, size = 1, prob = p)
  x1.tm1 <- x1.t
  x2.tm1 <- x2.t
  x3.tm1 <- x3.t
  a.tm1 <- a.t
  dat <- rbind(dat, cbind(1:n, rep(t,n), x1.tm1, x2.tm1, x3.tm1, a.tm1, Y.t, Rt))
}

colnames(dat) <- c("id", "time", "V1", "V2", "V3", "A", "Y", "Rt")
dat.long <- as.data.frame(dat); 

# names of time-dependent covariates
td.cov.names <- c("V1", "V2", "V3")

# ours
estimation.result <- 
  estimation.sample.splitting(dat.long, td.cov.names, n, T.max, delta.seq, 
                              nsplits=2, unisance.est.md = "RF", verbose=FALSE)
kvals <- estimation.result$kvals
ifvals <- estimation.result$ifvals
est.eff <- colMeans(kvals)

# compute 95% CI and CB
interval.list <- variance.bootstrap(ifvals)
eff.ll <- interval.list$eff.ll
eff.ul <- interval.list$eff.ul
eff.ll2 <- interval.list$eff.ll2
eff.ul2 <- interval.list$eff.ul2

# dropout rate
(n-sum(dat.long[dat.long$time==T.max,"Rt"]==1))/n

# ipsi with full data
dat.long_ <- dat.long[order(dat.long$id, dat.long$time), ]
x.td <- dat.long_[,td.cov.names]
Y <- dat.long_[dat.long_$time==T.max, "Y"]
A <- dat.long_$A
time <- dat.long_$time
id <- dat.long_$id
ipsi.res <- ipsi(Y, A, x.td, x.td, time, id, delta.seq, nsplits=nsplits, fit='rf', progress_bar=FALSE)

# ipsi with partially observed data
idc.t <- dat.long[dat.long$time==T.max,"Rt"]==1
idc <- rep(idc.t, each=2)
ipsi.res.1 <- ipsi(Y[idc.t], A[idc], x.td[idc,], x.td[idc,], time[idc], id[idc], 
                   delta.seq, nsplits=nsplits, fit='rf', progress_bar=FALSE)

# Draw curves
col.pal <- c("turquoise4", "dimgrey", "tomato4", "black")
max.y <- 3.5; min.y <- -0.1;
plot.delta.curve(delta.seq, est.eff, eff.ll, eff.ul, eff.ll2, eff.ul2, max.y=max.y, min.y=min.y,
                 tp=T.max, outcome.name = "Y", est.eff.col=col.pal[1], est.eff.lwd=2, par.new=T) # dropout adjusted (ours)
plot(delta.seq, ipsi.res$res$est, type="l", lwd=1.5, col=col.pal[2],
     xlab=NA, ylab=NA, ylim=c(min.y, max.y)) # ipsi with full data
par(new=T)
plot(delta.seq, ipsi.res.1$res$est, type="l", lwd=1.5, col=col.pal[3],
     xlab=NA, ylab=NA, ylim=c(min.y, max.y)) # ipsi with partially observed data
par(new=T)
plot(delta.seq, rep(1,length(delta.seq)), type="l", lwd=1.5, col=col.pal[4], lty=2,
     xlab=NA, ylab=NA, ylim=c(min.y, max.y)) # true
par(new=F)
legend(x = "topleft",          
       legend = c("ipsi DO adj", "ipsi full", "ipsi partial", "true"),  
       lty = c(1, 1, 1, 2),     
       col = col.pal,           
       lwd = 2)                 