# -----------------------------------------------------------------------------------
# Required functions & packages
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

# -----------------------------------------------------------------------------------
# Toy Data
# -----------------------------------------------------------------------------------

# data must be in long (not wide) form
n <- 500
T.max <- 4
p.td <- 5
time <- rep(1:T.max, n)
id <- rep(1:n, rep(T.max, n))
x.td <- matrix(rnorm(n * T.max * p.td), nrow = n * T.max)
x.nm <- NULL; for (i in 1:p.td) x.nm <- c(x.nm, paste("V",i,sep = ""));
colnames(x.td) <- x.nm
A <- rbinom(n * T.max, 1, .5)
Y <- rnorm(mean=1, n * T.max)
Rt <- rep(1, rep(T.max * n)) # R_2, ... , R_{T+1}

dat.long <- cbind(id, time, x.td, A, Y, Rt)

# add random dropout event
d.pr <- 0.5
id.dr <- rbinom(n,1,d.pr)
for (i in 1:n) {
  if (id.dr[i]){ # if dropout occurs
    dr.idx <- sample(1:(T.max-1), 1) + 1 # timepoint where to dropout
    dat.long[dat.long[,"id"]==i,"Rt"] <- c(rep(1, dr.idx - 1), rep(0, T.max - dr.idx + 1))
  }
}

# names of time-dependent covariates
td.cov.names <- x.nm


# -----------------------------------------------------------------------------------
#  Implementation of Algorithm 1
# -----------------------------------------------------------------------------------

# delta sequence
delta.seq <- c(seq(0.25, 1-0.05, by=0.05), seq(1, 3, by=0.1))

# sample splitting
nsplits = 2 # > 1

# mean value estimation
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

# Draw curve
plot.delta.curve(delta.seq, est.eff, eff.ll, eff.ul, eff.ll2, eff.ul2,
                 tp=T.max, outcome.name = "Y")


