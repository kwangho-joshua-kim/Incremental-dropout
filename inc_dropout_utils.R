####################################################################################
#   A set of utility functions for
#    "Incremental Intervention Effects in Studies 
#     with Dropout and Many Timepoints." 
####################################################################################


####################################################################################

# Implement Algorithm 1
estimation.sample.splitting <- 
  function(dat.long, td.cov.names, n, ntimes, delta.seq, 
           delta.omega=0.01, nsplits=2, unisance.est.md = "RF", verbose=TRUE){
    
    # ARGUMENTS:
    #   dat.long [dataframe]: dataframe that contains subject info, outcome of interest, time-varying covariates & treatment, etc.
    #   n [int]: numher of the subject
    #   ntimes [int]: the timepoint (t) at which the incremental effect is estimated
    #   delta.seq [vector]: a sequence of delta values
    #   nsplits [int]: the number of splits (K)
    #   delta.omega [float]: lower bound for omega
    
    # RETURNS:
    #   ifvals: influence function estimates
    #   kvals: mean estimate on group k units
    #   + some caches
    
    # setup storage
    n.delta <- length(delta.seq)
    ifvals <- matrix(nrow=n, ncol=n.delta)
    kvals <- matrix(nrow=nsplits,ncol=n.delta)
    est.eff <- rep(NA,n.delta)
    
    W.long.delta = list()
    cumW.long = list()
    V.long.delta = list()
    M.long.delta = list()
    cumW.long.delta = list()
    
    s <- sample(rep(1:nsplits,ceiling(n/nsplits))[1:n])
    uniq_ids <- 1:n
    times.vec <- ntimes:1
    
    dat.f <- as.data.frame(dat.long)
    # sort data in ascending order of time & id
    dat.f <- dat.f[order(dat.f$time, dat.f$id), ]
    # R_1,...,R_T (not R_2,...,R_{T+1})
    Rt.prev <- c(rep(1,n),dat.f$Rt[1:(length(dat.f$Rt)-n)])
    
    ## Big loop started
    for (split in 1:nsplits){ 
      if (verbose) print(paste("split",split)); flush.console()
      id.split <- uniq_ids[s == split]
      
      ## Step 1 & 2: estimate pi and omega functions
      if (verbose) print("  fitting pi and omega functions.."); flush.console()
      trt.mat <- data.frame(id = matrix(uniq_ids, nrow=n, ncol=1))
      omega.mat <- data.frame(id = matrix(uniq_ids, nrow=n, ncol=1))
      
      for (t in 1:ntimes) {
        id.t <- obsvbleIdSet(t-1, dat.f)
        x.trt.t <- buildCovMat(t, id.t, dat.f, ntimes, td.cov.names=td.cov.names)
        R.next.t <- dat.f[dat.f$time==t & dat.f$id %in% id.t, c("id", "Rt")]
        x.trt.R.t <- merge(x.trt.t, R.next.t, by="id")
        
        trt.t.formula <- as.formula(paste(paste("A_",t,sep=""),"~.",sep = ""))
        omega.t.formula <- as.formula(paste("Rt","~.",sep = ""))
        if (nsplits == 1) {
          warning("nsplits should be greater than 1: naive Z-estimator used", call. = FALSE)
          trt.mod <- ranger(trt.t.formula, dat=x.trt.t[,-1])
          omega.mod <- ranger(omega.t.formula, dat=x.trt.R.t[,-1])
        } else {
          if (unisance.est.md == "RF"){
            trt.mod <- ranger(trt.t.formula, dat=x.trt.t[!(x.trt.t$id %in% id.split), -1] , write.forest = TRUE)
            omega.mod <- ranger(omega.t.formula, dat=x.trt.R.t[!(x.trt.R.t$id %in% id.split), -1] , write.forest = TRUE)
          } else if (unisance.est.md == "SuperLearner") {
            trt.mod <- SuperLearner(Y = x.trt.t[!(x.trt.t$id %in% id.split), c(paste("A_",t,sep=""))],
                                    X = x.trt.t[!(x.trt.t$id %in% id.split), !names(x.trt.t) %in% c(paste("A_",t,sep=""), "id")], 
                                    family = gaussian(), verbose = FALSE, cvControl = list(V=2), 
                                    SL.library = c("SL.ksvm", "SL.ranger", "SL.kernelKnn"))
            omega.mod <- SuperLearner(Y = x.trt.R.t[!(x.trt.R.t$id %in% id.split), c("Rt")],
                                      X = x.trt.R.t[!(x.trt.R.t$id %in% id.split), !names(x.trt.R.t) %in% c("Rt", "id")], 
                                      family = gaussian(), verbose = FALSE, cvControl = list(V=2), 
                                      SL.library = c("SL.ranger", "SL.kernelKnn"))
            
          } else {
            stop("other methods to be implemented")
          }
        }
        
        ps.t <- if (unisance.est.md == "RF"){ 
          predict(trt.mod, data=x.trt.t[,!names(x.trt.t) %in% c(paste("A_",t,sep=""), "id")])$predictions
        } else if (unisance.est.md == "SuperLearner") {
          as.vector(predict(trt.mod, x.trt.t[,!names(x.trt.t) %in% c(paste("A_",t,sep=""), "id")], onlySL = TRUE)$pred)
        } else {
          stop("other methods to be implemented")
        }
        
        omega.t <- if (unisance.est.md == "RF"){ 
          predict(omega.mod, data=x.trt.R.t[,!names(x.trt.R.t) %in% c("Rt", "id")])$predictions
        } else if (unisance.est.md == "SuperLearner") {
          as.vector(predict(omega.mod, x.trt.R.t[,!names(x.trt.R.t) %in% c("Rt", "id")], onlySL = TRUE)$pred)
        } else {
          stop("other methods to be implemented")
        }
        
        # set the lower bound for omega.t
        omega.t[which(omega.t<delta.omega)] <- delta.omega
        
        id.ps.t <- cbind(id = id.t, ps.t); colnames(id.ps.t)[2] <- t; #colnames(id.ps.t)[2] <- paste("ps_",t,sep="")
        id.omega.t <- cbind(id = id.t, omega.t); colnames(id.omega.t)[2] <- t; #colnames(id.omega.t)[2] <- paste("omega_",t,sep="")
        trt.mat = merge(trt.mat, id.ps.t, by="id", all.x = TRUE)
        omega.mat = merge(omega.mat, id.omega.t, by="id", all.x = TRUE)
      }
      
      trt.long <- melt(trt.mat, id.vars = "id")
      omega.long <- melt(omega.mat, id.vars = "id")
      
      ## Step 3 - 5
      for (j in 1:n.delta){
        if (verbose) print(paste("  delta",j)); flush.console()
        delta <- delta.seq[j]
        
        # Step 3:
        W.long <- omega.long
        W.long$value <- 
          (delta*dat.f$A + 1-dat.f$A)/(delta*trt.long$value + 1-trt.long$value) * Rt.prev/omega.long$value
        
        cumW.mat <- as.data.frame(aggregate(W.long$value, by=list(W.long$id), cumprod)[[2]])
        cumW.mat$id <- uniq_ids; colnames(cumW.mat)[1:ntimes] <- 1:ntimes;
        cumW.long <- melt(cumW.mat, id.vars = "id")
        
        # Step 4: fit (pseudo) outcome models
        M.mat <- data.frame(id = matrix(uniq_ids, nrow=n, ncol=1))
        M.mod <- vector("list", ntimes) 
        id.tp1 <- obsvbleIdSet(ntimes, dat.f)
        id.M.tp1 <- dat.f[dat.f$time==ntimes & dat.f$id %in% id.tp1, c("id", "Y")] 
        if (verbose) print("    fitting psuedo regressions.."); flush.console()
        
        for (t in times.vec){
          colnames(id.M.tp1)[2] <- "Mt"
          x.trt.t <- buildCovMat(t, id.tp1, dat.f,  ntimes, td.cov.names=td.cov.names)
          x.M.t <- merge(x.trt.t, id.M.tp1, by="id")
          if (nsplits == 1) {
            M.mod[[t]] <- ranger(Mt ~ ., dat=x.M.t[,-1], write.forest = TRUE)
          } else {
            if (unisance.est.md == "RF"){
              M.mod[[t]] <- ranger(Mt ~ ., dat=x.M.t[!(x.M.t$id %in% id.split), -1], write.forest = TRUE)
            } else if (unisance.est.md == "SuperLearner") {
              M.mod[[t]] <- SuperLearner(Y = x.M.t[!(x.M.t$id %in% id.split), c("Mt")],
                                         X = x.M.t[!(x.M.t$id %in% id.split), !names(x.M.t) %in% c("Mt", "id")], 
                                         family = gaussian(), verbose = FALSE, cvControl = list(V=2), 
                                         SL.library = c("SL.ksvm", "SL.ranger", "SL.kernelKnn"))
            } else {
              stop("other methods to be implemented")
            }
          }
          
          # id set for prediction could be larger than the time ahead
          id.t <- obsvbleIdSet(t-1, dat.f)
          x.trt.t.pred <- buildCovMat(t, id.t, dat.f,  ntimes, td.cov.names=td.cov.names)
          a.idx = paste(c("A", t), collapse = "_")
          newx0 <- newx1 <- x.trt.t.pred[,-1]
          newx1[,a.idx] <- 1 # (H_t, 1)
          m1 <- if (unisance.est.md == "RF") {
            predict(M.mod[[t]], data=newx1)$predictions
          } else if (unisance.est.md == "SuperLearner") {
            as.vector(predict(M.mod[[t]], newx1, onlySL = TRUE)$pred)
          } else {
            stop("other methods to be implemented")
          }
          newx0[,a.idx] <- 0  # (H_t, 0)
          m0 <- if (unisance.est.md == "RF") {
            predict(M.mod[[t]], data=newx0)$predictions
          } else if (unisance.est.md == "SuperLearner") {
            as.vector(predict(M.mod[[t]], newx0, onlySL = TRUE)$pred)
          } else {
            stop("other methods to be implemented")
          }
          
          pi.t <- trt.mat[trt.mat$id %in% id.t,t+1]
          # recursive regression formula in E.1
          M.tp1 <- (delta*m1*pi.t + m0*(1-pi.t))/(delta*pi.t + 1-pi.t) 
          id.M.tp1 <- cbind(id = id.t, M.tp1); colnames(id.M.tp1)[2] <- t; 
          
          # step 4-b
          omega.t <- omega.mat[omega.mat$id %in% id.t,t+1]
          A.t <- dat.f[dat.f$id %in% id.t & dat.f$time == t, c("A")]
          R.tp1 <- dat.f[dat.f$id %in% id.t & dat.f$time==t, c("Rt")]
          M.t <- ((m1 - m0)*delta*(A.t - pi.t)*omega.t/(delta*pi.t + 1-pi.t) +
            delta*m1*(pi.t*omega.t - A.t*R.tp1) + m0*((1-pi.t)*omega.t - (1-A.t)*R.tp1)
          ) / (delta*A.t + 1-A.t)
          id.M.t <- cbind(id = id.t, M.t); colnames(id.M.t)[2] <- t;

          M.mat = merge(M.mat, id.M.t, by="id", all.x = TRUE)
          
          id.tp1 <- id.t
        }
        
        M.mat <- M.mat[,c(1,ncol(M.mat):2)]
        M.long <- melt(M.mat, id.vars = "id")
        
        #j indexes the deltas
        W.long.delta[[j]] <- W.long
        cumW.long.delta[[j]] <- cumW.long #cumulative incremental PS
        M.long.delta[[j]] <- M.long #outcome model results
        
        # Step 6:
        phi.val <- ((cumW.long$value*dat.f$Y*dat.f$Rt)[dat.f$time==ntimes] +
                      aggregate(cumW.long$value*M.long$value, by=list(dat.f$id), sum)[[2]])[s==split]
        
        ifvals[s==split,j] <- phi.val #influence function values
        kvals[split,j] <- mean(phi.val, na.rm = T) #mean in the k group
      } 
    }
    
    return(list(ifvals=ifvals, kvals=kvals, s=s,
                W.long.delta=W.long.delta, cumW.long.delta=cumW.long.delta,
                M.long.delta=M.long.delta))
  }

####################################################################################

####################################################################################

# Function to create (\bar{X}_t, \bar{A}_t) = (\bar{H}_t, A_t)
buildCovMat <- 
  function(t, id.t, dat, maxT, td.cov.names = NULL, append = FALSE, no.last.trt = FALSE)
  {
    # DESCRYPTION:
    #    Returns cbind(\bar{H}_t, A_t) with "id" column
    #    : i.e., historical values of covariate and exposure 
    #    for the specified id set (id.t) up to time t.
    #    If no.last.trt = TRUE, then we will just get \bar{H}_t
    
    # ARGUMENTS:
    #    dat [dataframe]
    #    t [int] timepoint  
    #    id.t =  observable "id"'s (not dropped-out) at time t
    
    # RETURNS:
    #   HA.t.bar [dataframe] cbind(\bar{H}_t, A_t) with id column
    
    if (t<1) return(NULL)
    
    # time-dependent covariates
    if (is.null(td.cov.names)) {
      warning('No time-dependent covariates used')
    }
    dat.id.t <- dat[dat$id %in% id.t,]
    dat.id.t <- dat.id.t[order(dat.id.t$time, dat.id.t$id),]
    
    # define a local function that returns td cov values at time t_ (<t)
    td.time.t <- function(t_, dat.id.t, td.var.names, keep.id = T) {
      x.td.t = dat.id.t[dat.id.t$time == t_, c("id", td.var.names), drop = F]
      td.names = NULL
      for (nm in td.var.names) td.names = c(td.names, paste(c(nm, t_), collapse = "_"))
      colnames(x.td.t)[-1] <- td.names
      if (keep.id == T) {
        return(x.td.t)
      } else {
        return(x.td.t[,-1, drop = F])  
      }
    }
    
    # build (\bar{H}_t, A_t)
    td.var.names <- c(td.cov.names, "A")
    HA.td.t.bar <- td.time.t(1, dat.id.t, td.var.names)
    if (t>1) {
      for (k in 2:t) {
        k.temp.x <- td.time.t(k, dat.id.t, td.var.names)
        HA.td.t.bar <- merge(HA.td.t.bar, k.temp.x, by="id")
      }
    }
    
    if (append == T & t < maxT) { # force to repeat the final value until maxT 
      for (k in (t+1):maxT) {
        k.temp <- td.time.t(t, dat.id.t, td.var.names)
        HA.td.t.bar <- cbind(HA.td.t.bar, k.temp)
      }
    }
    
    # drop the last column if needed (so we just have H_t)
    if (no.last.trt) HA.td.t.bar <- HA.td.t.bar[,-ncol(HA.td.t.bar), drop = F]
    
    # in case of using baseline cov... TBD
    # HA.t.bar <- if(is.null(x.cov.bs)) HA.td.t.bar else merge(x.bs.t, HA.td.t.bar, by="id")
    return(HA.td.t.bar)
    
  }

####################################################################################
####################################################################################

# Function to identify observable ids at time t
obsvbleIdSet <- 
  function(t, dat)
  {
    # DESCRYPTION:
    #    Returns the vector of id's that have not dropped out at t
    #    (i.e., {id_t: R_t = 1})
    #    (Note that we assume R_1 (corresponding to t=0) = 1)
    
    # ARGUMENTS:
    #    dat[dataframe]
    
    if (t==0) {
      id.t <- unique(dat$id) # 1:n
    } else{
      id.t <- dat[dat$time == t & dat$Rt == 1, "id"]  
    }
    
    return(id.t)
    
  }

####################################################################################

# compute asymptotic variance via bootstrapping
variance.bootstrap <- function(ifvals, nbs = 5000) {
  # DESCRYPTION:
  #   compute both 95% pointwise confidence interval (CI) and confidence band (CB)
  
  # ARGUMENTS:
  #   ifvals: value of estimated influence function 
  
  # RETURN(S):
  # eff.ll: lower bound of 95% pointwise CI
  # eff.ul: upper bound of 95% pointwise CI
  # eff.ll2: lower bound of 95% CB
  # eff.ul2: upper bound of 95% CB
  
  # 95% pointwise asymptotic variance
  sigma <- sqrt(apply(ifvals,2,var,na.rm=TRUE))
  ifvals.0 <- na.omit(ifvals)
  n0 <- dim(ifvals.0)[1]
  eff.ll <- est.eff-1.96*sigma/sqrt(n0); eff.ul <- est.eff+1.96*sigma/sqrt(n0)
  
  # confidence band based on multiplier bootstrapping
  eff.mat <- matrix(rep(est.eff,n0), nrow=n0, byrow=T)
  sig.mat <- matrix(rep(sigma,n0), nrow=n0, byrow=T)
  ifvals2 <- (ifvals.0 - eff.mat)/sig.mat
  mult <- matrix(2*rbinom(n0*nbs,1,.5) - 1, nrow=n0, ncol=nbs)
  
  maxvals <- sapply(1:nbs, function(col){
    max(abs(apply(mult[,col]*ifvals2,2,sum)/sqrt(n0))) 
  } 
  )
  calpha <- quantile(maxvals, 0.95)
  eff.ll2 <- est.eff - calpha*sigma/sqrt(n0); eff.ul2 <- est.eff + calpha*sigma/sqrt(n0)
  
  return(list(eff.ll=eff.ll, eff.ul=eff.ul, eff.ll2=eff.ll2, eff.ul2=eff.ul2))
}

####################################################################################

# draw curve for estimated incremental effects with 95% pointwise CI and CB
plot.delta.curve <- 
  function(delta.seq, est.eff, eff.ll, eff.ul, eff.ll2, eff.ul2, 
           tp="", outcome.name="Y", max.y=NULL, min.y=NULL, 
           est.eff.col="dimgrey", est.eff.lwd=1.5,  par.new=F) {
    
    outcome.name <- simpleCap(outcome.name)
    if (is.null(max.y) | is.null(min.y)) {
      max.y <- max(eff.ul2)+0.1; min.y <- min(eff.ll2) - 0.1;
    }
    plot(delta.seq, eff.ll2, type="l", col="firebrick", lty=3, lwd=2, xlab = NA, ylab = NA, ylim=c(min.y, max.y)); par(new=T)
    plot(delta.seq, eff.ul2, type="l", col="firebrick", lty=3, lwd=2, xlab = NA, ylab = NA, ylim=c(min.y, max.y)); par(new=T)
    plot(delta.seq, eff.ll, type="n", lty=3, lwd=2, xlab = NA, ylab = NA, ylim=c(min.y, max.y)); par(new=T)
    plot(delta.seq, eff.ul, type="n", lty=3, lwd=2, xlab = NA, ylab = NA, ylim=c(min.y, max.y)); par(new=T)
    polygon(c(delta.seq, rev(delta.seq)), c(eff.ul2, rev(eff.ll2)),
            col = "gainsboro", border = NA); par(new=T)
    polygon(c(delta.seq, rev(delta.seq)), c(eff.ul, rev(eff.ll)),
            col = "gray70", border = NA); par(new=T)
    plot(delta.seq, est.eff, type="l", lwd=1.5, col=est.eff.col, 
         xlab=NA, ylab=NA, ylim=c(min.y, max.y)) 
    abline(v = 1, col="darkblue", lwd=1.5, lty=3)
    title(main=paste("Estimated Probability of ",outcome.name," (T=", tp, ") ", sep=""),
          xlab=expression(paste("odds ratio ",delta,sep="")), ylab=expression(Psi(delta)))
    par(new=par.new)
  }


simpleCap <- function(x) {
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1,1)), substring(s, 2),
        sep="", collapse=" ")
}