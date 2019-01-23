#
#
# Variational Bayes for fitting a Gaussian Mixture Model
#
# a prior has to be specified which covers:
#   1) proportion of the components
#   2) a distribution over the mean vector of each component
#   3) a distribution over the covariance matrix of each component 
#


initGMM <- function(classes.numObs, means, means.numObs, Ws, Ws.numObs) {
  
  if (!all(length(classes.numObs)==
                  c(length(classes.numObs),ncol(means),length(means.numObs),length(Ws),length(Ws.numObs))))
    stop("length of arguments do not match")
  
  # initialization
  D <- nrow(means)
  chol.Ws <- lapply(Ws,chol)
  log.det.Ws <- sapply(chol.Ws,function(x) 2*sum(log(diag(x))))
  inv.Ws <- lapply(chol.Ws,chol2inv)
  
  train <- function(xobs,weights,nIter=50,randomize=FALSE,visualize="animation") {
    
    # effective number of observations
    numObs <- sum(weights^2) / weights^2
    
    # randomization assign prior means of not used component randomly to smpl points
    if (isTRUE(randomize)) {
      props <- classes.numObs / sum(classes.numObs)
      zeroComps <- which(props<0.001)
      means[,zeroComps] <<- xobs[,sample(ncol(xobs),length(zeroComps),replace=TRUE)]
      cat(paste0("randomize means of ",length(zeroComps)," components\n"))
    }
    
    # save the prior values because they are needed in each iteration
    classes.oldNumObs <- classes.numObs
    means.oldNumObs <- means.numObs
    oldMeans <- means
    oldWs <- Ws
    Ws.oldNumObs <- Ws.numObs
    # helpful quantities
    old.inv.Ws <- inv.Ws
    old.log.det.Ws <- log.det.Ws
    
    # start iterative procedure
    for (curIter in seq(nIter)) {
      
      ### E-step ###
      alpha.hat <- sum(classes.numObs)
      log.pi.tilde <- digamma(classes.numObs) - digamma(alpha.hat)
      log.Delta.tilde <- mapply(function(k) {
        sum(digamma((Ws.numObs[k]+1-(1:D))/2)) + D*log(2) + log.det.Ws[k]
      },k=seq_along(Ws))
      ## responsibilities
      log.resp <- mapply(function(k) {
        diff.xobs <- xobs - means[,k]
        log.pi.tilde[k] + 1/2*log.Delta.tilde[k] - D/(2*Ws.numObs[k]) - Ws.numObs[k]/2*colSums(diff.xobs*(Ws[[k]]%*%diff.xobs))  
      },k=seq_along(Ws))
      save.log.resp <- log.resp #debug
      log.resp <- apply(log.resp,1,function(x) x - max(x))
      resp <- exp(log.resp)
      resp <- apply(resp,2,function(x) x/sum(x))
      ## helpful quantities
      N <- rowSums(resp)
      if (any(is.na(N))) browser() #debug
      x.est <- mapply(function(k) {
        (if(N[k]>1e-10) 1/N[k] else 1)*colSums(resp[k,]*t(xobs))
      },k=seq_along(Ws))
      
      S.est <- mapply(function(k) {
        diff.obs <- xobs - x.est[,k]
        (if(N[k]>1e-10) 1/N[k] else 1)*diff.obs%*%(t(diff.obs)*resp[k,])
      },k=seq_along(Ws),SIMPLIFY=FALSE)
      
      ### M-step ###
      
      # update pseduobservations
      classes.numObs <<- classes.oldNumObs + N 
      means.numObs <<- means.oldNumObs + N
      Ws.numObs <<- Ws.oldNumObs + N
      
      # update components mean and Wishard-matrix (udpate in closure)
      means <<- t(1/means.numObs*(t(means)*means.oldNumObs + t(x.est)*N))
      
      inv.Ws <<- mapply(function(k) {
        diff.means <- x.est[,k] - oldMeans[,k]
        old.inv.Ws[[k]] + N[k]*S.est[[k]] + means.oldNumObs[k]*N[k]/means.numObs[k]*t(diff.means%*%t(diff.means)) #debug
      },k=seq_along(Ws),SIMPLIFY=FALSE)
      
      chol.inv.Ws <<- lapply(inv.Ws,chol)
      Ws <<- lapply(chol.inv.Ws,chol2inv)
      log.det.Ws <<- sapply(chol.inv.Ws,function(x) -2*sum(log(diag(x))))
      
      cat(paste0("Iteration: ",curIter," - Number of Components: ",sum(classes.numObs/sum(classes.numObs)>0.001),"\n"))
      # visualization
      if (visualize=="animation") {
        
        visualize(curIter, x.est,S.est, resp, xobs)
        Sys.sleep(1)
      }
    }
    if (visualize=="result")
      visualize(nIter, x.est,S.est, resp, xobs)
  }
  
  visualize <- function(curIter, x.est, S.est, resp, smpl) {
    # color according to class assignment
    classProp <- classes.numObs/sum(classes.numObs)
    plot(0,0,xlim=range(smpl[1,]),ylim=range(smpl[2,]),asp=1)
    colPal <- rep(c("red","green","blue","orange","magenta","cyan","brown","darkblue","black"),100)
    p.cidx <- apply(resp,2,which.max)
    for (i in seq_along(Ws)) {
      if (classProp[i]>0.001) {
        curSel <- which(p.cidx==i)
        points(smpl[1,curSel],smpl[2,curSel],col=colPal[i],cex=0.3)
        points(means[1,i],means[2,i],col=colPal[i],cex=1)
        points(x.est[1,i],x.est[2,i],col="black",cex=2,pch=7)
        lines(ellipse(S.est[[i]],centre=x.est[,i],level=0.68),col="black",lwd=2)
        lines(ellipse(inv.Ws[[i]]/Ws.numObs[i],centre=means[,i],level=0.68),col="blue",lwd=2)
      }
    }
  }
  
  printStatus <- function() {
    
    cat("### CLASS ASSIGNMENTS ###\n")
    print(classes.numObs/sum(classes.numObs))
    cat("### PSEUDO OBSERVATIONS oF COV ###\n")
    print(Ws.numObs)
    cat("### INVERSE WISHART MATRICES ###\n")
    print(mapply(function(k) inv.Ws[[k]]/Ws.numObs[k],k=seq_along(Ws),SIMPLIFY = FALSE))
  }
  
  # calculate predictive distribution
  getProbDens <- function(xmat) {
    
    require(mvtnorm)
    probDens <- matrix(0,ncol(xmat),length(Ws))
    alpha.hat <- sum(classes.numObs)
    
    for (k in seq_along(Ws))
    {
      scaleMat <- (1+means.numObs[k]) / ((Ws.numObs[k]+1-D)*means.numObs[k]) * inv.Ws[[k]]
      df <- Ws.numObs[k]+1-D
      probDens[,k] <- log(classes.numObs[k]) - log(alpha.hat) + dmvt(t(xmat),delta=means[,k],sigma = scaleMat, df=df,log=TRUE)
    }
    # rescale
    maxLog <- apply(probDens,1,max)
    probDens <- probDens - maxLog
    probDens <- log(rowSums(exp(probDens))) + maxLog
    return(probDens)
  }
  
  # sample probdens
  getSample <- function(num) {
    
    prop <- classes.numObs / sum(classes.numObs)
    compChain <- sample(seq_along(Ws),num,replace=TRUE,prob=prop)
    smpl <- matrix(0,num,D)
    
    for (k in seq_along(Ws)) {
      idx <- which(compChain==k)
      if (length(idx)>0) {
        # inefficient, could be moved into scope of closure
        scaleMat <- (1+means.numObs[k]) / ((Ws.numObs[k]+1-D)*means.numObs[k]) * inv.Ws[[k]]
        df <- Ws.numObs[k]+1-D
        # calculate
        smpl[idx,] <- rmvt(length(idx),delta=means[,k],sigma=scaleMat,df=df)
      }
    }
    # transpose the sample
    return(t(smpl))
  }
  
  # getter functions
  getProps <- function() { classes.numObs/sum(classes.numObs) }
  getMeans <- function() { means }
  getWs <- function() { Ws }
  getCovs <- function() { mapply(function(k) inv.Ws[[k]]/Ws.numObs[[k]],k=seq_along(Ws),SIMPLIFY=FALSE) }
  getWs.Obs <- function() { Ws.numObs }
  getMeans.Obs <- function() { means.numObs }
  getProps.Obs <- function() { classes.numObs }
  
  # export public functions
  publicFuns <- list(
    getProps=getProps, getMeans=getMeans, getWs=getWs,
    getCovs=getCovs, getWs.Obs=getWs.Obs, getMeans.Obs=getMeans.Obs,
    getProps.Obs=getProps.Obs,
    getProbDens=getProbDens,
    getSample=getSample,
    train=train
  )
  
  return(publicFuns)
}












