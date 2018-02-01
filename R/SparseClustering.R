# this file contains the hidden functions for the sparcl package
 
GetWCSS <- function(x, Cs, ws=NULL){
  wcss.perfeature <- numeric(ncol(x))
  for(k in unique(Cs)){
    whichers <- (Cs==k)
    if(sum(whichers)>1) wcss.perfeature <- wcss.perfeature + apply(scale(x[whichers,],center=TRUE, scale=FALSE)^2, 2, sum)
  }
  bcss.perfeature <- apply(scale(x, center=TRUE, scale=FALSE)^2, 2, sum)-wcss.perfeature
  if(!is.null(ws)) return(list(wcss.perfeature=wcss.perfeature, wcss=sum(wcss.perfeature), wcss.ws=sum(wcss.perfeature*ws),
                               bcss.perfeature=bcss.perfeature))
  if(is.null(ws)) return(list(wcss.perfeature=wcss.perfeature, wcss=sum(wcss.perfeature), bcss.perfeature=bcss.perfeature))
}
 
UpdateCs <- function(x, K, ws, Cs){
  x <- x[,ws!=0]
  z <- sweep(x, 2, sqrt(ws[ws!=0]), "*")
  nrowz <- nrow(z)
  mus <- NULL
  if(!is.null(Cs)){
    for(k in unique(Cs)){
      if(sum(Cs==k)>1) mus <- rbind(mus, apply(z[Cs==k,],2,mean))
      if(sum(Cs==k)==1) mus <- rbind(mus, z[Cs==k,])
    }
  }
  if(is.null(mus)){
    km <- kmeans(z, centers=K, nstart=10)
  } else {
    distmat <- as.matrix(dist(rbind(z, mus)))[1:nrowz, (nrowz+1):(nrowz+K)]
    nearest <- apply(distmat, 1, which.min)
    if(length(unique(nearest))==K){
      km <- kmeans(z, centers=mus)
    } else {
      km <- kmeans(z, centers=K, nstart=10)
    }
  }
  return(km$cluster)
}


BinarySearch <- function(numj, demj){
  if(l2n(numj)==0 || l2n(numj/demj)<=1) return(0)
  lam1 <- 0
  lam2 <- max(abs(numj)*sqrt(length(numj))) 
  iter <- 1
  while(iter<=32 && (lam2-lam1)>(1e-4)){
	amid <- (lam1+lam2)/2
    if(l2n(numj/(demj + amid))<=1){
      lam2 <- amid
    } else {
      lam1 <- amid
    }
    iter <- iter+1	
  }
  return((lam1+lam2)/2)
}



soft <- function(x,d){
  return(sign(x)*pmax(0, abs(x)-d))
}

l2n <- function(vec){
  return(sqrt(sum(vec^2)))
}


 
UpdateWs <- function(x, Cs, lambda1, lambda2, E, w0 = NULL){
  wcss.perfeature <- GetWCSS(x, Cs)$wcss.perfeature
  tss.perfeature <- GetWCSS(x, rep(1, nrow(x)))$wcss.perfeature
  Rj <- -wcss.perfeature+tss.perfeature
  nj <- numeric(length(Rj))
  ntapply <- tapply(E,E,length)
  nj[as.numeric(names(ntapply))] <- ntapply
  
  index1 <- E[,1]
  index2 <- E[,2]
  
  ## initialize all values
  if(is.null(w0)){
	  w0 <- rep(1/sqrt(length(Rj)), length(Rj))
  }
  w <- w0
  w_rand <- rnorm(length(w))
  z <- matrix(0,nrow=nrow(E), ncol=ncol(E))  
  u <- matrix(0,nrow=nrow(E), ncol=ncol(E))  
  z[] <- w[E]
  
  rho <- 1

  rt <- 1 ## primal residual
  st <- 1 ## dual residual

  ## ADMM loop
  ## ADMM varying rho parameter
  tau <- 2
  mu <- 10
  
  while(rt > 1e-4 | st > 1e-4){
	  ## update for z
	  b1 <- w[index1] + w[index2] + u[,1] + u[,2]
	  c1 <- w[index1] - w[index2] + u[,1] - u[,2]
	  S <- soft(c1, 2*lambda2/rho)
	  z_old <- z
	  z[,1] <- (b1 + S)/2
	  z[,2] <- (b1 - S)/2
  
	  ## update for w
	  aj <- numeric(length(Rj))
	  atapply <- tapply(rho*(z - u),E,sum)
	  aj[as.numeric(names(atapply))] <- atapply
  
	  numj <- pmax(aj + Rj - lambda1, 0)
	  demj <- rho * nj
	  lam <- BinarySearch(numj, demj)
	  w <- numj/(demj + lam)
  
	  ## update for u
	  u[] <- u + w[E] - z
    
	  ## update rho
  	  r <- w[E] - z
	  rt <- l2n(r)
	  deltaz <- z - z_old
	  ztapply <- tapply(deltaz,E,sum)
	  st <- rho * l2n(ztapply)
	  if(F){
		  if(rt > mu * st){
			  rho <- rho * tau
		  } else if(mu * rt < st){
			  rho <- rho / tau	  	
		  }	  	
	  
		  cat("rt:",rt,"\n")
		  cat("st:",st,"\n")
		  cat("rho:",rho,"\n")
	  }
  }
  
  return(w)
}

