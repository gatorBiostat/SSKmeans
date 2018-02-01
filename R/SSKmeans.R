##' Smoothed sparse Kmeans
##'
##' Perform sparse Kmeans to perform sample clustering and feature selection. In feature selection, we also want to incorporate spatio information such that adjacent voxels have similar coefficient.
##' @title SSKmeans
##' @param x a data matrix of dimension n * p, where n is number of samples to be clustered, p is number of features.
##' @param K pre-specified number of clusters
##' @param lambda1 tuning parameter for l1 norm lasso penalty. Large lambda1 will induce more feature weights to be 0.
##' @param lambda1 tuning parameter for the smoothness of feature selection. Large lambda will induce feature weights to be similar.
##' @param E network structure. E should be a m by 2 matrix, where m is total number of connections (edges) in the graph. 
##' For example, if feature 2 and feature 3 are connected in the graph, E[j,] <- c(2, 3). j = 1, ..., m.
##' @param nstart number of initialization for Kmeans.
##' @param silent if print out progress.
##' @param maxiter max number of iterations.
##' @return a list. The following items are included in the list.
##' \item{ws}{weight for each feature. Zero weight means the feature is not selected.}
##' \item{Cs}{Cluster Assignment}
##' \item{wcss}{within cluster sum of square}
##' \item{crit}{objective value}
##' \item{E}{network structure}
##' @author Caleb
##' @export
##' @examples
##' set.seed(11)
##' x <- matrix(rnorm(50*70),ncol=70)
##' x[1:25,1:20] <- x[1:25,1:20]+1
##' x <- scale(x, TRUE, TRUE)
##' # choose tuning parameter
##' lambda1 <- 1
##' lambda2 <- 1
##' E <- cbind(1:69,2:70)
##' K=2
##' nstart=20
##' silent=FALSE
##' maxiter=6
##' 
##' km.out <- SSKmeans(x,K=2,lambda1=1, lambda2=1, E=E)
##' 
SSKmeans <-function(x, K=NULL, lambda1 = 0, lambda2 = 0, E = NULL, nstart=20, silent=FALSE, maxiter=20){

  Cs <- kmeans(x, centers=K, nstart=nstart)$cluster
  ws <- rep(1/sqrt(ncol(x)), ncol(x)) # Start with equal weights on each feature
  ws.old <- rnorm(ncol(x))
  store.bcss.ws <- NULL
  niter <- 0

  while((sum(abs(ws-ws.old))/sum(abs(ws.old)))>1e-4 && niter<maxiter){
    if(!silent) cat(niter, fill=FALSE)
    niter <- niter+1
    ws.old <- ws
    if(niter>1) Cs <- UpdateCs(x, K, ws, Cs) # if niter=1, no need to update!!
	ws <- UpdateWs(x, Cs, lambda1, lambda2, E, ws)
    store.bcss.ws <- c(store.bcss.ws, sum(GetWCSS(x, Cs)$bcss.perfeature*ws))
  }
  out <- list(ws=ws, Cs=Cs, wcss=GetWCSS(x, Cs, ws), crit=store.bcss.ws, E=E)

  if(!silent) cat(fill=TRUE)
  class(out) <- "SSKmeans"
  return(out)
}

