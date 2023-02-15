rMSBM <- function(N, V,alpha_klq,pi_k,rho,sorted=T){
  #------------ Objectif ------------
  # Simuler des données à partir du modèle génératif MixtureSBM.

  #------------ Variables ------------
  # N : Nomber d'individus
  # K : Nombre de  clusters
  # V : Nombre de vues
  # Q : Nombre de composantes du mélange de vue

  # alpha_klq : array(K,K,Q) Tenseur de probabilité de lien
  # pi_k : vecteur(K) Probabilité d'appartenir aux clusters
  # rho : vecteur(Q) Probabilité d'appartenir au mélange de vue Q

  # Z : Variable latente individus -> appartenance aux clusters
  # W : Variable latente vues -> d'appartenance aux composantes de vues

  K = length(pi_k)
  Q = length(rho)
  A = rep(NA,N*N*V)
  attr(A,"dim") <- c(N,N,V)
  #------------ états cachés clustering ------------
  Z <- t(rmultinom(n=N, size=1, prob=pi_k))
  Z = apply(Z,1,which.max)
  if(sorted) {Z = Z[order(Z)]}

  #------------ mélange vues ------------
  W <- t(rmultinom(n=V, size=1, prob=rho))
  W <- apply(W,1,which.max)
  if(sorted) {W = W[order(W)]}

  #------------ Tenseur d'adjacence ------------

  for (v in 1:V){
    for(i in 1:(N-1)){
      for(j in (i+1):N){
        A[i,j,v] <- rbinom(1,1,prob=alpha_klq[Z[i],Z[j],W[v]])
        A[j,i,v] <- A[i,j,v]
      }
    }
    diag(A[,,v])<- 1
  }

  #------------ output ------------

  params <- list(N=N,K=K,V=V,Q=Q,pi_k=pi_k,rho=rho,alpha_klq=alpha_klq)
  simulation <- list(A=A,Z=Z,W=W)
  output <- list(params = params, simulation=simulation)
  return(output)
}
