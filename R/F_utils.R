one_hot_errormachine <- function(Z){
  #------------ Objectif ------------
    # Prend en vecteur de classfication en entrée pour le One-Hot-Encode
    # en prenant en compte l'erreur machine possible.

  #------------ Variables ------------
  n <- length(Z)
  K <- length(unique(Z))

  #------------ One Hot Encoding + erreur machine ------------
  mat <- matrix(.Machine$double.xmin,n,K)
  mat[cbind(1:n,Z)] <- 1-.Machine$double.xmin

  #------------ Output ------------
  return(mat)
}

diag_nulle<- function(A){
  # Permet de récupérer un tenseur avec un coef 0 sur la diagonale des deux premières dimensions
  V <- dim(A)[3]
  for(v in 1:V) {diag(A[,,v]) <- 0 }
  return(A)
}


trig_sup<-function(A,transp=F,diag=T){
  #------------ Objectif ------------
  # Permet de récupérer la matrice - où le tenseur - triangulaire supérieure

  #------------ Variables ------------
  # A : Matrice(.,.) ou array(.,.,V)

  #------------ Récupération triangulaire supérieure ------------
  if(length(dim(A))==3){
    V <- dim(A)[3]
    for(v in 1:V){
      tmp <- A[,,v]
      if(transp){tmp = t(tmp)}
      tmp[lower.tri(tmp,diag=diag)] <- 0
      A[,,v]<- tmp
    }
  } else {

  if(transp){A = t(A)}
  A[lower.tri(A,diag=T)] <- 0
}
  #------------ Output ------------
  return(A)
}

transpo<- function(A){
  #------------ Objectif ------------
  # Permet de transposer un tenseur sur les deux premières dimensions

  #------------ transposition ------------
  V <- dim(A)[3]
  for(v in 1:V) {A[,,v] <- t(A[,,v]) }
  #------------ Output ------------
  return(A)
}

sort_Z <- function(Z){
  #------------ Objectif ------------
  # Permet de réordonner un vecteur de labels (one-hot encoded ou non)

  #------------ Variables ------------
  # Z : vecteur de labels
  K = ncol(Z)

  #------------ Réordonnement ------------
  for(k in K:1){
    Z <- Z[order(Z[,k]),] # Réordonne les données
  }
  #------------ Output ------------
  return(Z)
}

CEM <- function(Z){
  #------------ Objectif ------------
  # transforme une matrice de label one-hot encoded en vecteur de labels

  #------------ Variable ------------
  n <- nrow(Z)

  #------------ Classfication EM ------------
  for(i in 1:n){
    tmp <- which.max(Z[i,])
    Z[i,tmp] = 1
    Z[i,-tmp] = 0
  }

  #------------ Output ------------
  return(Z)
}

log_Softmax <- function(log_X){
  K <- ncol(log_X)

  log_X <- log_X - apply(log_X,1,max)

  ## Now going back to exponential with the same normalization
  X <- exp(log_X) #(matrix(1,n,1) %*% pi) * exp(logX)
  X <- pmin(X,.Machine$double.xmax)
  X <- pmax(X,.Machine$double.xmin)
  X <- X / (rowSums(X) %*% matrix(1,1,K))
  X <- pmin(X,1-.Machine$double.xmin)
  X <- pmax(X,.Machine$double.xmin)

  return(X)
}

