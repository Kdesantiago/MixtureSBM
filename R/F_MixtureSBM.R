update_rho <- function(params){
  u<- params$u
  V = nrow(u)
  Q <- ncol(u)
  new_rho <- apply(u,2,sum)/V
  new_rho <- pmin(new_rho,1-.Machine$double.xmin)
  new_rho <- pmax(new_rho,.Machine$double.xmin)

  params$rho <-new_rho
  return(params)
}

update_pi <- function(params){
  tau<- params$tau
  N = nrow(tau)
  K = ncol(tau)
  new_pi <- apply(tau,2,sum)/N
  new_pi <- pmin(new_pi,1-.Machine$double.xmin)
  new_pi <- pmax(new_pi,.Machine$double.xmin)

  params$pi_k<-new_pi
  return(params)
}

update_u <- function(A,params){

  ##
  alpha <- params$alpha
  rho<- params$rho
  u<- params$u
  tau<- params$tau

  K = ncol(tau)
  V = dim(A)[3]
  Q <- ncol(u)

  log_u <- matrix(0,V,Q)

  X <- matrix(0,K,K)

  for (v in 1:V) {
    for (q in 1:Q) {
      X <- matrix(0,K,K)
      for(k in 1:K){
        for(l in 1:K){
          tmp <- A[,,v] * log(alpha[k,l,q] + .Machine$double.xmin) + (1-A[,,v]) * log(1-alpha[k,l,q] +.Machine$double.xmin)
          tmp[lower.tri(tmp,diag=T)] <- 0
          X[k,l] <- t(tau[,k]) %*% tmp %*% tau[,l]
        }
      }
      log_u[v,q] <- -log(rho[q]) + sum(X)
    }
  }



  u <- log_Softmax(log_u)

  params$u <- u
  return(params)
}

update_tau <- function(A,params,eps_conv=1e-3){

  alpha <- params$alpha
  pi_k<- params$pi_k
  u<- params$u
  Tau<- params$tau

  old_tau <- Tau + 1 # init random

  N = nrow(A)
  K = ncol(Tau)
  V = dim(A)[3]
  Q <- ncol(u)

  etp <- 0

  while( (sum(abs(old_tau - Tau)) > eps_conv) && etp<50 ){
    etp <- etp+1
    old_tau <- Tau
    log_tau <- matrix(0,N,K)

    for(i in 1:N){

      for (v in 1:V) {
        tmp <- A[,,v]
        diag(tmp) <- 0
        for (q in 1:Q) {
          for(k in 1:K){
            for(l in 1:K){
              log_tau[i,k] <- log_tau[i,k] + sum( u[v,q]* Tau[,l] *  ( tmp[i,] * log(alpha[k,l,q] + .Machine$double.xmin) + (1-tmp[i,]) * log(1-alpha[k,l,q] +.Machine$double.xmin)) )
            }
          }
        }
      }

      log_tau[i,k] <- log_tau[i,k] - log(pi_k[k])
    }



    Tau <- log_Softmax(log_tau)
  }

  params$tau <- Tau
  return(params)
}

update_alpha <-function(A,params){

  tau  <-   params$tau
  u <- params$u

  N <- nrow(A)
  V <- dim(A)[3]
  K <- ncol(tau)
  Q <- ncol(u)

  ones <- trig_sup(array(1,dim=c(N,N,V)))


  theta_ikjlvq <- tau %o% tau %o% u #  symétrie en k,l avec la transposée


  A <- trig_sup(A) + transpo(trig_sup(A,transp = T))

  tamp_alpha <- array(0, dim = c(K,K,Q))
  for(q in 1:Q){
    for(k in 1:K){
      for(l in k:K){
        tmp <- theta_ikjlvq[,k,,l,,q]
        tamp_alpha[k,l,q] <-  sum(tmp * A)/sum(tmp)
      }
    }

  }
  params$alpha <- tamp_alpha + apply(tamp_alpha,c(1,3),t) - array(apply(apply(tamp_alpha,3,diag),2,diag),dim=c(K,K,Q)) #estimation triang supérieur et on permute les coeffs
  return(params)
}


## Initialisation
initialisation_params <-function(A,K,Q,type_init="random"){
  N = nrow(A)
  V = dim(A)[3]
  params <- list()

  if(type_init=="Kmeans"){
    A_tmp <- apply(A,c(1,2),sum)
    tmp <- kmeans(A_tmp,centers = K,nstart = 50)

    params$pi_k <- tmp$size/sum(tmp$size)
    params$rho <- rep(1/Q, Q)
    params$tau <- one_hot_errormachine(tmp$cluster)
    params$u <- array(runif(N*K,0,1),dim=c(V,Q)) ; params$u <- params$u/apply(params$u,1,sum)
    params <- update_alpha(A,params)
    params <- update_u(A,params)
    params <- update_rho(params)
    params <- update_alpha(A,params)


  } else {
    params$pi_k <- rep(1/K, K)
    params$rho <- rep(1/Q, Q)
    params$tau <- array(runif(N*K,0,1),dim=c(N,K)) ; params$tau <- params$tau/apply(params$tau,1,sum)
    params$u <- array(runif(N*K,0,1),dim=c(V,Q)) ; params$u <- params$u/apply(params$u,1,sum)
    params <- update_alpha(A,params)
  }


  return(params)
}


## VEM

VEM_step<-function(A,params,alternate=T,eps_conv=1e-3){

  if(alternate){
    params <- update_u(A,params)
    params <- update_alpha(A,params)
    params <- update_rho(params)


    params <- update_tau(A,params,eps_conv)
    params <- update_alpha(A,params)
    params <- update_pi(params)


  } else{
    params <- update_u(A,params)
    params <- update_tau(A,params,eps_conv)
    params <- update_alpha(A,params)
    params <- update_pi(params)
    params <- update_rho(params)

  }


  return(params)
}


## ELBO

ELBO <- function(A,params){
  alpha <- params$alpha
  rho<- params$rho
  u<- params$u
  tau<- params$tau
  pi_k<- params$pi_k

  N = nrow(A)
  K = ncol(tau)
  V = dim(A)[3]
  Q <- ncol(u)

  for(v in 1:V){
    tmp <- A[,,v]
    tmp[lower.tri(tmp,diag=T)] <- 0
    A[,,v]<- tmp
  }
  theta_ikjlvq <- tau %o% tau %o% u

  res <- 0
  for(q in 1:Q){
    for(k in 1:K){
      for(l in k:K){
        tmp <- theta_ikjlvq[,k,,l,,q]
        tmp = trig_sup(tmp)
        res <- res + sum(tmp * (A*log(alpha[k,l,q] + .Machine$double.xmin) +(1-A) *log(1-alpha[k,l,q] + .Machine$double.xmin )))
      }
    }
  }





  #res <- res + sum(diag( t(tau) %*%(log(tau) +log(matrix(rep(pi_k,N),N,K,byrow=T))) ))


  #res + sum(diag( t(u) %*%(log(u * matrix(rep(rho,V),V,Q,byrow=T))) ))

  for(i in 1:(N)){
    for(k in 1:K){
      res <- res + tau[i,k] * ( log(pi_k[k]) +log(tau[i,k]) )
    }
  }
  for(v in 1:(V)){
    for(q in 1:Q){
      res <- res + u[v,q] * ( log(rho[q]) + log(u[v,q]) )
    }
  }

  return(res)
}

## Modèle



Mixture_SBM <-function(A,K,Q,tol,iter_max=10,n_init = 1,alternate=T, Verbose=T,eps_conv=1e-3,type_init="Kmeans"){

  #------------ Objectif ------------
  # Stochastic Block Model avec un mélange de vue, optimisé avec un Variation-EM.

  #------------ Variables ------------
  # A : array(N,N,V) - Tensor de matrices d'adjacence
  ## N : nombre d'individus
  ## V : nombre de vues
  # K : nombre de cluster
  # Q : nombre de composante de vues

  # tol : paramètre de convergence de l'elbo
  # iter_max : nombre d'itération maximale
  # n_init : nombre d'initialisation
  # alternate : Boolean, optimisation alternée ou non des paramètres de mélange de vues et de clustering
  # Verbose :
  # eps_conv : convergence partie point fixe de l'algorithme
  # type_init : Choix du type d'initialisation c("random","Kmeans")


  # param :
  ## param$tau
  ## param$u
  ## param$pi_k
  ## param$rho
  # param$alpha

  # elbo : Evidence Lower BOund du modèle
  #
  #
  #
  #
  #


  # Déclaration de variables
  #------------ Variables ------------
  N = nrow(A)
  V = dim(A)[3]
  output <- list()
  output$parametres <-list()
  output$n_iter <- rep(NA,n_init)
  output$elbo <- rep(NA,n_init)


  #------------ Boucle d'initialisation ------------
  for(init in 1:n_init)  {
    print(paste0("------------ Run ",init," ------------"))
    n_iter = 0
    params <- initialisation_params(A,K,Q,type_init)
    elbo <- ELBO(A,params)


    elbo_old <- -Inf
    params.old <- params

    params_best <- initialisation_params(A,K,Q,type_init)
    elbo_best <- -Inf

    #------------ Boucle VBEM - 1 run ------------
    while( (n_iter < iter_max) && (abs(elbo_old-elbo)>tol) ){
      elbo_old <- elbo
      params.old <- params

      n_iter = n_iter + 1
      if(Verbose){print(paste0("__ Intération n°", n_iter," __"))}

      params <- VEM_step(A,params,alternate,eps_conv)
      elbo <- ELBO(A,params)


      if(Verbose){print(paste0("L'Evidence Lower BOund vaut :",round(elbo,2)))}
      #if(elbo < elbo_old){ warning("l'ELBO diminue") ;}# params <- params.old}
      if(elbo>elbo_best){elbo_best <- elbo;params_best <- params}
    }

    output$parametres[[init]] <- params_best
    output$n_iter[init] <- n_iter
    output$elbo[init] <- elbo_best
  }


  #------------ Output ------------

  #output$best
  output$best <- output$parametres[[which.max(output$elbo)]]
  output$best$elbo <- output$elbo[which.max(output$elbo)]
  return(output)
}






ELBO <- function(A,params){
  alpha <- params$alpha
  rho<- params$rho
  u<- params$u
  tau<- params$tau
  pi_k<- params$pi_k

  N = nrow(A)
  K = ncol(tau)
  V = dim(A)[3]
  Q <- ncol(u)

  res <- 0
  for(v in 1:V){
    for(q in 1:Q){
      for(i in 1:(N-1)){
        for(j in (i+1):N){
          for(k in 1:K){
            for(l in k:K){
              #new_u[v,q] <- new_u[v,q] * (dbinom(A[i,j,v],1,alpha[k,l,q]) * rho[q])**(tau[i,k] * tau[j,l] )
              res <- res + tau[i,k] * tau[j,l] * u[v,q] * log(dbinom(A[i,j,v],1,alpha[k,l,q]) + .Machine$double.xmin)
            }
          }
        }
      }
    }
  }

  for(i in 1:(N)){
    for(k in 1:K){
      res <- res + tau[i,k] * ( log(pi_k[k] + .Machine$double.xmin) +log(tau[i,k]) + .Machine$double.xmin )
    }
  }
  for(v in 1:(V)){
    for(q in 1:Q){
      res <- res + u[v,q] * ( log(rho[q] + .Machine$double.xmin) + log(u[v,q]) + .Machine$double.xmin)
    }
  }
  return(res)
}
