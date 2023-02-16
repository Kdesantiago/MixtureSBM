#---------------- Update parameters ----------------

update_tau_bayesian<-function(A,params,eps_conv=1e-4){

  tau= params$tau
  u = params$u
  beta_k = params$beta_k
  eta = params$eta
  xi = params$xi

  N = dim(A)[1]
  K = ncol(tau)
  Q = ncol(u)
  V = nrow(u)
  log_tau <- matrix(0,N,K)

  A <- diag_nulle(A)
  old_tau <- tau + 1

  etp <- 0
  while( (sum(abs(old_tau - tau)) > eps_conv) && etp<50 ){
    etp <- etp+1
    old_tau <- tau
    for(i in 1:N){
      tau_tmp = tau
      tau_tmp[i,] = 0 # pour enlever le cas où i = j
      for(k in 1:K){
        for(q in 1:Q){
          log_tau[i,k] = log_tau[i,k] +  t(u[,q]) %*%(t(A[i,,]) %*% tau) %*%  ( digamma(eta[k,,q]) - digamma(xi[k,,q])) #Partie Aijv de la somme
        }
        log_tau[i,k] = log_tau[i,k] + sum( tau_tmp %*% (digamma(xi[k,,]) - digamma(eta[k,,] + xi[k,,])) %*% t(u)) # Deuxième partie (sans Aijv) de la somme

        log_tau[i,k] = log_tau[i,k] + digamma(beta_k[k]) - digamma(sum(beta_k)) # Partie esp(pi_k)
      }
    }
    tau <- log_Softmax(log_tau)

  }
  params$tau = tau
  return(params)
}

update_u_bayesian<-function(A,params){

  tau= params$tau
  theta = params$theta
  eta = params$eta
  xi = params$xi

  V = dim(A)[3]
  Q = length(theta)
  K = ncol(tau)
  log_u <- matrix(0,V,Q)

  A <- diag_nulle(A)
  A_trig_sup <- trig_sup(A)


  for(q in 1:Q){
    for(v in 1:V){
      for(k in 1:K){
        for(l in k:K){
          if(l == k){
            log_u[v,q] = log_u[v,q] + tau[,k] %*%  ( A_trig_sup[,,v]*( digamma(eta[k,l,q]) - digamma(xi[k,l,q])) + digamma(xi[k,l,q]) - digamma(eta[k,l,q] + xi[k,l,q]) ) %*% tau[,l]
          } else {
            log_u[v,q] = log_u[v,q] + tau[,k] %*%  (A[,,v]*(digamma(eta[k,l,q]) - digamma(xi[k,l,q])) + digamma(xi[k,l,q]) - digamma(eta[k,l,q] + xi[k,l,q]) ) %*% tau[,l]
          }
        }
      }
      log_u[v,q] = log_u[v,q] + digamma(theta[q]) - digamma(sum(theta))
    }
  }

  u <- log_Softmax(log_u)

  params$u <- u
  return(params)
}

update_beta_bayesian<-function(params){

  tau= params$tau
  beta_0 = params$beta_0

  beta_k <- beta_0 + apply(tau,2,sum)
  params$beta_k <- beta_k
  return(params)
}

update_theta_bayesian<-function(params){

  u= params$u
  theta_0 = params$theta_0

  theta <- theta_0 + apply(u,2,sum)
  params$theta <- theta
  return(params)
}

update_eta_bayesian<-function(A,params){

  eta_0 = params$eta_0
  tau= params$tau
  u= params$u

  K = ncol(tau)
  Q = ncol(u)
  eta <- array(0,c(K,K,Q))
  A = diag_nulle(A)
  A_trig_sup = trig_sup(A)

  for(q in 1:Q){
    for(k in 1:K){
      for(l in k:K){
        if(l == k){
          eta[k,l,q] = tau[,k] %*%  apply(A_trig_sup,c(1,2),function(x){x %*% u[,q]}) %*% tau[,l]
        } else {
          eta[k,l,q] = tau[,k] %*%  apply(A,c(1,2),function(x){x %*% u[,q]}) %*% tau[,l]
          eta[l,k,q] = eta[k,l,q]
        }
      }
    }
    eta[k,l,q] = eta[k,l,q] + eta_0[k,l,q]
  }
  params$eta <- pmax(eta,.Machine$double.eps)#eta
  return(params)
}

update_xi_bayesian<-function(A,params){

  xi_0 = params$xi_0
  tau= params$tau
  u= params$u

  K = ncol(tau)
  Q = ncol(u)
  xi <- array(0,c(K,K,Q))

  A = 1-A # Car xi dépend de 1 - Aijv
  A = diag_nulle(A)
  A_trig_sup = trig_sup(A)

  for(q in 1:Q){
    for(k in 1:K){
      for(l in k:K){
        if(l == k){
          xi[k,l,q] = t(tau[,k]) %*%  apply(A_trig_sup,c(1,2),function(x){x %*% u[,q]}) %*% tau[,l]
        } else {
          xi[k,l,q] = t(tau[,k]) %*%  apply(A,c(1,2),function(x){x %*% u[,q]}) %*% tau[,l]
          xi[l,k,q] = xi[k,l,q]
        }
      }
    }
    xi[k,l,q] = xi[k,l,q] + xi[k,l,q]
  }
  params$xi <-  pmax(xi,.Machine$double.eps) #xi
  return(params)
}

#---------------- Loss ----------------
Loss_BayesianMSBM <- function(params){

  tau = params$tau
  u = params$u
  beta_k = params$beta_k
  theta = params$theta
  eta = params$eta
  xi = params$xi
  beta_0 = params$beta_0
  theta_0 = params$theta_0
  eta_0 = params$eta_0
  xi_0 = params$xi_0

  K = ncol(tau)
  Q = ncol(u)

  res <- - sum(tau*log(tau)) - sum(u*log(u))

  res <- res + log(gamma(sum(beta_0)) * prod(gamma(beta_k)) /(gamma(sum(beta_k)) * prod(gamma(beta_0)))  )


  res <- res + log(gamma(sum(theta_0)) * prod(gamma(theta)) /(gamma(sum(theta)) * prod(gamma(theta_0)))  )

  tmp1 <- array(NA,c(K,K,Q))
  tmp2 <- array(NA,c(K,K,Q))
  for(k in 1:K){
    for(l in 1:K){
      for(q in 1:Q){
        tmp2[k,l,q] = beta(eta_0[k,l,q],xi_0[k,l,q])
        tmp1[k,l,q] = beta(eta[k,l,q],xi[k,l,q])
      }
    }
  }
  tmp2 <- pmax(tmp2,.Machine$double.xmin)
  tmp1 <- pmax(tmp1,.Machine$double.xmin)

  res <- res + sum(trig_sup(log(tmp1 / tmp2),diag = F))

  return(res)

}

#---------------- Variational Bayes EM ----------------

VBEM_step<-function(A,params,alternate=T,eps_conv=1e-3){

  if(alternate){
    params <- update_u_bayesian(A,params)
    params <- update_beta_bayesian(params)
    params <- update_theta_bayesian(params)
    params <- update_eta_bayesian(A,params)
    params <- update_xi_bayesian(A,params)

    params <- update_tau_bayesian(A,params,eps_conv)
    params <- update_beta_bayesian(params)
    params <- update_theta_bayesian(params)
    params <- update_eta_bayesian(A,params)
    params <- update_xi_bayesian(A,params)


  } else{
    params <- update_u_bayesian(A,params)
    params <- update_tau_bayesian(A,params,eps_conv)
    params <- update_beta_bayesian(params)
    params <- update_theta_bayesian(params)
    params <- update_eta_bayesian(A,params)
    params <- update_xi_bayesian(A,params)
  }


  return(params)
}

#---------------- Initialisation ----------------

initialisation_params_bayesian <-function(A,K,Q,beta_0=rep(1/2,K),theta_0=rep(1/2,Q),eta_0=array(rep(1/2,K*K*Q),c(K,K,Q)),xi_0=array(rep(1/2,K*K*Q),c(K,K,Q)),type_init="random"){
  N = nrow(A)
  V = dim(A)[3]
  params <- list()

  params$beta_0 = beta_0
  params$theta_0 = theta_0
  params$eta_0 = eta_0
  params$xi_0 = xi_0


  if(type_init=="Kmeans"){
    A_tmp <- apply(A,c(1,2),sum)
    tmp <- kmeans(A_tmp,centers = K,nstart = 50)

    params$tau <- one_hot_errormachine(tmp$cluster)
    params$u <- array(runif(N*K,0,1),dim=c(V,Q)) ; params$u <- params$u/apply(params$u,1,sum)

    params <- update_beta_bayesian(params)
    params <- update_theta_bayesian(params)
    params <- update_eta_bayesian(A,params)
    params <- update_xi_bayesian(A,params)

  } else {
    params$tau <- array(runif(N*K,0,1),dim=c(N,K)) ; params$tau <- params$tau/apply(params$tau,1,sum)
    params$u <- array(runif(N*K,0,1),dim=c(V,Q)) ; params$u <- params$u/apply(params$u,1,sum)

    params <- update_beta_bayesian(params)
    params <- update_theta_bayesian(params)
    params <- update_eta_bayesian(A,params)
    params <- update_xi_bayesian(A,params)
  }


  return(params)
}



#---------------- Model ----------------

BayesianMixture_SBM <-function(A,K,Q,beta_0=rep(1/2,K),theta_0=rep(1/2,Q),eta_0=array(rep(1/2,K*K*Q),c(K,K,Q)),xi_0=array(rep(1/2,K*K*Q),c(K,K,Q)),tol=1e-3,iter_max=10,n_init = 1,alternate=T, Verbose=T,eps_conv=1e-4,type_init="random"){


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
    params <- initialisation_params_bayesian(A,K=K,Q=Q,beta_0,theta_0,eta_0,xi_0,type_init)
    elbo = Loss_BayesianMSBM(params)


    elbo_old <- -Inf
    params.old <- params

    params_best <- initialisation_params_bayesian(A,K=K,Q=Q,beta_0,theta_0,eta_0,xi_0,type_init)
    elbo_best <- -Inf

    #------------ Boucle VBEM - 1 run ------------
    while( (n_iter < iter_max) && (abs(elbo_old-elbo)>tol) ){
      elbo_old <- elbo
      params.old <- params

      n_iter = n_iter + 1
      if(Verbose){print(paste0("__ Intération n°", n_iter," __"))}

      params <- VBEM_step(A,params,alternate,eps_conv)
      elbo <- Loss_BayesianMSBM(params)


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
