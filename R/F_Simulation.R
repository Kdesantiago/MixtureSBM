#---------------- Simulations basées sur le modèle ----------------

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


#---------------- Simulations Recouvrement ----------------

partition_K_barre <-function(K_barre,K){
  #---------------------- Variables ----------------------------
  # K : Nombre de clusters au sein de chaque composantes  (vecteur de taille Q)
  # K_barre : Nombre de vrais clusters
  ## Hypothèse : K[q] < K_barre pour tout q


  Q = length(K)
  clusters <- list()

  #---------------------- Partition clusters ----------------------------
  ## On commence par associer 1 classe à chaque cluster de K[q]
  for(q in 1:Q){
    tmp <- list()
    idx_K_barre = 1:K_barre
    idx <- sample(idx_K_barre,replace = F,size=K[q])
    idx_K_barre = idx_K_barre[-which(idx_K_barre %in% idx)]
    for(l in 1:K[q]){tmp[[l]] = idx[l]}

    ## On fait un tirage aléatoire pour chaque classe qui n'a pas été associée.
    for(idx in idx_K_barre){
      idx_K <- sample(1:K[q],size=1)
      tmp[[idx_K]] <- c(tmp[[idx_K]],idx)
    }

    clusters[[q]] <- tmp
  }

  # ----------------------  Output ----------------------------

  return(clusters) # Cluster[[Q]] [[K[Q]]]

}

Mat_lien_alpha <- function(clusters,K_barre,K){

  #---------------------- Variables ----------------------------
  # K : Nombre de clusters au sein de chaque composantes  (vecteur de taille Q)
  # K_barre : Nombre de vrais clusters
  # clusters : list indiquant les liens entre K_Barre et K ( Cluster[[Q]][[K[Q]]]  )

  Q = length(K)

  #---------------------- Simulation alpha ----------------------------
  alpha <- array(0,dim=c(K_barre,K_barre,Q))
  for(q in 1:Q){
    for(s in 1:K[q]){
      for(k in clusters[[q]][[s]]){
        for(l in clusters[[q]][[s]]){
          alpha[k,l,q] = 1
        }
      }
    }
  }
  return(alpha)
}

rSMB_partition <- function(N, V, K, pi_k,rho,sorted=T){

  #---------------------- Variables ----------------------------
  # N : Nombre d'observations à générer
  # V : Nombre de vues à générer
  # K : Nombre de clusters au sein de chaque composantes  (vecteur de taille Q)
  # K_barre : Nombre de vrais clusters
  # pi_k : Vecteur de probabilité d'appartenance aux clusters de K_barre
  # rho : Vecteur de probabilité d'appartenance aux clusters de W (mélange vue)
  # sorted : Boolean pour le réagencement des simulations.

  # clusters : list indiquant les liens entre K_Barre et K ( Cluster[[Q]][[K[Q]]]  )
  # alpha : array, probabilité de liaison entre 2 communautés. (dim = K_Barre,K_barre,Q)

  K_barre <- length(pi_k)

  #----------------------  K_Barre -> K ----------------------------

  clusters <- partition_K_barre(K_barre,K)

  # ----------------------  Alpha ----------------------------

  alpha <- Mat_lien_alpha(clusters,K_barre,K)

  # ----------------------  Z / W / A ----------------------------

  output <- rMSBM(N,V,alpha_klq = alpha,pi_k,rho,sorted=T)

  # ----------------------  Output ----------------------------

  output$params$clusters <- clusters
  output$params$K_barre <- K_barre
  output$params$K <- K
  return(output)
}




#---------------- Simulations Eclatement clusters ----------------

rKbarre_to_K <-function(params, size){
  #---------------------- Variables ----------------------------
  # idx_K : Indices des colonnes pour les clusters de chaque composante, afin de permettre le choix de quel cluster final influencera quels clusters au sein de chaque composantes
  # size : permet de choisir le nombre de tirage fait au sein des groupes pour composer les groupes éclatés.
  # clusters : Lien entre K_barre et K

  Q = params$Q
  K_barre = params$K_barre
  K = params$K

  #---------------------- Indice composantes ----------------------------
  idx_K <- list()
  for(q in 1:Q){
    idx_K[[q]] <- 1:K[q]
  }

  #---------------------- Eclatement des clusters ----------------------------
  clusters <- list()
  for(k in 1:(K_barre-1)){
    tmp2 <- list()
    for(q in 1:Q){
      if(identical(idx_K[[q]], integer(0)) ){next;} #vérification vecter indice pas vide

      #------ Partie du code modifiable pour ajuster les simulations ---------------

      ## Simulation qui fait le plus de sens :
      # si K multiple de K_barre : permet d'avoir un nombre de clusters qui font de bons groupes
      idx <- unique(sample(idx_K[[q]],size = size,replace = F))

      ## Simulation avec + ou - de sens dans la création des clusters :
      #Choisi un nombre size de clusters à chaque fois, avec possibilité d'en tirer les mêmes,
      #idx <- unique(sample(idx_K[[q]],size = size,replace = T))



      ## Simulation chaotique :
      # Nombre de clusters tirés aléatoirement au lieu de size
      #idx <- unique(sample(idx_K[[q]],size = sample(1:K[q],1),replace = T))

      #------------------------------------------------------------------
      idx_K[[q]] <- idx_K[[q]][-which(idx_K[[q]] %in% idx)]
      tmp2[[q]] <- idx
    }
    clusters[[k]] <- tmp2
  }


  #-----------------------------  Dernier cluster -----------------------------------

  # On se sert des clusters non utilisés pour former le dernier groupe
  k = K_barre
  tmp2 <- list()
  for(q in 1:Q){
    if(identical(idx_K[[q]], integer(0)) ){next;}
    tmp2[[q]] <-  idx_K[[q]]
  }
  clusters[[k]] <- tmp2

  #---------------------- Output ----------------------------
  return(clusters)
}


rZ <- function(params,Z_barre,clusters,sorted_by_Z=F){
  #---------------------- Variables ----------------------------
  # N : Nombre d'observations à générer
  # Q : Nombre de composantes du mélange de vues
  # K : Nombre de clusters au sein de chaque composantes  (vecteur de taille Q)
  # Z : Variables simulées d'appartenance à des clusters selon la composante du mélange de vue
  # Z_barre : vrai clustering des données.

  # sorted_by_Z : variable indiquant l'ordre de classement des variables
  # F -> dans l'ordre de Z_barre
  # T -> dans l'ordre des clusters générés dans chaque Z

  N <- params$N
  Q <- params$Q
  K <- params$K
  Z <- list()

  #---------------------- Génération des Z ----------------------------
  for(q in 1:Q){
    val <- matrix(0,N,K[q])
    for(i in 1:N){
      if (q > length(clusters[[Z_barre[i]]])){idx <- sample(1:K[q],1)}
      # Erreur cas où on a pas de groupe prédéfini pour le dernier élément des listes
      else if(is.null(clusters[[Z_barre[i]]][[q]])){idx <- sample(1:K[q],1)}
      # Erreur cas où on a pas de groupe prédéfini pour un élément (PAS DERNIER des listes).
      else if(length(clusters[[Z_barre[i]]][[q]]) ==1) {val[i,clusters[[Z_barre[i]]][[q]]] = 1  }
      # Erreur lorsqu'il y a qu'un seul cluster dans un le choix car sample considère ça comme un int
      # donc utilise la fonction sample.int au lieu de sample.
      else {idx <- sample(x=clusters[[Z_barre[i]]][[q]],size=1);val[i,idx] = 1  }

    }
    Z[[q]] <- max.col(val)
  }

  if(sorted_by_Z){
    for(q in Q:1){
      Z[[q]] <- Z[[q]][order(Z[[1]])]  #sort(max.col(Z2[[q]]))
    }
  }


  return(Z)
}


rmultimodal <- function(params,size=2,sorted_by_Z=F){
  #---------------------- Variables ----------------------------
  # params$N : Nombre d'observations à générer
  # params$V : Nombre de vues à générer
  # params$Q : Nombre de composantes du mélange de vues
  # params$K_barre : Nombre de vrais clusters
  # params$K : Nombre de clusters au sein de chaque composantes  (vecteur de taille Q)
  # clusters : Lien entre K_barre et K


  # params$p_Zbarre : vecteur de probabilité pour les clusters de Z_barre
  # params$p_W: vecteur de probabilité pour les clusters des vues

  # clusters : Lien entre K_barre et K
  # Z_barre : vrai clustering des données.
  # Z :Variables simulées d'appartenance à des clusters selon la composante du mélange de vue


  # size : Nombre d'éléments à mettre pour le lien entre K_barre et K
  # sorted_by_Z : variable indiquant l'ordre de classement des variables
  # F -> dans l'ordre de Z_barre
  # T -> dans l'ordre des clusters générés dans chaque Z
  # B : Matrices d'émissions
  # A : Matrices d'adjacences simulées

  #----------------------  K_Barre -> K ----------------------------

  clusters <- rKbarre_to_K(params, size=size)

  # ----------------------  Z_barre et Z ----------------------------
  Z_barre <- t(rmultinom(n=params$N,size =1,prob = params$p_Zbarre))
  Z_barre <- sort(max.col(Z_barre))

  Z <- rZ(params,Z_barre,clusters,sorted_by_Z=sorted_by_Z)

  # ----------------------  Vues et mélange de vues ----------------------------
  W <- t(rmultinom(n=params$V,size =1,prob = params$p_W ))
  W <- sort(max.col(W))


  # ----------------------  Matrices alpha ----------------------------

  # A mettre directement dans params ?

  B <- list()
  a = 1
  b=0.01
  for(q in 1:params$Q) {B[[q]] = matrix(b,params$K[q],params$K[q]); diag(B[[q]]) <- a}

  # ----------------------  Matrices d'adjacences ----------------------------
  A <- array(NA,dim=c(params$N,params$N,params$V))
  for(l in 1:params$V){
    for(i in 1:(params$N-1)){
      for(j in (i+1):params$N){
        q = W[l]
        Zi = Z[[q]][i]
        Zj =  Z[[q]][j]
        A[i,j,l] <- A[i,j,l] <- rbinom(1,1,prob=B[[q]][Zi,Zj])
        A[j,i,l] <- A[i,j,l]
      }
    }
    diag(A[,,l]) <- 1
  }

  # ----------------------  Output ----------------------------
  params$B <- B
  simulations <- list(A = A, Z = Z, Z_barre = Z_barre, clusters = clusters, W = W )
  output <- list(params=params,simulations=simulations)
  return(output)
}
