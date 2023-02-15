plot_adjency<- function(A){
  #------------ Objectif ------------
  # Permet de visualiser les différentes matrices d'adjacences du tensor A

  #------------ Variable ------------
  # A : array(.,.,V)
  V <- dim(A)[3]

  #------------ Visualisation ------------
  for(v in 1:V){
    image(A[,,v],axes=F)
  }
  image(apply(A,c(1,2),sum))
}
