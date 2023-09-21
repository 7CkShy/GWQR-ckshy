#' @export Kernel
Kernel <- function(Kname=FALSE, d, dstar=FALSE, h=FALSE, dN=FALSE){
  if (Kname == "binary" | Kname=="Binary"){
    W = ifelse(d<=dstar, yes=1, no=0)
  }

  if (Kname == "gaussian" | Kname=="Gaussian" | Kname == FALSE){
    W = exp(-0.5 * (d/h)^2)
  }

  if (Kname == "Exp" | Kname=="exp"){
    W = exp(-d/h)
  }

  if (Kname == "bi-square" | Kname=="Bi-square"){
    W = ifelse(d<h, yes=(1 - (d/h)^2)^2, no=0)
  }

  if (Kname == "rank" | Kname=="Rank"){
    W = exp(-rank(d)/h)
  }

  if (Kname == "adaptive" | Kname=="Adaptive"){
    W = ifelse(d<dN, yes=(1-(d/dN)^2)^2, no=0)
  }
  return(W)
}


# set.seed(1)
# d = matrix(rnorm(100,0,1),20,5)
# Kernel('rank', d, h=0.5)
