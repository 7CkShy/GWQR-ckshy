#' @title Local linear Quantile Regression
#'
#' @description fit a local linear quantile regression with rq function,
#' @param x The design matrix
#' @param y The response variable
#' @param h The bandwidth parameter of a kernel function
#' @param tau The tau quantile (0 to 1)
#' @param m The number of points at which the function is to be estimated
#' @param kernel The kernel function with several types. \code{norm} is to chose the Normal kernal
#'
#' @return parameter estimation at sample points.
#' @export llqr
#'
#' @examples
llqr <- function(x, y, h, tau, m=200, kernel=FALSE){
  points = seq(min(x), max(x), length.out=m)
  fv = points
  dv = FALSE
  for (i in 1:length(points)){
    z = x - points[i]
    w = dnorm(z/h)
    fit = rq(y~x, tau = tau, weights = w, method = 'br')
    fv[i] = fit$coef[1]
    dv[i] = fit$coef[2]
  }

  return(list(points=points, y=y))
}



