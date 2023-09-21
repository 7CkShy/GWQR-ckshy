#' @title Set a random quantile from 0 to 1
#' @describeIn It is to set a random quantile for a simulation if you have
#' no idea with setting quantile, return a tau quantile(numeric)
#' @param tau a input parameter with always TRUE
#' @import quantreg
#' @return \code{tau} the random quantile at round(x,2)
#' @export
#'
#' @examples
SetQuantile = function(tau=TRUE){
  c = runif(1)
  return(round(c,2))
}
