#' @export gwr.dist
gwr.dist <- function(dp.locat, rp.locat, focus = 0, p = 2, theta = 0, longlat = F)
{
  if (missing(dp.locat)||!is.numeric(dp.locat)||dim(dp.locat)[2]!=2)
    stop("Please input correct coordinates of data points")

  if (!missing(rp.locat)) {
    if (!is.numeric(rp.locat))
      stop("Please input correct coordinates of regression points")
    else
      rp.locat <- matrix(rp.locat, ncol = 2)
    if (focus < 0 || focus > length(rp.locat[,1]))
      stop("No regression point is fixed")
    dists <- gw_dist(dp.locat, rp.locat, focus - 1, ifelse(is.infinite(p), -1, p), theta, longlat, T)
  }
  else dists <- gw_dist(dp.locat, dp.locat, focus - 1, ifelse(is.infinite(p), -1, p), theta, longlat, F)
  dists
}
