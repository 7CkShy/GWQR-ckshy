#' geographical weighted quantile regression
#'
#' @param dp
#' @param d distance matrix
#' @param X input
#' @param y output
#' @param h bandwidth
#' @param fig logic(0 or 1) whether to plot the MSE by bandwidth
#' @param tau quantile
#' @param kernel kernel function to choose
#' @param alpha the step to update bandwidth
#' @param tolar
#' @param process
#' @param adaptive
#' @param iteration
#'
#' @return \code{beta} regression coef
#' @export
#'
#' @examples
GWQR <- function(dp=F, d, X, y, h, fig=F, tau, kernel=F, alpha=20,
                 tolar=1e-3, process=T, adaptive=F,iteration=200,method='br'){
  mse = c()
  H = c()
  n = nrow(X)
  p = ncol(X)
  U = matrix(rep(0, n*n), n, n)
  V = matrix(rep(0, n*n), n, n)
  X_inter = cbind(1,X)
  #W = matrix(rep(1),1,9)

  cv = matrix(rep(0, n*(p+1)), n, p+1)

  for (it in 1:iteration){  #h from 1 to 100
    if (kernel==FALSE){
      stop('choose a kernel function')}
    for (j in 1:n){
      # for (i in 1:n){
      #   U[i,i] = abs(dp[i, 1] - dp[j, 1])
      #   V[i,i] = abs(dp[i, 2] - dp[j, 2])
      #   }
      # Xnew = cbind(X, U%*%X, V%*%X)
      W = Kernel(Kname=kernel, d=d[j, ], h=h)
      fit = try(rq(Y ~ X, weights=W, tau=tau, method = method))
      if ('try-error' %in% class(fit)){
        #print('yes')
        break}
      else{
        #cv[j, ] = fit$coefficients[1:(dim(X)[2]+1)]
        cv[j, ] = fit$coefficients
        }
      }
    Y.fitted = rowSums(cv*X_inter)
    #Y.fitted = cv*X
    error = (Y - Y.fitted)^2
    #H = X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
    #sigma_hat = sum((error)^2)/(n - 2*tr(H) + tr(t(H) %*% H))
    MSE = round(mean(error), 4)
    R2 = 1 - sum(error) / sum((Y - mean(Y))^2)
    adjR2 = 1 - ((n-1)/(n-p-1))*(1-R2)
    h = h + alpha
    cat("h:", h, 'MSE:', MSE, 'iter:', it, '\n')
    mse[it] = MSE
    H[it] = h
    }
  if (fig){
    plot(H, mse, type='p', ylab='MSE', xlab='Bandwidth')
    fig = FALSE
    }

  return(list(beta=cv, MSE=MSE, h=h, R2 = R2, adjR2=adjR2))

  }


#fit = GWQR(dp=dp, d=d, X=X, y=Y, h=50, kernel='bi-square',
#           alpha=0.5, tau=0.5, iteration=20, fig=T)
