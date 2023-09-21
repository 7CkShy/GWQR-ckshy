#' @title Geographically Weighted Regression
#'
#' @param X Design matrix
#' @param Y Response variable
#' @param h Bandwidth of kernel function
#' @param d The distance matrix of sample points
#' @param kernal The name of kernel function. only \code{Norm},\code{Binary}, \code{Exp},
#' \code{bi-square}, \code{rank}, \code{adaptive} can chosen
#'
#' @return \code{beta} regression estimation
#' @return \code{Y.fitted}
#' @return \code{Weight}
#' @return \code{AICc}
#' @return \code{MSE}
#' @export
#'
#' @examples
GWR <- function(X, Y, h, d, maxit=500, kernel=FALSE, fig=F, parallel=FALSE, n_parallel=6,
                dstar=FALSE, alpha=20, tolar=1e-3, process=T, adaptive=F){
  n = nrow(X) #number of sample
  p = ncol(X) #number of variables(with intercept)
  iteration = 1
  error = 10000
  AICc = 0
  eps = 10000
  eps.old = 0
  AICc_plot = c()
  bandwidth_plot = c()
  mse_plot = c()
  #hmax = nrow(X)
  #uph = h + alpha
  #beta = matrix(rep(1, n*p), n, p) #initial beta matrix
  cat('Program start at:', as.character(Sys.time()), '\n')

  if (kernel==FALSE){
    #cat('kernel parameter error! please choose a kernel function')
    stop('kernel parameter error! please choose a kernel function')
  }
  if (h<0){
    #cat('h parameter error! bandwidth must be positive!')
    stop('h parameter error! bandwidth must be positive!')
  }
  if (parallel){
    cat('Setting parallel computing...','\n')
    cl.cores = detectCores(logical = F)
    cl <- makeCluster(cl.cores-1)
    registerDoParallel(cl)

    # cl <- makeCluster(n_parallel)
    # registerDoParallel(cl)
  }

  cat('-------------------Parameter Setting------------------', '\n')
  cat('Maxit:', maxit,';', 'Kernel function:', kernel, '\n')
  cat('fig:', fig,';', 'Tolar:', tolar, '\n')
  cat('alpha:', alpha,';', 'Process:', process, '\n')
  cat('parallel computing:', parallel,';', 'number of cores:', cl.cores-1, '\n')
  cat('------------------------------------------------------', '\n')
  cat('\n')


  H = matrix(1, n, n)
  while (abs(eps-eps.old) >= tolar & iteration < maxit){

    AICc.old = AICc
    eps.old = eps
    # beta <- foreach(i = 1:n, .combine = cbind,
    #                 .errorhandling = "remove",
    #                 .packages = c('LLQR')) %dopar%{
    #                   W = Kernel(Kname=kernel, d=d[i, ], h=h)
    #                   reg_gwr(X=X, Y=Y, W=W)
    #                 }

    if (parallel){
      beta <- foreach::foreach(
        i=1:n,
        .combine=cbind,
        .errorhandling = "pass", .export = 'GWR_reg',
        .packages = c("LLQR", "GWmodel"))%dopar% {
        W = diag(Kernel(Kname=kernel, d=d[i, ], h=h))
        GWR_reg(X, Y, W=W)

        }
      H <- foreach::foreach(
        i=1:n,
        .combine=rbind,
        .errorhandling = "pass", .export = 'GWR_hat',
        .packages = c("LLQR", "GWmodel"))%dopar% {
        W = diag(Kernel(Kname=kernel, d=d[i, ], h=h))
        GWR_hat(X=X, W=W, i=i)
        }
    }
    else{
      for (i in 1:n){
        W = diag(Kernel(Kname=kernel, d=d[i, ], h=h))
        #H[i, ] = X %*% ginv(t(X) %*% W %*% X) %*% t(X) %*% W
        XtWX_inv = try(solve(t(X) %*% W %*% X))
        if ('try-error' %in% class(XtWX_inv)){next}
        else{beta[i, ] = XtWX_inv %*% t(X) %*% W %*% Y}
      }
    XtWX_inv = solve(t(X) %*% W %*% X)
    H = X %*% XtWX_inv %*% t(X) %*% W
    }

    if(parallel){
    Y.fitted = rowSums(X * t(beta))}
    else{Y.fitted = rowSums(X * beta)}
    error = Y - Y.fitted #a vector

    sigma_hat = sum((error)^2)/(n - 2*tr(H) + tr(t(H) %*% H))
    AICc = 2 * n * log(sqrt(sigma_hat)) + n* log(2*pi) + n*((n + tr(H))/(n - 2 - tr(H)))
    BIC = n * log(sum(error^2) / n) + n * log(2.0 * pi) + log(n) * tr(H)
    MSE = mean((Y - Y.fitted)^2)
    mse_plot[iteration] = MSE
    bandwidth_plot[iteration] = h
    AICc_plot[iteration] = AICc
    h = h + alpha
    iteration = iteration + 1
    eps = abs(AICc - AICc.old)
    R2 = 1 - sum(error^2) / sum((Y - mean(Y))^2)
    adjR2 = 1 - ((n-1)/(n-p-1))*(1-R2)
    if (process){
      cat('h=', h,'\n')
      cat('AICc:', AICc, '\n')
      cat('BIC:', BIC, '\n')
      cat('MSE:', MSE, '\n')
      cat('epsilon:', eps, '\n')
    }

  }
###################plot######################
  if (fig){
    #plot(bandwidth_plot, mse_plot, type='p', xlab='bandwidth', ylab='MSE')
    plot(bandwidth_plot, AICc_plot, type='p', xlab='bandwidth', ylab='AICc')
    #print(H)
    #print(mse)
  }
###################plot######################
  cat('-------------------------Model Summary-----------------------', '\n')
  cat('R-squared', R2, '\n')
  cat('adjR-squared', adjR2, '\n')
  cat('Model Bandwidth:', h, '\n')
  cat('AICc:', AICc, '\n')
  cat('BIC:', BIC, '\n')
  cat('MSE:', MSE, '\n')
  if (!adaptive){
    cat('Bandwidth type: Fixed', '\n')
  }else{
    cat('Bandwidth type: Adaptive', '\n')
    }
  cat('---------------------------Paremeter Estimation------------------------------', '\n')
###############print summary of estimates############
  Min = c()
  Qu.25 = c()
  Median = c()
  Qu.75 = c()
  Max = c()
  Mean = c()
  for (i in 1:p){
    Min[i] = min(beta[, i])
    Qu.25[i] = quantile(beta[, i], 0.25)
    Median[i] = quantile(beta[, i], 0.5)
    Qu.75[i] =  quantile(beta[, i], 0.75)
    Max[i] = max(beta[, i])
    Mean[i] = mean(beta[, i])
  }

  cmat <- cbind(Min, Qu.25, Median, Qu.75, Max, Mean) #coef matrix
  colnames(cmat) <- c("Min.", "0.25 Qu.", "Median", "0.75 Qu.", 'Max.', 'Mean')
  #namelist = c()
  if (length(colnames(X))==0){
    if (sum(X[, 1])==length(X[, 1])){
      name = rep('x', p-1)
      index = 2:p
      namelist = paste(name, index, sep='')
      namelist = c('intercept', namelist)
      rownames(cmat) <- namelist
      printCoefmat(cmat)
    }
    else{
      name = rep('x', p)
      index = 1:p
      namelist = paste(name, index, sep='')
      rownames(cmat) <- namelist
      printCoefmat(cmat)
    }
  }
  else{
    rownames(cmat) <- colnames(X)
    printCoefmat(cmat)
  }

  cat('-----------------------------------------------------------------------------', '\n')
###############print summary of estimates############

  cat('Program end at', as.character(Sys.time()), '\n')

  stopImplicitCluster()
  return(list(beta=beta, Y.fitted = Y.fitted, Weight=W, AICc = AICc,
              MSE=MSE, bandwidth=h, Effective=n-2*tr(H)+tr(t(H)%*%H),
              bandwidth_plot = bandwidth_plot, mse_plot=mse_plot))
}

GWR_reg <- function(X,Y,W){
  #res = solve(t(X) %*% X) %*% t(X) %*% Y
  #W = Kernel(Kname=kernel, d=d[i, ], h=h)
  res = solve(t(X) %*% W %*% X) %*% t(X) %*% W %*% Y
  return(res)
}

GWR_hat <- function(X, W, i){
  #W = diag(Kernel(Kname='gaussian', d=d[i, ], h=40))
  res =  X[i,] %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
  return(res)
}


#fit = GWR(X, Y, h=40, d, maxit=200, fig=T, kernel='gaussian',
#          dstar=FALSE, alpha=0.1, tolar=1e-4, process = T)
# fit3 = GWR(X, Y, h=40, d=d, maxit=200, fig=T, kernel='gaussian',
#            dstar=FALSE, alpha=0.1, tolar=1e-4, process = T,
#            parallel = TRUE)

