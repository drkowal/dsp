#----------------------------------------------------------------------------
#' Simulate noisy observations from a function
#'
#' Builds upon the \code{make.signal()} function in the \code{wmtsa} package
#' to include Gaussian noise with a user-specificied root-signal-to-noise ratio.
#'
#' @param signalName string matching the "name" argument in the \code{make.signal()} function,
#' e.g. "bumps" or "doppler"
#' @param T number of points
#' @param RSNR root-signal-to-noise ratio
#' @param include_plot logical; if TRUE, include a plot of the simulated data and the true curve
#'
#' @return a list containing
#' \itemize{
#' \item the simulated function \code{y}
#' \item the true function \code{y_true}
#' \item the true observation standard devation \code{sigma_true}
#' }
#'
#' @note The root-signal-to-noise ratio is defined as RSNR = [sd of true function]/[sd of noise].
#'
#' @examples
#' sims = simUnivariate() # default simulations
#' names(sims) # variables included in the list
#'
#' @import wmtsa
#' @export
simUnivariate = function(signalName = "bumps", T = 200, RSNR = 10, include_plot = TRUE){

  # The true function:
  y_true = attr(make.signal(signalName, n=T), 'data')

  # Noise SD, based on RSNR (also put in a check for constant/zero functions)
  sigma_true = sd(y_true)/RSNR; if(sigma_true==0) sigma_true = sqrt(sum(y_true^2)/T)/RSNR + 10^-3

  # Simulate the data:
  y = y_true + sigma_true*rnorm(T)

  # Plot?
  if(include_plot) {t = seq(0, 1, length.out=T); plot(t, y, main = 'Simulated Data and True Curve'); lines(t, y_true, lwd=8, col='black') }

  # Return the raw data and the true values:
  list(y = y, y_true = y_true, sigma_true = sigma_true)
}
#----------------------------------------------------------------------------
#' Simulate noisy observations from a dynamic regression model
#'
#' Simulates data from a time series regression with dynamic regression coefficients.
#' The dynamic regression coefficients are selected using the options from the
#' \code{make.signal()} function in the \code{wmtsa} package.
#'
#' @param signalNames vector of strings matching the "name" argument in the \code{make.signal()} function,
#' e.g. "bumps" or "doppler"
#' @param T number of points
#' @param RSNR root-signal-to-noise ratio
#' @param p_0 number of true zero regression terms to include
#' @param include_intercept logical; if TRUE, the first column of X is 1's
#' @param scale_all logical; if TRUE, scale all regression coefficients to [0,1]
#' @param include_plot logical; if TRUE, include a plot of the simulated data and the true curve
#'
#' @return a list containing
#' \itemize{
#' \item the simulated function \code{y}
#' \item the simulated predictors \code{X}
#' \item the simulated dynamic regression coefficients \code{beta_true}
#' \item the true function \code{y_true}
#' \item the true observation standard devation \code{sigma_true}
#' }
#'
#' @note The number of predictors is \code{p = length(signalNames) + p_0}.
#'
#' @note The root-signal-to-noise ratio is defined as RSNR = [sd of true function]/[sd of noise].
#'
#' @examples
#' sims = simRegression() # default simulations
#' names(sims) # variables included in the list
#'
#' @import wmtsa
#' @export
simRegression = function(signalNames = c("bumps", "blocks"), T = 200, RSNR = 10, p_0 = 5, include_intercept = TRUE, scale_all = TRUE, include_plot = TRUE){

  # True number of signals
  p_true = length(signalNames)

  # Total number of covariates (non-intercept)
  p = p_true + p_0

  # Simulate the true regression signals
  beta_true = matrix(0, nr = T, nc = p)
  for(j in 1:p_true) beta_true[,j] = attr(make.signal(signalNames[j], n=T), 'data');
  if(scale_all) beta_true[,1:p_true] = apply(as.matrix(beta_true[,1:p_true]), 2, function(x) (x - min(x))/(max(x) - min(x)))

  # Simulate the predictors:
  X = matrix(rnorm(T*p), nr=T, nc = p)

  # If we want an intercept, simply replace the first column w/ 1s
  if(include_intercept) X[,1] = matrix(1, nr = nrow(X), nc = 1)

  # The true response function:
  y_true = rowSums(X*beta_true)

  # Noise SD, based on RSNR (also put in a check for constant/zero functions)
  sigma_true = sd(y_true)/RSNR; if(sigma_true==0) sigma_true = sqrt(sum(y_true^2)/T)/RSNR + 10^-3

  # Simulate the data:
  y = y_true + sigma_true*rnorm(T)

  # Plot?
  if(include_plot) {t = seq(0, 1, length.out=T); plot(t, y, main = 'Simulated Data and True Curve'); lines(t, y_true, lwd=8, col='black') }

  # Return the raw data and the true values:
  list(y = y, X = X, beta_true = beta_true, y_true = y_true, sigma_true = sigma_true)
}
#----------------------------------------------------------------------------
#' Initialize the evolution error variance parameters
#'
#' Compute initial values for evolution error variance parameters under the various options:
#' dynamic horseshoe prior ('DHS'), horseshoe prior ('HS'), or normal-inverse-gamma prior ('NIG')
#'
#' @param omega \code{T x p} matrix of evolution errors
#' @param evol_error the evolution error distribution; must be one of
#' 'DHS' (dynamic horseshoe prior), 'HS' (horseshoe prior), or 'NIG' (normal-inverse-gamma prior)
#' @return List of relevant components: \code{sigma_wt}, the \code{T x p} matrix of evolution standard deviations,
#' and additional parameters associated with the DHS and HS priors.
initEvolParams = function(omega, evol_error = "DHS"){

  # Check:
  if(!((evol_error == "DHS") || (evol_error == "HS") || (evol_error == "NIG"))) stop('Error type must be one of DHS, HS, or NIG')

  # Make sure omega is (n x p) matrix
  omega = as.matrix(omega); n = nrow(omega); p = ncol(omega)

  if(evol_error == "DHS") return(initDHS(omega))

  if(evol_error == "HS"){
    tauLambdaj = 1/omega^2;
    xiLambdaj = 1/(2*tauLambdaj); tauLambda = 1/(2*colMeans(xiLambdaj)); xiLambda = 1/(tauLambda + 1)

    # Parameters to store/return:
    return(list(sigma_wt = 1/sqrt(tauLambdaj), tauLambdaj = tauLambdaj, xiLambdaj = xiLambdaj, tauLambda = tauLambda, xiLambda = xiLambda))
  }
  if(evol_error == "NIG") return(list(sigma_wt = tcrossprod(rep(1,n), apply(omega, 2, function(x) sd(x, na.rm=TRUE)))))
}
#----------------------------------------------------------------------------
#' Initialize the evolution error variance parameters
#'
#' Compute initial values for evolution error variance parameters under the dynamic horseshoe prior
#'
#' @param omega \code{T x p} matrix of evolution errors
#' @return List of relevant components: the \code{T x p} evolution error SD \code{sigma_wt},
#' the \code{T x p} log-volatility \code{ht}, the \code{p x 1} log-vol unconditional mean(s) \code{dhs_mean},
#' the \code{p x 1} log-vol AR(1) coefficient(s) \code{dhs_phi},
#' the \code{T x p} log-vol innovation SD \code{sigma_eta_t} from the PG priors,
#' the \code{p x 1} initial log-vol SD \code{sigma_eta_0},
#' and the mean of log-vol means \code{dhs_mean0} (relevant when \code{p > 1})
initDHS = function(omega){

  # "Local" number of time points
  omega = as.matrix(omega)
  n = nrow(omega); p = ncol(omega)

  # Initialize the log-volatilities:
  ht = log(omega^2 + 0.0001)

  # Initialize the AR(1) model to obtain unconditional mean and AR(1) coefficient
  arCoefs = apply(ht, 2, function(x) arima(x, c(1,0,0))$coef)
  dhs_mean = arCoefs[2,]; dhs_phi = arCoefs[1,]; dhs_mean0 = mean(dhs_mean)

  # Initialize the SD of log-vol innovations simply using the expectation:
  sigma_eta_t = matrix(pi, nr = n-1, nc = p)
  sigma_eta_0 = rep(pi, p) # Initial value

  # Evolution error SD:
  sigma_wt = exp(ht/2)

  list(sigma_wt = sigma_wt, ht = ht, dhs_mean = dhs_mean, dhs_phi = dhs_phi, sigma_eta_t = sigma_eta_t, sigma_eta_0 = sigma_eta_0, dhs_mean0 = dhs_mean0)
}
#----------------------------------------------------------------------------
#' Initialize the parameters for the initial state variance
#'
#' The initial state SDs are assumed to follow standard half-Cauchy priors, C+(0,1).
#' This function initalizes the parameters for a PX-Gibbs sampler.
#'
#' @param mu0 \code{p x 1} vector of initial values (undifferenced)
#' @return List of relevant components: the \code{p x 1} evolution error SD \code{sigma_w0},
#' the \code{p x 1} precisions \code{tau_j0}, and the \code{p x 1} parameter-expanded RV's \code{xi_j0}
initEvol0 = function(mu0){
  p = length(mu0)

  # Local precisions:
  tau_j0 = 1/mu0^2; xi_j0 = 1/(tau_j0 + 1)

  list(sigma_w0 = 1/sqrt(tau_j0), tau_j0 = tau_j0, xi_j0 = xi_j0)
}
#----------------------------------------------------------------------------
#' Compute X'X
#'
#' Build the \code{Tp x Tp} matrix XtX using the Matrix() package
#' @param X \code{T x p} matrix of predictors
#' @return Block diagonal \code{Tp x Tp} Matrix (object) where each \code{p x p} block is \code{tcrossprod(matrix(X[t,]))}
#'
#' @note X'X is a one-time computing cost. Special cases may have more efficient computing options,
#' but the Matrix representation is important for efficient computations within the sampler.
#'
#' @import Matrix
#' @export
build_XtX = function(X){

  # Store the dimensions:
  T = nrow(X); p = ncol(X)

  # Store the matrix
  XtX = bandSparse(T*p, k = 0, diag = list(rep(1,T*p)), symm = TRUE)

  t.seq.p = seq(1, T*(p+1), by = p)

  for(t in 1:T){
    t.ind = t.seq.p[t]:(t.seq.p[t+1]-1)
    XtX[t.ind, t.ind] = tcrossprod(matrix(X[t,]))
  }
  XtX
}
#----------------------------------------------------------------------------
# Find (1-alpha)% credible BANDS for a function based on MCMC samples using Crainiceanu et al. (2007)
# sampFuns: (Nsims x m) matrix of Nsims MCMC samples
# alpha: confidence level
#' @export
credBands = function(sampFuns, alpha = .05){

  N = nrow(sampFuns); m = ncol(sampFuns)

  # Compute pointwise mean and SD of f(x):
  Efx = colMeans(sampFuns); SDfx = apply(sampFuns, 2, sd)

  # Compute standardized absolute deviation:
  Standfx = abs(sampFuns - tcrossprod(rep(1, N), Efx))/tcrossprod(rep(1, N), SDfx)

  # And the maximum:
  Maxfx = apply(Standfx, 1, max)

  # Compute the (1-alpha) sample quantile:
  Malpha = quantile(Maxfx, 1-alpha)

  # Finally, store the bands in a (m x 2) matrix of (lower, upper)
  cbind(Efx - Malpha*SDfx, Efx + Malpha*SDfx)
}
#----------------------------------------------------------------------------
#' Estimate the remaining time in the MCMC based on previous samples
#' @param nsi Current iteration
#' @param timer0 Initial timer value, returned from \code{proc.time()[3]}
#' @param nsims Total number of simulations
#' @param nrep Print the estimated time remaining every \code{nrep} iterations
#' @return Table of summary statistics using the function \code{summary}
computeTimeRemaining = function(nsi, timer0, nsims, nrep=1000){

  # Only print occasionally:
  if(nsi%%nrep == 0 || nsi==20) {
    # Current time:
    timer = proc.time()[3]

    # Simulations per second:
    simsPerSec = nsi/(timer - timer0)

    # Seconds remaining, based on extrapolation:
    secRemaining = (nsims - nsi -1)/simsPerSec

    # Print the results:
    if(secRemaining > 3600) {
      print(paste(round(secRemaining/3600, 1), "hours remaining"))
    } else {
      if(secRemaining > 60) {
        print(paste(round(secRemaining/60, 2), "minutes remaining"))
      } else print(paste(round(secRemaining, 2), "seconds remaining"))
    }
  }
}
#----------------------------------------------------------------------------
#' Summarize of effective sample size
#'
#' Compute the summary statistics for the effective sample size (ESS) across
#' posterior samples for possibly many variables
#'
#' @param postX An array of arbitrary dimension \code{(nsims x ... x ...)}, where \code{nsims} is the number of posterior samples
#' @return Table of summary statistics using the function \code{summary()}.
#'
#' @examples
#' # ESS for iid simulations:
#' rand_iid = rnorm(n = 10^4)
#' getEffSize(rand_iid)
#'
#' # ESS for several AR(1) simulations with coefficients 0.1, 0.2,...,0.9:
#' rand_ar1 = sapply(seq(0.1, 0.9, by = 0.1), function(x) arima.sim(n = 10^4, list(ar = x)))
#' getEffSize(rand_ar1)
#'
#' @import coda
#' @export
getEffSize = function(postX) {
  if(is.null(dim(postX))) return(effectiveSize(postX))
  summary(effectiveSize(as.mcmc(array(postX, c(dim(postX)[1], prod(dim(postX)[-1]))))))
}
#----------------------------------------------------------------------------
#' Compute the ergodic (running) mean.
#' @param x vector for which to compute the running mean
#' @return A vector \code{y} with each element defined by \code{y[i] = mean(x[1:i])}
#' @examples
#' # Compare:
#' ergMean(1:10)
#' mean(1:10)
#'
#'# Running mean for iid N(5, 1) samples:
#' x = rnorm(n = 10^4, mean = 5, sd = 1)
#' plot(ergMean(x))
#' abline(h=5)
#' @export
ergMean = function(x) {cumsum(x)/(1:length(x))}
#----------------------------------------------------------------------------
#' Compute the log-odds
#' @param x scalar or vector in (0,1) for which to compute the (componentwise) log-odds
#' @return A scalar or vector of log-odds
#' @examples
#' x = seq(0, 1, length.out = 10^3)
#' plot(x, logit(x))
#' @export
logit = function(x) {
  if(any(abs(x) > 1)) stop('x must be in (0,1)')
  log(x/(1-x))
}
#----------------------------------------------------------------------------
#' Compute the inverse log-odds
#' @param x scalar or vector for which to compute the (componentwise) inverse log-odds
#' @return A scalar or vector of values in (0,1)
#' @examples
#' x = seq(-5, 5, length.out = 10^3)
#' plot(x, invlogit(x))
#' @export
invlogit = function(x) exp(x - log(1+exp(x))) # exp(x)/(1+exp(x))

#----------------------------------------------------------------------------
#' Plot the Bayesian trend filtering fitted values
#'
#' Plot the BTF posterior means with posterior credible intervals (pointwise and joint),
#' the observed data, and true curves (if known)
#'
#' @param y the \code{T x 1} vector of time series observations
#' @param mu the \code{T x 1} vector of fitted values, i.e., posterior expectation of the state variables
#' @param postY the \code{nsims x T} matrix of posterior draws from which to compute intervals
#' @param y_true the \code{T x 1} vector of points along the true curve
#' @param t01 the observation points; if NULL, assume \code{T} equally spaced points from 0 to 1

#' @examples
#' simdata = simUnivariate(signalName = "doppler", T = 128, RSNR = 7, include_plot = FALSE)
#' y = simdata$y
#' mcmc_output = btf(y)
#' plot_fitted(y, mu = colMeans(mcmc_output$mu), postY = mcmc_output$yhat, y_true = simdata$y_true)
#' @import coda
#' @export
plot_fitted = function(y, mu, postY, y_true = NULL, t01 = NULL){
  T = length(y);
  if(is.null(t01)) t01 = seq(0, 1, length.out=T)

  dev.new(); par(mfrow=c(1,1), mai = c(1,1,1,1))
  dcip = HPDinterval(as.mcmc(postY)); dcib = credBands(postY)
  plot(t01, y, type='n', ylim=range(dcib, y, na.rm=TRUE), xlab = 't', ylab=expression(paste("y"[t])), main = 'Fitted Values: Conditional Expectation', cex.lab = 2, cex.main = 2, cex.axis = 2)
  polygon(c(t01, rev(t01)), c(dcib[,2], rev(dcib[,1])), col='gray50', border=NA)
  polygon(c(t01, rev(t01)), c(dcip[,2], rev(dcip[,1])), col='grey', border=NA)
  if(!is.null(y_true))  lines(t01, y_true, lwd=8, col='black', lty=6);
  lines(t01, y, type='p');
  lines(t01, mu, lwd=8, col = 'cyan');
}
#----------------------------------------------------------------------------
#' Univariate Slice Sampler from Neal (2008)
#'
#' Compute a draw from a univariate distribution using the code provided by
#' Radford M. Neal. The documentation below is also reproduced from Neal (2008).
#'
#' @param x0    Initial point
#' @param g     Function returning the log of the probability density (plus constant)
#' @param w     Size of the steps for creating interval (default 1)
#' @param m     Limit on steps (default infinite)
#' @param lower Lower bound on support of the distribution (default -Inf)
#' @param upper Upper bound on support of the distribution (default +Inf)
#' @param gx0   Value of g(x0), if known (default is not known)
#'
#' @return  The point sampled, with its log density attached as an attribute.
#'
#' @note The log density function may return -Inf for points outside the support
#' of the distribution.  If a lower and/or upper bound is specified for the
#' support, the log density function will not be called outside such limits.
uni.slice <- function (x0, g, w=1, m=Inf, lower=-Inf, upper=+Inf, gx0=NULL)
{
  # Check the validity of the arguments.

  if (!is.numeric(x0) || length(x0)!=1
      || !is.function(g)
      || !is.numeric(w) || length(w)!=1 || w<=0
      || !is.numeric(m) || !is.infinite(m) && (m<=0 || m>1e9 || floor(m)!=m)
      || !is.numeric(lower) || length(lower)!=1 || x0<lower
      || !is.numeric(upper) || length(upper)!=1 || x0>upper
      || upper<=lower
      || !is.null(gx0) && (!is.numeric(gx0) || length(gx0)!=1))
  {
    stop ("Invalid slice sampling argument")
  }

  # Keep track of the number of calls made to this function.
  #uni.slice.calls <<- uni.slice.calls + 1

  # Find the log density at the initial point, if not already known.

  if (is.null(gx0))
  { #uni.slice.evals <<- uni.slice.evals + 1
  gx0 <- g(x0)
  }

  # Determine the slice level, in log terms.

  logy <- gx0 - rexp(1)

  # Find the initial interval to sample from.

  u <- runif(1,0,w)
  L <- x0 - u
  R <- x0 + (w-u)  # should guarantee that x0 is in [L,R], even with roundoff

  # Expand the interval until its ends are outside the slice, or until
  # the limit on steps is reached.

  if (is.infinite(m))  # no limit on number of steps
  {
    repeat
    { if (L<=lower) break
      #uni.slice.evals <<- uni.slice.evals + 1
      if (g(L)<=logy) break
      L <- L - w
    }

    repeat
    { if (R>=upper) break
      #uni.slice.evals <<- uni.slice.evals + 1
      if (g(R)<=logy) break
      R <- R + w
    }
  }

  else if (m>1)  # limit on steps, bigger than one
  {
    J <- floor(runif(1,0,m))
    K <- (m-1) - J

    while (J>0)
    { if (L<=lower) break
      #uni.slice.evals <<- uni.slice.evals + 1
      if (g(L)<=logy) break
      L <- L - w
      J <- J - 1
    }

    while (K>0)
    { if (R>=upper) break
      #uni.slice.evals <<- uni.slice.evals + 1
      if (g(R)<=logy) break
      R <- R + w
      K <- K - 1
    }
  }

  # Shrink interval to lower and upper bounds.

  if (L<lower)
  { L <- lower
  }
  if (R>upper)
  { R <- upper
  }

  # Sample from the interval, shrinking it on each rejection.

  repeat
  {
    x1 <- runif(1,L,R)

    #uni.slice.evals <<- uni.slice.evals + 1
    gx1 <- g(x1)

    if (gx1>=logy) break

    if (x1>x0)
    { R <- x1
    }
    else
    { L <- x1
    }
  }

  # Return the point sampled, with its log density attached as an attribute.

  attr(x1,"log.density") <- gx1
  return (x1)

}
