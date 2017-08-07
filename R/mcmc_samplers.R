#' MCMC Sampler for Bayesian Trend Filtering
#'
#' Run the MCMC for Bayesian trend filtering with a penalty on first (D=1) or second (D=2)
#' differences of the conditional expectation.
#' The penalty is determined by the prior on the evolution errors, which include:
#' \itemize{
#' \item the dynamic horseshoe prior ('DHS');
#' \item the static horseshoe prior ('HS');
#' \item the normal-inverse-gamma prior ('NIG').
#' }
#' In each case, the evolution error is a scale mixture of Gaussians.
#' Sampling is accomplished with a (parameter-expanded) Gibbs sampler,
#' mostly relying on a dynamic linear model representation.

#' @param y the \code{T x 1} vector of time series observations
#' @param evol_error the evolution error distribution; must be one of
#' 'DHS' (dynamic horseshoe prior), 'HS' (horseshoe prior), or 'NIG' (normal-inverse-gamma prior)
#' @param D degree of differencing (D = 1 or D = 2)
#' @param nsims number of MCMC iterations to record
#' @param burnin length of the burnin, which is discarded
#' @param mcmc_params named list of parameters for which we store the MCMC output;
#' must be one or more of:
#' \itemize{
#' \item "mu" (conditional mean)
#' \item "yhat" (posterior predictive distribution)
#' \item "evol_sigma_t2" (evolution error variance)
#' \item "obs_sigma_t2" (observation error variance)
#' \item "dhs_phi" (DHS AR(1) coefficient)
#' \item "dhs_mean" (DHS AR(1) unconditional mean)
#' }
#'
#' @return A named list of the \code{nsims} MCMC samples for the parameters named in \code{mcmc_params}
#'
#' @note The data \code{y} may contain NAs, which will be treated with a simple imputation scheme
#' via an additional Gibbs sampling step. In general, rescaling \code{y} to have unit standard
#' deviation is recommended to avoid numerical issues.
#'
#' @examples
#' # Example 1: Bumps Data
#' simdata = simUnivariate(signalName = "bumps", T = 128, RSNR = 7, include_plot = TRUE)
#' y = simdata$y
#'
#' mcmc_output = btf(y)
#' plot_fitted(y, mu = colMeans(mcmc_output$mu), postY = mcmc_output$yhat, y_true = simdata$y_true)
#'
#' # Example 2: Doppler Data; longer series, more noise
#' simdata = simUnivariate(signalName = "doppler", T = 500, RSNR = 5, include_plot = TRUE)
#' y = simdata$y
#'
#' mcmc_output = btf(y, evol_error = 'DHS', D = 2,
#'                  mcmc_params = list('mu', 'yhat', 'dhs_phi', 'dhs_mean'))
#' plot_fitted(y, mu = colMeans(mcmc_output$mu), postY = mcmc_output$yhat, y_true = simdata$y_true)
#'
#'# And examine the AR(1) parameters for the log-volatility w/ traceplots:
#' plot(as.ts(mcmc_output$dhs_phi)) # AR(1) coefficient
#' plot(as.ts(mcmc_output$dhs_mean)) # Unconditional mean
#'
#'# Example 3: Blocks data (locally constant)
#' simdata = simUnivariate(signalName = "blocks", T = 1000, RSNR = 3, include_plot = TRUE)
#' y = simdata$y
#'
#' mcmc_output = btf(y, D = 1) # try D = 1 to approximate the locally constant behavior
#' plot_fitted(y, mu = colMeans(mcmc_output$mu), postY = mcmc_output$yhat, y_true = simdata$y_true)
#'
#' # Example 4: inductance plethsymography data
#' library(wavethresh); data(ipd);
#' # Subsample to better show the localized features:
#' y = as.numeric(ipd[round(seq(1, 4096, length.out = 512))])
#'
#' mcmc_output = btf(y)
#' plot_fitted(y, mu = colMeans(mcmc_output$mu), postY = mcmc_output$yhat)
#'
#' @export
btf = function(y, evol_error = 'DHS', D = 2, nsims = 10000, burnin = 5000, mcmc_params = list("mu", "yhat")){

  # Time points (in [0,1])
  T = length(y); t01 = seq(0, 1, length.out=T);

  # Begin by checking for missing values, then imputing (for initialization)
  is.missing = which(is.na(y)); any.missing = (length(is.missing) > 0)

  # Store the original (w/ missing) data, then impute the active "data"
  yna = y; y = approxfun(t01, y, rule = 2)(t01)

  # Initial SD (implicitly assumes a constant mean)
  sigma_e = sd(y, na.rm=TRUE); sigma_et = rep(sigma_e, T)

  # Initialize the conditional mean, mu, via sampling:
  mu = sampleBTF(y, obs_sigma_t2 = sigma_et^2, evol_sigma_t2 = 0.01*sigma_et^2, D = D)

  # Compute the evolution errors:
  omega = diff(mu, differences = D)

  # And the initial states:
  mu0 = as.matrix(mu[1:D,])

  # Initialize the evolution error variance paramters:
  evolParams = initEvolParams(omega, evol_error = evol_error)

  # Initial variance parameters:
  evolParams0 = initEvol0(mu0)

  # Store the MCMC output in separate arrays (better computation times)
  mcmc_output = vector('list', length(mcmc_params)); names(mcmc_output) = mcmc_params
  if(!is.na(match('mu', mcmc_params))) post_mu = array(NA, c(nsims, T))
  if(!is.na(match('yhat', mcmc_params))) post_yhat = array(NA, c(nsims, T))
  if(!is.na(match('obs_sigma_t2', mcmc_params))) post_obs_sigma_t2 = array(NA, c(nsims, T))
  if(!is.na(match('evol_sigma_t2', mcmc_params))) post_evol_sigma_t2 = array(NA, c(nsims, T))
  if(!is.na(match('dhs_phi', mcmc_params)) && evol_error == "DHS") post_dhs_phi = numeric(nsims)
  if(!is.na(match('dhs_mean', mcmc_params)) && evol_error == "DHS") post_dhs_mean = numeric(nsims)

  # Run the MCMC:
  timer0 = proc.time()[3] # For timing the sampler
  for(nsi in 1:(nsims + burnin)){

    # Impute missing values, if any:
    if(any.missing) y[is.missing] = mu[is.missing] + sigma_et[is.missing]*rnorm(length(is.missing))

    # Sample the states:
    mu = sampleBTF(y, obs_sigma_t2 = sigma_et^2, evol_sigma_t2 = c(evolParams0$sigma_w0^2, evolParams$sigma_wt^2), D = D)

    # Compute the evolution errors:
    omega = diff(mu, differences = D)

    # And the initial states:
    mu0 = as.matrix(mu[1:D,])

    # Sample the evolution error variance (and associated parameters):
    evolParams = sampleEvolParams(omega, evolParams, sigma_e/sqrt(T), evol_error)

    # Sample the observation error SD:
    if(evol_error == 'DHS') {
      sigma_e = uni.slice(sigma_e, g = function(x){
        -(T+2)*log(x) - 0.5*sum((y - mu)^2, na.rm=TRUE)/x^2 - log(1 + (sqrt(T)*exp(evolParams$dhs_mean0/2)/x)^2)
        }, lower = 0, upper = Inf)[1]
    }
    if(evol_error == 'HS') sigma_e = 1/sqrt(rgamma(n = 1, shape = T/2 + length(evolParams$xiLambda), rate = sum((y - mu)^2, na.rm=TRUE)/2 + T*sum(evolParams$xiLambda)))
    if(evol_error == 'IG') sigma_e = 1/sqrt(rgamma(n = 1, shape = T/2, rate = sum((y - mu)^2, na.rm=TRUE)/2))

    # Replicate for coding convenience:
    sigma_et = rep(sigma_e, T)

    # Sample the initial variance parameters:
    evolParams0 = sampleEvol0(mu0, evolParams0, A = 1)

    # Store the MCMC output:
    if(nsi > burnin){
      if(!is.na(match('mu', mcmc_params))) post_mu[nsi - burnin,] = mu
      if(!is.na(match('yhat', mcmc_params))) post_yhat[nsi - burnin,] = mu + sigma_et*rnorm(T)
      if(!is.na(match('obs_sigma_t2', mcmc_params))) post_obs_sigma_t2[nsi - burnin,] = sigma_et^2
      if(!is.na(match('evol_sigma_t2', mcmc_params))) post_evol_sigma_t2[nsi - burnin,] = c(evolParams0$sigma_w0^2, evolParams$sigma_wt^2)
      if(!is.na(match('dhs_phi', mcmc_params)) && evol_error == "DHS") post_dhs_phi[nsi - burnin] = evolParams$dhs_phi
      if(!is.na(match('dhs_mean', mcmc_params)) && evol_error == "DHS") post_dhs_mean[nsi - burnin] = evolParams$dhs_mean
    }
    computeTimeRemaining(nsi, timer0, (nsims + burnin), nrep = 1000)
  }

  if(!is.na(match('mu', mcmc_params))) mcmc_output$mu = post_mu
  if(!is.na(match('yhat', mcmc_params))) mcmc_output$yhat = post_yhat
  if(!is.na(match('obs_sigma_t2', mcmc_params))) mcmc_output$obs_sigma_t2 = post_obs_sigma_t2
  if(!is.na(match('evol_sigma_t2', mcmc_params))) mcmc_output$evol_sigma_t2 = post_evol_sigma_t2
  if(!is.na(match('dhs_phi', mcmc_params)) && evol_error == "DHS") mcmc_output$dhs_phi = post_dhs_phi
  if(!is.na(match('dhs_mean', mcmc_params)) && evol_error == "DHS") mcmc_output$dhs_mean = post_dhs_mean

  print(paste('Total time: ', round((proc.time()[3] - timer0)/60), 'minutes'))

  return (mcmc_output);
}
#----------------------------------------------------------------------------
#' MCMC Sampler for Bayesian Trend Filtering: Regression
#'
#' Run the MCMC for Bayesian trend filtering regression with a penalty on
#' first (D=1) or second (D=2) differences of each dynamic regression coefficient.
#' The penalty is determined by the prior on the evolution errors, which include:
#' \itemize{
#' \item the dynamic horseshoe prior ('DHS');
#' \item the static horseshoe prior ('HS');
#' \item the normal-inverse-gamma prior ('NIG').
#' }
#' In each case, the evolution error is a scale mixture of Gaussians.
#' Sampling is accomplished with a (parameter-expanded) Gibbs sampler,
#' mostly relying on a dynamic linear model representation.
#'
#' @param y the \code{T x 1} vector of time series observations
#' @param X the \code{T x p} matrix of time series predictors
#' @param evol_error the evolution error distribution; must be one of
#' 'DHS' (dynamic horseshoe prior), 'HS' (horseshoe prior), or 'NIG' (normal-inverse-gamma prior)
#' @param D degree of differencing (D = 1 or D = 2)
#' @param nsims number of MCMC iterations to record
#' @param burnin length of the burnin, which is discarded
#' @param mcmc_params named list of parameters for which we store the MCMC output;
#' must be one or more of:
#' \itemize{
#' \item "mu" (conditional mean)
#' \item "yhat" (posterior predictive distribution)
#' \item "beta" (dynamic regression coefficients)
#' \item "evol_sigma_t2" (evolution error variance)
#' \item "obs_sigma_t2" (observation error variance)
#' \item "dhs_phi" (DHS AR(1) coefficient)
#' \item "dhs_mean" (DHS AR(1) unconditional mean)
#' }
#'
#' @return A named list of the \code{nsims} MCMC samples for the parameters named in \code{mcmc_params}
#'
#' @note The data \code{y} may contain NAs, which will be treated with a simple imputation scheme
#' via an additional Gibbs sampling step. In general, rescaling \code{y} to have unit standard
#' deviation is recommended to avoid numerical issues.
#'
#' @examples
#' # Example 1: levelshift and doppler regression
#' simdata = simRegression(signalNames = c("levelshift", "doppler"), p_0 = 2)
#' y = simdata$y; X = simdata$X
#' mcmc_output = btf_reg(y, X)
#' for(j in 1:ncol(X))
#'  plot_fitted(rep(0, length(y)),
#'              mu = colMeans(mcmc_output$beta[,,j]),
#'              postY = mcmc_output$beta[,,j],
#'              y_true = simdata$beta_true[,j])
#'
#' # Example 2: jumpsine and blocks; longer time series, more zeros
#' simdata = simRegression(signalNames = c("jumpsine", "blocks"), p_0 = 5, T = 500)
#' y = simdata$y; X = simdata$X
#' mcmc_output = btf_reg(y, X, nsims = 1000, burnin = 1000) # Short MCMC run for a quick example
#' for(j in 1:ncol(X))
#'   plot_fitted(rep(0, length(y)),
#'               mu = colMeans(mcmc_output$beta[,,j]),
#'               postY = mcmc_output$beta[,,j],
#'               y_true = simdata$beta_true[,j])
#'
#' @export
btf_reg = function(y, X = NULL, evol_error = 'DHS', D = 2, nsims = 10000, burnin = 5000, mcmc_params = list("mu", "yhat", "beta")){

  # If no predictors are specified, just call the univariate trend filtering model: btf()
  if(is.null(X)) {
    mcmc_params[match("beta", mcmc_params)] = NULL # Remove "beta", since it does not apply for btf()
    return(btf(y, evol_error = evol_error, D = D, nsims = nsims, burnin = burnin, mcmc_params = mcmc_params))
  }

  # Time points (in [0,1])
  T = length(y); t01 = seq(0, 1, length.out=T);

  # Begin by checking for missing values, then imputing (for initialization)
  is.missing = which(is.na(y)); any.missing = (length(is.missing) > 0)

  # Store the original (w/ missing) data, then impute the active "data"
  yna = y; if(any.missing) y[is.missing] = mean(y, na.rm=TRUE)

  # Obtain XtX:
  XtX = build_XtX(X);

  # Number of predictors:
  p = ncol(X)

  # Initial SD (implicitly assumes a constant mean)
  sigma_e = sd(y, na.rm=TRUE); sigma_et = rep(sigma_e, T)

  # Initialize the dynamic regression coefficients, beta, via sampling:
  beta = sampleBTF_reg(y, X, obs_sigma_t2 = sigma_et^2, evol_sigma_t2 = matrix(0.01*sigma_et^2, nr = T, nc = p), XtX = XtX, D = D)

  # Conditional mean:
  mu = rowSums(X*beta)

  # Compute the evolution errors:
  omega = diff(beta, differences = D)

  # And the initial states:
  beta0 = matrix(beta[1:D,], nr = D)

  # Initialize the evolution error variance paramters:
  evolParams = initEvolParams(omega, evol_error = evol_error)

  # Initial variance parameters:
  evolParams0 = initEvol0(beta0)

  # Store the MCMC output in separate arrays (better computation times)
  mcmc_output = vector('list', length(mcmc_params)); names(mcmc_output) = mcmc_params
  if(!is.na(match('mu', mcmc_params))) post_mu = array(NA, c(nsims, T))
  if(!is.na(match('yhat', mcmc_params))) post_yhat = array(NA, c(nsims, T))
  if(!is.na(match('beta', mcmc_params))) post_beta = array(NA, c(nsims, T, p))
  if(!is.na(match('obs_sigma_t2', mcmc_params))) post_obs_sigma_t2 = array(NA, c(nsims, T))
  if(!is.na(match('evol_sigma_t2', mcmc_params))) post_evol_sigma_t2 = array(NA, c(nsims, T, p))
  if(!is.na(match('dhs_phi', mcmc_params)) && evol_error == "DHS") post_dhs_phi = array(NA, c(nsims, p))
  if(!is.na(match('dhs_mean', mcmc_params)) && evol_error == "DHS") post_dhs_mean = array(NA, c(nsims, p))

  # Run the MCMC:
  timer0 = proc.time()[3] # For timing the sampler
  for(nsi in 1:(nsims + burnin)){

    # Impute missing values, if any:
    if(any.missing) y[is.missing] = mu[is.missing] + sigma_et[is.missing]*rnorm(length(is.missing))

    # Sample the dynamic regression coefficients, beta:
    beta = sampleBTF_reg(y, X, obs_sigma_t2 = sigma_et^2, evol_sigma_t2 = rbind(matrix(evolParams0$sigma_w0^2, nr = D), evolParams$sigma_wt^2), XtX = XtX, D = D)

    # Conditional mean:
    mu = rowSums(X*beta)

    # Compute the evolution errors:
    omega = diff(beta, differences = D)

    # And the initial states:
    beta0 = matrix(beta[1:D,], nr = D)

    # Sample the evolution error variance (and associated parameters):
    evolParams = sampleEvolParams(omega, evolParams, sigma_e/sqrt(T*p), evol_error)

    # Sample the observation error SD:
    if(evol_error == 'DHS') {
      sigma_e = uni.slice(sigma_e, g = function(x){
        -(T+2)*log(x) - 0.5*sum((y - mu)^2, na.rm=TRUE)/x^2 - log(1 + (sqrt(T*p)*exp(evolParams$dhs_mean0/2)/x)^2)
      }, lower = 0, upper = Inf)[1]
    }
    if(evol_error == 'HS') sigma_e = 1/sqrt(rgamma(n = 1, shape = T/2 + length(evolParams$xiLambda), rate = sum((y - mu)^2, na.rm=TRUE)/2 + T*p*sum(evolParams$xiLambda)))
    if(evol_error == 'IG') sigma_e = 1/sqrt(rgamma(n = 1, shape = T/2, rate = sum((y - mu)^2, na.rm=TRUE)/2))

    # Replicate for coding convenience:
    sigma_et = rep(sigma_e, T)

    # Sample the initial variance parameters:
    evolParams0 = sampleEvol0(beta0, evolParams0, A = 1)

    # Store the MCMC output:
    if(nsi > burnin){
      if(!is.na(match('mu', mcmc_params))) post_mu[nsi - burnin,] = mu
      if(!is.na(match('yhat', mcmc_params))) post_yhat[nsi - burnin,] = mu + sigma_et*rnorm(T)
      if(!is.na(match('beta', mcmc_params))) post_beta[nsi - burnin,,] = beta
      if(!is.na(match('obs_sigma_t2', mcmc_params))) post_obs_sigma_t2[nsi - burnin,] = sigma_et^2
      if(!is.na(match('evol_sigma_t2', mcmc_params))) post_evol_sigma_t2[nsi - burnin,,] = rbind(matrix(evolParams0$sigma_w0^2, nr = D), evolParams$sigma_wt^2)
      if(!is.na(match('dhs_phi', mcmc_params)) && evol_error == "DHS") post_dhs_phi[nsi - burnin,] = evolParams$dhs_phi
      if(!is.na(match('dhs_mean', mcmc_params)) && evol_error == "DHS") post_dhs_mean[nsi - burnin,] = evolParams$dhs_mean
    }
    computeTimeRemaining(nsi, timer0, (nsims + burnin), nrep = 1000)
  }

  if(!is.na(match('mu', mcmc_params))) mcmc_output$mu = post_mu
  if(!is.na(match('yhat', mcmc_params))) mcmc_output$yhat = post_yhat
  if(!is.na(match('beta', mcmc_params))) mcmc_output$beta = post_beta
  if(!is.na(match('obs_sigma_t2', mcmc_params))) mcmc_output$obs_sigma_t2 = post_obs_sigma_t2
  if(!is.na(match('evol_sigma_t2', mcmc_params))) mcmc_output$evol_sigma_t2 = post_evol_sigma_t2
  if(!is.na(match('dhs_phi', mcmc_params)) && evol_error == "DHS") mcmc_output$dhs_phi = post_dhs_phi
  if(!is.na(match('dhs_mean', mcmc_params)) && evol_error == "DHS") mcmc_output$dhs_mean = post_dhs_mean

  print(paste('Total time: ', round((proc.time()[3] - timer0)/60), 'minutes'))

  return (mcmc_output);
}
#----------------------------------------------------------------------------
#' MCMC Sampler for B-spline Bayesian Trend Filtering
#'
#' Run the MCMC for B-spline fitting with a Bayesian trend filtering model on the
#' coefficients, i.e., a penalty on first (D=1) or second (D=2) differences of the
#' B-spline basis coefficients.
#' The penalty is determined by the prior on the evolution errors, which include:
#' \itemize{
#' \item the dynamic horseshoe prior ('DHS');
#' \item the static horseshoe prior ('HS');
#' \item the normal-inverse-gamma prior ('NIG').
#' }
#' In each case, the evolution error is a scale mixture of Gaussians.
#' Sampling is accomplished with a (parameter-expanded) Gibbs sampler,
#' mostly relying on a dynamic linear model representation.

#' @param y the \code{T x 1} vector of time series observations
#' @param x the \code{T x 1} vector of observation points; if NULL, assume equally spaced
#' @param num_knots the number of knots; if NULL, use the default of \code{max(20, min(ceiling(T/4), 150))}
#' @param evol_error the evolution error distribution; must be one of
#' 'DHS' (dynamic horseshoe prior), 'HS' (horseshoe prior), or 'NIG' (normal-inverse-gamma prior)
#' @param D degree of differencing (D = 1 or D = 2)
#' @param nsims number of MCMC iterations to record
#' @param burnin length of the burnin, which is discarded
#' @param mcmc_params named list of parameters for which we store the MCMC output;
#' must be one or more of:
#' \itemize{
#' \item "mu" (conditional mean)
#' \item "beta" (B-spline basis coefficients)
#' \item "yhat" (posterior predictive distribution)
#' \item "evol_sigma_t2" (evolution error variance)
#' \item "obs_sigma_t2" (observation error variance)
#' \item "dhs_phi" (DHS AR(1) coefficient)
#' \item "dhs_mean" (DHS AR(1) unconditional mean)
#' }
#'
#' @return A named list of the \code{nsims} MCMC samples for the parameters named in \code{mcmc_params}
#'
#' @note The data \code{y} may contain NAs, which will be treated with a simple imputation scheme
#' via an additional Gibbs sampling step. In general, rescaling \code{y} to have unit standard
#' deviation is recommended to avoid numerical issues.
#'
#' @note The primary advantages of \code{btf_bspline} over \code{btf} are
#' \enumerate{
#' \item Unequally-spaced points are handled automatically and
#' \item Computations are linear in the number of basis coefficients, which may be
#' substantially fewer than the number of time points.
#' }
#'
#'
#' @examples
#' # Example 1: Blocks data
#' simdata = simUnivariate(signalName = "blocks", T = 1000, RSNR = 3, include_plot = TRUE)
#' y = simdata$y
#' mcmc_output = btf_bspline(y, D = 1)
#' plot_fitted(y, mu = colMeans(mcmc_output$mu), postY = mcmc_output$yhat, y_true = simdata$y_true)
#'
#' # Example 2: motorcycle data (unequally-spaced points)
#' library(MASS)
#' y = scale(mcycle$accel) # Center and Scale for numerical stability
#' x = mcycle$times
#' plot(x, y, xlab = 'Time (ms)', ylab='Acceleration (g)', main = 'Motorcycle Crash Data')
#' mcmc_output = btf_bspline(y = y, x = x)
#' plot_fitted(y, mu = colMeans(mcmc_output$mu), postY = mcmc_output$yhat, t01 = x)
#'
#' # Example 3: inductance plethsymography data (w/o subsampling)
#' library(wavethresh); data(ipd);
#' y = as.numeric(ipd)
#' mcmc_output = btf_bspline(y, num_knots = 500,
#'                          nsims = 1000, burnin = 1000) # Short MCMC run for a quick example
#' plot_fitted(y, mu = colMeans(mcmc_output$mu), postY = mcmc_output$yhat)
#'
#' @import fda
#' @export
btf_bspline = function(y, x = NULL, num_knots = NULL,  evol_error = 'DHS', D = 2, nsims = 10000, burnin = 5000, mcmc_params = list("mu", "yhat")){

  # Length of time series
  T = length(y);

  # Observation points
  if(is.null(x)) x = seq(0, 1, length.out=T);
  # Rescale to (0,1):
  t01 = (x - min(x))/diff(range(x))

  # Compute B-spline basis matrix:
  if(is.null(num_knots)) num_knots = max(20, min(ceiling(T/4), 150))
  X = eval.basis(t01, create.bspline.basis(c(0,1), nbasis = num_knots))
  p = ncol(X)

  # In place of XtX, store the 4-bands of XtX:
  XtX = crossprod(X) #   XtX = Matrix(crossprod(X));
  XtX_bands = list(XtX_0 = diag(XtX),
                   XtX_1 = diag(XtX[-1,]),
                   XtX_2 = diag(XtX[-(1:2),]),
                   XtX_3 = diag(XtX[-(1:3),]))
  rm(XtX)

  # Begin by checking for missing values, then imputing (for initialization)
  is.missing = which(is.na(y)); any.missing = (length(is.missing) > 0)

  # Store the original (w/ missing) data, then impute the active "data"
  yna = y; y = approxfun(t01, y, rule = 2)(t01)

  # Initial SD (implicitly assumes a constant mean)
  sigma_e = sd(y, na.rm=TRUE)

  # Initialize the B-spline coefficients via sampling
  Xty = crossprod(X,y)
  beta = sampleBTF_bspline(y, X, obs_sigma2 = sigma_e^2, evol_sigma_t2 = rep(0.01*sigma_e^2, p), XtX_bands = XtX_bands, Xty = Xty, D = D)

  # And the conditional expectation:
  mu = as.numeric(X%*%beta)

  # Compute the evolution errors:
  omega = diff(beta, differences = D)

  # And the initial states:
  beta0 = matrix(beta[1:D,], nr = D)

  # Initialize the evolution error variance paramters:
  evolParams = initEvolParams(omega, evol_error = evol_error)

  # Initial variance parameters:
  evolParams0 = initEvol0(beta0)

  # Store the MCMC output in separate arrays (better computation times)
  mcmc_output = vector('list', length(mcmc_params)); names(mcmc_output) = mcmc_params
  if(!is.na(match('mu', mcmc_params))) post_mu = array(NA, c(nsims, T))
  if(!is.na(match('yhat', mcmc_params))) post_yhat = array(NA, c(nsims, T))
  if(!is.na(match('beta', mcmc_params))) post_beta = array(NA, c(nsims, p))
  if(!is.na(match('obs_sigma_t2', mcmc_params))) post_obs_sigma_t2 = array(NA, c(nsims, T))
  if(!is.na(match('evol_sigma_t2', mcmc_params))) post_evol_sigma_t2 = array(NA, c(nsims, p))
  if(!is.na(match('dhs_phi', mcmc_params)) && evol_error == "DHS") post_dhs_phi = numeric(nsims)
  if(!is.na(match('dhs_mean', mcmc_params)) && evol_error == "DHS") post_dhs_mean = numeric(nsims)

  # Run the MCMC:
  timer0 = proc.time()[3] # For timing the sampler
  for(nsi in 1:(nsims + burnin)){

    # Impute missing values, if any:
    if(any.missing) {
      y[is.missing] = mu[is.missing] + sigma_e*rnorm(length(is.missing))
      Xty = crossprod(X,y) # Recompute when missing values present
    }

    # Sample the dynamic regression coefficients, beta:
    beta = sampleBTF_bspline(y, X, obs_sigma2 = sigma_e^2, evol_sigma_t2 = rbind(matrix(evolParams0$sigma_w0^2, nr = D), evolParams$sigma_wt^2), XtX_bands = XtX_bands, Xty = Xty, D = D)

    # And the conditional expectation:
    mu = as.numeric(X%*%beta)

    # Compute the evolution errors:
    omega = diff(beta, differences = D)

    # And the initial states:
    beta0 = matrix(beta[1:D,], nr = D)

    # Sample the evolution error variance (and associated parameters):
    evolParams = sampleEvolParams(omega, evolParams, sigma_e/sqrt(p), evol_error)

    # Sample the observation error SD:
    if(evol_error == 'DHS') {
      sigma_e = uni.slice(sigma_e, g = function(x){
        -(T+2)*log(x) - 0.5*sum((y - mu)^2, na.rm=TRUE)/x^2 - log(1 + (sqrt(p)*exp(evolParams$dhs_mean0/2)/x)^2)
      }, lower = 0, upper = Inf)[1]
    }
    if(evol_error == 'HS') sigma_e = 1/sqrt(rgamma(n = 1, shape = T/2 + length(evolParams$xiLambda), rate = sum((y - mu)^2, na.rm=TRUE)/2 + p*sum(evolParams$xiLambda)))
    if(evol_error == 'IG') sigma_e = 1/sqrt(rgamma(n = 1, shape = T/2, rate = sum((y - mu)^2, na.rm=TRUE)/2))

    # Sample the initial variance parameters:
    evolParams0 = sampleEvol0(beta0, evolParams0, A = 1)

    # Store the MCMC output:
    if(nsi > burnin){
      if(!is.na(match('mu', mcmc_params))) post_mu[nsi - burnin,] = mu
      if(!is.na(match('yhat', mcmc_params))) post_yhat[nsi - burnin,] = mu + sigma_e*rnorm(T)
      if(!is.na(match('beta', mcmc_params))) post_beta[nsi - burnin,] = beta
      if(!is.na(match('obs_sigma_t2', mcmc_params))) post_obs_sigma_t2[nsi - burnin,] = sigma_et^2
      if(!is.na(match('evol_sigma_t2', mcmc_params))) post_evol_sigma_t2[nsi - burnin,] = evolParams$sigma_wt^2
      if(!is.na(match('dhs_phi', mcmc_params)) && evol_error == "DHS") post_dhs_phi[nsi - burnin] = evolParams$dhs_phi
      if(!is.na(match('dhs_mean', mcmc_params)) && evol_error == "DHS") post_dhs_mean[nsi - burnin] = evolParams$dhs_mean
    }
    computeTimeRemaining(nsi, timer0, (nsims + burnin), nrep = 1000)
  }

  if(!is.na(match('mu', mcmc_params))) mcmc_output$mu = post_mu
  if(!is.na(match('yhat', mcmc_params))) mcmc_output$yhat = post_yhat
  if(!is.na(match('beta', mcmc_params))) mcmc_output$beta = post_beta
  if(!is.na(match('obs_sigma_t2', mcmc_params))) mcmc_output$obs_sigma_t2 = post_obs_sigma_t2
  if(!is.na(match('evol_sigma_t2', mcmc_params))) mcmc_output$evol_sigma_t2 = post_evol_sigma_t2
  if(!is.na(match('dhs_phi', mcmc_params)) && evol_error == "DHS") mcmc_output$dhs_phi = post_dhs_phi
  if(!is.na(match('dhs_mean', mcmc_params)) && evol_error == "DHS") mcmc_output$dhs_mean = post_dhs_mean

  print(paste('Total time: ', round((proc.time()[3] - timer0)/60), 'minutes'))

  return (mcmc_output);
}


