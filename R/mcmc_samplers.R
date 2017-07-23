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
#' Sampling is accomplished with a (parameter-expanded) Gibbs sampler, mostly relying on a dynamic linear model representation.

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
#' plot_fitted(y, mu = colMeans(mcmc_output$mu), postY = mcmc_output$yhat, ytrue = simdata$ytrue)
#'
#' # Example 2: Doppler Data, longer series, more noise
#' simdata = simUnivariate(signalName = "doppler", T = 500, RSNR = 5, include_plot = TRUE)
#' y = simdata$y
#'
#' mcmc_output = btf(y, evol_error = 'DHS', D = 2, mcmc_params = list('mu', 'yhat', 'dhs_phi', 'dhs_mean'))
#' plot_fitted(y, mu = colMeans(mcmc_output$mu), postY = mcmc_output$yhat, ytrue = simdata$ytrue)
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
#' plot_fitted(y, mu = colMeans(mcmc_output$mu), postY = mcmc_output$yhat, ytrue = simdata$ytrue)
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
      if(!is.na(match('evol_sigma_t2', mcmc_params))) post_evol_sigma_t2[nsi - burnin,] = evolParams$sigma_wt^2
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

