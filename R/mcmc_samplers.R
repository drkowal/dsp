#' MCMC Sampler for Bayesian Trend Filtering
#'
#' Run the MCMC for Bayesian trend filtering with a penalty on zeroth (D = 0),
#' first (D = 1), or second (D = 2) differences of the conditional expectation.
#' The penalty is determined by the prior on the evolution errors, which include:
#' \itemize{
#' \item the dynamic horseshoe prior ('DHS');
#' \item the static horseshoe prior ('HS');
#' \item the Bayesian lasso ('BL');
#' \item the normal stochastic volatility model ('SV');
#' \item the normal-inverse-gamma prior ('NIG').
#' }
#' In each case, the evolution error is a scale mixture of Gaussians.
#' Sampling is accomplished with a (parameter-expanded) Gibbs sampler,
#' mostly relying on a dynamic linear model representation.

#' @param y the \code{T x 1} vector of time series observations
#' @param evol_error the evolution error distribution; must be one of
#' 'DHS' (dynamic horseshoe prior), 'HS' (horseshoe prior), 'BL' (Bayesian lasso), or 'NIG' (normal-inverse-gamma prior)
#' @param D degree of differencing (D = 0, D = 1, or D = 2)
#' @param useObsSV logical; if TRUE, include a (normal) stochastic volatility model
#' for the observation error variance
#' @param nsave number of MCMC iterations to record
#' @param nburn number of MCMC iterations to discard (burin-in)
#' @param nskip number of MCMC iterations to skip between saving iterations,
#' i.e., save every (nskip + 1)th draw
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
#' @param computeDIC logical; if TRUE, compute the deviance information criterion \code{DIC}
#' and the effective number of parameters \code{p_d}
#' @param verbose logical; should R report extra information on progress?
#'
#' @return A named list of the \code{nsave} MCMC samples for the parameters named in \code{mcmc_params}
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
#' out = btf(y)
#' plot_fitted(y, mu = colMeans(out$mu), postY = out$yhat, y_true = simdata$y_true)
#'
#' \dontrun{
#' # Example 2: Doppler Data; longer series, more noise
#' simdata = simUnivariate(signalName = "doppler", T = 500, RSNR = 5, include_plot = TRUE)
#' y = simdata$y
#'
#' out = btf(y)
#' plot_fitted(y, mu = colMeans(out$mu), postY = out$yhat, y_true = simdata$y_true)
#'
#' # And examine the AR(1) parameters for the log-volatility w/ traceplots:
#' plot(as.ts(out$dhs_phi)) # AR(1) coefficient
#' plot(as.ts(out$dhs_mean)) # Unconditional mean
#'
#'# Example 3: Blocks data (locally constant)
#' simdata = simUnivariate(signalName = "blocks", T = 1000, RSNR = 3, include_plot = TRUE)
#' y = simdata$y
#'
#' out = btf(y, D = 1) # try D = 1 to approximate the locally constant behavior
#' plot_fitted(y, mu = colMeans(out$mu), postY = out$yhat, y_true = simdata$y_true)
#' }
#'
#' @import spam
#' @export
btf = function(y, evol_error = 'DHS', D = 2, useObsSV = FALSE,
               nsave = 1000, nburn = 1000, nskip = 4,
               mcmc_params = list("mu", "yhat","evol_sigma_t2", "obs_sigma_t2", "dhs_phi", "dhs_mean"),
               computeDIC = TRUE,
               verbose = TRUE){

  # Convert to upper case:
  evol_error = toupper(evol_error)

  # For D = 0, return special case:
  if(D == 0){
    return(btf0(y = y, evol_error = evol_error, useObsSV = useObsSV,
                nsave = nsave, nburn = nburn, nskip = nskip,
                mcmc_params = mcmc_params,
                computeDIC = computeDIC,
                verbose = verbose))
  }

  # Time points (in [0,1])
  T = length(y); t01 = seq(0, 1, length.out=T);

  # Begin by checking for missing values, then imputing (for initialization)
  is.missing = which(is.na(y)); any.missing = (length(is.missing) > 0)

  # Impute the active "data"
  y = approxfun(t01, y, rule = 2)(t01)

  # Initial SD (implicitly assumes a constant mean)
  sigma_e = sd(y, na.rm=TRUE); sigma_et = rep(sigma_e, T)

  # Compute the Cholesky term (uses random variances for a more conservative sparsity pattern)
  chol0 = initChol.spam(T = T, D = D)

  # Initialize the conditional mean, mu, via sampling:
  mu = sampleBTF(y, obs_sigma_t2 = sigma_et^2, evol_sigma_t2 = 0.01*sigma_et^2, D = D, chol0 = chol0)

  # Compute the evolution errors:
  omega = diff(mu, differences = D)

  # And the initial states:
  mu0 = as.matrix(mu[1:D,])

  # Initialize the evolution error variance paramters:
  evolParams = initEvolParams(omega, evol_error = evol_error)

  # Initial variance parameters:
  evolParams0 = initEvol0(mu0)

  # SV parameters, if necessary:
  if(useObsSV) {svParams = initSV(y - mu); sigma_et = svParams$sigma_wt}

  # For HS MCMC comparisons:
  # evolParams$dhs_phi = 0

  # Store the MCMC output in separate arrays (better computation times)
  mcmc_output = vector('list', length(mcmc_params)); names(mcmc_output) = mcmc_params
  if(!is.na(match('mu', mcmc_params)) || computeDIC) post_mu = array(NA, c(nsave, T))
  if(!is.na(match('yhat', mcmc_params))) post_yhat = array(NA, c(nsave, T))
  if(!is.na(match('obs_sigma_t2', mcmc_params)) || computeDIC) post_obs_sigma_t2 = array(NA, c(nsave, T))
  if(!is.na(match('evol_sigma_t2', mcmc_params))) post_evol_sigma_t2 = array(NA, c(nsave, T))
  if(!is.na(match('dhs_phi', mcmc_params)) && evol_error == "DHS") post_dhs_phi = numeric(nsave)
  if(!is.na(match('dhs_mean', mcmc_params)) && evol_error == "DHS") post_dhs_mean = numeric(nsave)
  post_loglike = numeric(nsave)

  # Total number of MCMC simulations:
  nstot = nburn+(nskip+1)*(nsave)
  skipcount = 0; isave = 0 # For counting

  # Run the MCMC:
  if(verbose) timer0 = proc.time()[3] # For timing the sampler
  for(nsi in 1:nstot){

    # Impute missing values, if any:
    if(any.missing) y[is.missing] = mu[is.missing] + sigma_et[is.missing]*rnorm(length(is.missing))

    # Sample the states:
    mu = sampleBTF(y, obs_sigma_t2 = sigma_et^2, evol_sigma_t2 = c(evolParams0$sigma_w0^2, evolParams$sigma_wt^2), D = D, chol0 = chol0)

    # Compute the evolution errors:
    omega = diff(mu, differences = D)

    # And the initial states:
    mu0 = as.matrix(mu[1:D,])

    # Sample the initial variance parameters:
    evolParams0 = sampleEvol0(mu0, evolParams0, A = 1)

    # Sample the (observation and evolution) variances and associated parameters:
    if(useObsSV){
      # Evolution error variance + params:
      evolParams = sampleEvolParams(omega, evolParams, 1/sqrt(T), evol_error)

      # Observation error variance + params:
      svParams = sampleSVparams(omega = y - mu, svParams = svParams)
      sigma_et = svParams$sigma_wt

    } else {
      # Evolution error variance + params:
      evolParams = sampleEvolParams(omega, evolParams, sigma_e/sqrt(T), evol_error)

      if(evol_error == 'DHS') {
        sigma_e = uni.slice(sigma_e, g = function(x){
          -(T+2)*log(x) - 0.5*sum((y - mu)^2, na.rm=TRUE)/x^2 - log(1 + (sqrt(T)*exp(evolParams$dhs_mean0/2)/x)^2)
        }, lower = 0, upper = Inf)[1]
      }
      if(evol_error == 'HS') sigma_e = 1/sqrt(rgamma(n = 1, shape = T/2 + length(evolParams$xiLambda), rate = sum((y - mu)^2, na.rm=TRUE)/2 + T*sum(evolParams$xiLambda)))
      if(evol_error == 'BL') sigma_e = 1/sqrt(rgamma(n = 1, shape = T/2 + length(evolParams$tau_j)/2, rate = sum((y - mu)^2, na.rm=TRUE)/2 + T*sum((omega/evolParams$tau_j)^2)/2))
      if((evol_error == 'NIG') || (evol_error == 'SV')) sigma_e = 1/sqrt(rgamma(n = 1, shape = T/2, rate = sum((y - mu)^2, na.rm=TRUE)/2))

      # Replicate for coding convenience:
      sigma_et = rep(sigma_e, T)
    }

    # Store the MCMC output:
    if(nsi > nburn){
      # Increment the skip counter:
      skipcount = skipcount + 1

      # Save the iteration:
      if(skipcount > nskip){
        # Increment the save index
        isave = isave + 1

        # Save the MCMC samples:
        if(!is.na(match('mu', mcmc_params)) || computeDIC) post_mu[isave,] = mu
        if(!is.na(match('yhat', mcmc_params))) post_yhat[isave,] = mu + sigma_et*rnorm(T)
        if(!is.na(match('obs_sigma_t2', mcmc_params)) || computeDIC) post_obs_sigma_t2[isave,] = sigma_et^2
        if(!is.na(match('evol_sigma_t2', mcmc_params))) post_evol_sigma_t2[isave,] = c(evolParams0$sigma_w0^2, evolParams$sigma_wt^2)
        if(!is.na(match('dhs_phi', mcmc_params)) && evol_error == "DHS") post_dhs_phi[isave] = evolParams$dhs_phi
        if(!is.na(match('dhs_mean', mcmc_params)) && evol_error == "DHS") post_dhs_mean[isave] = evolParams$dhs_mean
        post_loglike[isave] = sum(dnorm(y, mean = mu, sd = sigma_et, log = TRUE))

        # And reset the skip counter:
        skipcount = 0
      }
    }
    if(verbose) computeTimeRemaining(nsi, timer0, nstot, nrep = 1000)
  }

  if(!is.na(match('mu', mcmc_params))) mcmc_output$mu = post_mu
  if(!is.na(match('yhat', mcmc_params))) mcmc_output$yhat = post_yhat
  if(!is.na(match('obs_sigma_t2', mcmc_params))) mcmc_output$obs_sigma_t2 = post_obs_sigma_t2
  if(!is.na(match('evol_sigma_t2', mcmc_params))) mcmc_output$evol_sigma_t2 = post_evol_sigma_t2
  if(!is.na(match('dhs_phi', mcmc_params)) && evol_error == "DHS") mcmc_output$dhs_phi = post_dhs_phi
  if(!is.na(match('dhs_mean', mcmc_params)) && evol_error == "DHS") mcmc_output$dhs_mean = post_dhs_mean

  # Also include the log-likelihood:
  mcmc_output$loglike = post_loglike

  if(computeDIC){
    # Log-likelihood evaluated at posterior means:
    loglike_hat = sum(dnorm(y,
                            mean = colMeans(post_mu),
                            sd = colMeans(sqrt(post_obs_sigma_t2)),
                            log = TRUE))

    # Effective number of parameters (Note: two options)
    p_d = c(2*(loglike_hat - mean(post_loglike)),
            2*var(post_loglike))
    # DIC:
    DIC = -2*loglike_hat + 2*p_d

    # Store the DIC and the effective number of parameters (p_d)
    mcmc_output$DIC = DIC; mcmc_output$p_d = p_d
  }

  if(verbose) print(paste('Total time: ', round((proc.time()[3] - timer0)), 'seconds'))

  return (mcmc_output);
}
#' MCMC Sampler for Bayesian Trend Filtering: D = 0
#'
#' Run the MCMC for Bayesian trend filtering with a penalty on the conditional expectation.
#' The penalty is determined by the prior on the evolution errors, which include:
#' \itemize{
#' \item the dynamic horseshoe prior ('DHS');
#' \item the static horseshoe prior ('HS');
#' \item the Bayesian lasso ('BL');
#' \item the normal stochastic volatility model ('SV');
#' \item the normal-inverse-gamma prior ('NIG').
#' }
#' In each case, the evolution error is a scale mixture of Gaussians.
#' Sampling is accomplished with a (parameter-expanded) Gibbs sampler,
#' mostly relying on a dynamic linear model representation.

#' @param y the \code{T x 1} vector of time series observations
#' @param evol_error the evolution error distribution; must be one of
#' 'DHS' (dynamic horseshoe prior), 'HS' (horseshoe prior), 'BL' (Bayesian lasso), or 'NIG' (normal-inverse-gamma prior)
#' @param useObsSV logical; if TRUE, include a (normal) stochastic volatility model
#' for the observation error variance
#' @param nsave number of MCMC iterations to record
#' @param nburn number of MCMC iterations to discard (burin-in)
#' @param nskip number of MCMC iterations to skip between saving iterations,
#' i.e., save every (nskip + 1)th draw
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
#' @param computeDIC logical; if TRUE, compute the deviance information criterion \code{DIC}
#' and the effective number of parameters \code{p_d}
#' @param verbose logical; should R report extra information on progress?
#'
#' @return A named list of the \code{nsave} MCMC samples for the parameters named in \code{mcmc_params}
#'
#' @note The data \code{y} may contain NAs, which will be treated with a simple imputation scheme
#' via an additional Gibbs sampling step. In general, rescaling \code{y} to have unit standard
#' deviation is recommended to avoid numerical issues.
btf0 = function(y, evol_error = 'DHS', useObsSV = FALSE,
                nsave = 1000, nburn = 1000, nskip = 4,
                mcmc_params = list("mu", "yhat","evol_sigma_t2", "obs_sigma_t2", "dhs_phi", "dhs_mean"),
                computeDIC = TRUE,
                verbose = TRUE){

  # Time points (in [0,1])
  T = length(y); t01 = seq(0, 1, length.out=T);

  # Begin by checking for missing values, then imputing (for initialization)
  is.missing = which(is.na(y)); any.missing = (length(is.missing) > 0)

  # Impute the active "data"
  y = approxfun(t01, y, rule = 2)(t01)

  # Initial SD (implicitly assumes a constant mean)
  sigma_e = sd(y, na.rm=TRUE); sigma_et = rep(sigma_e, T)

  # Initialize the conditional mean, mu, via sampling:
  mu = sampleBTF(y, obs_sigma_t2 = sigma_et^2, evol_sigma_t2 = 0.01*sigma_et^2, D = 0)

  # Initialize the evolution error variance paramters:
  evolParams = initEvolParams(omega = mu, evol_error = evol_error)

  # SV parameters, if necessary:
  if(useObsSV) {svParams = initSV(y - mu); sigma_et = svParams$sigma_wt}

  # Store the MCMC output in separate arrays (better computation times)
  mcmc_output = vector('list', length(mcmc_params)); names(mcmc_output) = mcmc_params
  if(!is.na(match('mu', mcmc_params)) || computeDIC) post_mu = array(NA, c(nsave, T))
  if(!is.na(match('yhat', mcmc_params))) post_yhat = array(NA, c(nsave, T))
  if(!is.na(match('obs_sigma_t2', mcmc_params)) || computeDIC) post_obs_sigma_t2 = array(NA, c(nsave, T))
  if(!is.na(match('evol_sigma_t2', mcmc_params))) post_evol_sigma_t2 = array(NA, c(nsave, T))
  if(!is.na(match('dhs_phi', mcmc_params)) && evol_error == "DHS") post_dhs_phi = numeric(nsave)
  if(!is.na(match('dhs_mean', mcmc_params)) && evol_error == "DHS") post_dhs_mean = numeric(nsave)
  post_loglike = numeric(nsave)

  # Total number of MCMC simulations:
  nstot = nburn+(nskip+1)*(nsave)
  skipcount = 0; isave = 0 # For counting

  # Run the MCMC:
  if(verbose) timer0 = proc.time()[3] # For timing the sampler
  for(nsi in 1:nstot){

    # Impute missing values, if any:
    if(any.missing) y[is.missing] = mu[is.missing] + sigma_et[is.missing]*rnorm(length(is.missing))

    # Sample the states:
    mu = sampleBTF(y, obs_sigma_t2 = sigma_et^2, evol_sigma_t2 = evolParams$sigma_wt^2, D = 0)

    # Sample the (observation and evolution) variances and associated parameters:
    if(useObsSV){
      # Evolution error variance + params:
      evolParams = sampleEvolParams(omega = mu, evolParams, 1/sqrt(T), evol_error)

      # Observation error variance + params:
      svParams = sampleSVparams(omega = y - mu, svParams = svParams)
      sigma_et = svParams$sigma_wt

    } else {
      # Evolution error variance + params:
      evolParams = sampleEvolParams(omega = mu, evolParams, sigma_e/sqrt(T), evol_error)

      # Sample the observation error SD:
      if(evol_error == 'DHS') {
        sigma_e = uni.slice(sigma_e, g = function(x){
          -(T+2)*log(x) - 0.5*sum((y - mu)^2, na.rm=TRUE)/x^2 - log(1 + (sqrt(T)*exp(evolParams$dhs_mean0/2)/x)^2)
        }, lower = 0, upper = Inf)[1]
      }
      if(evol_error == 'HS') sigma_e = 1/sqrt(rgamma(n = 1, shape = T/2 + length(evolParams$xiLambda), rate = sum((y - mu)^2, na.rm=TRUE)/2 + T*sum(evolParams$xiLambda)))
      if(evol_error == 'BL') sigma_e = 1/sqrt(rgamma(n = 1, shape = T/2 + length(evolParams$tau_j)/2, rate = sum((y - mu)^2, na.rm=TRUE)/2 + T*sum((mu/evolParams$tau_j)^2)/2))
      if((evol_error == 'NIG') || (evol_error == 'SV')) sigma_e = 1/sqrt(rgamma(n = 1, shape = T/2, rate = sum((y - mu)^2, na.rm=TRUE)/2))

      # Replicate for coding convenience:
      sigma_et = rep(sigma_e, T)
    }

    # Store the MCMC output:
    if(nsi > nburn){
      # Increment the skip counter:
      skipcount = skipcount + 1

      # Save the iteration:
      if(skipcount > nskip){
        # Increment the save index
        isave = isave + 1

        # Save the MCMC samples:
        if(!is.na(match('mu', mcmc_params)) || computeDIC) post_mu[isave,] = mu
        if(!is.na(match('yhat', mcmc_params))) post_yhat[isave,] = mu + sigma_et*rnorm(T)
        if(!is.na(match('obs_sigma_t2', mcmc_params)) || computeDIC) post_obs_sigma_t2[isave,] = sigma_et^2
        if(!is.na(match('evol_sigma_t2', mcmc_params))) post_evol_sigma_t2[isave,] =  evolParams$sigma_wt^2
        if(!is.na(match('dhs_phi', mcmc_params)) && evol_error == "DHS") post_dhs_phi[isave] = evolParams$dhs_phi
        if(!is.na(match('dhs_mean', mcmc_params)) && evol_error == "DHS") post_dhs_mean[isave] = evolParams$dhs_mean
        post_loglike[isave] = sum(dnorm(y, mean = mu, sd = sigma_et, log = TRUE))

        # And reset the skip counter:
        skipcount = 0
      }
    }
    if(verbose) computeTimeRemaining(nsi, timer0, nstot, nrep = 1000)
  }

  if(!is.na(match('mu', mcmc_params))) mcmc_output$mu = post_mu
  if(!is.na(match('yhat', mcmc_params))) mcmc_output$yhat = post_yhat
  if(!is.na(match('obs_sigma_t2', mcmc_params))) mcmc_output$obs_sigma_t2 = post_obs_sigma_t2
  if(!is.na(match('evol_sigma_t2', mcmc_params))) mcmc_output$evol_sigma_t2 = post_evol_sigma_t2
  if(!is.na(match('dhs_phi', mcmc_params)) && evol_error == "DHS") mcmc_output$dhs_phi = post_dhs_phi
  if(!is.na(match('dhs_mean', mcmc_params)) && evol_error == "DHS") mcmc_output$dhs_mean = post_dhs_mean

  # Also include the log-likelihood:
  mcmc_output$loglike = post_loglike

  if(computeDIC){
    # Log-likelihood evaluated at posterior means:
    loglike_hat = sum(dnorm(y,
                            mean = colMeans(post_mu),
                            sd = colMeans(sqrt(post_obs_sigma_t2)),
                            log = TRUE))

    # Effective number of parameters (Note: two options)
    p_d = c(2*(loglike_hat - mean(post_loglike)),
            2*var(post_loglike))
    # DIC:
    DIC = -2*loglike_hat + 2*p_d

    # Store the DIC and the effective number of parameters (p_d)
    mcmc_output$DIC = DIC; mcmc_output$p_d = p_d
  }

  if(verbose) print(paste('Total time: ', round((proc.time()[3] - timer0)), 'seconds'))

  return (mcmc_output);
}
#' Run the MCMC for sparse Bayesian trend filtering
#'
#' Sparse Bayesian trend filtering has two penalties:
#' (1) a penalty on the first (D = 1) or second (D = 2) differences of the conditional expectation and
#' (2) a penalty on the conditional expectation, i.e., shrinkage to zero.
#'
#' Each penalty is determined by a prior, which include:
#' \itemize{
#' \item the dynamic horseshoe prior ('DHS');
#' \item the static horseshoe prior ('HS');
#' \item the Bayesian lasso ('BL');
#' \item the normal stochastic volatility model ('SV');
#' \item the normal-inverse-gamma prior ('NIG').
#' }
#' In each case, the prior is a scale mixture of Gaussians.
#' Sampling is accomplished with a (parameter-expanded) Gibbs sampler,
#' mostly relying on a dynamic linear model representation.

#' @param y the \code{T x 1} vector of time series observations
#' @param evol_error the evolution error distribution; must be one of
#' 'DHS' (dynamic horseshoe prior), 'HS' (horseshoe prior), 'BL' (Bayesian lasso), or 'NIG' (normal-inverse-gamma prior)
#' @param zero_error the shrinkage-to-zero distribution; must be one of
#' 'DHS' (dynamic horseshoe prior), 'HS' (horseshoe prior), 'BL' (Bayesian lasso), or 'NIG' (normal-inverse-gamma prior)
#' @param D degree of differencing (D = 1, or D = 2)
#' @param nsave number of MCMC iterations to record
#' @param nburn number of MCMC iterations to discard (burin-in)
#' @param nskip number of MCMC iterations to skip between saving iterations,
#' i.e., save every (nskip + 1)th draw
#' @param mcmc_params named list of parameters for which we store the MCMC output;
#' must be one or more of:
#' \itemize{
#' \item "mu" (conditional mean)
#' \item "yhat" (posterior predictive distribution)
#' \item "evol_sigma_t2" (evolution error variance)
#' \item "zero_sigma_t2" (shrink-to-zero error variance)
#' \item "obs_sigma_t2" (observation error variance)
#' \item "dhs_phi" (DHS AR(1) coefficient for evolution error)
#' \item "dhs_mean" (DHS AR(1) unconditional mean for evolution error)
#' \item "dhs_phi_zero" (DHS AR(1) coefficient for shrink-to-zero error)
#' \item "dhs_mean_zero" (DHS AR(1) unconditional mean for shrink-to-zero error)
#' }
#' @param computeDIC logical; if TRUE, compute the deviance information criterion \code{DIC}
#' and the effective number of parameters \code{p_d}
#' @param verbose logical; should R report extra information on progress?
#'
#' @return A named list of the \code{nsave} MCMC samples for the parameters named in \code{mcmc_params}
#'
#' @note The data \code{y} may contain NAs, which will be treated with a simple imputation scheme
#' via an additional Gibbs sampling step. In general, rescaling \code{y} to have unit standard
#' deviation is recommended to avoid numerical issues.
#'
#' @examples
#' \dontrun{
#' # Example 1: Bumps Data
#' simdata = simUnivariate(signalName = "bumps", T = 128, RSNR = 7, include_plot = TRUE)
#' y = simdata$y
#'
#' out = btf_sparse(y)
#' plot_fitted(y, mu = colMeans(out$mu), postY = out$yhat, y_true = simdata$y_true)
#'
#' # Example 2: Doppler Data; longer series, more noise
#' simdata = simUnivariate(signalName = "doppler", T = 500, RSNR = 5, include_plot = TRUE)
#' y = simdata$y
#'
#' out = btf_sparse(y)
#' plot_fitted(y, mu = colMeans(out$mu), postY = out$yhat, y_true = simdata$y_true)
#'
#'# And examine the AR(1) parameters for the log-volatility w/ traceplots:
#' plot(as.ts(out$dhs_phi)) # AR(1) coefficient
#' plot(as.ts(out$dhs_mean)) # Unconditional mean
#'
#'# Example 3: Blocks data (locally constant)
#' simdata = simUnivariate(signalName = "blocks", T = 1000, RSNR = 3, include_plot = TRUE)
#' y = simdata$y
#'
#' out = btf_sparse(y, D = 1) # try D = 1 to approximate the locally constant behavior
#' plot_fitted(y, mu = colMeans(out$mu), postY = out$yhat, y_true = simdata$y_true)
#' }
#'
#' @import spam
#' @export
btf_sparse = function(y, evol_error = 'DHS', zero_error = 'DHS', D = 2,
                      nsave = 1000, nburn = 1000, nskip = 4,
                      mcmc_params = list("mu", "yhat","evol_sigma_t2", "obs_sigma_t2", "zero_sigma_t2", "dhs_phi", "dhs_mean","dhs_phi_zero", "dhs_mean_zero"),
                      computeDIC = TRUE,
                      verbose = TRUE){

  # Check differencing:
  if((D != 1) && (D != 2)) stop('btf_sparse() requires D = 1 or D = 2')

  # Convert to upper case:
  evol_error = toupper(evol_error)
  zero_error = toupper(zero_error)

  # Time points (in [0,1])
  T = length(y); t01 = seq(0, 1, length.out=T);

  # Begin by checking for missing values, then imputing (for initialization)
  is.missing = which(is.na(y)); any.missing = (length(is.missing) > 0)

  # Impute the active "data"
  y = approxfun(t01, y, rule = 2)(t01)

  # Initial SD (implicitly assumes a constant mean)
  sigma_e = sd(y, na.rm=TRUE); sigma_et = rep(sigma_e, T)

  # Compute the Cholesky term (uses random variances for a more conservative sparsity pattern)
  chol0 = initChol.spam(T = T, D = D)

  # Initialize the conditional mean, mu, via sampling:
  mu = sampleBTF(y, obs_sigma_t2 = sigma_et^2, evol_sigma_t2 = 0.01*sigma_et^2, D = D, chol0 = chol0)

  # Compute the evolution errors:
  omega = diff(mu, differences = D)

  # And the initial states:
  mu0 = as.matrix(mu[1:D,])

  # Initialize the evolution error variance paramters:
  evolParams = initEvolParams(omega, evol_error = evol_error)

  # Initial variance parameters:
  evolParams0 = initEvol0(mu0)

  # Initialize the shrink-to-zero parameters:
  zeroParams = initEvolParams(omega = mu, evol_error = zero_error)

  # For HS MCMC comparisons:
  # evolParams$dhs_phi = 0

  # Store the MCMC output in separate arrays (better computation times)
  mcmc_output = vector('list', length(mcmc_params)); names(mcmc_output) = mcmc_params
  if(!is.na(match('mu', mcmc_params)) || computeDIC) post_mu = array(NA, c(nsave, T))
  if(!is.na(match('yhat', mcmc_params))) post_yhat = array(NA, c(nsave, T))
  if(!is.na(match('obs_sigma_t2', mcmc_params)) || computeDIC) post_obs_sigma_t2 = array(NA, c(nsave, T))
  if(!is.na(match('evol_sigma_t2', mcmc_params))) post_evol_sigma_t2 = array(NA, c(nsave, T))
  if(!is.na(match('zero_sigma_t2', mcmc_params))) post_zero_sigma_t2 = array(NA, c(nsave, T))
  if(!is.na(match('dhs_phi', mcmc_params)) && evol_error == "DHS") post_dhs_phi = numeric(nsave)
  if(!is.na(match('dhs_mean', mcmc_params)) && evol_error == "DHS") post_dhs_mean = numeric(nsave)
  if(!is.na(match('dhs_phi_zero', mcmc_params)) && zero_error == "DHS") post_dhs_phi_zero = numeric(nsave)
  if(!is.na(match('dhs_mean_zero', mcmc_params)) && zero_error == "DHS") post_dhs_mean_zero = numeric(nsave)
  post_loglike = numeric(nsave)

  # Total number of MCMC simulations:
  nstot = nburn+(nskip+1)*(nsave)
  skipcount = 0; isave = 0 # For counting

  # Run the MCMC:
  if(verbose) timer0 = proc.time()[3] # For timing the sampler
  for(nsi in 1:nstot){

    # Impute missing values, if any:
    if(any.missing) y[is.missing] = mu[is.missing] + sigma_et[is.missing]*rnorm(length(is.missing))

    # Sample the observation SD:
    sigma_e = 1/sqrt(rgamma(n = 1, shape = T/2, rate = sum((y - mu)^2, na.rm=TRUE)/2))
    sigma_et = rep(sigma_e, T) # For coding convenience

    # Sample the states:
    #mu = sampleBTF(y, obs_sigma_t2 = sigma_et^2, evol_sigma_t2 = c(evolParams0$sigma_w0^2, evolParams$sigma_wt^2), D = D, chol0 = chol0)
    mu = sampleBTF_sparse(y,
                          obs_sigma_t2 = sigma_et^2,
                          evol_sigma_t2 = c(evolParams0$sigma_w0^2, evolParams$sigma_wt^2),
                          zero_sigma_t2 = zeroParams$sigma_wt^2,
                          D = D, chol0 = chol0)

    # Compute the evolution errors:
    omega = diff(mu, differences = D)

    # And the initial states:
    mu0 = as.matrix(mu[1:D,])

    # Sample the initial variance parameters:
    evolParams0 = sampleEvol0(mu0, evolParams0, A = 1)

    # Sample the evolution variance and associated parameters:
    evolParams = sampleEvolParams(omega = omega, evolParams = evolParams, sigma_e = 1/sqrt(T), evol_error = evol_error)

    # Sample the shrink-to-zero variance and associated parameters:
    zeroParams = sampleEvolParams(omega = mu, evolParams = zeroParams, sigma_e = 1/sqrt(T), evol_error = zero_error)

    # Store the MCMC output:
    if(nsi > nburn){
      # Increment the skip counter:
      skipcount = skipcount + 1

      # Save the iteration:
      if(skipcount > nskip){
        # Increment the save index
        isave = isave + 1

        # Save the MCMC samples:
        if(!is.na(match('mu', mcmc_params)) || computeDIC) post_mu[isave,] = mu
        if(!is.na(match('yhat', mcmc_params))) post_yhat[isave,] = mu + sigma_et*rnorm(T)
        if(!is.na(match('obs_sigma_t2', mcmc_params)) || computeDIC) post_obs_sigma_t2[isave,] = sigma_et^2
        if(!is.na(match('evol_sigma_t2', mcmc_params))) post_evol_sigma_t2[isave,] = c(evolParams0$sigma_w0^2, evolParams$sigma_wt^2)
        if(!is.na(match('zero_sigma_t2', mcmc_params))) post_zero_sigma_t2[isave,] = zeroParams$sigma_wt^2
        if(!is.na(match('dhs_phi', mcmc_params)) && evol_error == "DHS") post_dhs_phi[isave] = evolParams$dhs_phi
        if(!is.na(match('dhs_mean', mcmc_params)) && evol_error == "DHS") post_dhs_mean[isave] = evolParams$dhs_mean
        if(!is.na(match('dhs_phi_zero', mcmc_params)) && zero_error == "DHS") post_dhs_phi_zero[isave] = zeroParams$dhs_phi
        if(!is.na(match('dhs_mean_zero', mcmc_params)) && zero_error == "DHS") post_dhs_mean_zero[isave] = zeroParams$dhs_mean
        post_loglike[isave] = sum(dnorm(y, mean = mu, sd = sigma_et, log = TRUE))

        # And reset the skip counter:
        skipcount = 0
      }
    }
    if(verbose) computeTimeRemaining(nsi, timer0, nstot, nrep = 1000)
  }

  if(!is.na(match('mu', mcmc_params))) mcmc_output$mu = post_mu
  if(!is.na(match('yhat', mcmc_params))) mcmc_output$yhat = post_yhat
  if(!is.na(match('obs_sigma_t2', mcmc_params))) mcmc_output$obs_sigma_t2 = post_obs_sigma_t2
  if(!is.na(match('evol_sigma_t2', mcmc_params))) mcmc_output$evol_sigma_t2 = post_evol_sigma_t2
  if(!is.na(match('zero_sigma_t2', mcmc_params))) mcmc_output$zero_sigma_t2 = post_zero_sigma_t2
  if(!is.na(match('dhs_phi', mcmc_params)) && evol_error == "DHS") mcmc_output$dhs_phi = post_dhs_phi
  if(!is.na(match('dhs_mean', mcmc_params)) && evol_error == "DHS") mcmc_output$dhs_mean = post_dhs_mean
  if(!is.na(match('dhs_phi_zero', mcmc_params)) && zero_error == "DHS") mcmc_output$dhs_phi_zero = post_dhs_phi_zero
  if(!is.na(match('dhs_mean_zero', mcmc_params)) && zero_error == "DHS") mcmc_output$dhs_mean_zero = post_dhs_mean_zero

  # Also include the log-likelihood:
  mcmc_output$loglike = post_loglike

  if(computeDIC){
    # Log-likelihood evaluated at posterior means:
    loglike_hat = sum(dnorm(y,
                            mean = colMeans(post_mu),
                            sd = colMeans(sqrt(post_obs_sigma_t2)),
                            log = TRUE))

    # Effective number of parameters (Note: two options)
    p_d = c(2*(loglike_hat - mean(post_loglike)),
            2*var(post_loglike))
    # DIC:
    DIC = -2*loglike_hat + 2*p_d

    # Store the DIC and the effective number of parameters (p_d)
    mcmc_output$DIC = DIC; mcmc_output$p_d = p_d
  }

  if(verbose) print(paste('Total time: ', round((proc.time()[3] - timer0)), 'seconds'))

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
#' \item the Bayesian lasso ('BL');
#' \item the normal stochastic volatility model ('SV');
#' \item the normal-inverse-gamma prior ('NIG').
#' }
#' In each case, the evolution error is a scale mixture of Gaussians.
#' Sampling is accomplished with a (parameter-expanded) Gibbs sampler,
#' mostly relying on a dynamic linear model representation.
#'
#' @param y the \code{T x 1} vector of time series observations
#' @param X the \code{T x p} matrix of time series predictors
#' @param evol_error the evolution error distribution; must be one of
#' 'DHS' (dynamic horseshoe prior), 'HS' (horseshoe prior), 'BL' (Bayesian lasso), or 'NIG' (normal-inverse-gamma prior)
#' @param D degree of differencing (D = 1 or D = 2)
#' @param useObsSV logical; if TRUE, include a (normal) stochastic volatility model
#' for the observation error variance
#' @param nsave number of MCMC iterations to record
#' @param nburn number of MCMC iterations to discard (burin-in)
#' @param nskip number of MCMC iterations to skip between saving iterations,
#' i.e., save every (nskip + 1)th draw
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
#' @param use_backfitting logical; if TRUE, use backfitting to sample the predictors j=1,...,p
#' (faster, but usually less MCMC efficient)
#' @param computeDIC logical; if TRUE, compute the deviance information criterion \code{DIC}
#' and the effective number of parameters \code{p_d}
#' @param verbose logical; should R report extra information on progress?
#'
#' @return A named list of the \code{nsave} MCMC samples for the parameters named in \code{mcmc_params}
#'
#' @note The data \code{y} may contain NAs, which will be treated with a simple imputation scheme
#' via an additional Gibbs sampling step. In general, rescaling \code{y} to have unit standard
#' deviation is recommended to avoid numerical issues.
#'
#' @examples
#' # Example 1: all signals
#' simdata = simRegression(T = 200, p = 5, p_0 = 0)
#' y = simdata$y; X = simdata$X
#' out = btf_reg(y, X)
#' for(j in 1:ncol(X))
#'  plot_fitted(rep(0, length(y)),
#'              mu = colMeans(out$beta[,,j]),
#'              postY = out$beta[,,j],
#'              y_true = simdata$beta_true[,j])
#'
#' \dontrun{
#' # Example 2: some noise, longer series
#' simdata = simRegression(T = 500, p = 10, p_0 = 5)
#' y = simdata$y; X = simdata$X
#' out = btf_reg(y, X, nsave = 1000, nskip = 0) # Short MCMC run for a quick example
#' for(j in 1:ncol(X))
#'   plot_fitted(rep(0, length(y)),
#'               mu = colMeans(out$beta[,,j]),
#'               postY = out$beta[,,j],
#'               y_true = simdata$beta_true[,j])
#' }
#'
#' @import spam
#' @export
btf_reg = function(y, X = NULL, evol_error = 'DHS', D = 1, useObsSV = FALSE,
                   nsave = 1000, nburn = 1000, nskip = 4,
                   mcmc_params = list("mu", "yhat","beta","evol_sigma_t2", "obs_sigma_t2", "dhs_phi", "dhs_mean"),
                   use_backfitting = FALSE,
                   computeDIC = TRUE,
                   verbose = TRUE){

  # Convert to upper case:
  evol_error = toupper(evol_error)

  # If no predictors are specified, just call the univariate trend filtering model: btf()
  if(is.null(X)) {
    mcmc_params[match("beta", mcmc_params)] = NULL # Remove "beta", since it does not apply for btf()
    return(btf(y = y, evol_error = evol_error, D = D, nsave = nsave, nburn = nburn,
               nskip = nskip, mcmc_params = mcmc_params, computeDIC = computeDIC, verbose = verbose))
  }

  # Time points (in [0,1])
  T = length(y); t01 = seq(0, 1, length.out=T);

  # Begin by checking for missing values, then imputing (for initialization)
  is.missing = which(is.na(y)); any.missing = (length(is.missing) > 0)

  # Impute the active "data"
  if(any.missing) y[is.missing] = mean(y, na.rm=TRUE)

  # Obtain XtX
  XtX = build_XtX(X);

  # Number of predictors:
  p = ncol(X)

  # Initial SD (implicitly assumes a constant mean)
  sigma_e = sd(y, na.rm=TRUE); sigma_et = rep(sigma_e, T)

  # Compute the Cholesky term: use random variances for a more conservative sparsity pattern
  chol0 = initCholReg.spam(obs_sigma_t2 = abs(rnorm(T)),
                           evol_sigma_t2 = matrix(abs(rnorm(T*p)), nrow = T),
                           XtX = XtX, D = D)

  # Initialize the dynamic regression coefficients, beta, via sampling:
  beta = sampleBTF_reg(y, X, obs_sigma_t2 = sigma_et^2, evol_sigma_t2 = matrix(0.01*sigma_et^2, nr = T, nc = p), XtX = XtX, D = D, chol0 = chol0)

  # Conditional mean:
  mu = rowSums(X*beta)

  # Compute the evolution errors:
  omega = diff(beta, differences = D)

  # And the initial states:
  beta0 = matrix(beta[1:D,], nr = D)

  # Initialize the evolution error variance paramters:
  evolParams = initEvolParams(omega, evol_error = evol_error)

  # Initial variance parameters:
  evolParams0 = initEvol0(beta0, commonSD = FALSE)

  # SV parameters, if necessary:
  if(useObsSV) {svParams = initSV(y - mu); sigma_et = svParams$sigma_wt}

  # Store the MCMC output in separate arrays (better computation times)
  mcmc_output = vector('list', length(mcmc_params)); names(mcmc_output) = mcmc_params
  if(!is.na(match('mu', mcmc_params)) || computeDIC) post_mu = array(NA, c(nsave, T))
  if(!is.na(match('yhat', mcmc_params))) post_yhat = array(NA, c(nsave, T))
  if(!is.na(match('beta', mcmc_params))) post_beta = array(NA, c(nsave, T, p))
  if(!is.na(match('obs_sigma_t2', mcmc_params)) || computeDIC) post_obs_sigma_t2 = array(NA, c(nsave, T))
  if(!is.na(match('evol_sigma_t2', mcmc_params))) post_evol_sigma_t2 = array(NA, c(nsave, T, p))
  if(!is.na(match('dhs_phi', mcmc_params)) && evol_error == "DHS") post_dhs_phi = array(NA, c(nsave, p))
  if(!is.na(match('dhs_mean', mcmc_params)) && evol_error == "DHS") post_dhs_mean = array(NA, c(nsave, p))
  post_loglike = numeric(nsave)

  # Total number of MCMC simulations:
  nstot = nburn+(nskip+1)*(nsave)
  skipcount = 0; isave = 0 # For counting

  # Run the MCMC:
  if(verbose) timer0 = proc.time()[3] # For timing the sampler
  for(nsi in 1:nstot){

    # Impute missing values, if any:
    if(any.missing) y[is.missing] = mu[is.missing] + sigma_et[is.missing]*rnorm(length(is.missing))

    # Sample the dynamic regression coefficients, beta:
    # Backfitting is faster, but likely less efficient
    if(use_backfitting){
      beta = sampleBTF_reg_backfit(y, X, beta, obs_sigma_t2 = sigma_et^2, evol_sigma_t2 = rbind(matrix(evolParams0$sigma_w0^2, nr = D), evolParams$sigma_wt^2), D = D)
    } else beta = sampleBTF_reg(y, X, obs_sigma_t2 = sigma_et^2, evol_sigma_t2 = rbind(matrix(evolParams0$sigma_w0^2, nr = D), evolParams$sigma_wt^2), XtX = XtX, D = D, chol0 = chol0)

    # Conditional mean:
    mu = rowSums(X*beta)

    # Compute the evolution errors:
    omega = diff(beta, differences = D)

    # And the initial states:
    beta0 = matrix(beta[1:D,], nr = D)

    # Sample the initial variance parameters:
    evolParams0 = sampleEvol0(beta0, evolParams0, A = 1, commonSD = FALSE)

    # Sample the (observation and evolution) variances and associated parameters:
    if(useObsSV){
      # Evolution error variance + params:
      evolParams = sampleEvolParams(omega, evolParams, 1/sqrt(T*p), evol_error)

      # Observation error variance + params:
      svParams = sampleSVparams(omega = y - mu, svParams = svParams)
      sigma_et = svParams$sigma_wt

    } else {
      # Evolution error variance + params:
      evolParams = sampleEvolParams(omega, evolParams, sigma_e/sqrt(T*p), evol_error)

      # Sample the observation error SD:
      if(evol_error == 'DHS') {
        sigma_e = uni.slice(sigma_e, g = function(x){
          -(T+2)*log(x) - 0.5*sum((y - mu)^2, na.rm=TRUE)/x^2 - log(1 + (sqrt(T*p)*exp(evolParams$dhs_mean0/2)/x)^2)
        }, lower = 0, upper = Inf)[1]
      }
      #if(evol_error == 'HS') sigma_e = 1/sqrt(rgamma(n = 1, shape = T/2 + length(evolParams$xiLambda), rate = sum((y - mu)^2, na.rm=TRUE)/2 + T*p*sum(evolParams$xiLambda)))
      if(evol_error == 'HS') sigma_e = 1/sqrt(rgamma(n = 1, shape = T/2, rate = sum((y - mu)^2, na.rm=TRUE)/2))
      if(evol_error == 'BL') sigma_e = 1/sqrt(rgamma(n = 1, shape = T/2 + length(evolParams$tau_j)/2, rate = sum((y - mu)^2, na.rm=TRUE)/2 + T*p*sum((omega/evolParams$tau_j)^2)/2))
      if((evol_error == 'NIG') || (evol_error == 'SV')) sigma_e = 1/sqrt(rgamma(n = 1, shape = T/2, rate = sum((y - mu)^2, na.rm=TRUE)/2))

      # Replicate for coding convenience:
      sigma_et = rep(sigma_e, T)
    }


    # Store the MCMC output:
    if(nsi > nburn){
      # Increment the skip counter:
      skipcount = skipcount + 1

      # Save the iteration:
      if(skipcount > nskip){
        # Increment the save index
        isave = isave + 1

        # Save the MCMC samples:
        if(!is.na(match('mu', mcmc_params)) || computeDIC) post_mu[isave,] = mu
        if(!is.na(match('yhat', mcmc_params))) post_yhat[isave,] = mu + sigma_et*rnorm(T)
        if(!is.na(match('beta', mcmc_params))) post_beta[isave,,] = beta
        if(!is.na(match('obs_sigma_t2', mcmc_params)) || computeDIC) post_obs_sigma_t2[isave,] = sigma_et^2
        if(!is.na(match('evol_sigma_t2', mcmc_params))) {
          if(D == 0){
            post_evol_sigma_t2[isave,,] = evolParams$sigma_wt^2
          } else post_evol_sigma_t2[isave,,] = rbind(matrix(evolParams0$sigma_w0^2, nr = D), evolParams$sigma_wt^2)
        }
        if(!is.na(match('dhs_phi', mcmc_params)) && evol_error == "DHS") post_dhs_phi[isave,] = evolParams$dhs_phi
        if(!is.na(match('dhs_mean', mcmc_params)) && evol_error == "DHS") post_dhs_mean[isave,] = evolParams$dhs_mean
        post_loglike[isave] = sum(dnorm(y, mean = mu, sd = sigma_et, log = TRUE))

        # And reset the skip counter:
        skipcount = 0
      }
    }
    if(verbose) computeTimeRemaining(nsi, timer0, nstot, nrep = 1000)
  }

  if(!is.na(match('mu', mcmc_params))) mcmc_output$mu = post_mu
  if(!is.na(match('yhat', mcmc_params))) mcmc_output$yhat = post_yhat
  if(!is.na(match('beta', mcmc_params))) mcmc_output$beta = post_beta
  if(!is.na(match('obs_sigma_t2', mcmc_params))) mcmc_output$obs_sigma_t2 = post_obs_sigma_t2
  if(!is.na(match('evol_sigma_t2', mcmc_params))) mcmc_output$evol_sigma_t2 = post_evol_sigma_t2
  if(!is.na(match('dhs_phi', mcmc_params)) && evol_error == "DHS") mcmc_output$dhs_phi = post_dhs_phi
  if(!is.na(match('dhs_mean', mcmc_params)) && evol_error == "DHS") mcmc_output$dhs_mean = post_dhs_mean

  # Also include the log-likelihood:
  mcmc_output$loglike = post_loglike

  if(computeDIC){
    # Log-likelihood evaluated at posterior means:
    loglike_hat = sum(dnorm(y,
                            mean = colMeans(post_mu),
                            sd = colMeans(sqrt(post_obs_sigma_t2)),
                            log = TRUE))
    # Effective number of parameters (Note: two options)
    p_d = c(2*(loglike_hat - mean(post_loglike)),
            2*var(post_loglike))
    # DIC:
    DIC = -2*loglike_hat + 2*p_d

    # Store the DIC and the effective number of parameters (p_d)
    mcmc_output$DIC = DIC; mcmc_output$p_d = p_d
  }
  if(verbose) print(paste('Total time: ', round((proc.time()[3] - timer0)), 'seconds'))

  return (mcmc_output);
}
#----------------------------------------------------------------------------
#' MCMC Sampler for B-spline Bayesian Trend Filtering
#'
#' Run the MCMC for B-spline fitting with a Bayesian trend filtering model on the
#' coefficients, i.e., a penalty on zeroth (D=0), first (D=1), or second (D=2)
#' differences of the B-spline basis coefficients.
#' The penalty is determined by the prior on the evolution errors, which include:
#' \itemize{
#' \item the dynamic horseshoe prior ('DHS');
#' \item the static horseshoe prior ('HS');
#' \item the Bayesian lasso ('BL');
#' \item the normal stochastic volatility model ('SV');
#' \item the normal-inverse-gamma prior ('NIG').
#' }
#' In each case, the evolution error is a scale mixture of Gaussians.
#' Sampling is accomplished with a (parameter-expanded) Gibbs sampler,
#' mostly relying on a dynamic linear model representation.

#' @param y the \code{T x 1} vector of time series observations
#' @param x the \code{T x 1} vector of observation points; if NULL, assume equally spaced
#' @param num_knots the number of knots; if NULL, use the default of \code{max(20, min(ceiling(T/4), 150))}
#' @param evol_error the evolution error distribution; must be one of
#' 'DHS' (dynamic horseshoe prior), 'HS' (horseshoe prior), 'BL' (Bayesian lasso), or 'NIG' (normal-inverse-gamma prior)
#' @param D degree of differencing (D = 0, D = 1, or D = 2)
#' @param nsave number of MCMC iterations to record
#' @param nburn number of MCMC iterations to discard (burin-in)
#' @param nskip number of MCMC iterations to skip between saving iterations,
#' i.e., save every (nskip + 1)th draw
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
#' @param computeDIC logical; if TRUE, compute the deviance information criterion \code{DIC}
#' and the effective number of parameters \code{p_d}
#' @param verbose logical; should R report extra information on progress?
#'
#' @return A named list of the \code{nsave} MCMC samples for the parameters named in \code{mcmc_params}
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
#' @examples
#' # Example 1: Blocks data
#' simdata = simUnivariate(signalName = "blocks", T = 1000, RSNR = 3, include_plot = TRUE)
#' y = simdata$y
#' out = btf_bspline(y, D = 1)
#' plot_fitted(y, mu = colMeans(out$mu), postY = out$yhat, y_true = simdata$y_true)
#'
#' \dontrun{
#' # Example 2: motorcycle data (unequally-spaced points)
#' library(MASS)
#' y = scale(mcycle$accel) # Center and Scale for numerical stability
#' x = mcycle$times
#' plot(x, y, xlab = 'Time (ms)', ylab='Acceleration (g)', main = 'Motorcycle Crash Data')
#' out = btf_bspline(y = y, x = x)
#' plot_fitted(y, mu = colMeans(out$mu), postY = out$yhat, t01 = x)
#' }
#'
#' @import fda
#' @export
btf_bspline = function(y, x = NULL, num_knots = NULL, evol_error = 'DHS', D = 2,
                       nsave = 1000, nburn = 1000, nskip = 4,
                       mcmc_params = list("mu", "yhat", "beta", "evol_sigma_t2", "obs_sigma_t2", "dhs_phi", "dhs_mean"),
                       computeDIC = TRUE,
                       verbose = TRUE){

  # Convert to upper case:
  evol_error = toupper(evol_error)

  # For D = 0, return special case:
  if(D == 0){
    return(btf_bspline0(y = y, x = x, num_knots = num_knots, evol_error = evol_error,
                        nsave = nsave, nburn = nburn, nskip = nskip,
                        mcmc_params = mcmc_params,
                        computeDIC = computeDIC, verbose = verbose))
  }

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

  # Impute the active "data"
  y = approxfun(t01, y, rule = 2)(t01)

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
  if(!is.na(match('mu', mcmc_params)) || computeDIC) post_mu = array(NA, c(nsave, T))
  if(!is.na(match('yhat', mcmc_params))) post_yhat = array(NA, c(nsave, T))
  if(!is.na(match('beta', mcmc_params))) post_beta = array(NA, c(nsave, p))
  if(!is.na(match('obs_sigma_t2', mcmc_params)) || computeDIC) post_obs_sigma_t2 = array(NA, c(nsave, T))
  if(!is.na(match('evol_sigma_t2', mcmc_params))) post_evol_sigma_t2 = array(NA, c(nsave, p))
  if(!is.na(match('dhs_phi', mcmc_params)) && evol_error == "DHS") post_dhs_phi = numeric(nsave)
  if(!is.na(match('dhs_mean', mcmc_params)) && evol_error == "DHS") post_dhs_mean = numeric(nsave)
  post_loglike = numeric(nsave)

  # Total number of MCMC simulations:
  nstot = nburn+(nskip+1)*(nsave)
  skipcount = 0; isave = 0 # For counting

  # Run the MCMC:
  if(verbose) timer0 = proc.time()[3] # For timing the sampler
  for(nsi in 1:nstot){

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
    if(evol_error == 'BL') sigma_e = 1/sqrt(rgamma(n = 1, shape = T/2 + length(evolParams$tau_j)/2, rate = sum((y - mu)^2, na.rm=TRUE)/2 + p*sum((omega/evolParams$tau_j)^2)/2))
    if((evol_error == 'NIG') || (evol_error == 'SV')) sigma_e = 1/sqrt(rgamma(n = 1, shape = T/2, rate = sum((y - mu)^2, na.rm=TRUE)/2))

    # Sample the initial variance parameters:
    evolParams0 = sampleEvol0(beta0, evolParams0, A = 1)

    # Store the MCMC output:
    if(nsi > nburn){
      # Increment the skip counter:
      skipcount = skipcount + 1

      # Save the iteration:
      if(skipcount > nskip){
        # Increment the save index
        isave = isave + 1

        # Save the MCMC samples:
        if(!is.na(match('mu', mcmc_params)) || computeDIC) post_mu[isave,] = mu
        if(!is.na(match('yhat', mcmc_params))) post_yhat[isave,] = mu + sigma_e*rnorm(T)
        if(!is.na(match('beta', mcmc_params))) post_beta[isave,] = matrix(beta)
        if(!is.na(match('obs_sigma_t2', mcmc_params)) || computeDIC) post_obs_sigma_t2[isave,] = sigma_e^2
        if(!is.na(match('evol_sigma_t2', mcmc_params))) post_evol_sigma_t2[isave,] =  rbind(matrix(evolParams0$sigma_w0^2, nr = D), evolParams$sigma_wt^2)
        if(!is.na(match('dhs_phi', mcmc_params)) && evol_error == "DHS") post_dhs_phi[isave] = evolParams$dhs_phi
        if(!is.na(match('dhs_mean', mcmc_params)) && evol_error == "DHS") post_dhs_mean[isave] = evolParams$dhs_mean
        post_loglike[isave] = sum(dnorm(y, mean = mu, sd = sigma_e, log = TRUE))

        # And reset the skip counter:
        skipcount = 0
      }
    }
    if(verbose) computeTimeRemaining(nsi, timer0, nstot, nrep = 1000)
  }

  if(!is.na(match('mu', mcmc_params))) mcmc_output$mu = post_mu
  if(!is.na(match('yhat', mcmc_params))) mcmc_output$yhat = post_yhat
  if(!is.na(match('beta', mcmc_params))) mcmc_output$beta = post_beta
  if(!is.na(match('obs_sigma_t2', mcmc_params))) mcmc_output$obs_sigma_t2 = post_obs_sigma_t2
  if(!is.na(match('evol_sigma_t2', mcmc_params))) mcmc_output$evol_sigma_t2 = post_evol_sigma_t2
  if(!is.na(match('dhs_phi', mcmc_params)) && evol_error == "DHS") mcmc_output$dhs_phi = post_dhs_phi
  if(!is.na(match('dhs_mean', mcmc_params)) && evol_error == "DHS") mcmc_output$dhs_mean = post_dhs_mean

  # Also include the log-likelihood:
  mcmc_output$loglike = post_loglike

  if(computeDIC){
    # Log-likelihood evaluated at posterior means:
    loglike_hat = sum(dnorm(y,
                            mean = colMeans(post_mu),
                            sd = colMeans(sqrt(post_obs_sigma_t2)),
                            log = TRUE))
    # Effective number of parameters (Note: two options)
    p_d = c(2*(loglike_hat - mean(post_loglike)),
            2*var(post_loglike))
    # DIC:
    DIC = -2*loglike_hat + 2*p_d

    # Store the DIC and the effective number of parameters (p_d)
    mcmc_output$DIC = DIC; mcmc_output$p_d = p_d
  }
  if(verbose) print(paste('Total time: ', round((proc.time()[3] - timer0)), 'seconds'))

  return (mcmc_output);
}
#----------------------------------------------------------------------------
#' MCMC Sampler for Time-Varying Autoregression
#'
#' Run the MCMC for a time-varying autoregression with a penalty on
#' first (D=1) or second (D=2) differences of each dynamic autoregressive coefficient.
#' The penalty is determined by the prior on the evolution errors, which include:
#' \itemize{
#' \item the dynamic horseshoe prior ('DHS');
#' \item the static horseshoe prior ('HS');
#' \item the Bayesian lasso ('BL');
#' \item the normal stochastic volatility model ('SV');
#' \item the normal-inverse-gamma prior ('NIG').
#' }
#' In each case, the evolution error is a scale mixture of Gaussians.
#' Sampling is accomplished with a (parameter-expanded) Gibbs sampler,
#' mostly relying on a dynamic linear model representation.
#'
#' @param y the \code{T x 1} vector of time series observations
#' @param p_max the maximum order of lag to include
#' @param include_intercept logical; if TRUE, include a time-varying intercept
#' @param evol_error the evolution error distribution; must be one of
#' 'DHS' (dynamic horseshoe prior), 'HS' (horseshoe prior), 'BL' (Bayesian lasso), or 'NIG' (normal-inverse-gamma prior)
#' @param D degree of differencing (D = 1 or D = 2)
#' @param useObsSV logical; if TRUE, include a (normal) stochastic volatility model
#' for the observation error variance
#' @param nsave number of MCMC iterations to record
#' @param nburn number of MCMC iterations to discard (burin-in)
#' @param nskip number of MCMC iterations to skip between saving iterations,
#' i.e., save every (nskip + 1)th draw
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
#' @param computeDIC logical; if TRUE, compute the deviance information criterion \code{DIC}
#' and the effective number of parameters \code{p_d}
#' @param verbose logical; should R report extra information on progress?
#'
#' @return A named list of the \code{nsave} MCMC samples for the parameters named in \code{mcmc_params}
#'
#' @note The data \code{y} may NOT contain NAs. The data \code{y} will be used to construct
#' the predictor matrix (of lagged values), which is not permitted to include NAs.
#'
#' @examples
#' \dontrun{
#' # Example 1:
#' simdata = simUnivariate(signalName = "doppler", T = 128, RSNR = 7, include_plot = TRUE)
#' y = simdata$y
#'
#' # Note: in general should subtract off the sample mean
#' p = 6 # Lag
#' out = tvar(y, p_max = p, include_intercept = FALSE,
#'                    evol_error = 'DHS', D = 1,
#'                    mcmc_params = list('mu', 'yhat', 'beta', 'obs_sigma_t2'))
#'
#' for(j in 1:p)
#'   plot_fitted(rep(0, length(y) - p),
#'            mu = colMeans(out$beta[,,j]),
#'            postY = out$beta[,,j])
#'
#' plot_fitted(y[-(1:p)],
#'            mu = colMeans(out$mu),
#'            postY = out$yhat,
#'            y_true = simdata$y_true[-(1:p)])
#'
#' spec_TF = post_spec_dsp(post_ar_coefs = out$beta,
#'                    post_sigma_e = sqrt(out$obs_sigma_t2[,1]),
#'                    n.freq = 100)
#'
#' image(x = 1:(length(y)-p), y = spec_TF$freq, colMeans(log(spec_TF$post_spec)),
#'      xlab = 'Time', ylab = 'Freqency', main = 'Posterior Mean of Log-Spectrum')
#' }
#'
#' @export
tvar = function(y, p_max = 1, include_intercept = FALSE,
                evol_error = 'DHS', D = 2, useObsSV = FALSE,
                nsave = 1000, nburn = 1000, nskip = 4,
                mcmc_params = list("mu", "yhat", "beta", "evol_sigma_t2", "obs_sigma_t2", "dhs_phi", "dhs_mean"),
                computeDIC = TRUE,
                verbose = TRUE){

  # Convert to upper case:
  evol_error = toupper(evol_error)

  # Some quick checks:
  if((p_max < 0) || (p_max != round(p_max)))  stop('p_max must be a positive integer')
  if(any(is.na(y))) stop('NAs not permitted in y for tvar()')

  # Compute design matrix:
  X = getARpXmat(y, p_max, include_intercept)

  # Update the vector y to match the length of X:
  y = y[-(1:p_max)]

  # Once you have X, simply run the usual regression code:
  return(btf_reg(y = y, X = X, evol_error = evol_error, D = D, useObsSV = useObsSV,
                 nsave = nsave, nburn = nburn, nskip = nskip,
                 mcmc_params = mcmc_params, computeDIC = computeDIC, verbose = verbose))
}
#----------------------------------------------------------------------------
#' MCMC Sampler for B-spline Bayesian Trend Filtering: D = 0
#'
#' Run the MCMC for B-spline fitting with a penalty the B-spline basis coefficients.
#' The penalty is determined by the prior on the evolution errors, which include:
#' \itemize{
#' \item the dynamic horseshoe prior ('DHS');
#' \item the static horseshoe prior ('HS');
#' \item the Bayesian lasso ('BL');
#' \item the normal stochastic volatility model ('SV');
#' \item the normal-inverse-gamma prior ('NIG').
#' }
#' In each case, the evolution error is a scale mixture of Gaussians.
#' Sampling is accomplished with a (parameter-expanded) Gibbs sampler,
#' mostly relying on a dynamic linear model representation.

#' @param y the \code{T x 1} vector of time series observations
#' @param x the \code{T x 1} vector of observation points; if NULL, assume equally spaced
#' @param num_knots the number of knots; if NULL, use the default of \code{max(20, min(ceiling(T/4), 150))}
#' @param evol_error the evolution error distribution; must be one of
#' 'DHS' (dynamic horseshoe prior), 'HS' (horseshoe prior), 'BL' (Bayesian lasso), or 'NIG' (normal-inverse-gamma prior)
#' @param nsave number of MCMC iterations to record
#' @param nburn number of MCMC iterations to discard (burin-in)
#' @param nskip number of MCMC iterations to skip between saving iterations,
#' i.e., save every (nskip + 1)th draw
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
#' @param computeDIC logical; if TRUE, compute the deviance information criterion \code{DIC}
#' and the effective number of parameters \code{p_d}
#' @param verbose logical; should R report extra information on progress?
#'
#' @return A named list of the \code{nsave} MCMC samples for the parameters named in \code{mcmc_params}
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
#' @import fda
btf_bspline0 = function(y, x = NULL, num_knots = NULL, evol_error = 'DHS',
                        nsave = 1000, nburn = 1000, nskip = 4,
                        mcmc_params = list("mu", "yhat", "beta", "evol_sigma_t2", "obs_sigma_t2", "dhs_phi", "dhs_mean"),
                        computeDIC = TRUE,
                        verbose = TRUE){

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

  # Impute the active "data"
  y = approxfun(t01, y, rule = 2)(t01)

  # Initial SD (implicitly assumes a constant mean)
  sigma_e = sd(y, na.rm=TRUE)

  # Initialize the B-spline coefficients via sampling
  Xty = crossprod(X,y)
  beta = sampleBTF_bspline(y, X, obs_sigma2 = sigma_e^2, evol_sigma_t2 = rep(0.01*sigma_e^2, p), XtX_bands = XtX_bands, Xty = Xty, D = 0)

  # And the conditional expectation:
  mu = as.numeric(X%*%beta)

  # Initialize the evolution error variance paramters:
  evolParams = initEvolParams(omega = beta, evol_error = evol_error)

  # Store the MCMC output in separate arrays (better computation times)
  mcmc_output = vector('list', length(mcmc_params)); names(mcmc_output) = mcmc_params
  if(!is.na(match('mu', mcmc_params)) || computeDIC) post_mu = array(NA, c(nsave, T))
  if(!is.na(match('yhat', mcmc_params))) post_yhat = array(NA, c(nsave, T))
  if(!is.na(match('beta', mcmc_params))) post_beta = array(NA, c(nsave, p))
  if(!is.na(match('obs_sigma_t2', mcmc_params)) || computeDIC) post_obs_sigma_t2 = array(NA, c(nsave, T))
  if(!is.na(match('evol_sigma_t2', mcmc_params))) post_evol_sigma_t2 = array(NA, c(nsave, p))
  if(!is.na(match('dhs_phi', mcmc_params)) && evol_error == "DHS") post_dhs_phi = numeric(nsave)
  if(!is.na(match('dhs_mean', mcmc_params)) && evol_error == "DHS") post_dhs_mean = numeric(nsave)
  post_loglike = numeric(nsave)

  # Total number of MCMC simulations:
  nstot = nburn+(nskip+1)*(nsave)
  skipcount = 0; isave = 0 # For counting

  # Run the MCMC:
  if(verbose) timer0 = proc.time()[3] # For timing the sampler
  for(nsi in 1:nstot){

    # Impute missing values, if any:
    if(any.missing) {
      y[is.missing] = mu[is.missing] + sigma_e*rnorm(length(is.missing))
      Xty = crossprod(X,y) # Recompute when missing values present
    }

    # Sample the dynamic regression coefficients, beta:
    beta = sampleBTF_bspline(y, X, obs_sigma2 = sigma_e^2, evol_sigma_t2 = evolParams$sigma_wt^2, XtX_bands = XtX_bands, Xty = Xty, D = 0)

    # And the conditional expectation:
    mu = as.numeric(X%*%beta)

    # Sample the evolution error variance (and associated parameters):
    evolParams = sampleEvolParams(omega = beta, evolParams, sigma_e/sqrt(p), evol_error)

    # Sample the observation error SD:
    if(evol_error == 'DHS') {
      sigma_e = uni.slice(sigma_e, g = function(x){
        -(T+2)*log(x) - 0.5*sum((y - mu)^2, na.rm=TRUE)/x^2 - log(1 + (sqrt(p)*exp(evolParams$dhs_mean0/2)/x)^2)
      }, lower = 0, upper = Inf)[1]
    }
    if(evol_error == 'HS') sigma_e = 1/sqrt(rgamma(n = 1, shape = T/2 + length(evolParams$xiLambda), rate = sum((y - mu)^2, na.rm=TRUE)/2 + p*sum(evolParams$xiLambda)))
    if(evol_error == 'BL') sigma_e = 1/sqrt(rgamma(n = 1, shape = T/2 + length(evolParams$tau_j)/2, rate = sum((y - mu)^2, na.rm=TRUE)/2 + p*sum((beta/evolParams$tau_j)^2)/2))
    if((evol_error == 'NIG') || (evol_error == 'SV'))  sigma_e = 1/sqrt(rgamma(n = 1, shape = T/2, rate = sum((y - mu)^2, na.rm=TRUE)/2))

    # Store the MCMC output:
    if(nsi > nburn){
      # Increment the skip counter:
      skipcount = skipcount + 1

      # Save the iteration:
      if(skipcount > nskip){
        # Increment the save index
        isave = isave + 1

        # Save the MCMC samples:
        if(!is.na(match('mu', mcmc_params)) || computeDIC) post_mu[isave,] = mu
        if(!is.na(match('yhat', mcmc_params))) post_yhat[isave,] = mu + sigma_e*rnorm(T)
        if(!is.na(match('beta', mcmc_params))) post_beta[isave,] = matrix(beta)
        if(!is.na(match('obs_sigma_t2', mcmc_params)) || computeDIC) post_obs_sigma_t2[isave,] = sigma_e^2
        if(!is.na(match('evol_sigma_t2', mcmc_params))) post_evol_sigma_t2[isave,] = evolParams$sigma_wt^2
        if(!is.na(match('dhs_phi', mcmc_params)) && evol_error == "DHS") post_dhs_phi[isave] = evolParams$dhs_phi
        if(!is.na(match('dhs_mean', mcmc_params)) && evol_error == "DHS") post_dhs_mean[isave] = evolParams$dhs_mean
        post_loglike[isave] = sum(dnorm(y, mean = mu, sd = sigma_e, log = TRUE))

        # And reset the skip counter:
        skipcount = 0
      }
    }
    if(verbose) computeTimeRemaining(nsi, timer0, nstot, nrep = 1000)
  }

  if(!is.na(match('mu', mcmc_params))) mcmc_output$mu = post_mu
  if(!is.na(match('yhat', mcmc_params))) mcmc_output$yhat = post_yhat
  if(!is.na(match('beta', mcmc_params))) mcmc_output$beta = post_beta
  if(!is.na(match('obs_sigma_t2', mcmc_params))) mcmc_output$obs_sigma_t2 = post_obs_sigma_t2
  if(!is.na(match('evol_sigma_t2', mcmc_params))) mcmc_output$evol_sigma_t2 = post_evol_sigma_t2
  if(!is.na(match('dhs_phi', mcmc_params)) && evol_error == "DHS") mcmc_output$dhs_phi = post_dhs_phi
  if(!is.na(match('dhs_mean', mcmc_params)) && evol_error == "DHS") mcmc_output$dhs_mean = post_dhs_mean

  # Also include the log-likelihood:
  mcmc_output$loglike = post_loglike

  if(computeDIC){
    # Log-likelihood evaluated at posterior means:
    loglike_hat = sum(dnorm(y,
                            mean = colMeans(post_mu),
                            sd = colMeans(sqrt(post_obs_sigma_t2)),
                            log = TRUE))
    # Effective number of parameters (Note: two options)
    p_d = c(2*(loglike_hat - mean(post_loglike)),
            2*var(post_loglike))
    # DIC:
    DIC = -2*loglike_hat + 2*p_d

    # Store the DIC and the effective number of parameters (p_d)
    mcmc_output$DIC = DIC; mcmc_output$p_d = p_d
  }
  if(verbose) print(paste('Total time: ', round((proc.time()[3] - timer0)), 'seconds'))

  return (mcmc_output);
}
#----------------------------------------------------------------------------
#' MCMC Sampler for Linear Regression with Global-Local Priors
#'
#' Run the MCMC for linear regression with global-local priors on the regression
#' coefficients. The priors include:
#' \itemize{
#' \item the dynamic horseshoe prior ('DHS');
#' \item the static horseshoe prior ('HS');
#' \item the Bayesian lasso ('BL');
#' \item the normal stochastic volatility model ('SV');
#' \item the normal-inverse-gamma prior ('NIG').
#' }
#' In each case, the prior is a scale mixture of Gaussians.
#' Sampling is accomplished with a (parameter-expanded) Gibbs sampler,
#' mostly relying on a dynamic linear model representation.
#'
#' If p >= n, we use the fast sampling scheme of BHATTACHARYA et al (2016) for the
#' regression coefficients, which is O(n^2*p).
#' If p < n, we use the classical sampling scheme of RUE (2001), which is O(p^3).
#'
#' @param y the \code{n x 1} vector of response variables
#' @param X the \code{n x p} matrix of predictor variables
#' @param prior the prior distribution; must be one of
#' 'DHS' (dynamic horseshoe prior), 'HS' (horseshoe prior), 'BL' (Bayesian lasso), or 'NIG' (normal-inverse-gamma prior)
#' @param marginalSigma logical; if TRUE, marginalize over \code{beta} to sample \code{sigma}
#' @param nsave number of MCMC iterations to record
#' @param nburn number of MCMC iterations to discard (burin-in)
#' @param nskip number of MCMC iterations to skip between saving iterations,
#' i.e., save every (nskip + 1)th draw
#' @param mcmc_params named list of parameters for which we store the MCMC output;
#' must be one or more of:
#' \itemize{
#' \item "mu" (conditional mean)
#' \item "beta" (regression coefficients)
#' \item "yhat" (posterior predictive distribution)
#' \item "evol_sigma_t2" (evolution/prior error variance)
#' \item "obs_sigma_t2" (observation error variance)
#' \item "dhs_phi" (DHS AR(1) coefficient)
#' \item "dhs_mean" (DHS AR(1) unconditional mean)
#' }
#' @param computeDIC logical; if TRUE, compute the deviance information criterion \code{DIC}
#' and the effective number of parameters \code{p_d}
#' @param verbose logical; should R report extra information on progress?
#'
#' @return A named list of the \code{nsave} MCMC samples for the parameters named in \code{mcmc_params}
#'
#' @note The data \code{y} may contain NAs, which will be treated with a simple imputation scheme
#' via an additional Gibbs sampling step. In general, rescaling \code{y} to have unit standard
#' deviation is recommended to avoid numerical issues.
#'
#'
#' @examples
#' \dontrun{
#' # Case 1: p < n
#' n = 200; p = 50; RSNR = 10
#' X = matrix(rnorm(n*p), nr = n)
#' beta_true = simUnivariate(signalName = 'levelshift', p, include_plot = FALSE)$y_true
#' #beta_true = c(rep(1, 10), rep(0, p - 10))
#' y_true = X%*%beta_true; sigma_true = sd(y_true)/RSNR
#' y = y_true + sigma_true*rnorm(n)
#' out = bayesreg_gl(y, X)
#'
#' plot_fitted(lm(y ~ X - 1)$coef,
#'           mu = colMeans(out$beta),
#'           postY = out$beta,
#'           y_true = beta_true)
#'
#' # Case 2: p > n
#' n = 200; p = 1000; RSNR = 10
#' X = matrix(rnorm(n*p), nr = n)
#' #beta_true = simUnivariate(signalName = 'levelshift', p, include_plot = FALSE)$y_true
#' beta_true = c(rep(1, 20), rep(0, p - 20))
#' y_true = X%*%beta_true; sigma_true = sd(y_true)/RSNR
#' y = y_true + sigma_true*rnorm(n)
#' out = bayesreg_gl(y, X)
#'
#' plot_fitted(rep(0, p),
#'           mu = colMeans(out$beta),
#'           postY = out$beta,
#'           y_true = beta_true)
#' }
#' @export
bayesreg_gl = function(y, X, prior = 'DHS',
                       marginalSigma = TRUE,
                       nsave = 1000, nburn = 1000, nskip = 4,
                       mcmc_params = list("mu", "yhat", "beta", "evol_sigma_t2", "obs_sigma_t2", "dhs_phi", "dhs_mean"),
                       computeDIC = TRUE,
                       verbose = TRUE){

  # Convert to upper case:
  prior = toupper(prior)

  # Dimensions:
  n = nrow(X); p = ncol(X)

  # Redefine for package-consistent notation:
  evol_error = prior

  # Begin by checking for missing values, then imputing (for initialization)
  is.missing = which(is.na(y)); any.missing = (length(is.missing) > 0)

  # Impute the active "data"
  y[is.missing] = mean(y, na.rm=TRUE)

  # Initial SD (implicitly assumes a constant mean)
  sigma_e = sd(y, na.rm=TRUE); px_sigma_e = 1

  # Initialize the regression coefficients via sampling (p >= n) or OLS (p < n)
  if(p >= n){
    beta = sampleFastGaussian(Phi = X/sigma_e, Ddiag = rep(.01*sigma_e^2, p), alpha = y/sigma_e)
  } else {
    XtX = crossprod(X); Xty = crossprod(X, y)
    beta = lm(y ~ X - 1)$coef
  }

  # And the conditional expectation:
  mu = as.numeric(X%*%beta)

  # Initialize the evolution error variance paramters:
  evolParams = initEvolParams(omega = beta, evol_error = evol_error)

  # Store the MCMC output in separate arrays (better computation times)
  mcmc_output = vector('list', length(mcmc_params)); names(mcmc_output) = mcmc_params
  if(!is.na(match('mu', mcmc_params)) || computeDIC) post_mu = array(NA, c(nsave, n))
  if(!is.na(match('yhat', mcmc_params))) post_yhat = array(NA, c(nsave, n))
  if(!is.na(match('beta', mcmc_params))) post_beta = array(NA, c(nsave, p))
  if(!is.na(match('obs_sigma_t2', mcmc_params)) || computeDIC) post_obs_sigma_t2 = array(NA, c(nsave, n))
  if(!is.na(match('evol_sigma_t2', mcmc_params))) post_evol_sigma_t2 = array(NA, c(nsave, p))
  if(!is.na(match('dhs_phi', mcmc_params)) && evol_error == "DHS") post_dhs_phi = numeric(nsave)
  if(!is.na(match('dhs_mean', mcmc_params)) && evol_error == "DHS") post_dhs_mean = numeric(nsave)
  post_loglike = numeric(nsave)

  # Total number of MCMC simulations:
  nstot = nburn+(nskip+1)*(nsave)
  skipcount = 0; isave = 0 # For counting

  # Run the MCMC:
  if(verbose) timer0 = proc.time()[3] # For timing the sampler
  for(nsi in 1:nstot){

    # Impute missing values, if any:
    if(any.missing) {
      y[is.missing] = mu[is.missing] + sigma_e*rnorm(length(is.missing))
      if(p < n) Xty = crossprod(X, y)
    }
    # Sample the dynamic regression coefficients, beta:
    if(p >= n){
      beta = sampleFastGaussian(Phi = X/sigma_e,
                                Ddiag = sigma_e^2*as.numeric(evolParams$sigma_wt^2),
                                alpha = y/sigma_e)
    } else {
      # Sample the dynamic regression coefficients, beta:
      ch_Q = chol(sigma_e^-2*XtX + sigma_e^-2*diag(as.numeric(evolParams$sigma_wt^-2)))
      ell_beta = sigma_e^-2*Xty
      beta = backsolve(ch_Q,
                       forwardsolve(t(ch_Q), ell_beta) +
                         rnorm(p))
    }
    # And the observation error standard deviation:
    if(marginalSigma){
      # Marginalize over beta: more efficient, but less stable
      sigma_e = 1/sqrt(rgamma(n = 1,
                              shape = n/2 + 1/2,
                              #rate = px_sigma_e + crossprod(y, chol2inv(chol(crossprod(as.numeric(evolParams$sigma_wt)*t(X)) + diag(n))))%*%y))
                              rate = px_sigma_e + sum(forwardsolve(t(chol(crossprod(as.numeric(evolParams$sigma_wt)*t(X)) + diag(n))), y)^2)))

    } else {
      sigma_e = 1/sqrt(rgamma(n = 1,
                              shape = 1/2 + n/2 + p/2,
                              rate = px_sigma_e + sum((y - mu)^2, na.rm=TRUE)/2 +  sum((beta/evolParams$sigma_wt)^2)/2))
    }
    # Parameter-expanded piece
    px_sigma_e = rgamma(n = 1, shape = 1/2 + 1/2, rate = 1/sigma_e^2 + 1)

    # And the conditional expectation:
    mu = as.numeric(X%*%beta)

    # Sample the evolution error variance (and associated parameters):
    evolParams = sampleEvolParams(omega = beta/sigma_e, evolParams, 1/sqrt(n), evol_error)

    # Store the MCMC output:
    if(nsi > nburn){
      # Increment the skip counter:
      skipcount = skipcount + 1

      # Save the iteration:
      if(skipcount > nskip){
        # Increment the save index
        isave = isave + 1

        # Save the MCMC samples:
        if(!is.na(match('mu', mcmc_params)) || computeDIC) post_mu[isave,] = mu
        if(!is.na(match('yhat', mcmc_params))) post_yhat[isave,] = mu + sigma_e*rnorm(n)
        if(!is.na(match('beta', mcmc_params))) post_beta[isave,] = beta
        if(!is.na(match('obs_sigma_t2', mcmc_params)) || computeDIC) post_obs_sigma_t2[isave,] = sigma_e^2
        if(!is.na(match('evol_sigma_t2', mcmc_params))) post_evol_sigma_t2[isave,] = evolParams$sigma_wt^2
        if(!is.na(match('dhs_phi', mcmc_params)) && evol_error == "DHS") post_dhs_phi[isave] = evolParams$dhs_phi
        if(!is.na(match('dhs_mean', mcmc_params)) && evol_error == "DHS") post_dhs_mean[isave] = evolParams$dhs_mean
        post_loglike[isave] = sum(dnorm(y, mean = mu, sd = sigma_e, log = TRUE))

        # And reset the skip counter:
        skipcount = 0
      }
    }
    if(verbose) computeTimeRemaining(nsi, timer0, nstot, nrep = 1000)
  }

  if(!is.na(match('mu', mcmc_params))) mcmc_output$mu = post_mu
  if(!is.na(match('yhat', mcmc_params))) mcmc_output$yhat = post_yhat
  if(!is.na(match('beta', mcmc_params))) mcmc_output$beta = post_beta
  if(!is.na(match('obs_sigma_t2', mcmc_params))) mcmc_output$obs_sigma_t2 = post_obs_sigma_t2
  if(!is.na(match('evol_sigma_t2', mcmc_params))) mcmc_output$evol_sigma_t2 = post_evol_sigma_t2
  if(!is.na(match('dhs_phi', mcmc_params)) && evol_error == "DHS") mcmc_output$dhs_phi = post_dhs_phi
  if(!is.na(match('dhs_mean', mcmc_params)) && evol_error == "DHS") mcmc_output$dhs_mean = post_dhs_mean

  # Also include the log-likelihood:
  mcmc_output$loglike = post_loglike

  if(computeDIC){
    # Log-likelihood evaluated at posterior means:
    loglike_hat = sum(dnorm(y,
                            mean = colMeans(post_mu),
                            sd = colMeans(sqrt(post_obs_sigma_t2)),
                            log = TRUE))
    # Effective number of parameters (Note: two options)
    p_d = c(2*(loglike_hat - mean(post_loglike)),
            2*var(post_loglike))
    # DIC:
    DIC = -2*loglike_hat + 2*p_d

    # Store the DIC and the effective number of parameters (p_d)
    mcmc_output$DIC = DIC; mcmc_output$p_d = p_d
  }
  if(verbose) print(paste('Total time: ', round((proc.time()[3] - timer0)), 'seconds'))

  return (mcmc_output);
}
#' #----------------------------------------------------------------------------
#' #' MCMC Sampler for Bayesian Trend Filtering: Regression with VARIMA
#' #'
#' #' Run the MCMC for Bayesian trend filtering regression with a penalty on the
#' #' innovations from a VARIMA(1,1,0) model on the dynamic regression coefficients.
#' #' The penalty is determined by the prior on the evolution errors, which include:
#' #' \itemize{
#' #' \item the dynamic horseshoe prior ('DHS');
#' #' \item the static horseshoe prior ('HS');
#' #' \item the Bayesian lasso ('BL');
#' #' \item the normal stochastic volatility model ('SV');
#' #' \item the normal-inverse-gamma prior ('NIG').
#' #' }
#' #' In each case, the evolution error is a scale mixture of Gaussians.
#' #' Sampling is accomplished with a (parameter-expanded) Gibbs sampler,
#' #' mostly relying on a dynamic linear model representation.
#' #'
#' #' @param y the \code{T x 1} vector of time series observations
#' #' @param X the \code{T x p} matrix of time series predictors
#' #' @param evol_error the evolution error distribution; must be one of
#' #' 'DHS' (dynamic horseshoe prior), 'HS' (horseshoe prior), 'BL' (Bayesian lasso), or 'NIG' (normal-inverse-gamma prior)
#' #' @param useObsSV logical; if TRUE, include a (normal) stochastic volatility model
#' #' for the observation error variance
#' #' @param nsave number of MCMC iterations to record
#' #' @param nburn number of MCMC iterations to discard (burin-in)
#' #' @param nskip number of MCMC iterations to skip between saving iterations,
#' #' i.e., save every (nskip + 1)th draw
#' #' @param mcmc_params named list of parameters for which we store the MCMC output;
#' #' must be one or more of:
#' #' \itemize{
#' #' \item "mu" (conditional mean)
#' #' \item "yhat" (posterior predictive distribution)
#' #' \item "beta" (dynamic regression coefficients)
#' #' \item "Psi" (the VAR coefficient matrix)
#' #' \item "evol_sigma_t2" (evolution error variance)
#' #' \item "obs_sigma_t2" (observation error variance)
#' #' \item "dhs_phi" (DHS AR(1) coefficient)
#' #' \item "dhs_mean" (DHS AR(1) unconditional mean)
#' #' }
#' #' @param computeDIC logical; if TRUE, compute the deviance information criterion \code{DIC}
#' #' and the effective number of parameters \code{p_d}
#' #' @param verbose logical; should R report extra information on progress?
#' #'
#' #' @return A named list of the \code{nsave} MCMC samples for the parameters named in \code{mcmc_params}
#' #'
#' #' @note The data \code{y} may contain NAs, which will be treated with a simple imputation scheme
#' #' via an additional Gibbs sampling step. In general, rescaling \code{y} to have unit standard
#' #' deviation is recommended to avoid numerical issues.
#' #'
#' # @import KFAS vars
#' btf_reg_varima = function(y, X, evol_error = 'DHS', useObsSV = FALSE,
#'                           nsave = 1000, nburn = 1000, nskip = 4,
#'                           mcmc_params = list("mu", "yhat","beta","Psi", "evol_sigma_t2", "obs_sigma_t2", "dhs_phi", "dhs_mean"),
#'                           computeDIC = TRUE,
#'                           verbose = TRUE){
#'
#'   # Convert to upper case:
#'   evol_error = toupper(evol_error)
#'
#'   # Time points (in [0,1])
#'   T = length(y); t01 = seq(0, 1, length.out=T);
#'
#'   # Begin by checking for missing values, then imputing (for initialization)
#'   is.missing = which(is.na(y)); any.missing = (length(is.missing) > 0)
#'
#'   # Impute the active "data"
#'   if(any.missing) y[is.missing] = mean(y, na.rm=TRUE)
#'
#'   # Number of predictors:
#'   p = ncol(X)
#'
#'   # Useful construction of array:
#'   X.arr = array(0, c(1, 2*p, T)); for(i in 1:T) X.arr[,1:p,i] = X[i,]
#'
#'   # Initial SD (implicitly assumes a constant mean)
#'   sigma_e = sd(y, na.rm=TRUE); sigma_et = rep(sigma_e, T)
#'
#'   # Initialize the KFAS model object:
#'   kfas_model = SSModel(y~-1+SSMcustom(Z = array(X.arr, c(1, 2*p, T)),
#'                                       Q = array(diag(2*p), c(2*p, 2*p, T)),
#'                                       T = array(diag(2*p), c(2*p, 2*p, T)),
#'                                       R = array(diag(2*p), c(2*p, 2*p, T)),
#'                                       P1 = diag(10^4, 2*p)), H = array(1, c(1,1,T)))
#'   kfas_model$T[(p+1):(2*p), 1:p,] = diag(p); kfas_model$T[(p+1):(2*p), (p+1):(2*p),] = 0
#'   kfas_model$R[(p+1):(2*p), (p+1):(2*p), ] = 0
#'
#'   # Some initial values for prelim simulations:
#'   kfas_model$H = array(sigma_et^2, c(1,1,T)); #for(j in 1:p) kfas_model$Q[j,j,-T] = 0.01*sigma_e^2
#'
#'   # Simulate for initial values
#'   beta_samp = as.matrix(simulateSSM(kfas_model, "states", nsim = 1, antithetics=FALSE, filtered=FALSE)[,,1])
#'   beta = beta_samp[,1:p]; beta0 = matrix(beta_samp[1,], nr = 1);
#'
#'   # Conditional mean:
#'   mu = rowSums(X*beta)
#'
#'   # Initialize the VAR model:
#'   diff_beta = diff(beta)
#'   psi = unlist(lapply(VAR(diff_beta, p=1, "none")$varresult, coef))
#'   Psi = matrix(psi, nrow=p, byrow = FALSE)
#'
#'   # Prior mean and variance for psi = vec(Psi):
#'   mu_psi = matrix(diag(1, p))
#'   sigma_psi = sum((psi - mu_psi)^2)/p^2
#'   px_sigma_psi = 1
#'
#'   # And the residuals:
#'   omega = diff_beta[-1,] - t(tcrossprod(Psi, as.matrix(diff_beta[-(T-1),])))
#'   omega = rbind(as.numeric(beta[2,] - tcrossprod(kfas_model$T[1:p,,1], beta0)),
#'                 omega)
#'
#'   # Initialize the evolution error variance paramters:
#'   evolParams = initEvolParams(omega, evol_error = evol_error)
#'
#'   # Initial variance parameters:
#'   evolParams0 = initEvol0(beta0, commonSD = FALSE)
#'
#'   # SV parameters, if necessary:
#'   if(useObsSV) {svParams = initSV(y - mu); sigma_et = svParams$sigma_wt}
#'
#'   # Store the MCMC output in separate arrays (better computation times)
#'   mcmc_output = vector('list', length(mcmc_params)); names(mcmc_output) = mcmc_params
#'   if(!is.na(match('mu', mcmc_params)) || computeDIC) post_mu = array(NA, c(nsave, T))
#'   if(!is.na(match('yhat', mcmc_params))) post_yhat = array(NA, c(nsave, T))
#'   if(!is.na(match('beta', mcmc_params))) post_beta = array(NA, c(nsave, T, p))
#'   if(!is.na(match('obs_sigma_t2', mcmc_params)) || computeDIC) post_obs_sigma_t2 = array(NA, c(nsave, T))
#'   if(!is.na(match('evol_sigma_t2', mcmc_params))) post_evol_sigma_t2 = array(NA, c(nsave, T, p))
#'   if(!is.na(match('dhs_phi', mcmc_params)) && evol_error == "DHS") post_dhs_phi = array(NA, c(nsave, p))
#'   if(!is.na(match('dhs_mean', mcmc_params)) && evol_error == "DHS") post_dhs_mean = array(NA, c(nsave, p))
#'   if(!is.na(match('Psi', mcmc_params))) post_Psi = array(NA, c(nsave, p, p))
#'
#'   post_loglike = numeric(nsave)
#'
#'   # Total number of MCMC simulations:
#'   nstot = nburn+(nskip+1)*(nsave)
#'   skipcount = 0; isave = 0 # For counting
#'
#'   # Run the MCMC:
#'   if(verbose) timer0 = proc.time()[3] # For timing the sampler
#'   for(nsi in 1:nstot){
#'
#'     # Impute missing values, if any:
#'     if(any.missing) y[is.missing] = mu[is.missing] + sigma_et[is.missing]*rnorm(length(is.missing))
#'
#'     # Sample the dynamic regression coefficients, beta:
#'     evolParams$sigma_wt[which(evolParams$sigma_wt > 10^3, arr.ind = TRUE)] = 10^3
#'     kfas_model$T[1:p, 1:p,] = Psi + diag(p); kfas_model$T[1:p, (p+1):(2*p),] = -Psi
#'     kfas_model$H = array(sigma_et^2, c(1,1,T)); for(j in 1:p) kfas_model$Q[j,j,-T] = evolParams$sigma_wt[,j]^2
#'     diag(kfas_model$P1) = evolParams0$sigma_w0^2
#'     beta_samp = as.matrix(simulateSSM(kfas_model, "states", nsim = 1, antithetics=FALSE, filtered=FALSE)[,,1])
#'     beta = beta_samp[,1:p]; beta0 = matrix(beta_samp[1,], nr = 1);
#'
#'     # Conditional mean:
#'     mu = rowSums(X*beta)
#'
#'     # Sample the VAR matrix:
#'     #diff_beta = diff(beta)
#'     diff_beta = rbind(beta0[1,1:p] - beta0[1,-(1:p)],
#'                        diff(beta))
#'
#'     # This is going to be slow anyway:
#'     #Q_psi = diag(10^-6, p^2); ell_psi = 0
#'     Q_psi = diag(sigma_psi^-2, p^2); ell_psi = mu_psi/sigma_psi^2
#'     for(i in 1:(T-2)) {
#'       Q_psi = Q_psi + kronecker(tcrossprod(diff_beta[i,]),
#'                                 diag(as.numeric(1/evolParams$sigma_wt[i,]^2)))
#'       ell_psi = ell_psi + matrix(tcrossprod(1/evolParams$sigma_wt[i,]^2*diff_beta[i+1, ], diff_beta[i, ]))
#'     }
#'     chQ_psi = chol(Q_psi)
#'     psi = backsolve(chQ_psi,
#'                     forwardsolve(t(chQ_psi), ell_psi) +
#'                       rnorm(p^2))
#'     Psi = matrix(psi,nrow=p, byrow=FALSE)
#'
#'     # SD term:
#'     sigma_psi = 1/sqrt(rgamma(n = 1,
#'                               shape = p^2/2 + 1/2,
#'                               rate = sum((psi - mu_psi)^2)/2 + px_sigma_psi))
#'     px_sigma_psi = rgamma(n = 1, shape = 1/2 + 1/2, rate = 1/sigma_psi^2 + 1)
#'
#'     # And the residuals:
#'     omega = diff_beta[-1,] - t(tcrossprod(Psi, as.matrix(diff_beta[-T,])))
#'     #omega = rbind(as.numeric(beta[2,] - tcrossprod(kfas_model$T[1:p,,1], beta0)), omega)
#'
#'     # Sample the initial variance parameters:
#'     evolParams0 = sampleEvol0(beta0, evolParams0, A = 1, commonSD = FALSE)
#'
#'     # Sample the (observation and evolution) variances and associated parameters:
#'     if(useObsSV){
#'       # Evolution error variance + params:
#'       evolParams = sampleEvolParams(omega, evolParams, 1/sqrt(T*p), evol_error)
#'
#'       # Observation error variance + params:
#'       svParams = sampleSVparams(omega = y - mu, svParams = svParams)
#'       sigma_et = svParams$sigma_wt
#'
#'     } else {
#'       # Evolution error variance + params:
#'       evolParams = sampleEvolParams(omega, evolParams, sigma_e/sqrt(T*p), evol_error)
#'
#'       # Sample the observation error SD:
#'       if(evol_error == 'DHS') {
#'         sigma_e = uni.slice(sigma_e, g = function(x){
#'           -(T+2)*log(x) - 0.5*sum((y - mu)^2, na.rm=TRUE)/x^2 - log(1 + (sqrt(T*p)*exp(evolParams$dhs_mean0/2)/x)^2)
#'         }, lower = 0, upper = Inf)[1]
#'       }
#'       if(evol_error == 'HS') sigma_e = 1/sqrt(rgamma(n = 1, shape = T/2 + length(evolParams$xiLambda), rate = sum((y - mu)^2, na.rm=TRUE)/2 + T*p*sum(evolParams$xiLambda)))
#'       if(evol_error == 'BL') sigma_e = 1/sqrt(rgamma(n = 1, shape = T/2 + length(evolParams$tau_j)/2, rate = sum((y - mu)^2, na.rm=TRUE)/2 + T*p*sum((omega/evolParams$tau_j)^2)/2))
#'       if((evol_error == 'NIG') || (evol_error == 'SV')) sigma_e = 1/sqrt(rgamma(n = 1, shape = T/2, rate = sum((y - mu)^2, na.rm=TRUE)/2))
#'
#'       # Replicate for coding convenience:
#'       sigma_et = rep(sigma_e, T)
#'     }
#'
#'     # Store the MCMC output:
#'     if(nsi > nburn){
#'       # Increment the skip counter:
#'       skipcount = skipcount + 1
#'
#'       # Save the iteration:
#'       if(skipcount > nskip){
#'         # Increment the save index
#'         isave = isave + 1
#'
#'         # Save the MCMC samples:
#'         if(!is.na(match('mu', mcmc_params)) || computeDIC) post_mu[isave,] = mu
#'         if(!is.na(match('yhat', mcmc_params))) post_yhat[isave,] = mu + sigma_et*rnorm(T)
#'         if(!is.na(match('beta', mcmc_params))) post_beta[isave,,] = beta
#'         if(!is.na(match('obs_sigma_t2', mcmc_params)) || computeDIC) post_obs_sigma_t2[isave,] = sigma_et^2
#'         if(!is.na(match('evol_sigma_t2', mcmc_params))) post_evol_sigma_t2[isave,,] = rbind(matrix(evolParams0$sigma_w0[1:p]^2, nr = 1), evolParams$sigma_wt^2)
#'         if(!is.na(match('dhs_phi', mcmc_params)) && evol_error == "DHS") post_dhs_phi[isave,] = evolParams$dhs_phi
#'         if(!is.na(match('dhs_mean', mcmc_params)) && evol_error == "DHS") post_dhs_mean[isave,] = evolParams$dhs_mean
#'         if(!is.na(match('Psi', mcmc_params))) post_Psi[isave,,] = Psi
#'         post_loglike[isave] = sum(dnorm(y, mean = mu, sd = sigma_et, log = TRUE))
#'
#'         # And reset the skip counter:
#'         skipcount = 0
#'       }
#'     }
#'     if(verbose) computeTimeRemaining(nsi, timer0, nstot, nrep = 1000)
#'   }
#'
#'   if(!is.na(match('mu', mcmc_params))) mcmc_output$mu = post_mu
#'   if(!is.na(match('yhat', mcmc_params))) mcmc_output$yhat = post_yhat
#'   if(!is.na(match('beta', mcmc_params))) mcmc_output$beta = post_beta
#'   if(!is.na(match('obs_sigma_t2', mcmc_params))) mcmc_output$obs_sigma_t2 = post_obs_sigma_t2
#'   if(!is.na(match('evol_sigma_t2', mcmc_params))) mcmc_output$evol_sigma_t2 = post_evol_sigma_t2
#'   if(!is.na(match('dhs_phi', mcmc_params)) && evol_error == "DHS") mcmc_output$dhs_phi = post_dhs_phi
#'   if(!is.na(match('dhs_mean', mcmc_params)) && evol_error == "DHS") mcmc_output$dhs_mean = post_dhs_mean
#'   if(!is.na(match('Psi', mcmc_params))) mcmc_output$Psi = post_Psi
#'
#'   # Also include the log-likelihood:
#'   mcmc_output$loglike = post_loglike
#'
#'   if(computeDIC){
#'     # Log-likelihood evaluated at posterior means:
#'     loglike_hat = sum(dnorm(y,
#'                             mean = colMeans(post_mu),
#'                             sd = colMeans(sqrt(post_obs_sigma_t2)),
#'                             log = TRUE))
#'     # Effective number of parameters (Note: two options)
#'     p_d = c(2*(loglike_hat - mean(post_loglike)),
#'             2*var(post_loglike))
#'     # DIC:
#'     DIC = -2*loglike_hat + 2*p_d
#'
#'     # Store the DIC and the effective number of parameters (p_d)
#'     mcmc_output$DIC = DIC; mcmc_output$p_d = p_d
#'   }
#'   if(verbose) print(paste('Total time: ', round((proc.time()[3] - timer0)), 'seconds'))
#'
#'   return (mcmc_output);
#' }
