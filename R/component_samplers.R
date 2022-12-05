#----------------------------------------------------------------------------
#' Sampler for first or second order random walk (RW) Gaussian dynamic linear model (DLM)
#'
#' Compute one draw of the \code{T x 1} state variable \code{mu} in a DLM using back-band substitution methods.
#' This model is equivalent to the Bayesian trend filtering (BTF) model, assuming appropriate
#' (shrinkage/sparsity) priors for the evolution errors.
#'
#' @param y the \code{T x 1} vector of time series observations
#' @param obs_sigma_t2 the \code{T x 1} vector of observation error variances
#' @param evol_sigma_t2 the \code{T x 1} vector of evolution error variances
#' @param D the degree of differencing (one or two)
#' @param chol0 (optional) the \code{m x m} matrix of initial Cholesky factorization;
#' if NULL, use the \code{Matrix} package for sampling, otherwise use the \code{spam} package
#' @return \code{T x 1} vector of simulated states
#'
#' @note Missing entries (NAs) are not permitted in \code{y}. Imputation schemes are available.
#'
#' @examples
#' # Simulate some data:
#' T = 1000
#' y = seq(0, 10, length.out = T) + rnorm(T)
#' plot(y) # plot the data
#'
#' obs_sigma_t2 = rep(1, T)  # These could be static or dynamic
#' evol_sigma_t2 = rep(0.001, T)
#'
#' # Simulate one draw of the states:
#' mu = sampleBTF(y = y, obs_sigma_t2, evol_sigma_t2, D = 1)
#' lines(mu, lwd=8, col='blue') # add the states to the plot
#'
#' @import Matrix spam
#' @export
sampleBTF = function(y, obs_sigma_t2, evol_sigma_t2, D = 1, chol0 = NULL){

  # Some quick checks:
  if((D < 0) || (D != round(D)))  stop('D must be a positive integer')

  if(any(is.na(y))) stop('y cannot contain NAs')

  T = length(y)

  # Linear term:
  linht = y/obs_sigma_t2

  # Quadratic terms and solutions are computed differently, depending on D:

  if(D == 0){
    # Special case: no differencing

    # Posterior SDs and posterior means:
    postSD = 1/sqrt(1/obs_sigma_t2 + 1/evol_sigma_t2)
    postMean = (linht)*postSD^2

    # Sample the states:
    mu = rnorm(n = T, mean = postMean, sd = postSD)

  } else {
    # All other cases: positive integer differencing (D = 1 or D = 2)

    # Quadratic term (D = 1 or D = 2)
    QHt_Matrix = build_Q(obs_sigma_t2 = obs_sigma_t2, evol_sigma_t2 = evol_sigma_t2, D = D)

    if(!is.null(chol0)){
      # New sampler, based on spam package:

      # Sample the states:
      mu = matrix(rmvnorm.canonical(n = 1,
                                       b = linht,
                                       Q = as.spam.dgCMatrix(as(QHt_Matrix, "dgCMatrix")),
                                       Rstruct = chol0))
    } else {
      # Original sampler, based on Matrix package:

      # Cholesky of Quadratic term:
      chQht_Matrix = Matrix::chol(QHt_Matrix)

      # Sample the states:
      mu = as.matrix(Matrix::solve(chQht_Matrix,Matrix::solve(Matrix::t(chQht_Matrix), linht) + rnorm(T)))

    }
  }

  # And return the states:
  mu
}
#----------------------------------------------------------------------------
#' Sampler for first or second order random walk (RW) Gaussian dynamic linear model (DLM)
#' with additional shrinkage to zero
#'
#' Compute one draw of the \code{T x 1} state variable \code{mu} in a DLM using back-band substitution methods.
#' This model is equivalent to the Bayesian trend filtering (BTF) model, assuming appropriate
#' (shrinkage/sparsity) priors for the evolution errors, with an additional shrinkage-to-zero prior.
#'
#' @param y the \code{T x 1} vector of time series observations
#' @param obs_sigma_t2 the \code{T x 1} vector of observation error variances
#' @param evol_sigma_t2 the \code{T x 1} vector of evolution error variances
#' @param zero_sigma_t2 the \code{T x 1} vector of shrink-to-zero variances
#' @param D the degree of differencing (one or two)
#' @param chol0 (optional) the \code{m x m} matrix of initial Cholesky factorization;
#' if NULL, use the \code{Matrix} package for sampling, otherwise use the \code{spam} package
#' @return \code{T x 1} vector of simulated states
#'
#' @note Missing entries (NAs) are not permitted in \code{y}. Imputation schemes are available.
#'
#' @examples
#' # Simulate some data:
#' T = 1000
#' y = seq(0, 10, length.out = T) + rnorm(T)
#' plot(y) # plot the data
#'
#' obs_sigma_t2 = rep(1, T)  # These could be static or dynamic
#' evol_sigma_t2 = rep(0.001, T)
#' zero_sigma_t2 = rep(1, T)
#'
#' # Simulate one draw of the states:
#' mu = sampleBTF_sparse(y = y, obs_sigma_t2, evol_sigma_t2, zero_sigma_t2, D = 1)
#' lines(mu, lwd=8, col='blue') # add the states to the plot
#'
#' @import Matrix spam
#' @export
sampleBTF_sparse = function(y,
                            obs_sigma_t2,
                            evol_sigma_t2,
                            zero_sigma_t2,
                            D = 1, chol0 = NULL){

  # Some quick checks:
  if((D < 0) || (D != round(D)))  stop('D must be a positive integer')

  if(any(is.na(y))) stop('y cannot contain NAs')

  T = length(y)

  # Linear term:
  linht = y/obs_sigma_t2

  # Quadratic terms and solutions are computed differently, depending on D:

  if(D == 0){
    # Special case: no differencing

    # Posterior SDs and posterior means:
    postSD = 1/sqrt(1/obs_sigma_t2 + 1/evol_sigma_t2 + 1/zero_sigma_t2)
    postMean = (linht)*postSD^2

    # Sample the states:
    mu = rnorm(n = T, mean = postMean, sd = postSD)

  } else {
    # All other cases: positive integer differencing (D = 1 or D = 2)

    # Quadratic term (D = 1 or D = 2)
    #QHt_Matrix = build_Q(obs_sigma_t2 = obs_sigma_t2,
    QHt_Matrix = build_Q(obs_sigma_t2 = 1/(1/obs_sigma_t2 + 1/zero_sigma_t2),
                         evol_sigma_t2 = evol_sigma_t2,
                         D = D)

    if(!is.null(chol0)){
      # New sampler, based on spam package:

      # Sample the states:
      mu = matrix(rmvnorm.canonical(n = 1,
                                    b = linht,
                                    Q = as.spam.dgCMatrix(as(QHt_Matrix, "dgCMatrix")),
                                    Rstruct = chol0))
    } else {
      # Original sampler, based on Matrix package:

      # Cholesky of Quadratic term:
      chQht_Matrix = Matrix::chol(QHt_Matrix)

      # Sample the states:
      mu = as.matrix(Matrix::solve(chQht_Matrix,Matrix::solve(Matrix::t(chQht_Matrix), linht) + rnorm(T)))

    }
  }

  # And return the states:
  mu
}
#----------------------------------------------------------------------------
#' Sampler for first or second order random walk (RW) Gaussian dynamic linear model (DLM)
#'
#' Compute one draw of the \code{T x p} state variable \code{beta} in a DLM using back-band substitution methods.
#' This model is equivalent to the Bayesian trend filtering (BTF) model applied to \code{p}
#' dynamic regression coefficients corresponding to the design matrix \code{X},
#' assuming appropriate (shrinkage/sparsity) priors for the evolution errors.
#'
#' @param y the \code{T x 1} vector of time series observations
#' @param X the \code{T x p} matrix of time series predictors
#' @param obs_sigma_t2 the \code{T x 1} vector of observation error variances
#' @param evol_sigma_t2 the \code{T x p} matrix of evolution error variances
#' @param XtX the \code{Tp x Tp} matrix of X'X (one-time cost; see ?build_XtX)
#' @param D the degree of differencing (one or two)
#' @param chol0 (optional) the \code{m x m} matrix of initial Cholesky factorization;
#' if NULL, use the \code{Matrix} package for sampling, otherwise use the \code{spam} package
#' @return \code{T x p} matrix of simulated dynamic regression coefficients \code{beta}
#'
#' @note Missing entries (NAs) are not permitted in \code{y}. Imputation schemes are available.
#'
#' @import Matrix spam
#' @export
sampleBTF_reg = function(y, X, obs_sigma_t2, evol_sigma_t2, XtX, D = 1, chol0 = NULL){

  # Some quick checks:
  if((D < 0) || (D != round(D)))  stop('D must be a positive integer')

  if(any(is.na(y))) stop('y cannot contain NAs')

  # Dimensions of X:
  T = nrow(X); p = ncol(X)

  if(D == 1){
    # Lagged version of transposed precision matrix, with zeros as appropriate (needed below)
    t_evol_prec_lag_mat = matrix(0, nr = p, nc = T);
    t_evol_prec_lag_mat[,1:(T-1)] = t(1/evol_sigma_t2[-1,])

    # Diagonal of quadratic term:
    Q_diag = matrix(t(1/evol_sigma_t2) + t_evol_prec_lag_mat)

    # Off-diagonal of quadratic term:
    Q_off = matrix(-t_evol_prec_lag_mat)[-(T*p)]

    # Quadratic term:
    Qevol = bandSparse(T*p, k = c(0,p), diag = list(Q_diag, Q_off), symm = TRUE)

    # For checking via direct computation:
    # H1 = bandSparse(T, k = c(0,-1), diag = list(rep(1, T), rep(-1, T)), symm = FALSE)
    # IH = kronecker(as.matrix(H1), diag(p));
    # Q0 = t(IH)%*%diag(as.numeric(1/matrix(t(evol_sigma_t2))))%*%(IH)
    # print(sum((Qevol - Q0)^2))

  } else {
    if(D == 2){

      # Lagged x2 version of transposed precision matrix (recurring term)
      t_evol_prec_lag2 = t(1/evol_sigma_t2[-(1:2),])

      # Diagonal of quadratic term:
      Q_diag = t(1/evol_sigma_t2)
      Q_diag[,2:(T-1)] = Q_diag[,2:(T-1)] + 4*t_evol_prec_lag2
      Q_diag[,1:(T-2)] = Q_diag[,1:(T-2)] + t_evol_prec_lag2
      Q_diag = matrix(Q_diag)

      # Off-diagonal (1) of quadratic term:
      Q_off_1 = matrix(0, nr = p, nc = T);
      Q_off_1[,1] = -2/evol_sigma_t2[3,]
      Q_off_1[,2:(T-1)] = Q_off_1[,2:(T-1)] + -2*t_evol_prec_lag2
      Q_off_1[,2:(T-2)] = Q_off_1[,2:(T-2)] + -2*t_evol_prec_lag2[,-1]
      Q_off_1 = matrix(Q_off_1)

      # Off-diagonal (2) of quadratic term:
      Q_off_2 =  matrix(0, nr = p, nc = T); Q_off_2[,1:(T-2)] = t_evol_prec_lag2
      Q_off_2 = matrix(Q_off_2)

      # Quadratic term:
      Qevol = bandSparse(T*p, k = c(0, p, 2*p), diag = list(Q_diag, Q_off_1, Q_off_2), symm = TRUE)

      # For checking via direct computation:
      # H2 = bandSparse(T, k = c(0,-1, -2), diag = list(rep(1, T), c(0, rep(-2, T-1)), rep(1, T)), symm = FALSE)
      # IH = kronecker(as.matrix(H2), diag(p));
      # Q0 = t(IH)%*%diag(as.numeric(1/matrix(t(evol_sigma_t2))))%*%(IH)
      # print(sum((Qevol - Q0)^2))

    } else stop('sampleBTF_reg() requires D=1 or D=2')
  }

  # Quadratic term:
  Qobs = 1/rep(obs_sigma_t2, each = p)*XtX
  Qpost = Qobs + Qevol

  # Linear term:
  linht = matrix(t(X*as.numeric(y/obs_sigma_t2))) #matrix(t(X*tcrossprod(y/obs_sigma_t2, rep(1,p))))

  if(!is.null(chol0)){
    # Use spam sampler (new version)

    # Convert to spam object:
    QHt_Matrix = as.spam.dgCMatrix(as(Qpost, "dgCMatrix"))

    # NOTE: reorder (opposite of log-vol!)
    beta = matrix(rmvnorm.canonical(n = 1,
                                    b = linht,
                                    Q = QHt_Matrix,
                                    Rstruct = chol0),
                  nrow = T, byrow = TRUE)
  } else {
    # Use original sampler:

    # Cholesky:
    chQht_Matrix = Matrix::chol(Qpost)

    # NOTE: reorder (opposite of log-vol!)
    beta = matrix(Matrix::solve(chQht_Matrix,Matrix::solve(Matrix::t(chQht_Matrix), linht) + rnorm(T*p)), nr = T, byrow = TRUE)
  }
  beta
}
#----------------------------------------------------------------------------
#' (Backfitting) Sampler for first or second order random walk (RW) Gaussian dynamic linear model (DLM)
#'
#' Compute one draw of the \code{T x p} state variable \code{beta} in a DLM using back-band substitution methods.
#' This model is equivalent to the Bayesian trend filtering (BTF) model applied to \code{p}
#' dynamic regression coefficients corresponding to the design matrix \code{X},
#' assuming appropriate (shrinkage/sparsity) priors for the evolution errors. The sampler
#' here uses a backfitting method that draws each predictor j=1,...,p conditional on the
#' other predictors (and coefficients), which leads to a faster \code{O(Tp)} algorithm.
#' However, the MCMC may be less efficient.
#'
#' @param y the \code{T x 1} vector of time series observations
#' @param X the \code{T x p} matrix of time series predictors
#' @param beta the \code{T x p} matrix of previous dynamic regression coefficients
#' @param obs_sigma_t2 the \code{T x 1} vector of observation error variances
#' @param evol_sigma_t2 the \code{T x p} matrix of evolution error variances
#' @param D the degree of differencing (one or two)
#' @return \code{T x p} matrix of simulated dynamic regression coefficients \code{beta}
#'
#' @note Missing entries (NAs) are not permitted in \code{y}. Imputation schemes are available.
#'
#' @import Matrix
#' @export
sampleBTF_reg_backfit = function(y, X, beta, obs_sigma_t2, evol_sigma_t2, D = 1){

  # Some quick checks:
  if((D < 0) || (D != round(D)))  stop('D must be a positive integer')

  if(any(is.na(y))) stop('y cannot contain NAs')

  # Dimensions of X:
  T = nrow(X); p = ncol(X)

  # Sample each predictor curve via backfitting:
  for(j in sample(1:p,p)){
  #for(j in 1:p){
    # Subtract off non-j terms:
    y_nj = y - rowSums(X[,-j]*beta[,-j])

    # Linear term:
    linht = y_nj*X[,j]/obs_sigma_t2

    # Quadratic terms and solutions are computed differently, depending on D:

    if(D == 0){
      # Special case: no differencing

      # Posterior SDs and posterior means:
      postSD = 1/sqrt((X[,j]^2)/obs_sigma_t2 + 1/evol_sigma_t2)
      postMean = (linht)*postSD^2

      # Sample the states:
      beta[,j] = rnorm(n = T, mean = postMean, sd = postSD)

    } else {
      # Quadratic term (D = 1 or D = 2)
      # The likelihood precision term is simply X[,j]^2/obs_sigma_t2, so invert here:
      QHt_Matrix = build_Q(obs_sigma_t2 = obs_sigma_t2/X[,j]^2,
                           evol_sigma_t2 = evol_sigma_t2[,j],
                           D = D)

      # Cholesky of Quadratic term:
      chQht_Matrix = Matrix::chol(QHt_Matrix)

      # Sample the states:
      beta[,j] = as.matrix(Matrix::solve(chQht_Matrix,Matrix::solve(Matrix::t(chQht_Matrix), linht) + rnorm(T)))
    }
  }
  beta
}
#----------------------------------------------------------------------------
#' Sampler for first or second order random walk (RW) Gaussian dynamic linear model (DLM)
#'
#' Compute one draw of the \code{p x 1} B-spline basis coefficients \code{beta} in a DLM using
#' back-band substitution methods. The coefficients are penalized with a prior on the D = 0, D = 1, or
#' D = 2 differences. This model is equivalent to the Bayesian trend filtering (BTF) model
#' applied to \code{p x 1} vector of equally-spaced B-spline coefficients, with the basis matrix
#' serving as a design matrix in the observation equation.
#'
#' @param y the \code{T x 1} vector of time series observations
#' @param X the \code{T x p} basis matrix
#' @param obs_sigma2 the scalar observation error variance
#' @param evol_sigma_t2 the \code{p x 1} vector of evolution error variances
#' @param XtX_bands list with 4 vectors consistint of the 4-bands of XtX = crossprod(X) (one-time cost)
#' @param Xty the \code{p x 1} matrix crossprod(X,y), which is a one-time cost (assuming no missing entries in y)
#' @param D the degree of differencing (zero, one, or two)
#' @return \code{p x 1} vector of simulated basis coefficients \code{beta}
#'
#' @note Missing entries (NAs) are not permitted in \code{y}. Imputation schemes are available.
#'
#' @import Matrix
#' @export
sampleBTF_bspline = function(y, X, obs_sigma2, evol_sigma_t2, XtX_bands, Xty = NULL, D = 1){

  # Some quick checks:
  if((D < 0) || (D != round(D)))  stop('D must be a positive integer')

  if(any(is.na(y))) stop('y cannot contain NAs')

  # Dimensions of X:
  T = nrow(X); p = ncol(X)

  # Linear term:
  if(is.null(Xty)) Xty = crossprod(X, y)
  linht = 1/obs_sigma2*Xty

  # Quadratic terms and solutions are computed differently, depending on D:
  if(D == 0){
    # Special case: no differencing
    QHt_Matrix = bandSparse(p, k = c(0,1,2,3),
                            diag = list(XtX_bands$XtX_0/obs_sigma2,
                                        XtX_bands$XtX_1/obs_sigma2,
                                        XtX_bands$XtX_2/obs_sigma2,
                                        XtX_bands$XtX_3/obs_sigma2),
                            symm = TRUE)
  } else {
    # Prior/evoluation quadratic term: can construct directly for D = 1 or D = 2
    if(D == 1){
      QHt_Matrix = bandSparse(p, k = c(0,1,2,3),
                              diag = list(XtX_bands$XtX_0/obs_sigma2 + 1/evol_sigma_t2 + c(1/evol_sigma_t2[-1], 0),
                                          XtX_bands$XtX_1/obs_sigma2 + -1/evol_sigma_t2[-1],
                                          XtX_bands$XtX_2/obs_sigma2,
                                          XtX_bands$XtX_3/obs_sigma2),
                              symm = TRUE)
    } else {
      if(D == 2){
        QHt_Matrix = bandSparse(p, k = c(0,1,2,3),
                                diag = list(XtX_bands$XtX_0/obs_sigma2 + 1/evol_sigma_t2 + c(0, 4/evol_sigma_t2[-(1:2)], 0) + c(1/evol_sigma_t2[-(1:2)], 0, 0),
                                            XtX_bands$XtX_1/obs_sigma2 + c(-2/evol_sigma_t2[3], -2*(1/evol_sigma_t2[-(1:2)] + c(1/evol_sigma_t2[-(1:3)],0))),
                                            XtX_bands$XtX_2/obs_sigma2 + 1/evol_sigma_t2[-(1:2)],
                                            XtX_bands$XtX_3/obs_sigma2),
                                symm = TRUE)
      } else stop('sampleBTF_bspline() requires D=0, D=1, or D=2')
    }
  }

  # Cholesky:
  chQht_Matrix = Matrix::chol(QHt_Matrix)

  # And sample the basis coefficients:
  beta = Matrix::solve(chQht_Matrix,Matrix::solve(Matrix::t(chQht_Matrix), linht) + rnorm(p))

  # Return the sampled basis coefficients:
  beta
}
#----------------------------------------------------------------------------
#' Sample the latent log-volatilities
#'
#' Compute one draw of the log-volatilities using a discrete mixture of Gaussians
#' approximation to the likelihood (see Omori, Chib, Shephard, and Nakajima, 2007)
#' where the log-vols are assumed to follow an AR(1) model with time-dependent
#' innovation variances. More generally, the code operates for \code{p} independent
#' AR(1) log-vol processes to produce an efficient joint sampler in \code{O(Tp)} time.
#'
#' @param h_y the \code{T x p} matrix of data, which follow independent SV models
#' @param h_prev the \code{T x p} matrix of the previous log-vols
#' @param h_mu the \code{p x 1} vector of log-vol unconditional means
#' @param h_phi the \code{p x 1} vector of log-vol AR(1) coefficients
#' @param h_sigma_eta_t the \code{T x p} matrix of log-vol innovation standard deviations
#' @param h_sigma_eta_0 the \code{p x 1} vector of initial log-vol innovation standard deviations
#'
#' @return \code{T x p} matrix of simulated log-vols
#'
#' @note For Bayesian trend filtering, \code{p = 1}. More generally, the sampler allows for
#' \code{p > 1} but assumes (contemporaneous) independence across the log-vols for \code{j = 1,...,p}.
#'
#' @import Matrix
#' @import BayesLogit
sampleLogVols = function(h_y, h_prev, h_mu, h_phi, h_sigma_eta_t, h_sigma_eta_0){

  # Compute dimensions:
  h_prev = as.matrix(h_prev) # Just to be sure (T x p)
  n = nrow(h_prev); p = ncol(h_prev)

  # Mixture params: mean, variance, and weights
  # Kim, Shephard, Chib (1998) 7-component mixture:
  #m_st  = c(-11.40039, -5.24321, -9.83726, 1.50746,  -0.65098, 0.52478,  -2.35859)
  #v_st2 = c(5.795960,  2.613690, 5.179500, 0.167350, 0.640090, 0.340230, 1.262610)
  #q     = c(0.007300,  0.105560, 0.000020, 0.043950, 0.340010, 0.245660, 0.257500)

  # Omori, Chib, Shephard, Nakajima (2007) 10-component mixture:
  m_st  = c(1.92677, 1.34744, 0.73504, 0.02266, -0.85173, -1.97278, -3.46788, -5.55246, -8.68384, -14.65000)
  v_st2 = c(0.11265, 0.17788, 0.26768, 0.40611,  0.62699,  0.98583,  1.57469,  2.54498,  4.16591,   7.33342)
  q     = c(0.00609, 0.04775, 0.13057, 0.20674,  0.22715,  0.18842,  0.12047,  0.05591,  0.01575,   0.00115)

  # Add an offset: common for all times, but distict for each j=1,...,p
  yoffset = tcrossprod(rep(1,n),
                       apply(as.matrix(h_y), 2,
                             function(x) any(x^2 < 10^-16)*max(10^-8, mad(x)/10^6)))

  # This is the response in our DLM, log(y^2)
  ystar = log(h_y^2 + yoffset)

  # Sample the mixture components
  #z = draw.indicators(res = ystar-h_prev, nmix = list(m = m_st, v = v_st2, p = q))
  z = sapply(ystar-h_prev, ncind, m_st, sqrt(v_st2), q)

  # Subset mean and variances to the sampled mixture components; (n x p) matrices
  m_st_all = matrix(m_st[z], nr=n); v_st2_all = matrix(v_st2[z], nr=n)

  # Joint AWOL sampler for j=1,...,p:

  # Constant (but j-specific) mean
  h_mu_all = tcrossprod(rep(1,n), h_mu)

  # Constant (but j-specific) AR(1) coef
  h_phi_all = tcrossprod(rep(1,n), h_phi)

  # Linear term:
  linht = matrix((ystar - m_st_all - h_mu_all)/v_st2_all)

  # Evolution precision matrix (n x p)
  evol_prec_mat = matrix(0, nr = n, nc = p);
  evol_prec_mat[1,] = 1/h_sigma_eta_0^2;
  evol_prec_mat[-1,] = 1/h_sigma_eta_t^2;

  # Lagged version, with zeros as appropriate (needed below)
  evol_prec_lag_mat = matrix(0, nr = n, nc = p);
  evol_prec_lag_mat[1:(n-1),] = evol_prec_mat[-1,]

  # Diagonal of quadratic term:
  Q_diag = matrix(1/v_st2_all +  evol_prec_mat + h_phi_all^2*evol_prec_lag_mat)

  # Off-diagonal of quadratic term:
  Q_off = matrix(-h_phi_all*evol_prec_lag_mat)[-(n*p)]

  # Quadratic term:
  QHt_Matrix = bandSparse(n*p, k = c(0,1), diag = list(Q_diag, Q_off), symm = TRUE)
  #QHt_Matrix = as.spam.dgCMatrix(as(bandSparse(n*p, k = c(0,1), diag = list(Q_diag, Q_off), symm = TRUE),"dgCMatrix"))

  # Cholesky:
  chQht_Matrix = Matrix::chol(QHt_Matrix)

  # Sample the log-vols:
  hsamp = h_mu_all + matrix(Matrix::solve(chQht_Matrix,Matrix::solve(Matrix::t(chQht_Matrix), linht) + rnorm(length(linht))), nr = n)
  #hsamp = h_mu_all +matrix(rmvnorm.canonical(n = 1, b = linht, Q = QHt_Matrix, Rstruct = cholDSP0))


  # Return the (uncentered) log-vols
  hsamp
}
#----------------------------------------------------------------------------
#' Sampler evolution error variance parameters
#'
#' Compute one draw of evolution error variance parameters under the various options:
#' \itemize{
#' \item dynamic horseshoe prior ('DHS');
#' \item horseshoe prior ('HS');
#' \item normal-inverse-gamma prior ('NIG').
#' }
#'
#' @param omega \code{T x p} matrix of evolution errors
#' @param evolParams list of parameters pertaining to each \code{evol_error} type to be updated
#' @param sigma_e the observation error standard deviation; for (optional) scaling purposes
#' @param evol_error the evolution error distribution; must be one of
#' 'DHS' (dynamic horseshoe prior), 'HS' (horseshoe prior), or 'NIG' (normal-inverse-gamma prior)
#' @return List of relevant components in \code{evolParams}: \code{sigma_wt}, the \code{T x p} matrix of evolution standard deviations,
#' and additional parameters associated with the DHS and HS priors.
#'
#' @note The list \code{evolParams} is specific to each \code{evol_error} type,
#' but in each case contains the evolution error standard deviations \code{sigma_wt}.
#'
#' @note To avoid scaling by the observation standard deviation \code{sigma_e},
#' simply use \code{sigma_e = 1} in the functional call.
#'
#' @import stochvol
#' @export
sampleEvolParams = function(omega, evolParams,  sigma_e = 1, evol_error = "DHS"){

  # Check:
  if(!((evol_error == "DHS") || (evol_error == "HS") || (evol_error == "BL") || (evol_error == "SV") || (evol_error == "NIG"))) stop('Error type must be one of DHS, HS, BL, SV, or NIG')

  # Make sure omega is (n x p) matrix
  omega = as.matrix(omega); n = nrow(omega); p = ncol(omega)

  if(evol_error == "DHS") return(sampleDSP(omega, evolParams, sigma_e))

  if(evol_error == "HS"){

    # For numerical reasons, keep from getting too small
    hsOffset = tcrossprod(rep(1,n), apply(omega, 2, function(x) any(x^2 < 10^-16)*max(10^-8, mad(x)/10^6)))
    hsInput2 = omega^2 + hsOffset

    # Local scale params:
    evolParams$tauLambdaj = matrix(rgamma(n = n*p, shape = 1, rate = evolParams$xiLambdaj + hsInput2/2), nr = n)
    evolParams$xiLambdaj = matrix(rgamma(n = n*p, shape = 1, rate = evolParams$tauLambdaj + tcrossprod(rep(1,n), evolParams$tauLambda)), nr = n)

    # Global scale params:
    evolParams$tauLambda = rgamma(n = p, shape = 0.5 + n/2, colSums(evolParams$xiLambdaj) + evolParams$xiLambda)
    #evolParams$xiLambda = rgamma(n = p, shape = 1, rate = evolParams$tauLambda + 1/sigma_e^2)
    evolParams$xiLambda = rgamma(n = p, shape = 1, rate = evolParams$tauLambda + 1)

    evolParams$sigma_wt = 1/sqrt(evolParams$tauLambdaj)

    return(evolParams)
  }
  if(evol_error == "BL"){

    # For numerical reasons, keep from getting too small
    hsOffset = tcrossprod(rep(1,n), apply(omega, 2, function(x) any(x^2 < 10^-16)*max(10^-8, mad(x)/10^6)))
    hsInput2 = omega^2 + hsOffset

    # 1/tau_j^2 is inverse-gaussian (NOTE: this is very slow!)
    evolParams$tau_j = matrix(sapply(matrix(hsInput2), function(x){1/sqrt(rig(n = 1,
                                            mean = sqrt(evolParams$lambda2*sigma_e^2/x), # already square the input
                                            scale = 1/evolParams$lambda2))}), nr = n)
    # Note: should be better priors for lambda2
    evolParams$lambda2 = rgamma(n = 1,
                                shape = 1 + n*p,
                                rate = 2 + sum(evolParams$tau_j^2)/2)

    # For Bayesian lasso, scale by sigma_e:
    evolParams$sigma_wt = sigma_e*evolParams$tau_j

    return(evolParams)
  }
  if(evol_error == "SV") return(sampleSVparams(omega = omega, svParams = evolParams))
  #if(evol_error == "SV") return(sampleSVparams0(omega = omega, svParams = evolParams))
  if(evol_error == "NIG") {
    evolParams = list(sigma_wt = tcrossprod(rep(1,n),
                                            apply(omega, 2,
                                                  function(x) 1/sqrt(rgamma(n = 1, shape = n/2 + 0.01, rate = sum(x^2)/2 + 0.01)))))
    return(evolParams)
  }
}
#----------------------------------------------------------------------------
#' Sample the dynamic shrinkage process parameters
#'
#' Compute one draw for each of the parameters in the dynamic shrinkage process
#' for the special case in which the shrinkage parameter \code{kappa ~ Beta(alpha, beta)}
#' with \code{alpha = beta}. The primary example is the dynamic horseshoe process with
#' \code{alpha = beta = 1/2}.
#'
#' @param omega \code{T x p} matrix of evolution errors
#' @param evolParams list of parameters to be updated (see Value below)
#' @param sigma_e the observation error standard deviation; for (optional) scaling purposes
#' @param prior_dhs_phi the parameters of the prior for the log-volatilty AR(1) coefficient \code{dhs_phi};
#' either \code{NULL} for uniform on [-1,1] or a 2-dimensional vector of (shape1, shape2) for a Beta prior
#' on \code{[(dhs_phi + 1)/2]}
#' @param alphaPlusBeta For the symmetric prior kappa ~ Beta(alpha, beta) with alpha=beta,
#' specify the sum [alpha + beta]
#' @return List of relevant components:
#' \itemize{
#' \item the \code{T x p} evolution error standard deviations \code{sigma_wt},
#' \item the \code{T x p} log-volatility \code{ht}, the \code{p x 1} log-vol unconditional mean(s) \code{dhs_mean},
#' \item the \code{p x 1} log-vol AR(1) coefficient(s) \code{dhs_phi},
#' \item the \code{T x p} log-vol innovation standard deviations \code{sigma_eta_t} from the Polya-Gamma priors,
#' \item the \code{p x 1} initial log-vol SD \code{sigma_eta_0},
#' \item and the mean of log-vol means \code{dhs_mean0} (relevant when \code{p > 1})
#' }
#'
#' @note The priors induced by \code{prior_dhs_phi} all imply a stationary (log-) volatility process.
#'
#' @import BayesLogit
#' @export
sampleDSP = function(omega, evolParams, sigma_e = 1, prior_dhs_phi = c(10,2), alphaPlusBeta = 1){

  # Store the DSP parameters locally:
  ht = evolParams$ht; dhs_mean = evolParams$dhs_mean; dhs_phi = evolParams$dhs_phi; sigma_eta_t = evolParams$sigma_eta_t; sigma_eta_0 = evolParams$sigma_eta_0; dhs_mean0 = evolParams$dhs_mean0

  # "Local" number of time points
  ht = as.matrix(ht)
  n = nrow(ht); p = ncol(ht)

  # Sample the log-volatilities using AWOL sampler
  ht = sampleLogVols(h_y = omega, h_prev = ht, h_mu = dhs_mean, h_phi=dhs_phi, h_sigma_eta_t = sigma_eta_t, h_sigma_eta_0 = sigma_eta_0)

  # Compute centered log-vols for the samplers below:
  ht_tilde = ht - tcrossprod(rep(1,n), dhs_mean)

  # Sample AR(1) parameters
    # Note: dhs_phi = 0 means non-dynamic HS, while dhs_phi = 1 means RW, in which case we don't sample either
  if(!all(dhs_phi == 0) && !all(dhs_phi == 1)) dhs_phi = sampleAR1(h_yc = ht_tilde, h_phi = dhs_phi, h_sigma_eta_t = sigma_eta_t, prior_dhs_phi = prior_dhs_phi)

  # Sample the evolution error SD of log-vol (i.e., Polya-Gamma mixing weights)
  eta_t = ht_tilde[-1,] - tcrossprod(rep(1,n-1), dhs_phi)*ht_tilde[-n, ]       # Residuals
  sigma_eta_t = matrix(1/sqrt(rpg(num = (n-1)*p, h = alphaPlusBeta, z = eta_t)), nc = p) # Sample
  sigma_eta_0 = 1/sqrt(rpg(num = p, h = 1, z = ht_tilde[1,]))                # Sample the inital

  # Sample the unconditional mean(s), unless dhs_phi = 1 (not defined)
  if(!all(dhs_phi == 1)){
    if(p > 1){
      # Assume a hierarchy of the global shrinkage params across j=1,...,p
      muSample = sampleLogVolMu(h = ht, h_mu = dhs_mean, h_phi = dhs_phi, h_sigma_eta_t = sigma_eta_t, h_sigma_eta_0 = sigma_eta_0, h_log_scale = dhs_mean0);
      dhs_mean = muSample$dhs_mean
      dhs_mean0 = sampleLogVolMu0(h_mu = dhs_mean, h_mu0 = dhs_mean0, dhs_mean_prec_j = muSample$dhs_mean_prec_j, h_log_scale = log(sigma_e^2))
    } else {
      # p = 1
      muSample = sampleLogVolMu(h = ht, h_mu = dhs_mean, h_phi = dhs_phi, h_sigma_eta_t = sigma_eta_t, h_sigma_eta_0 = sigma_eta_0, h_log_scale = log(sigma_e^2));
      dhs_mean = dhs_mean0 = muSample$dhs_mean # save dhs_mean0 = dhs_mean for coding convenience later
    }
  } else {dhs_mean = rep(0, p); dhs_mean0 = 0} # When RW for log-vols, fix unconditional mean for identifiability

  # Evolution error SD:
  sigma_wt = exp(ht/2)

  # Return the same list, but with the new values
  list(sigma_wt = sigma_wt, ht = ht, dhs_mean = dhs_mean, dhs_phi = dhs_phi, sigma_eta_t = sigma_eta_t, sigma_eta_0 = sigma_eta_0, dhs_mean0 = dhs_mean0)
}
#----------------------------------------------------------------------------
#' Sampler for the stochastic volatility parameters
#'
#' Compute one draw of the normal stochastic volatility parameters.
#' The model assumes an AR(1) for the log-volatility.
#'
#' @param omega \code{T x p} matrix of errors
#' @param svParams list of parameters to be updated
#' @return List of relevant components in \code{svParams}: \code{sigma_wt}, the \code{T x p} matrix of standard deviations,
#' and additional parameters associated with SV model.
#'
#' @import stochvol
#' @export
sampleSVparams = function(omega, svParams){

  # Make sure omega is (n x p) matrix
  omega = as.matrix(omega); n = nrow(omega); p = ncol(omega)

  for(j in 1:p){
    # First, check for numerical issues:
    svInput = omega[,j]; #if(all(svInput==0)) {svInput = 10^-8} else svInput = svInput + sd(svInput)/10^8
    #hsOffset = tcrossprod(rep(1,n), apply(omega, 2, function(x) any(x^2 < 10^-16)*max(10^-8, mad(x)/10^6)))

    # Sample the SV parameters:
    svsamp = stochvol::svsample_fast_cpp(svInput,
                                         startpara = list(
                                           mu = svParams$svParams[1,j],
                                           phi = svParams$svParams[2,j],
                                           sigma = svParams$svParams[3,j]),
                                         startlatent = svParams$ht[,j])# ,priorphi = c(10^4, 10^4));
    # Update the parameters:
    svParams$svParams[,j] = svsamp$para[1:3];
    svParams$ht[,j] = svsamp$latent
  }
  # Finally, up the evolution error SD:
  svParams$sigma_wt = exp(svParams$ht/2)

  # Check for numerically large values:
  svParams$sigma_wt[which(svParams$sigma_wt > 10^3, arr.ind = TRUE)] = 10^3

  return(svParams)
}
#----------------------------------------------------------------------------
#' Sampler for the stochastic volatility parameters using same functions as DHS prior
#'
#' Compute one draw of the normal stochastic volatility parameters.
#' The model assumes an AR(1) for the log-volatility.
#'
#' @param omega \code{T x p} matrix of errors
#' @param svParams list of parameters to be updated
#' @return List of relevant components in \code{svParams}: \code{sigma_wt}, the \code{T x p} matrix of standard deviations,
#' and additional parameters associated with SV model.
#'
sampleSVparams0 = function(omega, svParams){

  # Make sure omega is (n x p) matrix
  omega = as.matrix(omega); n = nrow(omega); p = ncol(omega)

  # Store the parameters locally:
  ht = as.matrix(svParams$ht); sv_mean = svParams$svParams[1,]; sv_phi = svParams$svParams[2,]; sv_sigma = svParams$svParams[3,]

  # Sample the log-volatilities using AWOL sampler
  ht = sampleLogVols(h_y = omega, h_prev = ht, h_mu = sv_mean, h_phi = sv_phi,
                     h_sigma_eta_t = matrix(rep(sv_sigma, each = n-1), nrow = n-1), h_sigma_eta_0 = sv_sigma) # New part
  # Compute centered log-vols for the samplers below:
  ht_tilde = ht - tcrossprod(rep(1,n), sv_mean)

  # Sample the AR(1) parameters:
  sv_phi = sampleAR1(h_yc = ht_tilde, h_phi = sv_phi,
                     h_sigma_eta_t = matrix(rep(sv_sigma, each = n-1), nrow = n-1),
                     prior_dhs_phi = c(10, 2))
  # Sample the evolution error SD of log-vol
  eta_t = ht_tilde[-1,] - tcrossprod(rep(1,n-1), sv_phi)*ht_tilde[-n, ]       # Residuals
  sv_sigma = apply(eta_t, 2, function(x)
    1/sqrt(rgamma(n = 1, shape = length(x)/2 + 0.01, rate = sum(x^2)/2 + 0.01)))

  # Sample the mean parameters:
  y_mu = (ht[-1,] - tcrossprod(rep(1,n-1), sv_phi)*ht[-n,])/matrix(rep(sv_sigma, each = n-1), nrow = n-1);
  x_mu = tcrossprod(rep(1,n-1), 1 - sv_phi)/matrix(rep(sv_sigma, each = n-1), nrow = n-1)
  postSD = 1/sqrt(colSums(x_mu^2) + 1/10^2)
  postMean = (colSums(x_mu*y_mu))*postSD^2
  sv_mean = rnorm(n = p, mean = postMean, sd = postSD)

  # Evolution error SD:
  sigma_wt = exp(ht/2)

  # Update:
  svParams$sigma_wt = sigma_wt; svParams$ht = ht;
  svParams$svParams[1,] = sv_mean;
  svParams$svParams[2,] = sv_phi;
  svParams$svParams[3,] = sv_sigma

  # And return:
  return(svParams)
}
#----------------------------------------------------------------------------
#' Sample the AR(1) coefficient(s)
#'
#' Compute one draw of the AR(1) coefficient in a model with Gaussian innovations
#' and time-dependent innovation variances. In particular, we use the sampler for the
#' log-volatility AR(1) process with the parameter-expanded Polya-Gamma sampler. The sampler also applies
#' to a multivariate case with independent components.
#'
#' @param h_yc the \code{T x p} matrix of centered log-volatilities
#' (i.e., the log-vols minus the unconditional means \code{dhs_mean})
#' @param h_phi the \code{p x 1} vector of previous AR(1) coefficient(s)
#' @param h_sigma_eta_t the \code{T x p} matrix of log-vol innovation standard deviations
#' @param prior_dhs_phi the parameters of the prior for the log-volatilty AR(1) coefficient \code{dhs_phi};
#' either \code{NULL} for uniform on [-1,1] or a 2-dimensional vector of (shape1, shape2) for a Beta prior
#' on \code{[(dhs_phi + 1)/2]}
#'
#' @return \code{p x 1} vector of sampled AR(1) coefficient(s)
#'
#' @note For the standard AR(1) case, \code{p = 1}. However, the function applies more
#' generally for sampling \code{p > 1} independent AR(1) processes (jointly).
#'
#' @import truncdist
#' @export
sampleAR1 = function(h_yc, h_phi, h_sigma_eta_t, prior_dhs_phi = NULL){

  # Compute dimensions:
  n = nrow(h_yc); p = ncol(h_yc)

  # Loop over the j=1:p
  for(j in 1:p){

    # Compute "regression" terms for dhs_phi_j:
    y_ar = h_yc[-1,j]/h_sigma_eta_t[,j] # Standardized "response"
    x_ar = h_yc[-n,j]/h_sigma_eta_t[,j] # Standardized "predictor"

    # Using Beta distribution:
    if(!is.null(prior_dhs_phi)){

      # Check to make sure the prior params make sense
      if(length(prior_dhs_phi) != 2) stop('prior_dhs_phi must be a numeric vector of length 2')

      dhs_phi01 = (h_phi[j] + 1)/2 # ~ Beta(prior_dhs_phi[1], prior_dhs_phi[2])

      # Slice sampler when using Beta prior:
      dhs_phi01 = uni.slice(dhs_phi01, g = function(x){
        -0.5*sum((y_ar - (2*x - 1)*x_ar)^2) +
          dbeta(x, shape1 = prior_dhs_phi[1], shape2 = prior_dhs_phi[2], log = TRUE)
      }, lower = 0, upper = 1)[1]#}, lower = 0.005, upper = 0.995)[1] #

      h_phi[j] = 2*dhs_phi01 - 1

    } else {
      # For h_phi ~ Unif(-1, 1), the posterior is truncated normal
      h_phi[j] = rtrunc(n = 1, spec = 'norm',
                        a = -1, b = 1,
                        mean = sum(y_ar*x_ar)/sum(x_ar^2),
                        sd = 1/sqrt(sum(x_ar^2)))
    }
  }
  h_phi
}
#----------------------------------------------------------------------------
#' Sample the AR(1) unconditional means
#'
#' Compute one draw of the unconditional means in an AR(1) model with Gaussian innovations
#' and time-dependent innovation variances. In particular, we use the sampler for the
#' log-volatility AR(1) process with the parameter-expanded Polya-Gamma sampler. The sampler also applies
#' to a multivariate case with independent components.
#'
#' @param h the \code{T x p} matrix of log-volatilities
#' @param h_mu the \code{p x 1} vector of previous means
#' @param h_phi the \code{p x 1} vector of AR(1) coefficient(s)
#' @param h_sigma_eta_t the \code{T x p} matrix of log-vol innovation standard deviations
#' @param h_sigma_eta_0 the standard deviations of initial log-vols
#' @param h_log_scale prior mean from scale mixture of Gaussian (Polya-Gamma) prior, e.g. log(sigma_e^2) or dhs_mean0
#'
#' @return a list containing
#' \itemize{
#' \item the sampled mean(s) \code{dhs_mean} and
#' \item the sampled precision(s) \code{dhs_mean_prec_j} from the Polya-Gamma parameter expansion
#'}
#'
#' @import BayesLogit
#' @export
sampleLogVolMu = function(h, h_mu, h_phi, h_sigma_eta_t, h_sigma_eta_0, h_log_scale = 0){

  # Compute "local" dimensions:
  n = nrow(h); p = ncol(h)

  # Sample the precision term(s)
  dhs_mean_prec_j = rpg(num = p, h = 1, z = h_mu - h_log_scale)

  # Now, form the "y" and "x" terms in the (auto)regression
  y_mu = (h[-1,] - tcrossprod(rep(1,n-1), h_phi)*h[-n,])/h_sigma_eta_t;
  x_mu = tcrossprod(rep(1,n-1), 1 - h_phi)/h_sigma_eta_t

  # Include the initial sd?
  #if(!is.null(h_sigma_eta_0)){y_mu = rbind(h[1,]/h_sigma_eta_0, y_mu); x_mu = rbind(1/h_sigma_eta_0, x_mu)}
  y_mu = rbind(h[1,]/h_sigma_eta_0, y_mu);
  x_mu = rbind(1/h_sigma_eta_0, x_mu)

  # Posterior SD and mean:
  postSD = 1/sqrt(colSums(x_mu^2) + dhs_mean_prec_j)
  postMean = (colSums(x_mu*y_mu) + h_log_scale*dhs_mean_prec_j)*postSD^2
  dhs_mean = rnorm(n = p, mean = postMean, sd = postSD)

  list(dhs_mean = dhs_mean, dhs_mean_prec_j = dhs_mean_prec_j)
}
#----------------------------------------------------------------------------
#' Sample the mean of AR(1) unconditional means
#'
#' Compute one draw of the mean of unconditional means in an AR(1) model with Gaussian innovations
#' and time-dependent innovation variances (for p > 1). More generally, the sampler
#' applies to the "mean" parameter (on the log-scale) for a Polya-Gamma parameter expanded
#' hierarchical model.
#'
#' @param h_mu the \code{p x 1} vector of means
#' @param h_mu0 the previous mean of unconditional means
#' @param dhs_mean_prec_j the \code{p x 1} vector of precisions (from the Polya-Gamma parameter expansion)
#' @param h_log_scale prior mean from scale mixture of Gaussian (Polya-Gamma) prior, e.g. log(sigma_e^2)
#'
#' @return The sampled mean parameter \code{dhs_mean0}
#'
#' @note This sampler is particularly for \code{p > 1} and the setting in which we want hierarchical
#' shrinkage effects, e.g. predictor- and time-dependent shrinkage, predictor-dependent shrinkage,
#' and global shrinkage, with a natural hierarchical ordering.
#'
#' @import BayesLogit
#' @export
sampleLogVolMu0 = function(h_mu, h_mu0, dhs_mean_prec_j, h_log_scale = 0){

  dhs_mean_prec_0 = rpg(num = 1, h = 1, z = h_mu0 - h_log_scale)

  # Sample the common mean parameter:
  postSD = 1/sqrt(sum(dhs_mean_prec_j) + dhs_mean_prec_j)
  postMean = (sum(dhs_mean_prec_j*h_mu) + dhs_mean_prec_j*h_log_scale)*postSD^2
  rnorm(n = 1, mean = postMean, sd = postSD)
}
#----------------------------------------------------------------------------
#' Sample the parameters for the initial state variance
#'
#' The initial state SDs are assumed to follow half-Cauchy priors, C+(0,A),
#' where the SDs may be common or distinct among the states.
#'
#' This function samples the parameters for a PX-Gibbs sampler.
#'
#' @param mu0 \code{p x 1} vector of initial values (undifferenced)
#' @param evolParams0 list of relevant components (see below)
#' @param commonSD logical; if TRUE, use common SDs (otherwise distict)
#' @param A prior scale parameter from the half-Cauchy prior, C+(0,A)
#' @return List of relevant components:
#' the \code{p x 1} evolution error SD \code{sigma_w0}
#' and the \code{p x 1} parameter-expanded RV's \code{px_sigma_w0}
#' @export
sampleEvol0 = function(mu0, evolParams0, commonSD = FALSE, A = 1){

  # Store length locally:
  p = length(mu0)

  # For numerical stability:
  mu02offset = any(mu0^2 < 10^-16)*max(10^-8, mad(mu0)/10^6)
  mu02 = mu0^2 + mu02offset

  if(commonSD){
    # (Common) standard deviations:
    evolParams0$sigma_w0 = rep(1/sqrt(rgamma(n = 1, shape = p/2 + 1/2, rate = sum(mu02)/2 + evolParams0$px_sigma_w0[1])), p)

    # (Common) paramater expansion:
    evolParams0$px_sigma_w0 = rep(rgamma(n = 1, shape = 1/2 + 1/2, rate = 1/evolParams0$sigma_w0[1]^2 + 1/A^2), p)

  } else {
    # (Distinct) standard deviations:
    evolParams0$sigma_w0 = 1/sqrt(rgamma(n = p, shape = 1/2 + 1/2, rate = mu02/2 + evolParams0$px_sigma_w0))

    # (Distict) paramater expansion:
    #evolParams0$px_sigma_w0 = rgamma(n = p, shape = 1/2 + 1/2, rate = 1/evolParams0$sigma_w0^2 + 1/A^2)
    evolParams0$px_sigma_w0 = rgamma(n = p, shape = 1/2 + 1/2, rate = 1/evolParams0$sigma_w0^2 + 1/evolParams0$sigma_00^2)

    # Global standard deviations:
    evolParams0$sigma_00 = 1/sqrt(rgamma(n = 1, shape = p/2 + 1/2, rate = sum(evolParams0$px_sigma_w0) + evolParams0$px_sigma_00))

    # (Global) parameter expansion:
    evolParams0$px_sigma_00 = rgamma(n = 1, shape = 1/2 + 1/2, rate = 1/evolParams0$sigma_00^2 + 1/A^2)
  }

  # And return the list:
  evolParams0
}
#----------------------------------------------------------------------------
#' Sample a Gaussian vector using the fast sampler of BHATTACHARYA et al.
#'
#' Sample from N(mu, Sigma) where Sigma = solve(crossprod(Phi) + solve(D))
#' and mu = Sigma*crossprod(Phi, alpha):
#'
#' @param Phi \code{n x p} matrix (of predictors)
#' @param Ddiag \code{p x 1} vector of diagonal components (of prior variance)
#' @param alpha \code{n x 1} vector (of data, scaled by variance)
#' @return Draw from N(mu, Sigma), which is \code{p x 1}, and is computed in \code{O(n^2*p)}
#' @note Assumes D is diagonal, but extensions are available
#' @export
sampleFastGaussian = function(Phi, Ddiag, alpha){

  # Dimensions:
  Phi = as.matrix(Phi); n = nrow(Phi); p = ncol(Phi)

  # Step 1:
  u = rnorm(n = p, mean = 0, sd = sqrt(Ddiag))
  delta = rnorm(n = n, mean = 0, sd = 1)

  # Step 2:
  v = Phi%*%u + delta

  # Step 3:
  w = solve(crossprod(sqrt(Ddiag)*t(Phi)) + diag(n), #Phi%*%diag(Ddiag)%*%t(Phi) + diag(n)
            alpha - v)

  # Step 4:
  theta =  u + Ddiag*crossprod(Phi, w)

  # Return theta:
  theta
}
