% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/debiased_ivx_est.R
\name{debias_ivx}
\alias{debias_ivx}
\title{Debiased IVX for predictive regression}
\usage{
debias_ivx(
  w,
  y,
  d_ind,
  intercept = FALSE,
  standardize = TRUE,
  c_z = 5,
  a = 0.9,
  standardize_iv = TRUE,
  iid = TRUE,
  post_inference = TRUE,
  lambda_choice = vector("list", length(d_ind) + 1),
  lambda_seq = vector("list", length(d_ind) + 1),
  train_method = "timeslice",
  nlambda = 100,
  lambda_min_ratio = 1e-04,
  k = 10,
  initial_window = ceiling(nrow(w) * 0.7),
  horizon = 1,
  fixed_window = TRUE,
  skip = 0,
  zhang_zhang = TRUE
)
}
\arguments{
\item{w}{Matrix of all regressors}

\item{y}{Vector of dependent variable}

\item{d_ind}{Index for inference targets}

\item{intercept}{Whether to include intercept in Lasso regression}

\item{c_z}{Parameter in constructing IV (Phillips and Lee, 2016) \deqn{z = \sum_{j=0}^{n-1} (1 - c_z / n^a)^j \Delta d_{t-j}}}

\item{a}{Parameter in constructing IV (Phillips and Lee, 2016)}

\item{standardize_iv}{Whether to standardize the IV}

\item{iid}{boolean indicating whether we want to adjust the long-run variance}

\item{post_inference}{boolean indicating whether to conduct post-lasso inference}

\item{lambda_choice}{Choice of lambda for Lasso regression; List of length length(d_ind) + 1: each element = NULL or = a number if user has a specific choice of tuning parameter}

\item{lambda_seq}{pre-specified sequence of tuning parameter for parameter tuning; Useful in calibration of tuning parameter based on the rate conditions in the asymptotic theory; List of length length(d_ind) + 1: Each element = NULL or a vector of tuning parameters}

\item{train_method}{The parameter tuning method
\itemize{
 \item{"timeslice"}{https://topepo.github.io/caret/data-splitting.html#time}:
     By combining initial window, horizon, fixed window and skip, we can control the sample splitting.
     Roll_block: Setting initial_window = horizon = floor(nrow(x) / k), fixed_window = False, and skip = floor(nrow(x) / k) - 1
     Period-by-period rolling: skip = 0. 
 \item "cv": Cross-validation based on block splits.
 \item "cv_random": Cross-validition based on random splits.
 \item "aic", "bic", "aicc", "hqc": based on information criterion.
 }}

\item{nlambda}{number of candidate lambdas}

\item{lambda_min_ratio}{# lambda_min_ratio * lambda_max = lambda_min (default: 0.0001): Determines the search range of lambda}

\item{k}{k-fold cv if "cv" is chosen (default: 10)}

\item{initial_window}{length of initial window for "timeslice" method}

\item{horizon}{length of horizon for "timeslice" method}

\item{fixed_window}{whether to use fixed window for "timeslice" method}

\item{skip}{length of skip for "timeslice" method}

\item{zhangzhang}{boolean indicating whether to conduct Zhang and Zhang (2014) debiased IVX}
}
\value{
A list contains
\item{theta_hat_las}{Estimate of delta from Lasso regression}
\item{theta_hat_ivx}{Estimate of delta from debiased IVX}
\item{sigma_hat}{Estimate of standard error of theta_hat_ivx}
\item{lambda_hat}{Chosen tuning parameter for Lasso regression}
\item{phi_hat}{Estimate of the frequency of 0s and L1 norm of std/nonstd coefficients in the second stage}
}
\description{
\deqn{y_t = w_{t-1} theta + u_t} where data is already aligned to incorporate the lagged regressor
}
