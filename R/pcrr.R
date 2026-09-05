#' Fit a Parametric Regression Model for the Cumulative Incidence Function
#'
#' Fits parametric regression models for the cumulative incidence function 
#' in competing risks data using the generalized odds rate transformation 
#' model by using various parametric baseline functions.
#' 

#' @details
#' The cumulative incidence function of cause k is modelled through the generalized
#' odds-rate link (Dabrowska and Doksum, 1988)
#' \deqn{F_k(t; \mathbf{Z}) = 1 - \{1 + \alpha_k \exp(\mathbf{Z}^{\top}\boldsymbol{\beta}_k) u_k(t)\}^{-1/\alpha_k},}
#' where \eqn{u_k(t)} is a parametric baseline function that is
#' nondecreasing in \eqn{t}. With
#' \code{distribution = "gompertz2"} the baseline (Jeong and Fine, 2007) is the two-parameter Gompertz
#' cumulative hazard
#' \deqn{u_k(t) = \dfrac{\tau_k}{\rho_k}(e^{\rho_k t} - 1),}
#' while \code{distribution = "gompertz3"} uses its three-parameter extension
#' (Haile et al., 2016)
#' \deqn{u_k(t) = \dfrac{\tau_k}{\rho_k \eta_k}
#'   (e^{(\eta_k e^{(\rho_k t)})} - e^{\eta_k}).}
#' The extra parameter \eqn{\eta_k} bends the baseline hazard so that unimodal
#' shapes can be represented; the two-parameter baseline only allows monotone
#' ones. As \eqn{\eta_k \to 0} the three-parameter baseline converges to the
#' two-parameter baseline, so \code{"gompertz2"} is nested within
#' \code{"gompertz3"} and the hypothesis \eqn{H_0: \eta_k = 0} can be tested;
#' this test is reported by \code{\link{summary.pcrr}}.
#'
#' The hazard is unimodal when \eqn{\rho_k > 0} and \eqn{-1 < \eta_k < 0}, or
#' when \eqn{\rho_k < 0} and \eqn{\eta_k < -1}, in which case it peaks at
#' \eqn{\frac{1}{\rho_k} \log(-\frac{1}{\eta_k})}. Fitting \code{"gompertz3"} to data with no
#' unimodal hazard usually makes the extra parameter unidentifiable, which shows
#' up as extreme estimates and a singular information matrix.
#' 
#' With \code{distribution = "logistic"}, the baseline is the three-parameter
#' modified logistic model (Cheng, 2009), with cumulative baseline function
#' \deqn{u_k(t) = -\log\left\{1 - \dfrac{p_k e^{b_k(t-c_k)} - p_k e^{-b_k c_k}}
#' {1 + e^{b_k(t-c_k)}}\right\}.}
#' The parameter \eqn{p_k} determines the long-term level of the cumulative
#' incidence, with \eqn{u_k(t) \to -\log(1-p_k)} as \eqn{t\to\infty},
#' while \eqn{b_k} controls the rate of increase and \eqn{c_k} determines the
#' location of the rise. The baseline hazard can be unimodal under
#' \eqn{\log\left(\frac{1+p_k e^{-b_k c_k}}{1-p_k}\right)\geq -2b_k c_k,}
#' in which case its peak occurs at \eqn{x_{\mathrm{mh}} = \frac{1}{2b_k}
#' \log\left(\frac{1+p_k e^{-b_k c_k}}{1-p_k}\right)+c_k.}
#' 
#' \emph{Parameter space.} Following Haile et al. (2016), the baseline requires
#' \eqn{\tau_k > 0}, while \eqn{\rho_k} and \eqn{\eta_k} can take any finite
#' real values. Accordingly, only \eqn{\tau_k} is constrained during estimation.
#' Following Cheng (2009), the three-parameter modified logistic baseline
#' requires \eqn{b_k > 0}, \eqn{0 < p_k < 1}, while \eqn{c_k} can take any finite
#' real value. Accordingly, \eqn{b_k} and \eqn{p_k} are constrained during
#' estimation, whereas \eqn{c_k} is unconstrained.
#'  The link parameter \eqn{\alpha_k} is likewise unconstrained, since the generalized
#' odds-rate transformation of Dabrowska and Doksum (1988) is defined for all
#' real \eqn{\alpha_k}. Note, however, that when \eqn{\alpha_k < 0}, the CIF is defined for all
#' \eqn{t \geq 0} such that
#' \eqn{1 + \alpha_k \exp(\mathbf{Z}^{\top}\boldsymbol{\beta}_k) u_k(t) \geq 0}.
#' Once this quantity becomes negative, the CIF is no longer defined for
#' subsequent values of \eqn{t}.
#' In \code{\link{predict.pcrr}}, the CIF is set to 1 from the first time
#' point at which this quantity becomes negative and remains at 1 thereafter,
#' whereas \code{\link{cure.pcrr}} returns \code{NA} for the corresponding
#' cure fraction.
#'
#' \emph{Notation.} Haile et al. (2016) write the three baseline parameters as
#' \eqn{(\alpha, \beta, \eta)}. They appear here as \code{rho}, \code{tau} and
#' \code{eta}, because \code{alpha} and \code{beta} already denote the link
#' parameter and the regression coefficients of Jeong and Fine (2007).
#'
#' Model parameters are estimated by maximum likelihood using \code{\link[stats]{nlminb}}.
#' When \code{variance = TRUE}, standard errors are obtained from the inverse of the observed information matrix. 
#' The returned object is of class \code{"pcrr"} and supports the S3 methods 
#' \code{print()}, \code{summary()}, \code{predict()}, and \code{cure()}.
#'
#' @seealso
#' \code{\link{print.pcrr}}, 
#' \code{\link{summary.pcrr}}, 
#' \code{\link{predict.pcrr}}, 
#' \code{\link{plot.predict.pcrr}}, 
#' \code{\link{cure.pcrr}}
#' 
#'
#' @param ftime numeric vector of failure or censoring times.
#' @param fstatus numeric vector indicating the event status, with a unique
#'  code for each failure type and a separate code for censored observations.
#' @param cov numeric matrix (\code{nobs x ncovs}) of fixed covariates.
#' @param distribution a character string specifying the baseline distribution.
#'  Available options are \code{"gompertz2"} (two-parameter Gompertz),
#'  \code{"gompertz3"} (three-parameter Gompertz), and
#'  \code{"logistic"} (three-parameter modified logistic).
#'  The abbreviations \code{"gom2"}, \code{"gom3"}, and \code{"logi"}
#'  are also accepted.
#' @param dist an optional alias for \code{distribution}. 
#'  If both \code{dist} and \code{distribution} are specified, 
#'  \code{dist} takes precedence.
#' @param failcode integer code in \code{fstatus} indicating the failure type of interest.
#' @param cencode integer code in \code{fstatus} indicating censored observations.
#' @param na.action function specifying how missing values in \code{ftime},
#'  \code{fstatus}, or \code{cov} are handled. The default is
#'  \code{na.omit}; \code{na.fail} is also supported.
#' @param gtol relative convergence tolerance (\code{rel.tol}) used by 
#' \code{nlminb()}. The default value is \code{1e-6}.
#' @param maxiter maximum number of optimization iterations. internally, \code{eval.max} is set to \code{3 * maxiter}.
#' @param init initial values of the model parameters. For 
#'  \code{distribution = "gompertz2"}, these are given as 
#'  \code{c(alpha1, rho1, tau1, beta11, beta12, ..., alpha2, rho2, tau2, 
#'  beta21, beta22, ...)}, and for \code{distribution = "gompertz3"} as 
#'  \code{c(alpha1, rho1, tau1, eta1, beta11, beta12, ..., alpha2, rho2, 
#'  tau2, eta2, beta21, beta22, ...)}. For 
#'  \code{distribution = "logistic"}, these are given as
#'  \code{c(alpha1, b1, c1, p1, beta11, beta12, ..., alpha2, b2, c2, p2, 
#'  beta21, beta22, ...)}. The expected length is \code{k * (3 + p)}, 
#'  \code{k * (4 + p)}, and \code{k * (4 + p)}, respectively. The first 
#'  block corresponds to the event of interest (\code{failcode}), followed 
#'  by the competing events in ascending order of their event codes. The 
#'  censoring code (\code{cencode}) is excluded.
#' @param variance logical value indicating whether the variance-covariance
#'  matrix is computed. The default is \code{TRUE}. If \code{FALSE}, only the
#'  maximum likelihood estimates and the score vector are computed.
#'
#' @return
#' An object of class \code{"pcrr"}, which is a list containing the following components:
#' \item{\code{coef}}{estimated model parameters. For \code{distribution = "gompertz2"},
#'  these are given as \code{c(alpha1, rho1, tau1, beta11, beta12, ...,
#'  alpha2, rho2, tau2, beta21, beta22, ...)}; for
#'  \code{distribution = "gompertz3"}, as
#'  \code{c(alpha1, rho1, tau1, eta1, beta11, beta12, ...,
#'  alpha2, rho2, tau2, eta2, beta21, beta22, ...)}; and for
#'  \code{distribution = "logistic"}, as
#'  \code{c(alpha1, b1, c1, p1, beta11, beta12, ..., alpha2, b2, c2, p2,
#'  beta21, beta22, ...)}.}
#' \item{\code{loglik}}{maximized log-likelihood value.}
#' \item{\code{init}}{initial parameter values used to start the optimization. The parameters are stored in the same order as \code{coef}.}
#' \item{\code{score}}{score vector evaluated at the maximum likelihood estimates.}
#' \item{\code{inf}}{observed information matrix (the negative Hessian matrix); \code{NULL} if \code{variance = FALSE}.}
#' \item{\code{invinf}}{inverse of the observed information matrix; \code{NULL} if \code{variance = FALSE}.}
#' \item{\code{converged}}{logical value indicating whether the optimization converged.}
#' \item{\code{iter}}{number of iterations performed by the optimization algorithm.}
#' \item{\code{message}}{message returned by the optimization routine.}
#' \item{\code{call}}{matched function call.}
#' \item{\code{n}}{total number of observations in the original data set.}
#' \item{\code{n_missing}}{number of observations removed due to missing values.}
#' \item{\code{k}}{number of event categories, excluding the censoring category.}
#' \item{\code{p}}{number of covariates.}
#' \item{\code{distribution}}{baseline distribution used for the model.}
#' \item{\code{cov_names}}{names of the covariates.}
#' \item{\code{mapping}}{mapping between the original event labels and the internal event codes.}
#' \item{\code{maxtime}}{maximum observed follow-up time.}
#' 
#' 
#' @importFrom survival coxph Surv
#' @importFrom stats D model.matrix na.fail na.omit nlminb pnorm qnorm printCoefmat setNames uniroot approx
#' @importFrom graphics lines legend abline axis par points
#' @importFrom rlang .data
#' @importFrom utils head
#' 
#' @examples
#' ## Example 1: Fit the model with default initial values
#' ## Two competing events (including the event of interest) and two covariates
#' 
#' set.seed(2026)
#'
#' # Covariates
#' z1 <- rbinom(500, 1, 0.5)
#' z2 <- rnorm(500)
#'
#' # Simulate event and censoring times
#' t1 <- rexp(500, rate = 0.15 * exp(-0.5 * z1 + 0.2 * z2))
#' t2 <- rexp(500, rate = 0.10 * exp(0.3 * z1))
#' tc <- runif(500, 0, 3)
#'
#' # Observed time and event indicator
#' time  <- pmin(t1, t2, tc)
#' event <- ifelse(time == t1, 1, ifelse(time == t2, 2, 0))
#'
#' fit1 <- pcrr(ftime = time, fstatus = event, cov = cbind(z1 = z1, z2 = z2))
#' summary(fit1)
#'
#' pred1 <- predict(fit1, cov = rbind(c(0, 0.13), c(1, -0.15), c(0, 0.40)))
#' print(pred1)
#' 
#' plot(pred1)
#'
#' cure(fit1, cov = rbind(c(0, 0.13), c(1, -0.15), c(0, 0.40)))
#'
#' @examples
#' ## Example 2: Fit the model with user-specified initial values
#' ## Three competing events (including the event of interest) and two covariates
#' 
#' set.seed(10)
#'
#' # Random competing risks data
#' time  <- rexp(500)
#' event <- sample(c(2, 4, 7, 8), 500, replace = TRUE)
#' z <- cbind(z1 = rbinom(500, 1, 0.3),
#'            z2 = rnorm(500))
#'
#' # Initial values:
#' # (alpha, rho, tau, beta1, beta2) for each event
#' user_init <- c(1, -1, 0.5,  1,  1,      # event 7 : event of interest (failcode)
#'                0.01, -1, 0.25, -1, 0,   # event 4 : competing event
#'                1, -1, 0.25, -1, -1)     # event 8 : competing event
#' 
#' fit2 <- pcrr(ftime = time, fstatus = event, cov = z, 
#'              failcode = 7, cencode = 2, init = user_init)
#' summary(fit2)
#'
#' pred2 <- predict(fit2, cov = rbind(c(0, 0.15), c(1, -0.30), c(0, 0.70)))
#' print(pred2)
#' 
#' plot(pred2)
#'
#' cure(fit2, cov = rbind(c(0, 0.15), c(1, -0.30), c(0, 0.70)))
#' 
#' @references
#' Jeong, J.-H. and Fine, J. P. (2007). Parametric regression on the
#' cumulative incidence function. \emph{Biostatistics}, 8(2), 184--196.
#' 
#' Dabrowska, D. M. and Doksum, K. A. (1988). Estimation and testing in a two-sample
#' generalized odds-rate model. \emph{Journal of the American Statistical Association},
#' 83(403), 744--749.
#'
#' Haile, S. R., Jeong, J.-H., Chen, X. and Cheng, Y. (2016). A 3-parameter Gompertz
#' distribution for survival data with competing risks, with an application to
#' breast cancer data. \emph{Journal of Applied Statistics}, 43(12), 2239--2253.
#' 
#' Cheng, Y. (2009). Modeling Cumulative Incidences of Dementia and Dementia-Free Death
#'  Using a Novel Three-Parameter Logistic Function.
#' \emph{The International Journal of Biostatistics}, 5(1), Article 29.
#'  
#' @export
pcrr <- function(ftime, fstatus, cov, distribution="gompertz2", dist=NULL, failcode=1, cencode=0,
                 na.action=na.omit, gtol=1e-6, maxiter=300, init, variance=TRUE) {

  if (!is.null(dist)) distribution <- dist
  
  if(distribution != "gompertz2" && distribution != "gom2" 
     && distribution != "gompertz3" && distribution != "gom3"
     && distribution != "logistic" && distribution != "logi"){
    warning("Unknown distribution. It runs by default. (gompertz2)")
    distribution <- "gompertz2"
  }
  if (distribution == "gom2") distribution <- "gompertz2"
  if (distribution == "gom3") distribution <- "gompertz3"
  if (distribution == "logi") distribution <- "logistic"
  
  if (!identical(na.action, na.omit) && !identical(na.action, na.fail)) {
    warning("Wrong na.action. It runs by default. (na.omit)")
    na.action <- na.omit
  }
  if (!identical(variance, TRUE) && !identical(variance, FALSE)) {
    warning("Wrong variance notaion. It runs by default. (variance = TRUE)")
    variance <- TRUE
  }
  
  
  call <- match.call()
  
  cov_name <- deparse(substitute(cov))
  cov_vars <- colnames(as.matrix(cov))
  
  
  # Handling missing values
  user_data <- data.frame(ftime = ftime, fstatus = fstatus, cov)
  
  na_others <- c("", ".", "-", "NA", "na", "N/A", "n/a", "NULL", "null")
  user_data[] <- lapply(user_data, function(x) {
    x[trimws(as.character(x)) %in% na_others] <- NA
    x
  })
  user_data <- na.action(user_data)
  
  
  N <- length(ftime)
  n <- nrow(user_data)
  N_mis <- 0
  
  if (N != n) {
    N_mis <- N - n
    cat(format(N_mis),'cases omitted due to missing values\n')
  }
  
  ftime   <- as.numeric(user_data$ftime)
  fstatus <- user_data$fstatus
  cov <- as.matrix(user_data[, c(-1, -2), drop = FALSE])
  
  
  # Dummy coding
  if (!(cencode  %in% fstatus)) stop("cencode is not in fstatus.")
  if (!(failcode %in% fstatus)) stop("failcode is not in fstatus.")
  
  other_events_codes <- sort(setdiff(unique(fstatus), c(cencode, failcode)))
  codes_ordered <- c(cencode, failcode, other_events_codes)
  mapping <- as.numeric(codes_ordered[-1])
  fevent  <- match(fstatus, codes_ordered) - 1
  
  dummy <- model.matrix(~ factor(fevent))
  
  delta <- as.matrix(dummy[ , -1, drop = FALSE])
  
  K <- ncol(delta)
  P <- ncol(cov)
  
  x <- ftime
  z <- matrix(as.numeric(cov), nrow = n, ncol = P)
  
  
  # initial value check
  init_ok <- FALSE
  if (!missing(init)) {
    n_param <- switch(distribution, gompertz2 = 3 + P, gompertz3 = 4 + P, logistic = 4 + P)
    if (length(init) == K * n_param) {
      init_ok <- TRUE
    } else {
      warning("Invalid length of 'init'. Default initial values will be used.")
    }
  }
  

  
  # assign covariate names
  if (!is.null(cov_vars) && all(!is.na(cov_vars) & nzchar(cov_vars))) {
    cov_names <- cov_vars               
  } else if (P == 1) {
    cov_names <- cov_name             
  } else {
    cov_names <- paste0(cov_name, seq_len(P))
  }
  colnames(z) <- cov_names
  


  
  # Kernel Operations
  if (distribution == "gompertz2"){
    if (init_ok) {
      theta_init <- init
    } else {
      theta_init <- .init_values_gom2(x, delta, z)
    }
    
    val_mle <- suppressWarnings(
      tryCatch(.log_lik_gom2(x = x, delta = delta, z = z, theta = theta_init),
        error = function(e) NaN)
    )
    
    if (!is.finite(val_mle) || abs(val_mle) >= 1e+100) {
      warning("The log-likelihood evaluated at the initial values returned NaN. ",
              "Optimization will proceed, but may converge to an incorrect or degenerate solution.")
    }
    
    theta_mle <- .estimate_mle_gom2(x, delta, z, theta_init, gtol, maxiter)
    score_hessian <- .score_hessian_gom2(x, delta, z, theta_mle$par, variance)
  } else if (distribution == "gompertz3"){
    if (init_ok) {
      theta_init <- init
    } else {
      theta_init <- .init_values_gom3(x, delta, z)
    }
    
    val_mle <- suppressWarnings(
      tryCatch(.log_lik_gom3(x = x, delta = delta, z = z, theta = theta_init),
               error = function(e) NaN)
    )
    
    if (!is.finite(val_mle) || abs(val_mle) >= 1e+100) {
      warning("The log-likelihood evaluated at the initial values returned NaN. ",
              "Optimization will proceed, but may converge to an incorrect or degenerate solution.")
    }
    
    theta_mle <- .estimate_mle_gom3(x, delta, z, theta_init, gtol, maxiter)
    score_hessian <- .score_hessian_gom3(x, delta, z, theta_mle$par, variance)
  } else if (distribution == "logistic"){
    if (init_ok) {
      theta_init <- init
    } else {
      theta_init <- .init_values_logi(x, delta, z)
    }
    
    val_mle <- suppressWarnings(
      tryCatch(.log_lik_logi(x = x, delta = delta, z = z, theta = theta_init),
               error = function(e) NaN)
    )
    
    if (!is.finite(val_mle) || abs(val_mle) >= 1e+100) {
      warning("The log-likelihood evaluated at the initial values returned NaN. ",
              "Optimization will proceed, but may converge to an incorrect or degenerate solution.")
    }
    
    theta_mle <- .estimate_mle_logi(x, delta, z, theta_init, gtol, maxiter)
    score_hessian <- .score_hessian_logi(x, delta, z, theta_mle$par, variance)
  } 
  
  
  # Display param name
  if (distribution == "gompertz2"){
    display_names <- character(K * (3 + P))
    for (k in 1:K){
      tag <- paste0("event", mapping[k], " : ")
      display_names[(k - 1) * (3 + P) + 1] <- paste0(tag, "alpha")
      display_names[(k - 1) * (3 + P) + 2] <- paste0(tag, "rho")
      display_names[(k - 1) * (3 + P) + 3] <- paste0(tag, "tau")
      for (p in 1:P){
        display_names[(k - 1) * (3 + P) + p + 3] <- paste0(tag, cov_names[p])
      }
    }
  } else if (distribution == "gompertz3"){
    display_names <- character(K * (4 + P))
    for (k in 1:K){
      tag <- paste0("event", mapping[k], " : ")
      display_names[(k - 1) * (4 + P) + 1] <- paste0(tag, "alpha")
      display_names[(k - 1) * (4 + P) + 2] <- paste0(tag, "rho")
      display_names[(k - 1) * (4 + P) + 3] <- paste0(tag, "tau")
      display_names[(k - 1) * (4 + P) + 4] <- paste0(tag, "eta")
      for (p in 1:P){
        display_names[(k - 1) * (4 + P) + p + 4] <- paste0(tag, cov_names[p])
      }
    }
  } else if (distribution == "logistic"){
    display_names <- character(K * (4 + P))
    for (k in 1:K){
      tag <- paste0("event", mapping[k], " : ")
      display_names[(k - 1) * (4 + P) + 1] <- paste0(tag, "alpha")
      display_names[(k - 1) * (4 + P) + 2] <- paste0(tag, "b")
      display_names[(k - 1) * (4 + P) + 3] <- paste0(tag, "c")
      display_names[(k - 1) * (4 + P) + 4] <- paste0(tag, "p")
      for (p in 1:P){
        display_names[(k - 1) * (4 + P) + p + 4] <- paste0(tag, cov_names[p])
      }
    }
  }
  
  par <- theta_mle$par
  sco <- score_hessian$score
  hess <- score_hessian$hessian
  
  names(par) <- display_names
  if (!is.null(sco)) names(sco) <- display_names
  if (variance && !is.null(hess)) dimnames(hess) <- list(display_names, display_names)
  
  inv_hess <- NULL 
  if (!is.null(hess)){
    inv_hess <- tryCatch(solve(-hess),
                         error = function(e) {
                           warning("Hessian matrix is singular or non-invertible.")
                           return(NULL)
                         })
  }

  
  # Define Class 'pcrr'
  cls <- list(coef      = par,
              loglik    = if(!is.finite(theta_mle$objective) || theta_mle$objective == 1e+100) NaN else -theta_mle$objective,
              init      = theta_init,
              score     = score_hessian$score,
              inf       = if (variance && !is.null(hess)) -hess else NULL,
              invinf    = if (variance) inv_hess else NULL,
              converged = theta_mle$convergence == 0,
              iter      = theta_mle$iterations,
              message   = theta_mle$message,
              call      = call,
              n         = N,
              n_missing = N_mis,
              k = K,
              p = P,
              distribution = distribution,
              cov_names = cov_names,
              mapping = mapping,
              maxtime = max(x)
  )
  class(cls) <- "pcrr"
  cls
}


#' Print a Fitted Parametric Competing Risks Regression Model
#'
#' Prints the estimated regression coefficients. If the model was fitted with
#' \code{variance = TRUE}, standard errors and two-sided p-values are also shown.
#'
#' @param x object of class \code{"pcrr"}.
#' @param ... further arguments passed to or from other methods.
#'
#' @return
#' The fitted model is printed. The function returns the input object invisibly.
#'
#' @seealso
#' \code{\link{pcrr}}, 
#' \code{\link{summary.pcrr}}
#'
#' @export
print.pcrr <- function(x, ...){
  P <- x$p
  K <- x$k

  # per-event block : gompertz2 -> 3 + P, gompertz3, logistic -> 4 + P
  block  <- if (x$distribution == "gompertz2") 3 + P else 4 + P
  bstart <- block - P + 1        # first beta slot inside a block
  
  cat("convergence:", x$converged, "\n")
  
  if (!x$converged) cat("[", x$message, "]\n")
  
  if (is.null(x$cov_names)) {
    cov_names <- paste0("beta", seq_len(P))
  } else {
    cov_names <- x$cov_names
  }
  
  if (is.null(x$invinf)) {
    cat("\nRegression coefficients:\n")
    
    for (k in seq_len(K)) {
      idx <- ((k - 1) * block + bstart):(k * block)
      
      out <- data.frame(
        Estimate = x$coef[idx],
        row.names = cov_names
      )
      
      cat("\nEvent", x$mapping[k], ":\n")
      print(out, ...)
    }
    
    return(invisible())
  }
  
  se <- sqrt(diag(x$invinf))
  pv <- 2 * (1 - pnorm(abs(x$coef) / se))
  
  cat("\nRegression coefficients:\n")
  
  for (k in seq_len(K)) {
    
    idx <- ((k - 1) * block + bstart):(k * block)
    
    out <- data.frame(
      Estimate   = x$coef[idx],
      Std.Error  = se[idx],
      P.value    = signif(pv[idx], 4),
      row.names  = cov_names
    )
    
    cat("\nEvent", x$mapping[k], ":\n")
    print(out, ...)
  }
  
  invisible()
}

#' Summarize a Fitted Parametric Competing Risks Regression Model
#'
#' Produces a summary of a fitted \code{"pcrr"} object, including
#' regression coefficients, confidence intervals, baseline parameter
#' estimates, and link function tests.
#'
#' @details
#' Four tables are produced.
#'
#' \emph{Regression coefficients.} For each covariate and cause the estimate
#' \eqn{\hat{\beta}}, its standard error, the Wald statistic
#' \eqn{z = \hat{\beta}/\mathrm{se}(\hat{\beta})} and the two-sided p-value are
#' reported. The Wald test assesses \eqn{H_0\!: \beta = 0}, that is, whether the
#' covariate affects the cumulative incidence of that cause. Standard errors come
#' from the inverse observed information matrix, so this table requires
#' \code{variance = TRUE}.
#'
#' \emph{Confidence intervals.} The same coefficients are shown on the
#' exponentiated scale, with intervals formed as
#' \eqn{\exp\{\hat{\beta} \pm z_{1-\alpha/2}\,\mathrm{se}(\hat{\beta})\}}.
#' Building the interval on the log scale and exponentiating keeps it positive and
#' asymmetric, as is appropriate for a ratio. The column \code{exp(-coef)} gives
#' the same effect with the reference group reversed.
#'
#' \emph{Parameters.} The distribution-specific baseline parameters are reported
#' without exponentiation, since they are not ratios. For
#' \code{distribution = "gompertz2"}, the parameters \eqn{\rho_k} and
#' \eqn{\tau_k} are reported. For \code{distribution = "gompertz3"}, the
#' additional shape parameter \eqn{\eta_k} is also reported. For
#' \code{distribution = "logistic"}, the baseline parameters \eqn{b_k},
#' \eqn{c_k}, and \eqn{p_k} are reported. The link parameter \eqn{\alpha_k} is
#' also reported for all three distributions.
#'
#' \emph{Link function tests.} Two Wald tests are reported for each
#' \eqn{\alpha_k}: \eqn{H_0\!: \alpha_k = 0}, under which the model reduces to a
#' proportional hazards model, and \eqn{H_0\!: \alpha_k = 1}, under which it
#' reduces to a proportional odds model. They are computed as
#' \eqn{(\hat{\alpha}_k - 0)/\mathrm{se}(\hat{\alpha}_k)} and
#' \eqn{(\hat{\alpha}_k - 1)/\mathrm{se}(\hat{\alpha}_k)}, and indicate whether a
#' simpler, more interpretable link is compatible with the data.
#'
#' \emph{Baseline shape test.} For \code{distribution = "gompertz3"} an
#' additional Wald test of \eqn{H_0: \eta_k = 0} is reported. Under this null the
#' three-parameter baseline collapses to the two-parameter one, so the test asks
#' whether the extra shape parameter is needed at all: a small p-value supports
#' \code{"gompertz3"}, while a large one indicates that \code{"gompertz2"} is
#' adequate. The test is omitted for \code{"gompertz2"} and \code{"logistic"} fits.
#' 
#' @param object object of class \code{"pcrr"}.
#' @param conf_int confidence level for confidence intervals. The default is \code{0.95}.
#' @param digits number of digits to print.
#' @param ... further arguments passed to or from other methods.
#'
#' @return
#' An object of class \code{"summary.pcrr"}, a list whose components include
#' \code{coef} (regression coefficients), \code{conf_int} (confidence intervals
#' on the exponentiated scale), \code{baseline} (the estimated link parameter
#' \code{alpha} and distribution-specific baseline parameters: \code{rho} and
#' \code{tau} for \code{"gompertz2"}, \code{rho}, \code{tau}, and \code{eta} for
#' \code{"gompertz3"}, and \code{b}, \code{c}, and \code{p} for
#' \code{"logistic"}), \code{link_ph} and \code{link_po} (link function tests),
#' and \code{shape} (the \eqn{\eta = 0} test for \code{"gompertz3"}, and
#' \code{NULL} for \code{"gompertz2"} and \code{"logistic"}).
#'
#' @seealso
#' \code{\link{pcrr}}, 
#' \code{\link{print.summary.pcrr}}
#'
#' @export
summary.pcrr <- function(object, conf_int = 0.95, digits = max(options()$digits - 4, 3), ...){
  est <- object$coef
  se  <- sqrt(diag(object$invinf))
  zst <- est / se
  pv  <- 2 * (1 - pnorm(abs(zst)))
  
  K <- object$k
  P <- object$p
  
  if (object$distribution == "gompertz2"){
    idx_alpha <- (seq_len(K) - 1) * (3 + P) + 1
    idx_rho   <- (seq_len(K) - 1) * (3 + P) + 2
    idx_tau   <- (seq_len(K) - 1) * (3 + P) + 3
    idx_nonbeta <- c(idx_alpha, idx_rho, idx_tau)
  } else if (object$distribution == "gompertz3"){
    idx_alpha <- (seq_len(K) - 1) * (4 + P) + 1
    idx_rho   <- (seq_len(K) - 1) * (4 + P) + 2
    idx_tau   <- (seq_len(K) - 1) * (4 + P) + 3
    idx_eta   <- (seq_len(K) - 1) * (4 + P) + 4
    idx_nonbeta <- c(idx_alpha, idx_rho, idx_tau, idx_eta)
  } else if (object$distribution == "logistic"){
    idx_alpha <- (seq_len(K) - 1) * (4 + P) + 1
    idx_b   <- (seq_len(K) - 1) * (4 + P) + 2
    idx_c   <- (seq_len(K) - 1) * (4 + P) + 3
    idx_p   <- (seq_len(K) - 1) * (4 + P) + 4
    idx_nonbeta <- c(idx_alpha, idx_b, idx_c, idx_p)
  }
  
  is_alpha <- logical(length(est))
  is_alpha[idx_alpha] <- TRUE
  is_beta  <- rep(TRUE, length(est))
  is_beta[idx_nonbeta] <- FALSE
  
  # Regression coefficient table
  coef_tab <- cbind(est[is_beta], exp(est[is_beta]), se[is_beta], zst[is_beta], pv[is_beta])
  dimnames(coef_tab) <- list(names(est)[is_beta], c("coef", "exp(coef)", "se(coef)", "z value", "Pr(>|z|)"))
  
  # Confidence interval
  a  <- (1 - conf_int) / 2
  a  <- c(a, 1 - a)
  zq <- qnorm(a)
  ci_tab <- cbind(exp(est[is_beta]), exp(-est[is_beta]),
                  exp(est[is_beta] + zq[1] * se[is_beta]),
                  exp(est[is_beta] + zq[2] * se[is_beta]))
  dimnames(ci_tab) <- list(names(est)[is_beta],
                           c("exp(coef)", "exp(-coef)",
                             paste0(format(100 * a, trim = TRUE, digits = 4), "%")))
  
  # Basis parameters (alpha, rho, tau)
  base_tab <- cbind(est[!is_beta], se[!is_beta], zst[!is_beta], pv[!is_beta])
  dimnames(base_tab) <- list(names(est)[!is_beta], c("est", "se", "z value", "Pr(>|z|)"))
  
  # Link function testing : H0: alpha=0 -> PH,  alpha=1 -> PO.
    z_ph <- (est[is_alpha] - 0) / se[is_alpha]
    z_po <- (est[is_alpha] - 1) / se[is_alpha]
    link_ph <- cbind(est = est[is_alpha], `z value` = z_ph, `Pr(>|z|)` = 2 * (1 - pnorm(abs(z_ph))))
    link_po <- cbind(est = est[is_alpha], `z value` = z_po, `Pr(>|z|)` = 2 * (1 - pnorm(abs(z_po))))


  # Baseline shape test : H0 : eta = 0, under which gompertz3 reduces to gompertz2
  shape_tab <- NULL
  if (object$distribution == "gompertz3"){
    z_eta <- est[idx_eta] / se[idx_eta]
    shape_tab <- cbind(est = est[idx_eta], `z value` = z_eta,
                       `Pr(>|z|)` = 2 * (1 - pnorm(abs(z_eta))))
  }
  
  out <- list(call = object$call,
              converged = object$converged,
              message = object$message,
              n = object$n,
              n_missing = object$n_missing,
              loglik = object$loglik,
              distribution = object$distribution,
              mapping = object$mapping,
              digits = digits,
              coef = coef_tab,
              conf_int = ci_tab,
              baseline = base_tab,
              shape = shape_tab,
              link_ph = link_ph,
              link_po = link_po
  )
  class(out) <- "summary.pcrr"
  out
}

#' Print a Summary of a Parametric Competing Risks Regression Model
#'
#' Prints a summary object produced by \code{\link{summary.pcrr}},
#' including regression coefficients, confidence intervals, baseline
#' parameter estimates, and link function tests.
#'
#' @param x object of class \code{"summary.pcrr"}.
#' @param digits number of significant digits to display.
#' @param ... further arguments passed to or from other methods.
#'
#' @return
#' The summary is printed. The function returns the input object invisibly.
#'
#' @seealso
#' \code{\link{pcrr}}, 
#' \code{\link{summary.pcrr}}
#'
#' @export
print.summary.pcrr <- function(x, digits = x$digits, ...){
  cat("Parametric Competing Risks Regression\n\n")
  if (!is.null(x$call)){ cat("Call:\n"); dput(x$call); cat("\n") }
  if (!x$converged){ cat("pcrr converged:", x$converged, "\n", "[", x$message, "]\n"); return(invisible()) }
  
  savedig <- options(digits = digits); on.exit(options(savedig))
  
  cat("Baseline distribution:", x$distribution, "\n")
  cat("Event codes: ", x$mapping[1], " = failure of interest",
      if (length(x$mapping) > 1)
        paste0(",    ", paste(x$mapping[-1], collapse = ", "), " = competing"),
      "\n\n", sep = "")
  cat("Regression coefficients:\n")
  printCoefmat(x$coef, digits = digits, signif.stars = TRUE,
               has.Pvalue = TRUE, P.values = TRUE, cs.ind = 1:3, tst.ind = 4)
  cat("\n")
  print(x$conf_int);                                    cat("\n")
  cat("Parameters:\n")
  printCoefmat(x$baseline, digits = digits, signif.stars = TRUE,
               has.Pvalue = TRUE, P.values = TRUE, cs.ind = 1:2, tst.ind = 3)
  cat("\n")
  
  # Link function tests (alpha = 0, alpha = 1)
  if (!is.null(x$link_ph)){
    cat("Link function tests\n\n")
    cat("H0: alpha = 0 (Proportional Hazards model)\n")
    printCoefmat(x$link_ph, digits = digits, signif.stars = TRUE,
                 has.Pvalue = TRUE, P.values = TRUE, cs.ind = 1, tst.ind = 2)
    cat("\n")
    cat("H0: alpha = 1 (Proportional Odds model)\n")
    printCoefmat(x$link_po, digits = digits, signif.stars = TRUE,
                 has.Pvalue = TRUE, P.values = TRUE, cs.ind = 1, tst.ind = 2)
    cat("\n")
  }
  
  # Baseline shape test (eta = 0)
  if (!is.null(x$shape)){
    cat("Baseline shape test\n\n")
    cat("H0: eta = 0 (two-parameter Gompertz baseline)\n")
    printCoefmat(x$shape, digits = digits, signif.stars = TRUE,
                 has.Pvalue = TRUE, P.values = TRUE, cs.ind = 1, tst.ind = 2)
    cat("\n")
  }
  
  cat("Num. cases =", x$n - x$n_missing)
  if (x$n_missing > 0)
    cat(" (", x$n_missing, " cases omitted due to missing values)", sep = "")
  cat("\n")
  cat("Log-likelihood =", format(x$loglik, nsmall = 4), "\n")
  invisible()
}



#' Predict Cumulative Incidence Functions
#'
#' Computes predicted cumulative incidence functions and predicted 
#' subdistribution hazards for covariate profiles from a fitted 
#' \code{"pcrr"} model.
#'
#' @details
#' The cumulative incidence function is obtained from the fitted
#' generalized odds-rate (GOR) transformation model:
#' \deqn{\widehat{F}_k(t;\mathbf{z}) = 1 - \left[1+ \widehat{\alpha}_k 
#' \exp(\mathbf{z}^{\top}\widehat{\boldsymbol{\beta}}_k) \widehat{u}_k(t)
#' \right]^{-1/\widehat{\alpha}_k}}
#' where \eqn{\widehat{u}_k(t)} is the fitted cumulative
#' baseline hazard corresponding to the selected distribution.
#'
#' Predictions are obtained separately for each row of \code{cov},
#' with one cumulative incidence curve returned for each covariate
#' profile.
#'
#' The fitted subdistribution hazard is evaluated as
#' \deqn{\lambda_k^{CI}(t;\mathbf{z}) = \frac{
#' \exp(\mathbf{z}^{\top}\widehat{\boldsymbol{\beta}}_k)\widehat{\lambda}_{0k}(t)}
#' {1+\widehat{\alpha}_k\exp(\mathbf{z}^{\top}\widehat{\boldsymbol{\beta}}_k)
#' \widehat{u}_k(t)}.}
#' where \eqn{\widehat{\lambda}_{0k}(t)}
#' denotes the fitted baseline hazard function and \eqn{\widehat{u}_k(t)}
#' is its cumulative form.
#'
#' For the three-parameter Gompertz distribution, a maximum
#' of the baseline hazard exists when
#' \eqn{\rho_k > 0, -1 < \eta_k < 0} or
#' \eqn{\rho_k < 0, \eta_k < -1}.
#' In this case, the baseline maximum hazard time is
#' \deqn{x_{mh}^{base} = \dfrac{1}{\rho_k}\log(-\dfrac{1}{\eta_k}).}
#'
#' The three-parameter Gompertz baseline hazard has a unimodal shape
#' when the corresponding parameter conditions are satisfied.
#' However, after applying the GOR transformation, the resulting
#' subdistribution hazard may exhibit either a unimodal or a U-shaped
#' pattern, depending on the transformation parameter and covariate profile.
#'
#' A sufficient condition for the subdistribution hazard to be unimodal is
#' \eqn{\eta_k<0,\quad \rho_k\neq0,\quad
#' \rho_k(1+\eta_k) -
#' \alpha_k \exp({\mathbf Z^\top\boldsymbol\beta_k})
#' \tau_k e^{\eta_k}>0,\quad
#' 1+\alpha_k \exp({\mathbf Z^\top\boldsymbol\beta_k})u_k(\infty)>0.}
#' However, when a boundary point \eqn{t_{\mathrm{boundary}}} exists,
#' the subdistribution hazard may still be unimodal even if these
#' conditions are not satisfied. Thus, the above conditions are sufficient
#' but not necessary for unimodality.
#'
#' A U-shaped subdistribution hazard occurs when
#' \eqn{\eta_k>0,\quad \rho_k>0,\quad \alpha_k>0,\quad
#' \rho_k(1+\eta_k) -
#' \alpha_k \exp({\mathbf Z^\top\boldsymbol\beta_k})
#' \tau_k e^{\eta_k}<0.}
#'
#' For the three-parameter Modified Logistic distribution, a maximum
#' of the baseline hazard exists when
#' \eqn{\log\left(\frac{1+p_ke^{-b_kc_k}}{1-p_k}\right)\leq-2b_kc_k.}
#' In this case, the baseline maximum hazard time is
#' \deqn{x_{mh}^{base} =
#' \dfrac{1}{2b_k}\log\left(\dfrac{1+p_ke^{-b_kc_k}}{1-p_k}\right)+c_k.}
#'
#' The three-parameter Modified Logistic baseline hazard has a unimodal
#' shape when the corresponding parameter conditions are satisfied.
#' Unlike the three-parameter Gompertz distribution, however, the
#' subdistribution hazard after applying the GOR transformation cannot
#' exhibit a U-shaped pattern. The condition for the subdistribution
#' hazard to be unimodal is
#' \eqn{e^{b_kc_k}+p_k\left\{1-\alpha_k\exp(\mathbf Z^\top\boldsymbol\beta_k)\right\}-1>0.}
#'
#' Because covariates modify the subdistribution hazard through the
#' GOR transformation, the hazard extrema for an individual
#' covariate profile generally differ from the corresponding baseline
#' values.
#' These profile-specific maximum or minimum hazard times are obtained
#' numerically by solving
#' \eqn{\frac{d}{dt}\lambda_k^{CI}(t;\mathbf{z})=0.}
#'
#' @param object object of class \code{"pcrr"}.
#' @param cov numeric matrix of covariate values. Rows correspond to
#'  covariate profiles and columns correspond to covariates included
#'  in the fitted model.
#' @param times optional vector of time points at which predictions are
#'  evaluated. If omitted, 200 equally spaced time points over the
#'  observed time range are used.
#' @param event integer code identifying the event type for prediction.
#'  If omitted, the first event type (the event of interest) in the fitted model is used.
#' @param ... further arguments passed to or from other methods.
#'
#' @return
#' An object of class \code{"predict.pcrr"}, containing:
#' \item{pred}{matrix containing time points in the first column and
#'  predicted cumulative incidence functions in subsequent columns.}
#' \item{baseline_hazard}{matrix containing the baseline hazard evaluated
#'  at the specified time points.}
#' \item{subdistribution_hazard}{matrix containing predicted
#'  subdistribution hazards for each covariate profile.}
#' \item{xmh_base}{baseline maximum hazard time for the
#'  \code{"gompertz3"} and \code{"logistic"} distribution when a maximum exists;
#'  otherwise \code{NULL}.}
#' \item{xmh_obs}{vector of covariate-specific maximum or minimum times of the
#'  subdistribution hazard obtained numerically for
#'  \code{"gompertz3"} and \code{"logistic"} models; otherwise \code{NULL}.}
#' \item{t_boundary}{vector of covariate-specific time points at which
#'  the GOR transformation becomes undefined, i.e.,
#'  \eqn{1+\widehat{\alpha}_k\exp(\mathbf{z}^{\top}\boldsymbol{\widehat{\beta}}_k)
#'  \widehat{u}_k(t)=0}. If this quantity remains positive for all
#'  \eqn{t\geq0}, \code{Inf} is returned. These values define the
#'  valid time domain for cumulative incidence function and
#'  subdistribution hazard calculations.}
#' \item{distribution}{the baseline distribution used in the fitted model.}
#' \item{eta}{estimated shape parameter of the three-parameter
#'  Gompertz distribution. This parameter influences the shape of the
#'  baseline hazard function, including unimodal and U-shaped
#'  patterns. Returned only for
#'  \code{distribution = "gompertz3"}; otherwise \code{NULL}.}
#' \item{labels}{labels corresponding to each covariate profile.}
#' \item{event}{event type for which predictions were obtained.}
#'
#' @seealso
#' \code{\link{pcrr}},
#' \code{\link{print.predict.pcrr}},
#' \code{\link{plot.predict.pcrr}}
#'
#' @export
predict.pcrr <- function(object, cov, times = NULL, event = NULL, ...){
  P <- object$p
  if (is.null(event)) event <- object$mapping[1]      # default : failcode
  k <- match(event, object$mapping)
  if (is.na(k)) stop("event must be one of: ", paste(object$mapping, collapse = ", "))
  if (is.null(times)) times <- seq(0, object$maxtime, length.out = 200)
  
  if (!is.matrix(cov)) {
    if (is.vector(cov)){
      if (length(cov) %% P == 0) cov <- matrix(cov, ncol = P, byrow = TRUE)
      else stop("cov must have ", P, " column(s).")
    } else cov <- as.matrix(cov)
  }
  if (!is.numeric(cov)) stop("cov must be numeric.")
  if (ncol(cov) != P) stop("cov must have ", P, " column(s).")
  

  
  if (object$distribution == "gompertz2" || object$distribution == "gompertz3"){
    u_k <- function(t, rho, tau, eta, tol){
      if (abs(rho) < tol && abs(eta) < tol){
        u <- tau * t
      } else if (abs(rho) < tol){
        u <- tau * exp(eta) * t
      } else if (abs(eta) < tol){
        u <- tau * expm1(rho * t) / rho
      } else {
        u <- tau * exp(eta) * expm1(eta * expm1(rho * t)) / (rho * eta)
      }
      u
    }
    du_k <- function(t, rho, tau, eta, tol){
      if (abs(rho) < tol && abs(eta) < tol){
        du <- tau + 0 * t
      } else if (abs(rho) < tol){
        du <- tau * exp(eta) + 0 * t
      } else if (abs(eta) < tol){
        du <- tau * exp(rho * t)
      } else {
        du <- tau * exp(rho * t) * exp(eta * exp(rho * t))
      }
      du
    }
  } else if (object$distribution == "logistic"){
    u_k <- function(t, b, c, p, tol){
      inner <- - p + p*(1+exp(-b*c))/(1+exp(b*(t-c)))
      u <- rep(NaN, length(t))
      valid <- 1 + inner > 0
      u[valid] <- -log1p(inner[valid])
      u
    }
    du_k <- function(t, b, c, p, tol){
      du <- b * p * (1 + exp(-b * c)) * exp(b * (t - c))/((1 + exp(b * (t - c)))*(1 + p * exp(-b * c)+(1 - p)*exp(b * (t - c))))
      du
    }
  }
  
  
  if (object$distribution == "gompertz2"){
    alpha <- object$coef[(k - 1) * (3 + P) + 1]
    rho   <- object$coef[(k - 1) * (3 + P) + 2]
    tau   <- object$coef[(k - 1) * (3 + P) + 3]
    eta   <- 0
    beta  <- object$coef[((k - 1) * (3 + P) + 4):(k * (3 + P))]
    
    u <- u_k(times, rho, tau, eta, 1e-12)
    base_haz <- du_k(times, rho, tau, eta, 1e-12)
  } else if (object$distribution == "gompertz3"){
    alpha <- object$coef[(k - 1) * (4 + P) + 1]
    rho   <- object$coef[(k - 1) * (4 + P) + 2]
    tau   <- object$coef[(k - 1) * (4 + P) + 3]
    eta   <- object$coef[(k - 1) * (4 + P) + 4]
    beta  <- object$coef[((k - 1) * (4 + P) + 5):(k * (4 + P))]
    
    u <- u_k(times, rho, tau, eta, 1e-12)
    base_haz <- du_k(times, rho, tau, eta, 1e-12)
  } else if (object$distribution == "logistic"){
    alpha <- object$coef[(k - 1) * (4 + P) + 1]
    b   <- object$coef[(k - 1) * (4 + P) + 2]
    c   <- object$coef[(k - 1) * (4 + P) + 3]
    p   <- object$coef[(k - 1) * (4 + P) + 4]
    beta  <- object$coef[((k - 1) * (4 + P) + 5):(k * (4 + P))]
    
    u <- u_k(times, b, c, p, 1e-12)
    base_haz <- du_k(times, b, c, p, 1e-12)
  }
  
  
  
  # CIF
  alpha_tol <- 1e-8
  cif <- matrix(0, nrow = length(times), ncol = nrow(cov))
  for (j in seq_len(nrow(cov))){
    ezb <- exp(sum(cov[j, ] * beta))
    if (abs(alpha) < alpha_tol){
      cif[, j] <- 1.0 - exp(-ezb * u)
      next
    }
    base <- 1.0 + alpha * ezb * u
    valid <- (!is.nan(u)) & (base > 1e-12)
    cif[valid, j] <- 1.0 - base[valid]^(-1.0 / alpha)
    cif[!valid, j] <- 1.0 # treat undefined region as CIF = 1
    if (any(!valid)){
      undef_t <- times[which(!valid)[1]]
      warning("CIF cannot be evaluated after time ", undef_t, " for obs ", j, ".")
    }
  }
  
 
  # maximum baseline hazard rate
  xmh_base <- NULL
  if (object$distribution == "gompertz3"){
    if ((rho > 0 && eta > -1 && eta < 0) || (rho < 0 && eta < -1)){
      calc_xmh <- (1 / rho) * log(-1.0 / eta)
      if (calc_xmh > 0) xmh_base <- unname(calc_xmh)
    }
  } else if (object$distribution == "logistic"){
    if (p > 0 && p < 1){
      calc_xmh <- (log1p(p*exp(-b * c)) - log1p(-p)) / (2 * b) + c
      if (calc_xmh > 0) xmh_base <- unname(calc_xmh)
    }
  }
  
  
  # subdistribution hazard
  sub_haz  <- matrix(0, nrow = length(times), ncol = nrow(cov))
  for (j in seq_len(nrow(cov))){
    ezb <- exp(sum(cov[j, ] * beta))
    base <- 1.0 + alpha * ezb * u
    valid <- (!is.nan(u)) & (base > 1e-12)
    sub_haz[valid, j] <- ezb * base_haz[valid] / base[valid]
    sub_haz[!valid, j] <- NA_real_ # treat undefined region as NA for prediction
    if (any(!valid)){
      undef_t <- times[which(!valid)[1]]
      warning("Subdistribution Hazard cannot be evaluated after time ", undef_t, " for obs ", j, ".")
    }
  }
  
  
  
  # boundary time 
  t_boundary <- rep(Inf, nrow(cov))
  
  if (object$distribution == "gompertz2" || object$distribution == "gompertz3") {
    # the time when 1 + alpha * ezb * u(t) = 0
    for (j in seq_len(nrow(cov))){
      ezb <- exp(sum(cov[j, ] * beta))
      if (alpha < 0){
        if ((rho == 0 || (rho > 0 && eta >= 0)) # u(Inf) -> Inf
            || (rho < 0 && 1 + alpha * ezb * u_k(Inf, rho, tau, eta, 1e-12) < 0)
            || (rho > 0 && eta < 0 && 1 + alpha * ezb * u_k(Inf, rho, tau, eta, 1e-12) < 0)){
          inner <- function(t, alpha, rho, tau, eta, ezb){
            u <- u_k(t, rho, tau, eta, 1e-12)
            1 + alpha * ezb * u
          }
          time_inner <- max(times)
          upper <- inner(time_inner, alpha, rho, tau, eta, ezb)
          iter <- 0
          while (upper > 0){
            time_inner <- time_inner * 10
            upper <- inner(time_inner, alpha, rho, tau, eta, ezb)
            iter <- iter + 1
            if(iter > 20) break
          }
          if (iter < 20){
            res <- tryCatch(uniroot(inner, interval = c(0, time_inner), tol = 1e-12, alpha = alpha, rho = rho, tau = tau, eta = eta, ezb = ezb), 
                            error = function(e) NULL)
            if (!is.null(res)) t_boundary[j] <- res$root
          } 
        }
      }
    }
  } else if (object$distribution == "logistic"){
    for (j in seq_len(nrow(cov))){
      ezb <- exp(sum(cov[j, ] * beta))
      if (alpha > 0 && p > 1) {
        t_boundary[j] <- (log1p(p*exp(-b * c)) - log(p - 1)) / b + c
      } else if (alpha < 0 && (p > 1 - exp(1/(alpha * ezb)))){
        t_boundary[j] <- (log1p(p*exp(-b * c) - exp(1/(alpha * ezb))) - log(exp(1/(alpha * ezb)) + p - 1)) / b + c
      } 
    }
  }
  
  
  
  # subdistribution hazard peak
  xmh_obs <- NULL
  if (object$distribution == "gompertz3"){
    xmh_obs <- rep(NA_real_, nrow(cov))
    
    lambda_prime <- function(t, alpha, rho, tau, eta, ezb){
      u <- u_k(t, rho, tau, eta, 1e-12)
      du <- du_k(t, rho, tau, eta, 1e-12)
      rho * (1 + eta * exp(rho * t)) * (1 + alpha * ezb * u) - alpha * ezb * du
    }
    
    for (j in seq_len(nrow(cov))){
      ezb <- exp(sum(cov[j, ] * beta))

      lambda_prime_t0 <- lambda_prime(0, alpha, rho, tau, eta, ezb)
      if (abs(lambda_prime_t0) < 1e-12) {xmh_obs[j] <- 0; next}
      if (is.finite(t_boundary[j])) {
        lambda_prime_t_max <- lambda_prime(t_boundary[j], alpha, rho, tau, eta, ezb)
        if (lambda_prime_t0 * lambda_prime_t_max > 0) next
      } else {
        if (!((eta > 0 && lambda_prime_t0 < 0) || (eta < 0 && lambda_prime_t0 > 0))) next
      }
      
      if (is.finite(t_boundary[j])) {
        t_max <- max(0, t_boundary[j]*(1-1e-8))
      } else {
        t_max <- max(times)
        lambda_prime_t_max <- lambda_prime(t_max, alpha, rho, tau, eta, ezb)
        iter <- 0
        while (lambda_prime_t0 * lambda_prime_t_max > 0) {
          t_max <- t_max * 10
          lambda_prime_t_max <- lambda_prime(t_max, alpha, rho, tau, eta, ezb)
          iter <- iter + 1
          if(iter > 20){
            t_max <- NA_real_
            break
          }
        }
      }
      if (!is.finite(t_max)) next
      res <- tryCatch(uniroot(lambda_prime, c(0, t_max), tol = 1e-12, alpha = alpha, rho = rho, tau = tau, eta = eta, ezb = ezb),
                      error = function(e) {warning("Failed to find root: ", conditionMessage(e)); NULL})
      if (is.null(res)) next
      xmh_obs[j] <- res$root
    }
    
    names(xmh_obs) <- paste0("obs ", seq_len(nrow(cov)))
    
  } else if (object$distribution == "logistic"){
    xmh_obs <- rep(NA_real_, nrow(cov))
    
    lambda_prime <- function(t, alpha, b, c, p, ezb){
      u <- u_k(t, b, c, p, 1e-12)
      du <- du_k(t, b, c, p, 1e-12)
      b*(1+p*exp(-b*c)-(1-p)*exp(2*b*(t-c)))/((1+exp(b*(t-c)))*(1+p*exp(-b*c)+(1-p)*exp(b*(t-c)))) * (1 + alpha * ezb * u) - alpha * ezb * du
    }
    
    for (j in seq_len(nrow(cov))){
      ezb <- exp(sum(cov[j, ] * beta))

      # condition of unimodal
      if (alpha * ezb * p < exp(b*c) + p - 1){
        lambda_prime_t0 <- lambda_prime(0, alpha, b, c, p, ezb)
        if (is.finite(t_boundary[j])) {
          t_max <- max(0, t_boundary[j]*(1-1e-8))
        } else {
          t_max <- max(times)
          lambda_prime_t_max <- lambda_prime(t_max, alpha, b, c, p, ezb)
          iter <- 0
          while (lambda_prime_t0 * lambda_prime_t_max > 0) {
            t_max <- t_max * 10
            lambda_prime_t_max <- lambda_prime(t_max, alpha, b, c, p, ezb)
            iter <- iter + 1
            if(iter > 20){
              t_max <- NA_real_
              break
            }
          }
        }
        if (!is.finite(t_max)) next
        res <- tryCatch(uniroot(lambda_prime, c(0, t_max), tol = 1e-12, alpha = alpha, b = b, c = c, p = p, ezb = ezb),
                        error = function(e) {warning("Failed to find root: ", conditionMessage(e)); NULL})
        if (is.null(res)) next
        xmh_obs[j] <- res$root
      } else if (abs(exp(b*c) + p *(1 - alpha * ezb) -1) < 1e-12) {
        xmh_obs[j] <- 0
      }
    }
    names(xmh_obs) <- paste0("obs ", seq_len(nrow(cov)))
  }
  
  colnames(sub_haz) <- paste0("obs ", seq_len(nrow(cov)))
  names(t_boundary) <- paste0("obs ", seq_len(nrow(cov)))
  
  
  
  # curve labels for the legend, e.g. "obs 1", "obs 2", ...
  labs <- paste0("obs ", seq_len(nrow(cov)))
  pred <- cbind(times, cif)
  colnames(pred) <- c("time", labs)
  base_haz <- cbind(time = times, `baseline hazard` = base_haz)
  sub_haz <- cbind(time = times, sub_haz)

  out <- list(pred = pred,
              baseline_hazard = base_haz,
              subdistribution_hazard = sub_haz,
              xmh_base = xmh_base,
              xmh_obs = xmh_obs,
              t_boundary = t_boundary,
              distribution = object$distribution,
              eta = if(object$distribution=="gompertz3") eta else NULL,
              labels = labs,
              event = event
              )
  
  class(out) <- "predict.pcrr"
  
  out
}


#' Print Predicted Cumulative Incidence Functions
#'
#' Prints a summary of predicted cumulative incidence functions produced by
#' \code{\link{predict.pcrr}}, including event type, maximum hazard rate time
#' (for \code{"gompertz3"} and \code{"logistic"}), and predicted values.
#'
#' @param x object of class \code{"predict.pcrr"}.
#' @param digits number of decimal places to format the maximum hazard rate time.
#' @param ... further arguments passed to \code{\link[base]{print}}.
#'
#' @return
#' The function prints the prediction information and returns the input object invisibly.
#'
#' @seealso
#' \code{\link{pcrr}}, 
#' \code{\link{predict.pcrr}}
#'
#' @export
print.predict.pcrr <- function(x, digits = 4, ...){
  cat("Predicted Cumulative Incidence (Event: ", x$event, ")\n", sep = "")
  
  if (x$distribution == "gompertz3"){
    if (is.null(x$xmh_base)) {
      cat("\nMaximum hazard rate time (baseline): None\n\n")
    } else {
      cat("\nMaximum hazard rate time (baseline)\n")
      cat("  baseline :", format(round(x$xmh_base, digits)), "\n\n")
    }
    if (!is.null(x$xmh_obs)){
      if (x$eta < 0){
        cat("\nStationary point of subdistribution hazard\n")
        cat("  Type: Maximum (unimodal)\n")
      } else if (x$eta > 0){
        cat("\nStationary point of subdistribution hazard\n")
        cat("  Type: Minimum (U-shape)\n")
      } else {
        cat("\nStationary point of subdistribution hazard\n")
        cat("  Type: Gompertz2 case (no stationary point)\n")
      }
      
      for (j in seq_along(x$xmh_obs)){
        cat(sprintf("  %-8s : %s\n",
                    names(x$xmh_obs)[j],
                    ifelse(is.finite(x$xmh_obs[j]),
                           format(round(x$xmh_obs[j], digits)),
                           "None")))
      }
    }
  } else if (x$distribution == "logistic"){
    if (is.null(x$xmh_base)) {
      cat("\nMaximum hazard rate time (baseline): None\n\n")
    } else {
      cat("\nMaximum hazard rate time (baseline)\n")
      cat("  baseline :", format(round(x$xmh_base, digits)), "\n\n")
    }
    if (!is.null(x$xmh_obs)){
        cat("\nStationary point of subdistribution hazard\n")
        cat("  Type: Maximum (unimodal)\n")
      
      
      for (j in seq_along(x$xmh_obs)){
        cat(sprintf("  %-8s : %s\n",
                    names(x$xmh_obs)[j],
                    ifelse(is.finite(x$xmh_obs[j]),
                           format(round(x$xmh_obs[j], digits)),
                           "None")))
      }
    }
  }
  cat("\n\n")
  print(x$pred, ...)
  invisible(x)
}

#' Plot Predicted Cumulative Incidence Functions
#'
#' Plots cumulative incidence functions produced by \code{\link{predict.pcrr}}.
#'
#' @param x object of class \code{"predict.pcrr"}.
#' @param lty line types for curves.
#' @param color line colors for curves.
#' @param ylim limits for the y-axis.
#' @param xmin minimum value of the x-axis.
#' @param xmax maximum value of the x-axis.
#' @param xlab label for the x-axis.
#' @param ylab label for the y-axis.
#' @param legend logical value indicating whether a legend is drawn. 
#'  If \code{TRUE} (default), a legend identifying each curve is displayed.
#' @param legend.pos position of the legend, passed to \code{\link[graphics]{legend}}. 
#'  The default is \code{"topleft"}.
#' @param legend.title optional title for the legend.
#' @param lwd line width for the curves.
#' @param main main title. If \code{NULL}, the event being plotted is used.
#' @param hazard logical value. If \code{TRUE}, the subdistribution hazard is
#'  drawn in a left-hand panel next to the cumulative incidence.
#'  Default is \code{TRUE}.
#' @param mark.xmh logical value. If \code{TRUE} (default) and profile-specific
#'  hazard turning points (\code{xmh_obs}) exist, peak or minimum hazard
#'  points and vertical reference lines are added to the plots for
#'  \code{"gompertz3"}, while peak hazard points and vertical reference lines
#'  are added for \code{"logistic"}.
#' @param max.marks maximum number of curves for which profile-specific peak hazard
#'  vertical reference lines are drawn. If there are more curves than this threshold,
#'  vertical lines are omitted to avoid overlap, keeping only the peak points on each curve.
#'  Default is 4.
#' @param xmh.lty line type for the vertical reference lines at peak hazard times. Default is 3 (dotted).
#' @param xmh.col color of the baseline \eqn{x_{mh}} reference line. Default is
#'  \code{"black"}, so that the covariate-free baseline peak stands out from the
#'  profile-specific lines, which follow the color of their own curve.
#' @param ... additional graphical parameters passed to \code{\link[graphics]{plot}}.
#'
#' @return
#' Produces a plot of predicted cumulative incidence functions and returns the input object invisibly.
#'
#' @seealso
#' \code{\link{pcrr}}, 
#' \code{\link{predict.pcrr}}
#'
#' @export
plot.predict.pcrr <- function(x, lty = seq_along(x$labels),
                              color = seq_along(x$labels) + 1,
                              ylim = NULL, xmin = 0,
                              xmax = max(x$pred[, 1], na.rm = TRUE),
                              xlab = "Time", ylab = "Cumulative Incidence",
                              legend = TRUE, legend.pos = "topleft",
                              legend.title = NULL, lwd = 2, main = NULL,
                              hazard = TRUE, mark.xmh = TRUE, max.marks = 4,
                              xmh.lty = 3, xmh.col = "black", ...){
  ncurve <- length(x$labels)
  
  lty   <- rep(lty,   length.out = ncurve)
  color <- rep(color, length.out = ncurve)
  
  labs <- x$labels
  if (is.null(labs) || length(labs) != ncurve){
    labs <- colnames(x$pred)[-1]
    if (is.null(labs)) labs <- paste("curve", seq_len(ncurve))
  }
  
  ev <- x$event
  if (is.null(main) && !is.null(ev)) main <- paste("Event", ev)
  
  show_haz <- isTRUE(hazard) && !is.null(x$subdistribution_hazard)
  

  marks <- NULL
  if (isTRUE(mark.xmh) && !is.null(x$xmh_obs)) marks <- x$xmh_obs
  show_lines <- !is.null(marks) && ncurve <= max.marks
  
  base_x <- NULL
  if (isTRUE(mark.xmh) && !is.null(x$xmh_base))
    if (x$xmh_base >= xmin && x$xmh_base <= xmax) base_x <- x$xmh_base
  
  if (show_haz){
    oldpar <- graphics::par(mfrow = c(1, 2))
    on.exit(graphics::par(oldpar), add = TRUE)
    
    hz <- x$subdistribution_hazard[, -1, drop = FALSE]
    plot(c(xmin, xmax), c(0, max(hz, na.rm = TRUE) * 1.05), type = "n",
         xlab = xlab, ylab = "Subdistribution Hazard", main = main, ...)
    if (!is.null(base_x)) abline(v = base_x, lty = xmh.lty, col = xmh.col, lwd = 1.6)
    if (show_lines)
      for (j in seq_len(ncurve))
        if (is.finite(marks[j]) && marks[j] >= xmin && marks[j] <= xmax)
          abline(v = marks[j], lty = xmh.lty, col = color[j])
    for (j in seq_len(ncurve))
      lines(x$subdistribution_hazard[, 1], hz[, j], lty = lty[j], col = color[j], lwd = lwd)
    if (!is.null(marks))
      for (j in seq_len(ncurve))
        if (is.finite(marks[j]) && marks[j] >= xmin && marks[j] <= xmax)
          points(marks[j], approx(x$subdistribution_hazard[, 1], hz[, j], xout = marks[j])$y,
                 pch = 19, col = color[j], cex = 0.9)
    if (isTRUE(legend)){
      lg <- labs
      if (!is.null(marks) && any(is.finite(marks)))
        lg <- ifelse(is.finite(marks), sprintf("%s   x_mh = %.2f", labs, marks),labs)
      lc <- color; ll <- lty; lw <- rep(lwd, ncurve)
      if (!is.null(base_x)){
        lg <- c(lg, sprintf("baseline x_mh = %.2f", base_x))
        lc <- c(lc, xmh.col); ll <- c(ll, xmh.lty); lw <- c(lw, 1.6)
      }
      legend(legend.pos, legend = lg, lty = ll, col = lc, lwd = lw,
             title = legend.title, bty = "n")
    }
  }
  
  if (is.null(ylim)) ylim <- c(0, max(x$pred[, -1], na.rm = TRUE))
  plot(c(xmin, xmax), ylim, type = "n", xlab = xlab, ylab = ylab, main = main, ...)
  
  if (!is.null(base_x)) abline(v = base_x, lty = xmh.lty, col = xmh.col, lwd = 1.6)
  if (show_lines)
    for (j in seq_len(ncurve))
      if (is.finite(marks[j]) && marks[j] >= xmin && marks[j] <= xmax)
        abline(v = marks[j], lty = xmh.lty, col = color[j])
  
  for (j in seq_len(ncurve))
    lines(x$pred[, 1], x$pred[, j + 1], lty = lty[j], col = color[j], lwd = lwd)
  
  if (!is.null(marks))
    for (j in seq_len(ncurve))
      if (is.finite(marks[j]) && marks[j] >= xmin && marks[j] <= xmax)
        points(marks[j], approx(x$pred[, 1], x$pred[, j + 1], xout = marks[j])$y,
               pch = 19, col = color[j], cex = 0.9)
  
  if (isTRUE(legend)){
    leg <- labs
    if (show_lines && any(is.finite(marks)))
      leg <- ifelse(is.finite(marks), sprintf("%s   x_mh = %.2f", labs, marks), labs)
    lc <- color; ll <- lty; lw <- rep(lwd, ncurve)
    if (!is.null(base_x)){
      leg <- c(leg, sprintf("baseline x_mh = %.2f", base_x))
      lc <- c(lc, xmh.col); ll <- c(ll, xmh.lty); lw <- c(lw, 1.6)
    }
    legend(legend.pos, legend = leg, lty = ll, col = lc, lwd = lw,
           title = legend.title, bty = "n")
  }
  
  invisible(x)
}

#' Estimate Cure Fractions
#'
#' Generic function for estimating cure fractions from fitted
#' competing risks regression models.
#'
#' @details
#' The cure fraction is the limiting probability of never experiencing the
#' event, \eqn{1 - F_k(\infty; \mathbf{z})}. Writing
#' \deqn{A_k(t; \mathbf{z}) =
#' 1 + \widehat{\alpha}_k
#' \exp(\mathbf{z}^{\top}\widehat{\boldsymbol{\beta}}_k)
#' \widehat{u}_k(t),}
#' the fitted cumulative incidence function is
#' \eqn{\widehat{F}_k(t; \mathbf{z}) =
#' 1 - A_k(t; \mathbf{z})^{-1/\widehat{\alpha}_k}}, so a cure fraction exists only
#' if \eqn{A_k(t; \mathbf{z})} converges to a finite positive limit. Two conditions are
#' therefore involved: whether the fitted baseline cumulative hazard
#' \eqn{\widehat{u}_k(t)} converges as \eqn{t \to \infty}, and whether
#' \eqn{A_k(t; \mathbf{z})} stays positive over \eqn{(0, \infty)}. The second one
#' depends on the covariates through
#' \eqn{\exp(\mathbf{z}^{\top}\widehat{\boldsymbol{\beta}}_k)}, so different
#' covariate profiles from a single fit can fall into different cases.
#'
#' \emph{Baseline limit.} The fitted baseline cumulative hazard
#' \eqn{\widehat{u}_k(t)} increases from \eqn{\widehat{u}_k(0) = 0}. Its limiting
#' behavior depends on the baseline distribution:
#'
#' \itemize{
#'   \item For \code{"gompertz2"}, \eqn{\widehat{u}_k(t)} converges to
#'   \eqn{-\frac{\widehat{\tau}_k}{\widehat{\rho}_k}} when
#'   \eqn{\widehat{\rho}_k < 0}. It diverges when
#'   \eqn{\widehat{\rho}_k \geq 0}.
#'
#'   \item For \code{"gompertz3"}, \eqn{\widehat{u}_k(t)} converges when
#'   \eqn{\widehat{\rho}_k < 0}, with limit
#'   \eqn{-\frac{\widehat{\tau}_k}{\widehat{\rho}_k\widehat{\eta}_k}
#'   (e^{\widehat{\eta}_k}-1)}, and also when
#'   \eqn{\widehat{\rho}_k > 0} together with
#'   \eqn{\widehat{\eta}_k < 0}, with limit
#'   \eqn{-\frac{\widehat{\tau}_k}{\widehat{\rho}_k\widehat{\eta}_k}
#'   e^{\widehat{\eta}_k}}.
#'
#'   \item For \code{"logistic"}, under the parameter space
#'   \eqn{\widehat{b}_k > 0}, \eqn{-\infty < \widehat{c}_k < \infty}, and
#'   \eqn{0 < \widehat{p}_k < 1}, \eqn{\widehat{u}_k(t)} always converges to
#'   \eqn{-\log(1-\widehat{p}_k)} as \eqn{t \to \infty}. Thus,
#'   \eqn{\widehat{u}_k(t)} does not diverge within the parameter space.
#' }
#'
#' Thus, a positive \eqn{\widehat{\rho}_k} does not by itself rule out a cure
#' fraction under \code{"gompertz3"}, unlike under \code{"gompertz2"}.
#'
#' \emph{Cure classification.}
#' The existence of a cure fraction is determined by the limiting behavior of
#' \eqn{\widehat{u}_k(t)} and the sign of \eqn{A_k(t;\mathbf{z})}.
#'
#' \describe{
#'   \item{\code{"cure"}}{
#'     \eqn{u_k(\infty)<\infty} and
#'     \eqn{A_k(\infty;\mathbf{z})>0}.
#'     The cumulative incidence function levels off below one, yielding a
#'     positive cure fraction.
#'   }
#'   \item{\code{"no cure (asymptotic)"}}{
#'     (\eqn{u_k(\infty)=\infty}, \eqn{A_k(\infty;\mathbf{z})>0}) or
#'     (\eqn{u_k(\infty)<\infty}, \eqn{A_k(\infty;\mathbf{z})=0}).
#'     The cumulative incidence function is defined for all finite
#'     \eqn{t} and converges to one, so the cure fraction is zero.
#'   }
#'   \item{\code{"no cure (finite support)"}}{
#'     \eqn{A_k(\infty;\mathbf{z})<0}.
#'     \eqn{A_k(t;\mathbf{z})} reaches zero at a finite time
#'     \eqn{t_{\mathrm{boundary}}}, at which point the cumulative incidence
#'     function reaches one and is undefined thereafter. Thus, a cure fraction
#'     is not defined.
#'   }
#' }
#'
#'
#' @param object a fitted model object.
#' @param ... further arguments passed to or from methods.
#'
#' @return
#' A numeric vector containing the estimated cure fractions.
#'
#'
#' @export
cure <- function(object, ...) UseMethod("cure")

#' @rdname cure
#'
#' @param cov numeric matrix of covariate values.
#' @param event event type for which the cure fraction is computed.
#'  The default is the failure type of interest.
#'
#' @return
#' An object of class \code{"cure.pcrr"} containing the following components:
#' \describe{
#'   \item{\code{cure}}{a numeric vector of estimated cure fractions, with one
#'   entry per row of \code{cov}. the value is positive for \code{"cure"}
#'   profiles, zero for \code{"no cure (asymptotic)"} profiles, and
#'   \code{NA} for \code{"no cure (finite support)"} profiles.}
#'   \item{\code{status}}{a character vector classifying each covariate profile
#'   as \code{"cure"}, \code{"no cure (asymptotic)"}, or
#'   \code{"no cure (finite support)"}.}
#'   \item{\code{t_boundary}}{a numeric vector giving the time at which
#'   \eqn{A_k(t; \mathbf{z})} reaches zero, if such a finite time exists, and
#'   \code{Inf} otherwise. It is finite only for
#'   \code{"no cure (finite support)"} profiles.}
#' }
#'
#' For profiles classified as \code{"no cure (finite support)"}, the
#' cumulative incidence function reaches the boundary of the GOR
#' transformation domain, which is 1, at \code{t_boundary} and is undefined afterwards.
#'
#' @seealso
#' \code{\link{pcrr}}
#' \code{\link{print.cure.pcrr}}
#'
#' @export
cure.pcrr <- function(object, cov, event = NULL, ...){
  P <- object$p
  if (is.null(event)) event <- object$mapping[1]      # default : failcode
  k <- match(event, object$mapping)
  if (is.na(k)) stop("event must be one of: ", paste(object$mapping, collapse = ", "))
  
  if (!is.matrix(cov)) {
    if (is.vector(cov)){
      if (length(cov) %% P == 0) cov <- matrix(cov, ncol = P, byrow = TRUE)
      else stop("cov must have ", P, " column(s).")
    } else cov <- as.matrix(cov)
  }
  if (!is.numeric(cov)) stop("cov must be numeric.")
  if (ncol(cov) != P) stop("cov must have ", P, " column(s).")
  
  
  if (object$distribution == "gompertz2" || object$distribution == "gompertz3"){
    u_k <- function(t, rho, tau, eta, tol){
      if (abs(rho) < tol && abs(eta) < tol){
        u <- tau * t
      } else if (abs(rho) < tol){
        u <- tau * exp(eta) * t
      } else if (abs(eta) < tol){
        u <- tau * expm1(rho * t) / rho
      } else {
        u <- tau * exp(eta) * expm1(eta * expm1(rho * t)) / (rho * eta)
      }
      u
    }
  } else if (object$distribution == "logistic"){
    u_k <- function(t, b, c, p, tol){
      inner <- - p + p*(1+exp(-b*c))/(1+exp(b*(t-c)))
      
      if (1 + inner > 0) u <- -log1p(inner)
      else u <- NaN
      u
    }
  }
  
  if (object$distribution == "gompertz2"){
    alpha <- object$coef[(k - 1) * (3 + P) + 1]
    rho   <- object$coef[(k - 1) * (3 + P) + 2]
    tau   <- object$coef[(k - 1) * (3 + P) + 3]
    eta   <- 0
    beta  <- object$coef[((k - 1) * (3 + P) + 4):(k * (3 + P))]
    
    u_Inf <- u_k(Inf, rho, tau, eta, 1e-12)
  } else if (object$distribution == "gompertz3"){
    alpha <- object$coef[(k - 1) * (4 + P) + 1]
    rho   <- object$coef[(k - 1) * (4 + P) + 2]
    tau   <- object$coef[(k - 1) * (4 + P) + 3]
    eta   <- object$coef[(k - 1) * (4 + P) + 4]
    beta  <- object$coef[((k - 1) * (4 + P) + 5):(k * (4 + P))]
    
    u_Inf <- u_k(Inf, rho, tau, eta, 1e-12)
  } else if (object$distribution == "logistic"){
    alpha <- object$coef[(k - 1) * (4 + P) + 1]
    b   <- object$coef[(k - 1) * (4 + P) + 2]
    c   <- object$coef[(k - 1) * (4 + P) + 3]
    p   <- object$coef[(k - 1) * (4 + P) + 4]
    beta  <- object$coef[((k - 1) * (4 + P) + 5):(k * (4 + P))]
    
    u_Inf <- u_k(Inf, b, c, p, 1e-12)
  }
  
  t_boundary <- suppressWarnings(predict.pcrr(object, cov, event = event)$t_boundary)
  
  cure <- rep(NA_real_, nrow(cov))
  status <- character(nrow(cov))
  
  for (j in seq_len(nrow(cov))){
    ezb <- exp(sum(cov[j, ] * beta))

    if (is.finite(t_boundary[j])){
      status[j] <- "no cure (finite support)"
      cure[j] <- NA_real_
    } else{
      if (is.finite(u_Inf)){
        if (abs(alpha) < 1e-12) cure[j] <- exp(-ezb * u_Inf)
        else cure[j] <- (1.0 + alpha * ezb * u_Inf)^(-1.0 / alpha)
        status[j] <- "cure"
      } else {
        cure[j] <- 0
        status[j] <- "no cure (asymptotic)"
      }
    }
  }

  result <- list(
    cure = unname(cure),
    status = unname(status),
    t_boundary = unname(t_boundary)
  )
  
  class(result) <- "cure.pcrr"
  
  result
}


#' Print Estimated Cure Fractions
#'
#' Prints a summary of estimated cure fractions produced by
#' \code{\link{cure.pcrr}}, including the estimated cure fraction,
#' its classification, and the boundary time of the cumulative incidence
#' function when a finite boundary exists.
#'
#' @param x object of class \code{"cure.pcrr"}.
#' @param digits number of decimal places to format the cure fractions
#'   and finite boundary times.
#' @param ... further arguments passed to \code{\link[base]{print}}.
#'
#' @return
#' The function prints the cure fraction information and returns the input
#' object invisibly.
#'
#' @seealso
#' \code{\link{pcrr}},
#' \code{\link{cure.pcrr}}
#'
#' @export
print.cure.pcrr <- function(x, digits = 8, ...) {
  
  for (j in seq_along(x$cure)) {
    
    width <- digits + 4
    cure_txt <- if (is.na(x$cure[j])) {
      sprintf(paste0("%-", width, "s"), "NA")
    } else {
      sprintf(paste0("%-", width, ".", digits, "f"), x$cure[j])
    }
    
    if (is.finite(x$t_boundary[j])) interval_txt <- sprintf("[0, %s]", x$t_boundary[j])
    else interval_txt <- sprintf("[0, %s)", x$t_boundary[j])
    
    cat("obs", j, ":", cure_txt,"<", x$status[j], ">  ",interval_txt, "\n")
  }
  
  invisible(x)
}


# Kernel Function
.init_values_gom2 <- function(x, delta, z){
  K <- ncol(delta)
  P <- ncol(z)
  
  # theta_init : vector
  # Order : c(alpha1, rho1, tau1, beta11, beta12, ..., alpha2, rho2, tau2, beta21, beta22, ..., alpha3, rho3, tau3...)
  theta_init <- numeric(K * (3 + P))
  
  # [theta Index]
  # alpha_k : theta[(k - 1) * (3 + P) + 1]
  # rho_k : theta[(k - 1) * (3 + P) + 2]
  # tau_k : theta[(k - 1) * (3 + P) + 3]
  # beta_kp : theta[(k - 1) * (3 + P) + p + 3]
  
  
  # [alpha] : fit 1
  for (k in 1:K){
    theta_init[(k - 1) * (3 + P) + 1] <- 1.0
  }
  
  
  # [rho, tau] : baseline MLE
  # baseline function
  baseline <- function(theta, x, delta, k_num)
  {
    k <- k_num
    sum( -delta[, k] * (theta[1] * x + log(theta[2])) +
           theta[2] * expm1(theta[1] * x) / theta[1] )
  }
  
  # init tau, rho
  for (k in 1:K){
    fit <- suppressWarnings(nlminb(start=c(-0.5,0.5), objective=baseline, x=x, delta=delta, k_num=k))
    if (fit$convergence == 0){
      theta_init[(k - 1) * (3 + P) + 2] <- fit$par[1]
      theta_init[(k - 1) * (3 + P) + 3] <- fit$par[2]
    }else{
      stop("Baseline Convergence Failed: ", fit$message, call. = FALSE)
    }
  }
  
  
  # [beta] : Cox Model
  for (k in 1:K){
    fit <- suppressWarnings(coxph(Surv(x, delta[, k]) ~ z))
    for (p in 1:P){
      theta_init[(k - 1) * (3 + P) + p + 3] <- fit$coef[p]
    }
  }
  
  theta_init
}
.log_lik_gom2 <- function(x, delta, z, theta){
  
  if (any(!is.finite(theta))) return(1e+100)
  
  N <- length(x)
  K <- ncol(delta)
  P <- ncol(z)
  
  tol <- 1e-8
  
  log_l <- numeric(N) # return value : -sum(log_l)
  
  # event case
  for (k in 1:K) {
    alpha <- theta[(k - 1) * (3 + P) + 1]
    rho   <- theta[(k - 1) * (3 + P) + 2]
    tau   <- theta[(k - 1) * (3 + P) + 3]
    beta  <- theta[((k - 1) * (3 + P) + 4):(k * (3 + P))]
    
    zb <- drop(z %*% beta)
    
    if (abs(rho) < tol) {
      u <- tau * x
    } else {
      u <- tau * expm1(rho * x) / rho
    }
    
    if (abs(alpha) < tol){
      log_GOR_fk_term <- -exp(zb) * u
    } else {
      log_GOR_fk_term <- -(1.0 / alpha + 1.0) * log1p(alpha * exp(zb) * u)
    }
    
    log_l <- log_l + delta[, k] * (log(tau) + zb + rho * x + log_GOR_fk_term)
    
  }
  # censoring case
  Fk_sum_term <- numeric(N)
  for (k in 1:K) {
    alpha <- theta[(k - 1) * (3 + P) + 1]
    rho   <- theta[(k - 1) * (3 + P) + 2]
    tau   <- theta[(k - 1) * (3 + P) + 3]
    beta  <- theta[((k - 1) * (3 + P) + 4):(k * (3 + P))]
    
    zb <- drop(z %*% beta)
    
    if (abs(rho) < tol) {
      u <- tau * x
    } else {
      u <- tau * expm1(rho * x) / rho
    }
    
    if (abs(alpha) < tol){
      GOR_Fk <- 1.0 - exp(-exp(zb) * u)
    } else {
      GOR_Fk <- 1.0 - (1.0 + alpha * exp(zb) * u)^(-1.0 / alpha)
    }
    
    Fk_sum_term <- Fk_sum_term + GOR_Fk
  }
  
  censor <- (rowSums(delta) == 0)
  log_l[censor] <- log_l[censor] + log1p(-Fk_sum_term[censor])
  
  sum_log_l <- sum(log_l)
  if (!is.finite(sum_log_l)) return(1e+100)
  
  return(-sum_log_l)
}
.estimate_mle_gom2 <- function(x, delta, z, theta_init, gtol, maxiter){
  
  N <- length(x)
  K <- ncol(delta)
  P <- ncol(z)
  
  # Estimate
  lower <- rep(-Inf, K * (3 + P))
  for (k in 1:K){
    lower[(k - 1) * (3 + P) + 3] <- 1e-6   # tau > 0
  }
  control <- list(rel.tol = gtol, iter.max  = maxiter, eval.max = 3 * maxiter)
  mle <- suppressWarnings(nlminb(start = theta_init, objective = .log_lik_gom2,
                                 x = x, delta = delta, z = z, control = control, lower = lower))
  
  mle
}
.score_hessian_gom2 <- function(x, delta, z, theta_mle, variance){
  
  N <- length(x)
  K <- ncol(delta)
  P <- ncol(z)
  
  # symbol name : alpha1, rho1, tau1, beta11, beta12, ..., alpha2, ...
  # symbol name : delta1, delta2, ...
  # symbol name : z1, z2, ...
  
  
  # Make Expression Text
  txt_zb <- character(K)
  txt_u <- character(K)
  for (k in 1:K){
    txt_zb[k] <- ""
    for (p in 1:P){
      txt_zb[k] <- paste0(txt_zb[k], "z", p, " * beta", k, p)
      if (p != P){
        txt_zb[k] <- paste0(txt_zb[k], " + ")
      }
    }
    txt_u[k] <- paste0(txt_u[k], "tau", k, " * expm1(rho", k, " * x) / rho", k)
    
  }
  txt_log_term <- character(K)
  txt_F <- character(K)
  txt_sum_F <- ""
  for (k in 1:K){
    txt_log_term[k] <- paste0(txt_log_term[k], "- (1 / alpha", k, " + 1) * ", "log1p(alpha", k, " * exp(", txt_zb[k], ") * ", txt_u[k], ")")
    txt_F[k] <- paste0(txt_F[k], "1 - (1 + alpha", k, " * exp(", txt_zb[k], ") * ", txt_u[k], ")^(-1 / alpha", k, ")")
    
    txt_sum_F <- paste0(txt_sum_F, "(", txt_F[k], ")")
    if (k != K) txt_sum_F <- paste0(txt_sum_F, " + ")
  }
  txt_main <- ""
  for (k in 1:K){
    # event term
    txt_main <- paste0(txt_main, "delta", k, " * (log(tau", k, ") + ", txt_zb[k], " + rho", k, " * x ", txt_log_term[k], ") + ")
  }
  for (k in 1:K){
    # censoring term
    if (k == 1) txt_main <- paste0(txt_main, "(1")
    txt_main <- paste0(txt_main, " - delta", k)
    if (k == K) txt_main <- paste0(txt_main, ") * ")
  }
  txt_main <- paste0(txt_main, "log1p(-(", txt_sum_F, "))")
  
  # Save as Expression
  log_lik_exprs <- parse(text = txt_main)
  
  
  
  # Parameter Name Definition
  n_param <- K * (3 + P)
  param_names <- character(n_param)
  for (k in 1:K){
    # alpha
    param_names[(k - 1) * (3 + P) + 1] <- paste0("alpha", k)
    # rho
    param_names[(k - 1) * (3 + P) + 2] <- paste0("rho", k)
    # tau
    param_names[(k - 1) * (3 + P) + 3] <- paste0("tau", k)
    # beta
    for (p in 1:P){
      param_names[(k - 1) * (3 + P) + p + 3] <- paste0("beta", k, p)
    }
  }
  
  # 0 prevention
  for (k in 1:K) {
    if (abs(theta_mle[(k - 1) * (3 + P) + 1]) == 0) {
      theta_mle[(k - 1) * (3 + P) + 1] <- 1e-100
      warning("The parameter value was very small and replaced with 1e-100. (alpha)")
    }
    if (abs(theta_mle[(k - 1) * (3 + P) + 2]) == 0) {
      theta_mle[(k - 1) * (3 + P) + 2] <- 1e-100
      warning("The parameter value was very small and replaced with 1e-100. (rho)")
    }
  }
  
  names(theta_mle) <- param_names
  # Dynamic Variable Allocation (Caution!)
  for(name in param_names){ # parameter symbol : alpha1, rho1, ...
    assign(name, theta_mle[name])
  }
  for (k in 1:K){ # event symbol : delta1, delta2, ...
    assign(paste0("delta", k), delta[, k])
  }
  for (p in 1:P){ # cov symbol : z1, z2, ...
    assign(paste0("z", p), z[, p])
  }
  
  
  
  # Score Expression Vector
  score_exprs <- setNames(vector("list", length(param_names)), param_names)
  for (p in param_names) {
    score_exprs[[p]] <- D(log_lik_exprs, p)
  }
  # Score
  score <- numeric(n_param)
  #names(score) <- param_names
  for (p in 1:n_param){
    score[p] <- sum( eval( score_exprs[[p]] ) )
  }
  
  if (!variance) return(list(score = score, hessian = NULL))
  
  
  # Hessian Expression Matrix
  hessian_exprs <- matrix(vector("list", n_param * n_param), nrow = n_param, ncol = n_param)
  dimnames(hessian_exprs) <- list(param_names, param_names)
  for (i in 1:n_param) {
    for (j in 1:n_param) {
      hessian_exprs[[i, j]] <- D(score_exprs[[i]], param_names[j])
    }
  }
  # Hessian
  hessian <- matrix(0, nrow = n_param, ncol = n_param)
  #dimnames(hessian) <- list(param_names, param_names)
  for (p1 in 1:n_param){
    for (p2 in 1:n_param){
      hessian[p1, p2] <- sum( eval( hessian_exprs[p1, p2][[1]] ) )
    }
  }
  
  list(score = score, hessian = hessian)
}

.init_values_gom3 <- function(x, delta, z){
  K <- ncol(delta)
  P <- ncol(z)
  
  # theta_init : vector
  # Order : c(alpha1, rho1, tau1, eta1, beta11, ..., alpha2, rho2, tau2, eta2, beta21, ...)
  theta_init <- numeric(K * (4 + P))
  
  # Warm start from the 2-parameter Gompertz fit : the gompertz2 MLE gives good
  # starting values for (alpha, rho, tau, beta); eta then starts near 0, where
  # gompertz3 coincides with gompertz2. This avoids the optimizer getting stuck.
  g2 <- tryCatch({
    t2 <- .init_values_gom2(x, delta, z)
    m2 <- .estimate_mle_gom2(x, delta, z, t2, 1e-6, 200)
    if (m2$convergence == 0) m2$par else NULL
  }, error = function(e) NULL)
  
  if (!is.null(g2)){
    for (k in 1:K){
      a2 <- g2[(k - 1) * (3 + P) + 1]
      r2 <- g2[(k - 1) * (3 + P) + 2]
      t2 <- g2[(k - 1) * (3 + P) + 3]
      b2 <- g2[((k - 1) * (3 + P) + 4):(k * (3 + P))]
      theta_init[(k - 1) * (4 + P) + 1] <- max(a2, 0.1)      # alpha
      theta_init[(k - 1) * (4 + P) + 2] <- r2                # rho
      theta_init[(k - 1) * (4 + P) + 3] <- max(t2, 1e-4)     # tau
      theta_init[(k - 1) * (4 + P) + 4] <- -0.1              # eta (near gompertz2)
      for (p in 1:P)
        theta_init[(k - 1) * (4 + P) + p + 4] <- b2[p]       # beta
    }
    return(theta_init)
  }
  
  # ---- fallback : independent initialization (used only if the gompertz2 warm start fails) ----
  
  # [theta Index]
  # alpha_k : theta[(k - 1) * (4 + P) + 1]
  # rho_k   : theta[(k - 1) * (4 + P) + 2]
  # tau_k   : theta[(k - 1) * (4 + P) + 3]
  # eta_k   : theta[(k - 1) * (4 + P) + 4]
  # beta_kp : theta[(k - 1) * (4 + P) + p + 4]
  
  
  # [alpha] : fit 1
  for (k in 1:K){
    theta_init[(k - 1) * (4 + P) + 1] <- 1.0
  }
  
  
  # [rho, tau] : baseline MLE (2-parameter Gompertz, same as gompertz2)
  # baseline function
  baseline <- function(theta, x, delta, k_num)
  {
    k <- k_num
    sum( -delta[, k] * (theta[1] * x + log(theta[2])) +
           theta[2] * expm1(theta[1] * x) / theta[1] )
  }
  
  # init tau, rho : try several starting points and lower bounds for stability
  for (k in 1:K){
    starts <- list(c(-0.5, 0.5), c(-0.1, 0.1), c(-1.0, 0.3), c(0.05, 0.2))
    ok <- FALSE
    for (st in starts){
      fit <- suppressWarnings(nlminb(start = st, objective = baseline,
                                     lower = c(-Inf, 1e-6), x = x, delta = delta, k_num = k))
      if (fit$convergence == 0 && is.finite(fit$objective)){
        theta_init[(k - 1) * (4 + P) + 2] <- fit$par[1]
        theta_init[(k - 1) * (4 + P) + 3] <- fit$par[2]
        ok <- TRUE
        break
      }
    }
    if (!ok){
      # last resort : rough moment-based guess
      theta_init[(k - 1) * (4 + P) + 2] <- -0.1
      theta_init[(k - 1) * (4 + P) + 3] <- max(mean(delta[, k]) / mean(x), 1e-3)
    }
  }
  
  
  # [eta] : start near 0 so that gompertz3 begins at the gompertz2 fit
  for (k in 1:K){
    theta_init[(k - 1) * (4 + P) + 4] <- -0.5
  }
  
  
  # [beta] : Cox Model
  for (k in 1:K){
    fit <- suppressWarnings(coxph(Surv(x, delta[, k]) ~ z))
    for (p in 1:P){
      theta_init[(k - 1) * (4 + P) + p + 4] <- fit$coef[p]
    }
  }
  
  theta_init
}
.log_lik_gom3 <- function(x, delta, z, theta){
  
  if (any(!is.finite(theta))) return(1e+100)
  
  N <- length(x)
  K <- ncol(delta)
  P <- ncol(z)
  
  tol <- 1e-8
  
  log_l <- numeric(N) # return value : -sum(log_l)
  
  # event case
  for (k in 1:K){
    alpha  <- theta[(k - 1) * (4 + P) + 1]
    rho  <- theta[(k - 1) * (4 + P) + 2]
    tau <- theta[(k - 1) * (4 + P) + 3]
    eta <- theta[(k - 1) * (4 + P) + 4]
    beta  <- theta[((k - 1) * (4 + P) + 5):(k * (4 + P))]
    
    zb <- drop(z %*% beta)
    
    if (abs(rho) < tol && abs(eta) < tol){
      u <- tau * x
    } else if (abs(rho) < tol){
      u <- tau * exp(eta) * x
    } else if (abs(eta) < tol){
      u <- tau * expm1(rho * x) / rho
    } else {
      u <- tau * exp(eta) * expm1(eta * expm1(rho * x)) / (rho * eta)
    }
    
    if (abs(alpha) < tol){
      log_GOR_fk_term <- -exp(zb) * u
    } else {
      log_GOR_fk_term <- -(1.0 / alpha + 1.0) * log1p(alpha * exp(zb) * u)
    }
    
    log_l <- log_l + delta[, k] * (log(tau) + zb + rho * x + eta * exp(rho * x) + log_GOR_fk_term)
    
  }
  # censoring case
  Fk_sum_term <- numeric(N)
  for (k in 1:K){
    alpha  <- theta[(k - 1) * (4 + P) + 1]
    rho  <- theta[(k - 1) * (4 + P) + 2]
    tau <- theta[(k - 1) * (4 + P) + 3]
    eta <- theta[(k - 1) * (4 + P) + 4]
    beta  <- theta[((k - 1) * (4 + P) + 5):(k * (4 + P))]
    
    zb <- drop(z %*% beta)
    
    if (abs(rho) < tol && abs(eta) < tol){
      u <- tau * x
    } else if (abs(rho) < tol){
      u <- tau * exp(eta) * x
    } else if (abs(eta) < tol){
      u <- tau * expm1(rho * x) / rho
    } else {
      u <- tau * exp(eta) * expm1(eta * expm1(rho * x)) / (rho * eta)
    }
    
    if (abs(alpha) < tol){
      GOR_Fk <- 1.0 - exp(-exp(zb) * u)
    } else {
      GOR_Fk <- 1.0 - (1.0 + alpha * exp(zb) * u)^(-1.0 / alpha)
    }
    
    Fk_sum_term <- Fk_sum_term + GOR_Fk
  }
  censor <- (rowSums(delta) == 0)
  log_l[censor] <- log_l[censor] + log1p(-Fk_sum_term[censor])
  
  sum_log_l <- sum(log_l)
  if (!is.finite(sum_log_l)) return(1e+100)
  
  return(-sum_log_l)
}
.estimate_mle_gom3 <- function(x, delta, z, theta_init, gtol, maxiter){
  N <- length(x)
  K <- ncol(delta)
  P <- ncol(z)
  
  lower <- rep(-Inf, K * (4 + P))
  for (k in 1:K){
    lower[(k - 1) * (4 + P) + 3] <- 1e-6   # tau > 0
  }
  
  # Coarse control
  control_coarse <- list(rel.tol = 1e-4, iter.max = 50, eval.max = 150)
  # Fine control
  control_fine <- list(rel.tol = gtol, iter.max = maxiter, eval.max = 3 * maxiter)
  
  fit_coarse <- function(st){
    suppressWarnings(nlminb(start = st, objective = .log_lik_gom3,
                            x = x, delta = delta, z = z,
                            lower = lower, control = control_coarse))
  }
  fit_fine <- function(st){
    suppressWarnings(nlminb(start = st, objective = .log_lik_gom3,
                            x = x, delta = delta, z = z,
                            lower = lower, control = control_fine))
  }
  
  # initial parameter
  starts <- list(theta_init)
  rho_start <- c(-0.3, 0.3)
  eta_start <- c(-1, -0.5, -0.1)
  for (r in rho_start){
    for (e in eta_start){
      st <- theta_init
      for (k in 1:K){
        st[(k - 1) * (4 + P) + 2] <- r
        st[(k - 1) * (4 + P) + 4] <- e
      }
      starts[[length(starts) + 1]] <- st
    }
  }
  
  # step 1
  fits_coarse <- lapply(starts, function(st){
    fit <- fit_coarse(st)
    list(par = fit$par, objective = fit$objective)
  })
  
  # select top 2
  finite_idx <- which(sapply(fits_coarse, function(f) is.finite(f$objective)))
  
  if (length(finite_idx) == 0){
    return(fit_fine(theta_init))
  }
  
  top_idx <- finite_idx[order(sapply(fits_coarse[finite_idx], `[[`, "objective"))]
  top_idx <- head(top_idx, 2)
  
  # step 2
  fits_fine <- lapply(top_idx, function(i){
    fit_fine(fits_coarse[[i]]$par)
  })
  
  fits_fine <- Filter(function(f) is.finite(f$objective), fits_fine)
  
  if (length(fits_fine) == 0){
    return(fit_fine(theta_init))
  }
  
  mle <- fits_fine[[which.min(sapply(fits_fine, `[[`, "objective"))]]
  mle
}
.score_hessian_gom3 <- function(x, delta, z, theta_mle, variance){
  
  N <- length(x)
  K <- ncol(delta)
  P <- ncol(z)
  
  # symbol name : alpha1, rho1, tau1, eta1, beta11, beta12, ..., alpha2, ...
  # symbol name : delta1, delta2, ...
  # symbol name : z1, z2, ...
  
  # Make Expression Text
  txt_zb <- character(K)
  txt_u <- character(K)
  for (k in 1:K){
    txt_zb[k] <- ""
    for (p in 1:P){
      txt_zb[k] <- paste0(txt_zb[k], "z", p, " * beta", k, p)
      if (p != P){
        txt_zb[k] <- paste0(txt_zb[k], " + ")
      }
    }
    txt_u[k] <- paste0(txt_u[k], "tau", k, " * exp(eta", k, ") * expm1(eta", k, " * expm1(rho", k, " * x)) / (rho", k, " * eta", k, ")")
  }
  txt_log_term <- character(K)
  txt_F <- character(K)
  txt_sum_F <- ""
  for (k in 1:K){
    txt_log_term[k] <- paste0(txt_log_term[k], "- (1 / alpha", k, " + 1) * ", "log1p(alpha", k, " * exp(", txt_zb[k], ") * ", txt_u[k], ")")
    txt_F[k] <- paste0(txt_F[k], "1 - (1 + alpha", k, " * exp(", txt_zb[k], ") * ", txt_u[k], ")^(-1 / alpha", k, ")")
    
    txt_sum_F <- paste0(txt_sum_F, "(", txt_F[k], ")")
    if (k != K) txt_sum_F <- paste0(txt_sum_F, " + ")
  }
  txt_main <- ""
  for (k in 1:K){
    # event term
    txt_main <- paste0(txt_main, "delta", k, " * (log(tau", k, ") + ", txt_zb[k], " + rho", k, " * x + ",
                       "eta", k, " * exp(rho", k, " * x) ", txt_log_term[k], ") + ")
  }
  for (k in 1:K){
    # censoring term
    if (k == 1) txt_main <- paste0(txt_main, "(1")
    txt_main <- paste0(txt_main, " - delta", k)
    if (k == K) txt_main <- paste0(txt_main, ") * ")
  }
  txt_main <- paste0(txt_main, "log1p(-(", txt_sum_F, "))")
  
  
  
  # Save as Expression
  log_lik_exprs <- parse(text = txt_main)
  
  
  
  # Parameter Name Definition
  n_param <- K * (4 + P)
  param_names <- character(n_param)
  for (k in 1:K){
    # alpha
    param_names[(k - 1) * (4 + P) + 1] <- paste0("alpha", k)
    # rho
    param_names[(k - 1) * (4 + P) + 2] <- paste0("rho", k)
    # tau
    param_names[(k - 1) * (4 + P) + 3] <- paste0("tau", k)
    # eta
    param_names[(k - 1) * (4 + P) + 4] <- paste0("eta", k)
    # beta
    for (p in 1:P){
      param_names[(k - 1) * (4 + P) + p + 4] <- paste0("beta", k, p)
    }
  }
  
  
  # 0 prevention
  for (k in 1:K) {
    if (abs(theta_mle[(k - 1) * (4 + P) + 1]) == 0) {
      theta_mle[(k - 1) * (4 + P) + 1] <- 1e-100
      warning("The parameter value was very small and replaced with 1e-100. (alpha)")
    }
    if (abs(theta_mle[(k - 1) * (4 + P) + 2]) == 0) {
      theta_mle[(k - 1) * (4 + P) + 2] <- 1e-100
      warning("The parameter value was very small and replaced with 1e-100. (rho)")
    }
    if (abs(theta_mle[(k - 1) * (4 + P) + 4]) == 0) {
      theta_mle[(k - 1) * (4 + P) + 4] <- 1e-100
      warning("The parameter value was very small and replaced with 1e-100. (eta)")
    }
  }
  
  names(theta_mle) <- param_names
  # Dynamic Variable Allocation (Caution!)
  for(name in param_names){ # parameter symbol : alpha1, rho1, ...
    assign(name, theta_mle[name])
  }
  for (k in 1:K){ # event symbol : delta1, delta2, ...
    assign(paste0("delta", k), delta[, k])
  }
  for (p in 1:P){ # cov symbol : z1, z2, ...
    assign(paste0("z", p), z[, p])
  }
  
  
  
  # Score Expression Vector
  score_exprs <- setNames(vector("list", length(param_names)), param_names)
  for (p in param_names) {
    score_exprs[[p]] <- D(log_lik_exprs, p)
  }
  # Score
  score <- numeric(n_param)
  #names(score) <- param_names
  for (p in 1:n_param){
    score[p] <- sum( eval( score_exprs[[p]] ) )
  }
  
  if (!variance) return(list(score = score, hessian = NULL))
  
  
  # Hessian Expression Matrix
  hessian_exprs <- matrix(vector("list", n_param * n_param), nrow = n_param, ncol = n_param)
  dimnames(hessian_exprs) <- list(param_names, param_names)
  for (i in 1:n_param) {
    for (j in 1:n_param) {
      hessian_exprs[[i, j]] <- D(score_exprs[[i]], param_names[j])
    }
  }
  # Hessian
  hessian <- matrix(0, nrow = n_param, ncol = n_param)
  #dimnames(hessian) <- list(param_names, param_names)
  for (p1 in 1:n_param){
    for (p2 in 1:n_param){
      hessian[p1, p2] <- sum( eval( hessian_exprs[p1, p2][[1]] ) )
    }
  }
  
  list(score = score, hessian = hessian)
}

.init_values_logi <- function(x, delta, z){
  K <- ncol(delta)
  P <- ncol(z)
  
  # theta_init : vector
  # Order : c(alpha1, b1, c1, p1, beta11, beta12, ..., alpha2, b2, c2, p2, beta21, beta22, ..., alpha3, b3, c3...)
  theta_init <- numeric(K * (4 + P))
  
  # [theta Index]
  # alpha_k : theta[(k - 1) * (4 + P) + 1]
  # b_k : theta[(k - 1) * (4 + P) + 2]
  # c_k : theta[(k - 1) * (4 + P) + 3]
  # p_k : theta[(k - 1) * (4 + P) + 4]
  # beta_kp : theta[(k - 1) * (4 + P) + p + 4]
  
  
  # [alpha] : fit 1
  for (k in 1:K){
    theta_init[(k - 1) * (4 + P) + 1] <- 1.0
  }
  
  
  # [b, c, p] : baseline MLE
  # baseline function
  baseline <- function(theta, x, delta, k_num)
  {
    k <- k_num
    
    b <- theta[1]
    c <- theta[2]
    p <- theta[3]
    
    bx <- b * (x - c)
    
    log_hazard <- log(b) + log(p) + log1p(exp(-b * c)) + bx - log1p(exp(bx)) -
                  log1p(p * exp(-b * c) +(1 - p) * exp(bx))
    
    u <- log1p(exp(bx)) - log1p(p * exp(-b * c) + (1 - p) * exp(bx))
    
    sum(-delta[, k] * log_hazard + u)
  }
  
  # init tau, rho
  for (k in 1:K){
    fit <- suppressWarnings(nlminb(start=c(0.5,0,0.5), objective=baseline, x=x, delta=delta, k_num=k))
    if (fit$convergence == 0){
      theta_init[(k - 1) * (4 + P) + 2] <- fit$par[1]
      theta_init[(k - 1) * (4 + P) + 3] <- fit$par[2]
      theta_init[(k - 1) * (4 + P) + 4] <- fit$par[3]
    }else{
      stop("Baseline Convergence Failed: ", fit$message, call. = FALSE)
    }
  }
  
  
  # [beta] : Cox Model
  for (k in 1:K){
    fit <- suppressWarnings(coxph(Surv(x, delta[, k]) ~ z))
    for (p in 1:P){
      theta_init[(k - 1) * (4 + P) + p + 4] <- fit$coef[p]
    }
  }
  
  theta_init
}
.log_lik_logi <- function(x, delta, z, theta){
  
  if (any(!is.finite(theta))) return(1e+100)
  
  N <- length(x)
  K <- ncol(delta)
  P <- ncol(z)
  
  tol <- 1e-8
  
  log_l <- numeric(N) # return value : -sum(log_l)
  
  # event case
  for (k in 1:K) {
    alpha  <- theta[(k - 1) * (4 + P) + 1]
    b  <- theta[(k - 1) * (4 + P) + 2]
    c <- theta[(k - 1) * (4 + P) + 3]
    p <- theta[(k - 1) * (4 + P) + 4]
    beta  <- theta[((k - 1) * (4 + P) + 5):(k * (4 + P))]
    
    zb <- drop(z %*% beta)
    
    u <- -log1p(- p + p * (1 + exp(-b * c)) / (1 + exp(b * (x - c))))
    
    if (abs(alpha) < tol){
      log_GOR_fk_term <- -exp(zb) * u
    } else {
      log_GOR_fk_term <- -(1.0 / alpha + 1.0) * log1p(alpha * exp(zb) * u)
    }
    
    log_l <- log_l + delta[, k] * (zb + log(b) + log(p) + log1p(exp(-b * c)) + b * (x - c) -
                                     log1p(exp(b * (x - c))) - log1p(p * exp(-b * c) + (1 - p) * exp(b * (x - c))) + log_GOR_fk_term)
    
  }
  # censoring case
  Fk_sum_term <- numeric(N)
  for (k in 1:K) {
    alpha  <- theta[(k - 1) * (4 + P) + 1]
    b  <- theta[(k - 1) * (4 + P) + 2]
    c <- theta[(k - 1) * (4 + P) + 3]
    p <- theta[(k - 1) * (4 + P) + 4]
    beta  <- theta[((k - 1) * (4 + P) + 5):(k * (4 + P))]
    
    zb <- drop(z %*% beta)
    
    u <- -log1p(- p + p * (1 + exp(-b * c)) / (1 + exp(b * (x - c))))
    
    if (abs(alpha) < tol){
      GOR_Fk <- 1.0 - exp(-exp(zb) * u)
    } else {
      GOR_Fk <- 1.0 - (1.0 + alpha * exp(zb) * u)^(-1.0 / alpha)
    }
    
    Fk_sum_term <- Fk_sum_term + GOR_Fk
    
  }
  censor <- (rowSums(delta) == 0)
  log_l[censor] <- log_l[censor] + log1p(-Fk_sum_term[censor])
  
  sum_log_l <- sum(log_l)
  if (!is.finite(sum_log_l)) return(1e+100)
  
  return(-sum_log_l)
}
.estimate_mle_logi <- function(x, delta, z, theta_init, gtol, maxiter){
  
  N <- length(x)
  K <- ncol(delta)
  P <- ncol(z)
  
  # Estimate
  lower <- rep(-Inf, K * (4 + P))
  upper <- rep(Inf, K * (4 + P))
  for (k in 1:K){
    lower[(k - 1) * (4 + P) + 2] <- 1e-6   # b > 0
    lower[(k - 1) * (4 + P) + 4] <- 1e-6   # p > 0
    upper[(k - 1) * (4 + P) + 4] <- 1-1e-6 # p < 1
  }
  control <- list(rel.tol = gtol, iter.max  = maxiter, eval.max = 3 * maxiter)
  mle <- suppressWarnings(nlminb(start = theta_init, objective = .log_lik_logi,
                                 x = x, delta = delta, z = z, control = control, lower = lower, upper = upper))
  
  mle
}
.score_hessian_logi <- function(x, delta, z, theta_mle, variance){
  
  N <- length(x)
  K <- ncol(delta)
  P <- ncol(z)
  
  # symbol name : alpha1, b1, c1, p1, beta11, beta12, ..., alpha2, ...
  # symbol name : delta1, delta2, ...
  # symbol name : z1, z2, ...
  
  # Make Expression Text
  txt_zb <- character(K)
  txt_u <- character(K)
  for (k in 1:K){
    txt_zb[k] <- ""
    for (p in 1:P){
      txt_zb[k] <- paste0(txt_zb[k], "z", p, " * beta", k, p)
      if (p != P){
        txt_zb[k] <- paste0(txt_zb[k], " + ")
      }
    }
    txt_u[k] <- paste0(txt_u[k], "log1p(-p", k, " + p", k, " * (exp(-b", k, " * c", k, ") + 1) / (1 + exp(b", k ," * (x - c", k, "))))")
  }
  txt_log_term <- character(K)
  txt_F <- character(K)
  txt_sum_F <- ""
  for (k in 1:K){
    txt_log_term[k] <- paste0(txt_log_term[k], "- (1 / alpha", k, " + 1) * ", "log1p(-alpha", k, " * exp(", txt_zb[k], ") * ", txt_u[k], ")")
    txt_F[k] <- paste0(txt_F[k], "1 - (1 - alpha", k, " * exp(", txt_zb[k], ") * ", txt_u[k], ")^(-1 / alpha", k, ")")
    
    txt_sum_F <- paste0(txt_sum_F, "(", txt_F[k], ")")
    if (k != K) txt_sum_F <- paste0(txt_sum_F, " + ")
  }
  txt_main <- ""
  for (k in 1:K){
    # event term
    txt_main <- paste0(txt_main, "delta", k, " * (", txt_zb[k], " + log(b", k, ")", " + log(p", k, ")",
                       " + log1p(exp(-b", k, " * c", k, ")) + b", k," * (x - c", k, ")",
                       " - log1p(exp(b", k, " * (x - c", k, ")))",
                       " - log1p(p", k, " * exp(-b", k, " * c", k, ") + (1 - p", k, ") * exp(b", k, " * (x - c", k, "))) ",
                       txt_log_term[k], ") + ")
  }
  for (k in 1:K){
    # censoring term
    if (k == 1) txt_main <- paste0(txt_main, "(1")
    txt_main <- paste0(txt_main, " - delta", k)
    if (k == K) txt_main <- paste0(txt_main, ") * ")
  }
  txt_main <- paste0(txt_main, "log1p(-(", txt_sum_F, "))")
  
  
  
  # Save as Expression
  log_lik_exprs <- parse(text = txt_main)
  
  
  
  # Parameter Name Definition
  n_param <- K * (4 + P)
  param_names <- character(n_param)
  for (k in 1:K){
    # alpha
    param_names[(k - 1) * (4 + P) + 1] <- paste0("alpha", k)
    # b
    param_names[(k - 1) * (4 + P) + 2] <- paste0("b", k)
    # c
    param_names[(k - 1) * (4 + P) + 3] <- paste0("c", k)
    # p
    param_names[(k - 1) * (4 + P) + 4] <- paste0("p", k)
    # beta
    for (p in 1:P){
      param_names[(k - 1) * (4 + P) + p + 4] <- paste0("beta", k, p)
    }
  }
  
  
  # 0 prevention
  for (k in 1:K) {
    if (abs(theta_mle[(k - 1) * (4 + P) + 1]) == 0) {
      theta_mle[(k - 1) * (4 + P) + 1] <- 1e-100
      warning("The parameter value was very small and replaced with 1e-100. (alpha)")
    }
  }
  
  names(theta_mle) <- param_names
  # Dynamic Variable Allocation (Caution!)
  for(name in param_names){ # parameter symbol : alpha1, b1, ...
    assign(name, theta_mle[name])
  }
  for (k in 1:K){ # event symbol : delta1, delta2, ...
    assign(paste0("delta", k), delta[, k])
  }
  for (p in 1:P){ # cov symbol : z1, z2, ...
    assign(paste0("z", p), z[, p])
  }
  
  
  
  # Score Expression Vector
  score_exprs <- setNames(vector("list", length(param_names)), param_names)
  for (p in param_names) {
    score_exprs[[p]] <- D(log_lik_exprs, p)
  }
  # Score
  score <- numeric(n_param)
  #names(score) <- param_names
  for (p in 1:n_param){
    score[p] <- sum( eval( score_exprs[[p]] ) )
  }
  
  if (!variance) return(list(score = score, hessian = NULL))
  
  
  # Hessian Expression Matrix
  hessian_exprs <- matrix(vector("list", n_param * n_param), nrow = n_param, ncol = n_param)
  dimnames(hessian_exprs) <- list(param_names, param_names)
  for (i in 1:n_param) {
    for (j in 1:n_param) {
      hessian_exprs[[i, j]] <- D(score_exprs[[i]], param_names[j])
    }
  }
  # Hessian
  hessian <- matrix(0, nrow = n_param, ncol = n_param)
  #dimnames(hessian) <- list(param_names, param_names)
  for (p1 in 1:n_param){
    for (p2 in 1:n_param){
      hessian[p1, p2] <- sum( eval( hessian_exprs[p1, p2][[1]] ) )
    }
  }
  
  list(score = score, hessian = hessian)
  
}
