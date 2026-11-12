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
#' The baseline hazard is unimodal when \eqn{\rho_k > 0} and \eqn{-1 < \eta_k < 0}, or
#' when \eqn{\rho_k < 0} and \eqn{\eta_k < -1}, in which case it peaks at
#' \eqn{\frac{1}{\rho_k} \log(-\frac{1}{\eta_k})}. Fitting \code{"gompertz3"} to data with no
#' unimodal hazard usually makes the extra parameter unidentifiable, which shows
#' up as extreme estimates and a singular information matrix.
#' 
#' \emph{Notation.} Haile et al. (2016) write the three baseline parameters as
#' \eqn{(\alpha, \beta, \eta)}. They appear here as \code{rho}, \code{tau} and
#' \code{eta}, because \code{alpha} and \code{beta} already denote the link
#' parameter and the regression coefficients of Jeong and Fine (2007).
#' 
#' With \code{distribution = "logistic"}, the baseline is the three-parameter
#' modified logistic model (Cheng, 2009), with cumulative baseline function
#' \deqn{u_k(t) = -\log\left\{1 - \dfrac{p_k e^{b_k(t-c_k)} - p_k e^{-b_k c_k}}
#' {1 + e^{b_k(t-c_k)}}\right\}.}
#' Here \eqn{b_k} controls the rate of increase and \eqn{c_k} determines the
#' location of the rise, while \eqn{u_k(t) \to -\log(1-p_k)} as
#' \eqn{t\to\infty}, so the baseline cumulative hazard is finite for every
#' \eqn{p_k < 1}.
#'
#' The parameter \eqn{p_k} is the asymptote of the baseline cumulative
#' incidence \eqn{1-\exp\{-u_k(t)\}}, not of the fitted curve. The asymptote
#' implied by the fitted model is
#' \eqn{1 - \{1 + \alpha_k \exp(\mathbf{Z}^{\top}\boldsymbol{\beta}_k)
#' u_k(\infty)\}^{-1/\alpha_k}}, a monotone function of \eqn{p_k} that equals
#' it only as \eqn{\alpha_k \to 0}. Under the proportional odds link
#' (\eqn{\alpha_k = 1}) the two differ materially, so the long-term event
#' probability should be taken from \code{\link{cure.pcrr}} rather than read
#' off \eqn{\widehat{p}_k}. The baseline hazard can be unimodal under
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
#'
#' Note that this situation cannot arise in the fits that are reported. Every
#' model case fixes \eqn{\alpha_k} at 0 or 1 (see \emph{Model cases} below),
#' so \eqn{1 + \alpha_k \exp(\mathbf{Z}^{\top}\boldsymbol{\beta}_k) u_k(t)
#' \geq 1} throughout and no finite boundary exists. The discussion above
#' therefore applies to the unconstrained fit used for the assumption tests,
#' whose coefficients are returned in \code{coef}.
#'
#'
#' \emph{Model cases.} Only two values of \eqn{\alpha_k} admit a direct
#' reading: at \eqn{\alpha_k = 0} the quantity
#' \eqn{\exp(\mathbf{Z}^{\top}\boldsymbol{\beta}_k)} is a subdistribution
#' hazard ratio (proportional hazards, PH), and at \eqn{\alpha_k = 1} it is an
#' odds ratio (proportional odds, PO). At any other value it is neither.
#' \code{pcrr} therefore treats the unconstrained \eqn{\widehat{\alpha}_k} as
#' a diagnostic rather than as a result, and proceeds in three stages.
#'
#' \enumerate{
#'   \item The model is fitted once with every \eqn{\alpha_k} free.
#'   \item For each cause, \eqn{H_0\!: \alpha_k = 0} and
#'   \eqn{H_0\!: \alpha_k = 1} are tested by Wald statistics at
#'   \code{sig.level}. A cause for which neither hypothesis is rejected admits
#'   both links; one for which exactly one is rejected admits the other.
#'   \item The admissible combinations across causes are enumerated as
#'   numbered \emph{model cases}, and the model is refitted at each of them
#'   with the corresponding \eqn{\alpha_k} held fixed.
#' }
#'
#' The case numbering is shown by \code{\link{print.pcrr}} and is the
#' \code{case} argument of \code{\link{summary.pcrr}},
#' \code{\link{predict.pcrr}}, \code{\link{cure.pcrr}},
#' \code{\link{plot.pcrr}} and \code{\link{plot.predict.pcrr}}. The number of
#' cases is the product of the number of links admitted by each cause, so it
#' lies between 1 and \eqn{2^K}; the first cause varies fastest in the
#' enumeration. There is no argument for specifying a link by hand, since a
#' link outside the admissible set is one the data have rejected.
#'
#' \code{pcrr} stops with an error, and returns no object, when no model
#' case survives: for some cause both hypotheses are rejected, one cannot be
#' tested and the other is rejected, or neither can be tested because the
#' standard error of \eqn{\widehat{\alpha}_k} is not finite. The error names
#' the causes responsible. It also stops, before any test, if the
#' unconstrained fit does not converge or its Hessian is singular.
#'
#' \emph{Diagnostics.} \code{\link{plot.pcrr}} draws two screens for each
#' model case. The first overlays the fitted baseline \eqn{u_k(t)} on a
#' semiparametric estimate that assumes the same link, Fine--Gray under PH and
#' a proportional odds subdistribution model under PO, so that a gap between
#' the curves points at the baseline alone; this is the primary check. The
#' second shows Cox--Snell residuals, which are informative only below the
#' shaded limit \eqn{\min_i L_i}. When the fitted cumulative incidence is
#' improper, each residual follows \eqn{\mathrm{Exp}(1)} only up to a
#' subject-specific limit \eqn{L_i}, so beyond the smallest one the
#' 45-degree line no longer applies. Fitted cumulative incidence and
#' subdistribution hazard curves at chosen covariate profiles are produced by
#' \code{\link{predict.pcrr}} and drawn by \code{\link{plot.predict.pcrr}}.
#'
#' Model parameters are estimated by maximum likelihood using
#' \code{\link[stats]{nlminb}}, with standard errors obtained from the inverse
#' of the observed information matrix. The returned object is of class
#' \code{"pcrr"} and supports the S3 methods \code{print()}, \code{summary()},
#' \code{plot()}, \code{predict()}, and \code{cure()}.
#'
#' @seealso
#' \code{\link{print.pcrr}}, 
#' \code{\link{plot.pcrr}}, 
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
#'  \code{nlminb()}. The default value is \code{1e-10}.
#' @param maxiter maximum number of optimization iterations. Internally,
#'  \code{eval.max} is set to \code{3 * maxiter}.
#' @param init a user-specified initial parameter vector.
#'  See \code{\link{pcrr-parameter-order}} for the parameter ordering.
#' @param variance logical value indicating whether variance estimates and
#'  Wald tests are computed for the possible model assumption cases.
#'  Variance estimation for the unconstrained GOR model is always performed
#'  to conduct the PH and PO model assumption tests. If \code{TRUE}, the
#'  inverse information matrix is additionally computed for each possible
#'  model assumption case, allowing Wald inference for those models. If
#'  \code{FALSE}, the optimization for each possible model assumption case
#'  is still performed, but their variance estimates and Wald tests are not
#'  computed. The default is \code{TRUE}.
#' @param sig.level numeric value specifying the significance level used for
#'  the model assumption tests. The default is \code{0.05}.
#'
#' @return
#' An object of class \code{"pcrr"}, which is a list containing the following components:
#' \item{\code{coef}}{estimated model parameters, stored according to the parameter ordering described in \code{\link{pcrr-parameter-order}}.}
#' \item{\code{loglik}}{maximized log-likelihood value.}
#' \item{\code{init}}{initial parameter values used to start the optimization,
#'  stored according to the parameter ordering described in \code{\link{pcrr-parameter-order}}.}
#' \item{\code{score}}{score vector evaluated at the maximum likelihood estimates.}
#' \item{\code{inf}}{observed information matrix (the negative Hessian matrix); \code{NULL} if \code{variance = FALSE}.}
#' \item{\code{invinf}}{inverse of the observed information matrix; \code{NULL} if \code{variance = FALSE}.}
#' \item{\code{converged}}{logical value indicating whether the optimization converged.}
#' \item{\code{iter}}{number of iterations performed by the optimization algorithm.}
#' \item{\code{message}}{message returned by the optimization routine.}
#' \item{\code{call}}{matched function call.}
#' \item{\code{x}}{The processed \code{ftime} vector containing the observed event or censoring times.}
#' \item{\code{event}}{The event of interest used to identify the cause-specific event observations.}
#' \item{\code{z}}{The processed \code{cov} matrix containing the covariate values.}
#' \item{\code{cencode}}{The code used to identify censored observations in the \code{fstatus} vector.}
#' \item{\code{n}}{total number of observations in the original data set.}
#' \item{\code{n_missing}}{number of observations removed due to missing values.}
#' \item{\code{k}}{number of event categories, excluding the censoring category.}
#' \item{\code{p}}{number of covariates.}
#' \item{\code{distribution}}{baseline distribution used for the model.}
#' \item{\code{cov_names}}{names of the covariates.}
#' \item{\code{mapping}}{mapping between the original event labels and the internal event codes.}
#' \item{\code{maxtime}}{maximum observed follow-up time.}
#' \item{\code{case_all}}{matrix with one row per model case and one column
#'  per cause, holding the fixed \eqn{\alpha_k}: 0 for PH and 1 for PO.}
#' \item{\code{case_model}}{labels describing all possible model assumption cases.}
#' \item{\code{mle_case_all}}{the results of \code{nlminb()} optimization for all possible model assumption cases.}
#' \item{\code{sco_case_all}}{score vectors evaluated at the MLEs for all possible model assumption cases.}
#' \item{\code{hess_case_all}}{hessian matrices evaluated at the MLEs for all possible model assumption cases.} 
#' 
#' 
#' @importFrom cmprsk crr
#' @importFrom timereg prop.odds.subdist Event
#' @importFrom survival coxph Surv survfit
#' @importFrom stats D model.matrix na.fail na.omit nlminb pnorm qnorm printCoefmat setNames uniroot approx as.formula
#' @importFrom graphics abline box legend lines mtext par plot.new points rect title
#' @importFrom grDevices dev.interactive devAskNewPage recordGraphics
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
#' print(fit1)
#' 
#' plot(fit1)
#' 
#' summary(fit1, case = c(1, 3))
#'
#' pred1 <- predict(fit1, cov = rbind(c(0, 0.13), c(1, -0.15), c(0, 0.40)), case = c(1, 2))
#' print(pred1)
#' 
#' plot(pred1, case = 1)
#'
#' cure(fit1, cov = rbind(c(0, 0.13), c(1, -0.15), c(0, 0.40)))
#'
#'
#'
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
#' print(fit2)
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
#' Cheng, Y. (2009). Modeling Cumulative Incidences of Dementia and Dementia-Free Death
#'  Using a Novel Three-Parameter Logistic Function.
#' \emph{The International Journal of Biostatistics}, 5(1), Article 29.
#'
#' Cox, D. R. and Snell, E. J. (1968). A general definition of residuals.
#' \emph{Journal of the Royal Statistical Society: Series B}, 30(2), 248--265.
#'
#' Dabrowska, D. M. and Doksum, K. A. (1988). Estimation and testing in a two-sample
#' generalized odds-rate model. \emph{Journal of the American Statistical Association},
#' 83(403), 744--749.
#'
#' Eriksson, F., Li, J., Scheike, T., and Zhang, M.-J. (2015).
#' The proportional odds cumulative incidence model for competing risks.
#' \emph{Biometrics}, 71(3), 687--695.
#' 
#' Fine, J. P. and Gray, R. J. (1999). A proportional hazards model for the
#' subdistribution of a competing risk. \emph{Journal of the American
#' Statistical Association}, 94(446), 496--509.
#'
#' Haile, S. R., Jeong, J.-H., Chen, X. and Cheng, Y. (2016). A 3-parameter Gompertz
#' distribution for survival data with competing risks, with an application to
#' breast cancer data. \emph{Journal of Applied Statistics}, 43(12), 2239--2253.
#'
#' Jeong, J.-H. and Fine, J. P. (2007). Parametric regression on the
#' cumulative incidence function. \emph{Biostatistics}, 8(2), 184--196.
#'
#' Loynes, R. M. (1969). On Cox and Snell's general definition of residuals.
#' \emph{Journal of the Royal Statistical Society: Series B}, 31(1), 103--106.
#' 
#' @export
pcrr <- function(ftime, fstatus, cov, distribution="gompertz2", dist=NULL, failcode=1, cencode=0,
                 na.action=na.omit, gtol=1e-10, maxiter=300, init, variance=TRUE, sig.level=0.05) {

  if (!is.null(dist)) distribution <- dist
  distribution <- tolower(distribution)
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
  
  if (!is.numeric(sig.level) || length(sig.level) != 1 || !is.finite(sig.level) ||
      sig.level <= 0 || sig.level >= 1) {
    warning("Invalid significance level. Using the default value (sig.level = 0.05).")
    sig.level <- 0.05
  }
  
  call <- match.call()
  
  cov_name <- deparse(substitute(cov))
  cov_vars <- colnames(as.matrix(cov))
  
  if (length(ftime) != length(fstatus) || length(ftime) != NROW(cov)) {
    stop("The lengths of ftime, fstatus, and cov must be equal.")
  }
  
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
    message(format(N_mis), " cases omitted due to missing values")
  }
  
  ftime   <- as.numeric(user_data$ftime)
  fstatus <- user_data$fstatus
  cov <- as.matrix(user_data[, c(-1, -2), drop = FALSE])
  
  
  # Dummy coding
  if (!(cencode  %in% fstatus)) stop("cencode is not in fstatus.")
  if (!(failcode %in% fstatus)) stop("failcode is not in fstatus.")
  if (cencode == failcode) stop("cencode and failcode must be different.")
  
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

  
  # Kernel Operations
  if (distribution == "gompertz2"){
    .init_values <- .init_values_gom2
    .log_lik <- .log_lik_gom2
    .estimate_mle <- .estimate_mle_gom2
    .score_hessian <- .score_hessian_gom2
  } else if (distribution == "gompertz3"){
    .init_values <- .init_values_gom3
    .log_lik <- .log_lik_gom3
    .estimate_mle <- .estimate_mle_gom3
    .score_hessian <- .score_hessian_gom3
  } else if (distribution == "logistic"){
    .init_values <- .init_values_logi
    .log_lik <- .log_lik_logi
    .estimate_mle <- .estimate_mle_logi
    .score_hessian <- .score_hessian_logi
  } 
  
  if (init_ok) theta_init <- init
  else theta_init <- .init_values(x, delta, z)
  
  val_mle <- suppressWarnings(tryCatch(.log_lik(x = x, delta = delta, z = z, theta = theta_init),
                                       error = function(e) NaN))
  if (!is.finite(val_mle) || abs(val_mle) >= 1e+100) {
    warning("The log-likelihood evaluated at the initial values returned NaN. ",
            "Optimization will proceed, but may converge to an incorrect or degenerate solution.")
  }
  
  theta_mle <- .estimate_mle(x, delta, z, theta_init, gtol, maxiter, NULL)
  score_hessian <- .score_hessian(x, delta, z, theta_mle$par, TRUE, NULL)
  
  if(!is.finite(theta_mle$objective) || theta_mle$objective >= 1e+100) {
    theta_mle$convergence <- 1L
    theta_mle$message <- "did not move from the starting values (infeasible region)"
    theta_mle$objective <- NaN
  }
  
  if (theta_mle$convergence != 0) {
    stop("The unconstrained (GOR) fit did not converge",
         if (nzchar(theta_mle$message)) paste0(": ", theta_mle$message), 
         ". Model assumption tests cannot proceed.")
  }
  
  par <- theta_mle$par
  sco <- score_hessian$score
  hess <- score_hessian$hessian
  
  if (is.null(sco)) {
    stop("The Score vector at the unconstrained (GOR) fit could not be computed. The function cannot proceed.")
  }
  if (is.null(hess)) {
    stop("The Hessian matrix at the unconstrained (GOR) fit could not be computed. The function cannot proceed.")
  }
  
  names(par) <- display_names  
  names(sco) <- display_names  
  dimnames(hess) <- list(display_names, display_names)
  
  inv_hess <- tryCatch(
    solve(-hess), error = function(e) {
      stop("The Hessian matrix for the unconstrained (GOR) fit is singular or non-invertible. The function cannot proceed.")
    }
  )

  
  # Model assumption testing:
  # H0: alpha = 0 corresponds to PH; H0: alpha = 1 corresponds to PO.

  fixed_model <- integer(K)
  signif_level <- sig.level
  
  est <- par
  
  d <- diag(inv_hess)
  if (any(!is.finite(d)) || any(d <= 0)) {
    warning("Some variance estimates are non-positive; ",
            "the PH/PO tests may be unreliable.")
  }
  se <- sqrt(d)

  
  if (distribution == "gompertz2"){
    idx_alpha <- (seq_len(K) - 1) * (3 + P) + 1
  } else if (distribution == "gompertz3"){
    idx_alpha <- (seq_len(K) - 1) * (4 + P) + 1
  } else if (distribution == "logistic"){
    idx_alpha <- (seq_len(K) - 1) * (4 + P) + 1
  }
  
  is_alpha <- logical(length(par))
  is_alpha[idx_alpha] <- TRUE
  
  
  z_ph <- (est[is_alpha] - 0) / se[is_alpha]
  z_po <- (est[is_alpha] - 1) / se[is_alpha]
  p_ph <- 2 * (1 - pnorm(abs(z_ph)))
  p_po <- 2 * (1 - pnorm(abs(z_po)))
  
  for (k in 1:K) {
    if (!is.finite(p_ph[k]) && !is.finite(p_po[k])) {
      message("Neither the PH nor the PO hypothesis could be tested for event ", mapping[k], ".")
      fixed_model[k] <- -1
      
    } else if (!is.finite(p_ph[k])) {
      if (p_po[k] > signif_level) {
        message("The PH hypothesis could not be tested, but the PO hypothesis ",
                "was not rejected at the ", signif_level, " significance level for event ", mapping[k], ".")
        fixed_model[k] <- 1
      } else {
        message("The PH hypothesis could not be tested, and the PO hypothesis ",
                "was rejected at the ", signif_level, " significance level for event ", mapping[k], ".")
        fixed_model[k] <- -1
      }
      
    } else if (!is.finite(p_po[k])) {
      if (p_ph[k] > signif_level) {
        message("The PO hypothesis could not be tested, but the PH hypothesis ",
                "was not rejected at the ", signif_level, " significance level for event ", mapping[k], ".")
        fixed_model[k] <- 0
      } else {
        message("The PO hypothesis could not be tested, and the PH hypothesis ",
                "was rejected at the ", signif_level, " significance level for event ", mapping[k], ".")
        fixed_model[k] <- -1
      }
      
    } else if (p_ph[k] <= signif_level && p_po[k] <= signif_level) {
      message("At the ", signif_level, " significance level, both the PH and PO hypotheses are rejected ",
              "for event ", mapping[k], ".")
      fixed_model[k] <- -1
      
    } else if (p_ph[k] <= signif_level && p_po[k] > signif_level) {
      message("At the ", signif_level, " significance level, the PH hypothesis is rejected, but the PO ",
              "hypothesis is not rejected for event ", mapping[k], ".")
      fixed_model[k] <- 1
      
    } else if (p_ph[k] > signif_level && p_po[k] <= signif_level) {
      message("At the ", signif_level, " significance level, the PO hypothesis is rejected, but the PH ",
              "hypothesis is not rejected for event ", mapping[k], ".")
      fixed_model[k] <- 0
      
    } else if (p_ph[k] > signif_level && p_po[k] > signif_level) {
      message("At the ", signif_level, " significance level, neither the PH nor the PO hypothesis is ",
              "rejected for event ", mapping[k], ".")
      fixed_model[k] <- 2
      
    } else {
      warning("An unexpected result occurred during the model assumption tests ",
              "for event ", mapping[k], ".")
      fixed_model[k] <- -1
    }
  }

  
  
  if (any(fixed_model == -1)) {
    stop("Neither the PH nor the PO transformation model is appropriate for the data for event(s) ",
         paste(mapping[which(fixed_model == -1)], collapse = ", "), ". The function cannot proceed.")
  }
  
  # all cases
  case_all <- as.matrix(expand.grid(rep(list(0:1), K)))
  colnames(case_all) <- paste0("event ", mapping)
  for (k in 1:K){
    if (fixed_model[k] == 2) next
    idx <- case_all[, k] == fixed_model[k]
    case_all <- case_all[idx, , drop = FALSE]
  }
  
  case_model <- character(nrow(case_all))
  for (i in 1:nrow(case_all)) {
    case_model[i] <- paste0("[case ", i, "] ")
    for (k in 1:K) {
      case_model[i] <- paste0(case_model[i], if (k > 1) ", " else "","event ", mapping[k], " : ", ifelse(case_all[i, k] == 0, "PH", "PO"))
    }
  }
  

  # fitting MLE for all case
  mle_case_all <- vector("list", nrow(case_all))
  for (i in 1:nrow(case_all)) {
    mle <- .estimate_mle(x, delta, z, theta_init, gtol, maxiter, case_all[i, ])
    if (!is.finite(mle$objective) || mle$objective >= 1e+100) {
      mle$convergence <- 1L
      mle$message <- "did not move from the starting values (infeasible region)"
      mle$objective <- NaN
    }
    mle_case_all[[i]] <- mle
    names(mle_case_all[[i]]$par) <- display_names
  }
  names(mle_case_all) <- case_model
  
  # calculate score, hessian for all case
  display_names2 <- display_names[!is_alpha]
  
  sco_case_all <- vector("list", nrow(case_all))    
  if (variance) {    
    hess_case_all <- vector("list", nrow(case_all))  
  } else {    
    hess_case_all <- NULL  
  }
 
  for (i in 1:nrow(case_all)) {
    sco_hess <- .score_hessian(x, delta, z, mle_case_all[[i]]$par, variance, case_all[i, ])
    sco_case_all[[i]] <- sco_hess$score
    names(sco_case_all[[i]]) <- display_names2
    if (variance) {    
      hess_case_all[[i]] <- sco_hess$hessian    
      dimnames(hess_case_all[[i]]) <- list(display_names2, display_names2)    
    }
  }
  names(sco_case_all) <- case_model
  if (variance) {    
    names(hess_case_all) <- case_model  
  }
  

  
  # Define Class 'pcrr'
  cls <- list(coef      = par,
              loglik    = -theta_mle$objective,
              init      = theta_init,
              score     = sco,
              inf       = -hess,
              invinf    = inv_hess,
              converged = theta_mle$convergence == 0,
              iter      = theta_mle$iterations,
              message   = theta_mle$message,
              call      = call,
              x         = x,
              event     = fstatus,
              z         = z,
              cencode   = cencode,
              n         = N,
              n_missing = N_mis,
              k = K,
              p = P,
              distribution = distribution,
              cov_names = cov_names,
              mapping = mapping,
              maxtime = max(x),
              case_all = case_all,
              case_model = case_model,
              mle_case_all = mle_case_all,
              sco_case_all = sco_case_all,
              hess_case_all = hess_case_all
  )
  class(cls) <- "pcrr"
  cls
}


#' Print a Fitted Parametric Competing Risks Regression Model
#'
#' Prints the convergence status and model assumption tests for a fitted
#' parametric competing risks regression model. The model assumption tests
#' are based on the estimated shape parameter \eqn{\alpha} from the
#' generalized odds-rate (GOR) model, in which \eqn{\alpha} is estimated
#' without imposing a fixed value.
#'
#' For each event, the function tests the proportional hazards (PH) and
#' proportional odds (PO) model assumptions based on the estimated
#' \eqn{\alpha}. The PH model corresponds to \eqn{H_0: \alpha = 0}, while
#' the PO model corresponds to \eqn{H_0: \alpha = 1}. For each hypothesis,
#' the estimated value of \eqn{\alpha}, its standard error, Wald z-statistic,
#' and two-sided p-value are displayed.
#'
#' The function also displays all possible model assumption cases considered
#' in the fitting procedure. Each case specifies whether the PH or PO model
#' is assumed for each event.
#'
#'
#' @param x an object of class \code{"pcrr"}, representing a fitted
#'   parametric competing risks regression model.
#' @param digits the number of significant digits to use when printing the
#'   results. Defaults to \code{max(options()$digits - 4, 3)}.
#' @param ... further arguments passed to the printing methods.
#'
#' @return
#' The input object \code{x}, returned invisibly.
#'
#' @seealso
#' \code{\link{pcrr}},
#' \code{\link{summary.pcrr}}
#'
#' @export
print.pcrr <- function(x, digits = max(options()$digits - 4, 3), ...) {
  P <- x$p
  K <- x$k
  
  
  if (x$distribution == "gompertz2") {
    block <- 3 + P
  } else {
    block <- 4 + P
  }
  idx_alpha <- (seq_len(K) - 1) * block + 1
  
  
  
  cat("convergence : ", x$converged, "  (iteration : ", x$iter, ")\n", sep = "")
  
  if (!isTRUE(x$converged)) {
    if (!is.null(x$message) && nzchar(x$message)) {
      cat("[", x$message, "]\n", sep = "")
    }
    
    return(invisible(x))
  }
  
  
  
  est <- x$coef
  se <- sqrt(diag(x$invinf))
  

  
  # Model assumption tests

  cat("\n")
  cat("========================================\n")
  cat("Model assumption tests\n")
  cat("========================================\n\n")

  
  alpha_est <- est[idx_alpha]
  alpha_se  <- se[idx_alpha]
  
  z_ph <- (alpha_est - 0) / alpha_se
  p_ph <- 2 * (1 - pnorm(abs(z_ph)))
  
  z_po <- (alpha_est - 1) / alpha_se
  p_po <- 2 * (1 - pnorm(abs(z_po)))
  
  link_ph <- data.frame(est = alpha_est, se = alpha_se, `z value` = z_ph, `Pr(>|z|)` = p_ph, check.names = FALSE)
  
  link_po <- data.frame(est = alpha_est, se = alpha_se, `z value` = z_po, `Pr(>|z|)` = p_po, check.names = FALSE)
  
  
  cat("[Proportional Hazards]\n")
  cat("H0 : alpha = 0\n\n")
  
  printCoefmat(link_ph, digits = digits, signif.stars = FALSE, has.Pvalue = TRUE,
               P.values = TRUE, cs.ind = 1:2, tst.ind = 3)
  
  cat("\n----------------------------------------\n\n")

  

  cat("[Proportional Odds]\n")
  cat("H0 : alpha = 1\n\n")
  
  printCoefmat(link_po, digits = digits, signif.stars = FALSE, has.Pvalue = TRUE,
               P.values = TRUE, cs.ind = 1:2, tst.ind = 3)
  
  cat("\n========================================\n")

  
  cat("\nThe possible model cases are as follows :\n")
  for (i in seq_along(x$case_model)) {
    cat(x$case_model[i], "\n")
  }
    
  
  
  invisible(x)
}


#' Diagnostic Plots for a Parametric Competing Risks Regression Model
#'
#' @description
#' Draws diagnostic plots for selected model cases and events of a fitted
#' \code{"pcrr"} object. Two types of diagnostics are available:
#' a baseline comparison and a Cox--Snell residual plot.
#' Either plot can be omitted using \code{baseline} or \code{coxsnell}.
#' The underlying curves are returned invisibly.
#'
#' @details
#' \strong{Baseline diagnostic.}
#' Compares the fitted parametric baseline function with a semiparametric
#' estimate for each selected event. Both are evaluated at zero covariates.
#' The semiparametric reference is obtained using \code{\link[cmprsk]{crr}}
#' for the proportional hazards (PH) link and
#' \code{\link[timereg]{prop.odds.subdist}} for the proportional odds (PO)
#' link. Both reference estimates are computed, but only the one corresponding
#' to the link specified by each model case is displayed.
#'
#' The comparison helps assess whether the chosen baseline distribution
#' adequately describes the data. Under the PO link, the baseline function
#' represents the baseline odds rather than the baseline cumulative hazard.
#'
#' \strong{Cox--Snell residual diagnostic.}
#' For cause \eqn{k}, the Cox--Snell residual of subject \eqn{i} is calculated
#' by evaluating the fitted cumulative subdistribution hazard at the observed
#' follow-up time:
#' \deqn{
#'   r_{ki} = \widehat{\Lambda}_{ki}(X_i),
#'   \qquad X_i = \min(T_i, C_i).
#' }
#' Residuals are calculated for all subjects, including those who experience
#' the event of interest, a competing event, or censoring.
#'
#' The diagnostic plot compares the weighted Nelson--Aalen estimate of the
#' cumulative hazard of the residuals with the 45-degree line. At each
#' distinct residual value \eqn{u_j} corresponding to an event of interest,
#' the estimate is updated using the number of events and the weighted
#' risk-set size.
#'
#' Subjects contribute to the risk set while their residuals have not yet
#' been exceeded. Subjects experiencing a competing event can remain in the
#' risk set beyond their observed residual, with an inverse-probability-of-
#' censoring weight
#' \deqn{
#'   q_{ki}(u) =
#'   \frac{\widehat{G}\{\widehat{\Lambda}_{ki}^{-1}(u)\}}
#'        {\widehat{G}(X_i)}.
#' }
#' Here, \eqn{\widehat{G}(t)} is the Kaplan--Meier estimate of the censoring
#' survival function, evaluated at the left limit and estimated without
#' covariates.
#'
#' To avoid extrapolating the censoring survival function beyond the observed
#' follow-up, the contribution of each subject is restricted to residual
#' values satisfying
#' \eqn{u \le L_{ki}}, where
#' \eqn{L_{ki} = \widehat{\Lambda}_{ki}(X_{(\max)})} and
#' \eqn{X_{(\max)}} is the largest observed follow-up time.
#'
#' Under a correctly specified model, the Cox--Snell residuals should behave
#' as a sample from the unit exponential distribution, accounting for
#' competing events and right censoring. Consequently, the estimated
#' cumulative hazard should be close to the 45-degree line.
#' The final jump is omitted from the plot, unless it is the only jump,
#' because it is based on the smallest risk set and can be unstable.
#' It is nevertheless retained in the returned results.
#'
#' No goodness-of-fit p-value is provided.
#'
#' \strong{Reproducibility.}
#' The \code{prop.odds.subdist()} function uses random jitter to handle
#' tied event times. A fixed seed is therefore used for each fit, and the
#' random number state is restored after the function finishes.
#'
#' \strong{Plot layout.}
#' For each selected model case, the baseline diagnostic is drawn first,
#' followed by the Cox--Snell diagnostic. Each screen contains one panel
#' per selected event. If the selected events do not fit on one page,
#' additional pages are used. The case label and the quantity being plotted
#' are displayed in the header of each page.
#'
#' \strong{Warnings and missing plots.}
#' A model case that did not converge is still plotted, but a warning is
#' issued and its diagnostics should be interpreted with caution.
#' If a semiparametric fit fails or returns non-finite values, a warning
#' is issued and the corresponding baseline panel displays only the
#' parametric curve. A non-converged \code{crr()} fit is retained, with
#' a warning that its estimate may be unreliable.
#' If no usable events are available for a Cox--Snell panel, a note is
#' displayed instead of the plot.
#' Setting \code{baseline = FALSE} skips the semiparametric fits.
#'
#' @param x an object of class \code{"pcrr"}. The original data must be
#'   stored in the object as \code{x}, \code{event}, and \code{z}.
#' @param case integer vector specifying the model cases to plot, as numbered
#'   in \code{\link{print.pcrr}}. By default, all available cases are selected.
#'   Duplicates are removed, and plots follow the order supplied.
#' @param event numeric vector specifying the event codes to plot, as stored
#'   in \code{x$mapping}. By default, all events are selected. Duplicates
#'   are removed, and plots follow the order supplied.
#' @param baseline,coxsnell logical. Whether to draw the baseline diagnostic
#'   and Cox--Snell residual diagnostic, respectively. Both default to
#'   \code{TRUE}. At least one must be enabled. A disabled diagnostic is
#'   not computed, and its corresponding returned component is \code{NULL}.
#' @param col.par colour of the parametric baseline curve and the
#'   Nelson--Aalen estimate in the Cox--Snell plot.
#' @param col.semi colour of the semiparametric baseline curve.
#' @param lty.par,lty.semi line types of the parametric and semiparametric
#'   baseline curves, respectively. The Nelson--Aalen curve is always solid.
#' @param ref.col,ref.lty colour and line type of the 45-degree reference line.
#' @param ylim y-axis limits for all baseline panels. If \code{NULL},
#'   the limits are determined separately for each panel unless
#'   \code{common.ylim = TRUE}.
#' @param lim.coxsnell limits applied to both axes of every Cox--Snell panel.
#'   If \code{NULL}, the limits are determined from the plotted values.
#' @param xmin,xmax range of the time axis for the baseline plots.
#'   \code{xmin} must be a single finite non-negative number.
#'   If \code{xmax} is \code{NULL}, the axis extends to the largest observed
#'   follow-up time, while each parametric curve ends at the largest observed
#'   event time for its corresponding cause. Otherwise, \code{xmax} must be
#'   a finite number greater than \code{xmin}, and the parametric curves
#'   are evaluated up to that time.
#' @param xlab,ylab axis labels for the baseline plots. The \code{ylab}
#'   is also used in the page header. Under the PO link, the plotted
#'   quantity is baseline odds rather than cumulative hazard.
#' @param xlab.coxsnell,ylab.coxsnell axis labels for the Cox--Snell plots.
#'   The \code{xlab.coxsnell} is also used in the page header.
#' @param legend logical. Whether to display a legend in each panel.
#' @param legend.pos legend positions. The first value specifies the position
#'   in baseline panels and the second in Cox--Snell panels. A single value
#'   is used for both.
#' @param legend.title optional title for the legends. Defaults to \code{NULL}.
#' @param lwd line width of the curves. The 45-degree reference line always
#'   has width 1.
#' @param main optional title or titles replacing the model case labels
#'   in the page headers. Values are recycled over the selected cases.
#' @param n.grid number of evaluation points for the parametric baseline
#'   curves. Must be a single integer greater than or equal to 2.
#'   Defaults to 300.
#' @param mfrow panel layout for each screen, specified as
#'   \code{c(nrow, ncol)}. If \code{NULL}, the current \code{par("mfrow")}
#'   setting is used when it has been modified; otherwise, a layout is
#'   determined automatically.
#' @param common.ylim logical. Whether to use common axis limits across
#'   all panels within each screen. For the Cox--Snell screen, the same
#'   limits are applied to both axes. Ignored when \code{ylim} or
#'   \code{lim.coxsnell} is specified.
#' @param ask logical. Whether to prompt before drawing each new page.
#'   If \code{NULL}, prompting is enabled when multiple pages are drawn
#'   on an interactive device.
#' @param ... additional graphical parameters passed to
#'   \code{\link[graphics]{plot.default}} when setting up each panel.
#'
#' @return
#' Invisibly returns a list containing the following components.
#' \itemize{
#'   \item \code{semi_base_cum_haz_list}: Semiparametric baseline estimates
#'   for each selected event. Each element contains \code{PH}, the
#'   Fine--Gray baseline cumulative hazard estimate, and \code{PO},
#'   the proportional odds baseline estimate. A component is \code{NULL}
#'   if the corresponding fit fails or returns non-finite values.
#'
#'   \item \code{para_base_cum_haz_list}: Parametric baseline estimates
#'   for each selected model case and event. Each matrix contains
#'   \code{time} and \code{base_cum}.
#'
#'   \item \code{cox_snell_residual_list}: Cox--Snell residual diagnostic
#'   results for each selected model case and event. Each matrix contains
#'   \code{residual}, the distinct residual values at event times, and
#'   \code{NA_estimated_cumhaz}, the corresponding weighted Nelson--Aalen
#'   estimates. The final row is retained even when it is omitted from
#'   the plot.
#' }
#'
#' The components corresponding to disabled diagnostics are \code{NULL}.
#'
#' @seealso
#' \code{\link{pcrr}},
#' \code{\link{print.pcrr}},
#' \code{\link{summary.pcrr}},
#' \code{\link{predict.pcrr}},
#' \code{\link[cmprsk]{crr}},
#' \code{\link[timereg]{prop.odds.subdist}}
#'
#' @export
plot.pcrr <- function(x, case = NULL, event = NULL,
                      baseline = TRUE, coxsnell = TRUE,
                      col.par = "black", col.semi = "blue",
                      lty.par = 1, lty.semi = 2,
                      ref.col = "red", ref.lty = 2,
                      ylim = NULL, lim.coxsnell = NULL,
                      xmin = 0, xmax = NULL,
                      xlab = "Time",
                      ylab = "Baseline Cumulative Hazard",
                      xlab.coxsnell = "Cox-Snell Residual",
                      ylab.coxsnell = "Estimated Cumulative Hazard",
                      legend = TRUE,
                      legend.pos = c("bottomright", "topleft"),
                      legend.title = NULL, lwd = 2, main = NULL,
                      n.grid = 300, mfrow = NULL, common.ylim = FALSE,
                      ask = NULL, ...) {
  
  if (!is.logical(baseline) || length(baseline) != 1L || is.na(baseline))
    stop("`baseline` must be TRUE or FALSE.", call. = FALSE)
  
  if (!is.logical(coxsnell) || length(coxsnell) != 1L || is.na(coxsnell))
    stop("`coxsnell` must be TRUE or FALSE.", call. = FALSE)
  
  if (!baseline && !coxsnell)
    stop("At least one of `baseline` and `coxsnell` must be TRUE.", call. = FALSE)
  
  if (is.null(x$x) || is.null(x$event) || is.null(x$z))
    stop("The original data are not stored in the fitted object, so the ",
         "diagnostics cannot be computed.", call. = FALSE)
  
  if (!is.numeric(xmin) || length(xmin) != 1L || !is.finite(xmin) || xmin < 0) 
    stop("`xmin` must be a single finite non-negative number.", call. = FALSE)
  
  if (!is.null(xmax)) {
    if (!is.numeric(xmax) || length(xmax) != 1L || !is.finite(xmax)) 
      stop("`xmax` must be NULL or a single finite number.", call. = FALSE)
    
    if (xmax <= xmin) stop("`xmax` must be greater than `xmin`.", call. = FALSE)
  }
  
  if (!is.numeric(n.grid) || length(n.grid) != 1L ||
      !is.finite(n.grid) || n.grid < 2 || n.grid != floor(n.grid)) 
    stop("`n.grid` must be a single integer greater than or equal to 2.", call. = FALSE)
  n.grid <- as.integer(n.grid)
  
  if (is.null(xmax) && xmin >= max(x$x))
    stop("`xmin` must be smaller than the end of follow-up (", format(max(x$x)), ").", call. = FALSE)
  
  
  case_all <- x$case_all
  n_case_all <- nrow(case_all)
  
  avail_case <- seq_len(n_case_all)
  
  if (is.null(case)) {
    case <- avail_case
  } else {
    if (!is.numeric(case) || length(case) == 0 || any(!is.finite(case)) ||
        any(case != floor(case)) || !all(case %in% avail_case))
      stop("`case` must be one or more of: ",
           paste(avail_case, collapse = ", "),
           "\n\nThe possible model cases are as follows :\n",
           paste(x$case_model, collapse = "\n"), call. = FALSE)
    case <- as.integer(unique(case))
  }
  
  n_case <- length(case)
  
  
  mapping <- x$mapping
  if (is.null(event)) {
    event <- mapping
  } else {
    if (!is.numeric(event) || length(event) == 0 || !all(event %in% mapping))
      stop("`event` must be one or more of: ", paste(mapping, collapse = ", "), call. = FALSE)
    event <- unique(event)
  }
  
  
  bad <- case[vapply(case, function(j) { x$mle_case_all[[j]]$convergence != 0 }, logical(1))]
  
  if (length(bad) > 0) {
    warning("Case(s) ", paste(bad, collapse = ", "),
            " did not converge; their diagnostics should not be interpreted.",
            call. = FALSE)
  }
  
  
  model <- x$case_model
  
  
  distribution <- x$distribution
  cencode <- x$cencode
  
  K <- x$k
  P <- x$p
  ftime <- x$x
  fevent <- x$event
  cov <- x$z
  
  
  # re indexing
  fevent_re <- ifelse(fevent == cencode, 0, match(fevent, mapping)) 
  data_re <- data.frame(ftime, fevent_re)
  
  
  tol <- 1e-12
  # 0. define function
  if (distribution == "gompertz2" || distribution == "gompertz3"){
    u_k <- function(t, rho, tau, eta, tol) {
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
    
    res_PH <- function(t, z, rho, tau, eta, beta, tol) {
      ezb <- exp(z %*% beta)
      u <- u_k(t, rho, tau, eta, tol)
      
      ezb * u
    }
    res_PO <- function(t, z, rho, tau, eta, beta, tol) {
      ezb <- exp(z %*% beta)
      u <- u_k(t, rho, tau, eta, tol)
      
      log1p(ezb * u)
    }
    
    inv_res_PH <- function(u, z, rho, tau, eta, beta, tol) {
      ezb <- exp(z %*% beta)
      
      if (abs(rho) < tol && abs(eta) < tol){
        
        return(u / (tau * ezb))
        
      } else if (abs(rho) < tol){
        
        return(u / (tau * exp(eta) * ezb))
        
      } else if (abs(eta) < tol){
        inner <- rho * u / (tau * ezb)
        result <- rep(Inf, length(inner))
        valid <- inner > -1
        
        result[valid] <- log1p(inner[valid]) / rho
        return(result)
        
      } else {
        inner <- exp(eta) + (rho * eta / tau) * u / ezb
        result <- rep(Inf, length(inner))
        valid <- inner > 0 & log(inner) / eta > 0
        
        result[valid] <- log(log(inner[valid]) / eta) / rho
        return(result)
      }
      
    }
    inv_res_PO <- function(u, z, rho, tau, eta, beta, tol) {
      ezb <- exp(z %*% beta)
      
      if (abs(rho) < tol && abs(eta) < tol){
        
        return(expm1(u) / (tau * ezb))
        
      } else if (abs(rho) < tol){
        
        return(expm1(u) / (tau * exp(eta) * ezb))
        
      } else if (abs(eta) < tol){
        inner <- rho * expm1(u) / (tau * ezb)
        result <- rep(Inf, length(inner))
        valid <- inner > -1
        
        result[valid] <- log1p(inner[valid]) / rho
        return(result)
        
      } else {
        inner <- exp(eta) + (rho * eta / tau) * expm1(u) / ezb
        result <- rep(Inf, length(inner))
        valid <- inner > 0 & log(inner) / eta > 0
        
        result[valid] <- log(log(inner[valid]) / eta) / rho
        return(result)
      }
    }
    
  } else if (distribution == "logistic"){
    u_k <- function(t, b, c, p, tol){
      inner <- - p + p*(1+exp(-b*c))/(1+exp(b*(t-c)))
      # u <- rep(NaN, length(t))
      # valid <- 1 + inner > 0
      # u[valid] <- -log1p(inner[valid])
      u <- -log1p(inner)
      u
    }
    
    res_PH <- function(t, z, b, c, p, beta, tol) {
      ezb <- exp(z %*% beta)
      u <- u_k(t, b, c, p, tol)
      
      ezb * u
    }
    res_PO <- function(t, z, b, c, p, beta, tol) {
      ezb <- exp(z %*% beta)
      u <- u_k(t, b, c, p, tol)
      
      log1p(ezb * u)
    }
    
    inv_res_PH <- function(u, z, b, c, p, beta, tol) {
      ezb <- exp(z %*% beta)
      
      numerator <- p * (1 + exp(- b * c))
      denominator <- p - 1 + exp(- u / ezb)
      
      inner <- numerator / denominator
      
      result <- rep(Inf, length(inner))
      
      valid <- abs(denominator) > tol & inner > 1
      
      result[valid] <- c + log(inner[valid] - 1) / b
      return(result)
    }
    inv_res_PO <- function(u, z, b, c, p, beta, tol) {
      ezb <- exp(z %*% beta)
      
      numerator <- p * (1 + exp(- b * c))
      denominator <- p - 1 + exp(- expm1(u) / ezb)
      
      inner <- numerator / denominator
      
      result <- rep(Inf, length(inner))
      
      valid <- abs(denominator) > tol & inner > 1
      
      result[valid] <- c + log(inner[valid] - 1) / b
      return(result)
    }
    
  }
  
  
  
  
  semi_base_cum_haz_list <- NULL
  para_base_cum_haz_list <- NULL
  if (baseline){
    # 1. parametric vs. semi-parametric
    
    
    
    had_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
    
    if (had_seed) old_seed <- get(".Random.seed", envir = globalenv())
    
    on.exit({
      if (had_seed) {
        assign(".Random.seed", old_seed, envir = globalenv())
      } else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
        rm(".Random.seed", envir = globalenv())
      }
    }, add = TRUE)
    
    
    # i) semi-parametric
    semi_base_cum_haz_list <- vector("list", length(event))
    names(semi_base_cum_haz_list) <- paste0("event ", event)
    
    for (e in 1:length(event)) {
      semi_base_cum_haz_list[[e]] <- list(PH = NULL, PO = NULL)
      
      fit_PH <- tryCatch(
        crr(ftime = ftime, fstatus = fevent, cov1 = cov, failcode = event[e], cencode = x$cencode, variance = FALSE),
        error = function(err) {
          warning("crr() failed for event ", event[e], ": ", conditionMessage(err), call. = FALSE)
          NULL
        })

      if (!is.null(fit_PH)) {
        if (!isTRUE(fit_PH$converged))
          warning("crr() did not converge for event ", event[e], "; its curve may be unreliable.", call. = FALSE)
        if (all(is.finite(fit_PH$bfitj))) {
          semi_base_cum_haz_list[[e]][[1]] <- cbind(time = fit_PH$uftime, base_cum_haz = cumsum(fit_PH$bfitj))
        } else {
          warning("crr() returned non-finite values for event ", event[e], "; its curve is omitted.", call. = FALSE)
        }
      }
 
      set.seed(1)
      fit_PO <- tryCatch(
        prop.odds.subdist(Event(ftime, fevent_re) ~ cov, data = data_re, cause = match(event[e], mapping), cens.model = "KM", n.sim = 0),
        error = function(err) {
          warning("prop.odds.subdist() failed for event ", event[e], ": ", conditionMessage(err), call. = FALSE)
          NULL
        })
      
      
      if (!is.null(fit_PO)) {
        if (all(is.finite(fit_PO$cum[, 1:2]))) {
          semi_base_cum_haz_list[[e]][[2]] <- cbind(time = fit_PO$cum[, 1], base_cum_odds = fit_PO$cum[, 2])
        } else {
          warning("prop.odds.subdist() returned non-finite values for event ", event[e], "; its curve is omitted.", call. = FALSE)
        }
      }
    }
    
    # ii) parametric
    para_base_cum_haz_list <- vector("list", n_case)
    names(para_base_cum_haz_list) <- model[case]
    
    for (l in 1:n_case) {
      j <- case[l]
      coef <- x$mle_case_all[[j]]$par
      
      para_base_cum_haz_list[[l]] <- vector("list", length(event))
      names(para_base_cum_haz_list[[l]]) <- paste0("event ", event)
      
      
      for (e in 1:length(event)) {
        k <- match(event[e], mapping)
        
        if (is.null(xmax)) max_time <- max(ftime[fevent == event[e]])
        else max_time <- xmax
        times <- seq(xmin, max_time, length.out = n.grid)
        
        if (distribution == "gompertz2"){
          alpha <- coef[(k - 1) * (3 + P) + 1]
          rho   <- coef[(k - 1) * (3 + P) + 2]
          tau   <- coef[(k - 1) * (3 + P) + 3]
          eta   <- 0
          beta  <- coef[((k - 1) * (3 + P) + 4):(k * (3 + P))]
          
          uk <- u_k(times, rho, tau, eta, tol)
        } else if (distribution == "gompertz3"){
          alpha <- coef[(k - 1) * (4 + P) + 1]
          rho   <- coef[(k - 1) * (4 + P) + 2]
          tau   <- coef[(k - 1) * (4 + P) + 3]
          eta   <- coef[(k - 1) * (4 + P) + 4]
          beta  <- coef[((k - 1) * (4 + P) + 5):(k * (4 + P))]
          
          uk <- u_k(times, rho, tau, eta, tol)
        } else if (distribution == "logistic"){
          alpha <- coef[(k - 1) * (4 + P) + 1]
          b   <- coef[(k - 1) * (4 + P) + 2]
          c   <- coef[(k - 1) * (4 + P) + 3]
          p   <- coef[(k - 1) * (4 + P) + 4]
          beta  <- coef[((k - 1) * (4 + P) + 5):(k * (4 + P))]
          
          uk <- u_k(times, b, c, p, tol)
        }
        
        para_base_cum_haz_list[[l]][[e]] <- cbind(time = times, base_cum = uk)
        
      }
    }
  }
  
  
  
  
  cox_snell_residual_list <- NULL
  
  if (coxsnell) {
    # 2. cox-snell residual
    cox_snell_residual_list <- vector("list", n_case)
    names(cox_snell_residual_list) <- model[case]
    
    
    # KM estimate
    fit_SC <- survfit(Surv(ftime, fevent == cencode) ~ 1)
    summ <- summary(fit_SC)
    SC_KM <- cbind(time = summ$time, surv = summ$surv)
    
    G_KM <- function(t) {
      c(1, SC_KM[, 2])[findInterval(t, SC_KM[, 1], left.open = TRUE) + 1L]
    }
    
    
    for (l in 1:n_case) {
      j <- case[l]
      coef <- x$mle_case_all[[j]]$par
      
      cox_snell_residual_list[[l]] <- vector("list", length(event))
      names(cox_snell_residual_list[[l]]) <- paste0("event ", event)
      
      
      for (e in 1:length(event)) {
        k <- match(event[e], mapping)
        
        if (distribution == "gompertz2"){
          alpha <- coef[(k - 1) * (3 + P) + 1]
          rho   <- coef[(k - 1) * (3 + P) + 2]
          tau   <- coef[(k - 1) * (3 + P) + 3]
          eta   <- 0
          beta  <- coef[((k - 1) * (3 + P) + 4):(k * (3 + P))]
          
        } else if (distribution == "gompertz3"){
          alpha <- coef[(k - 1) * (4 + P) + 1]
          rho   <- coef[(k - 1) * (4 + P) + 2]
          tau   <- coef[(k - 1) * (4 + P) + 3]
          eta   <- coef[(k - 1) * (4 + P) + 4]
          beta  <- coef[((k - 1) * (4 + P) + 5):(k * (4 + P))]
          
        } else if (distribution == "logistic"){
          alpha <- coef[(k - 1) * (4 + P) + 1]
          b   <- coef[(k - 1) * (4 + P) + 2]
          c   <- coef[(k - 1) * (4 + P) + 3]
          p   <- coef[(k - 1) * (4 + P) + 4]
          beta  <- coef[((k - 1) * (4 + P) + 5):(k * (4 + P))]
          
        }
        
        
        # Residual
        t_max <- max(ftime)
        
        if (case_all[j,k] == 0) {  # PH
          if (distribution != "logistic") {
            residual <- res_PH(ftime, cov, rho, tau, eta, beta, tol)
            L <- res_PH(t_max, cov, rho, tau, eta, beta, tol)
          } else {
            residual <- res_PH(ftime, cov, b, c, p, beta, tol)
            L <- res_PH(t_max, cov, b, c, p, beta, tol)
          }
          
          
        } else if (case_all[j,k] == 1) {  # PO
          if (distribution != "logistic") {
            residual <- res_PO(ftime, cov, rho, tau, eta, beta, tol)
            L <- res_PO(t_max, cov, rho, tau, eta, beta, tol)
          } else {
            residual <- res_PO(ftime, cov, b, c, p, beta, tol)
            L <- res_PO(t_max, cov, b, c, p, beta, tol)
          }
          
        }
        
        
        # NA estimate
        u_j <- sort(unique(residual[fevent == event[e]]))
        
        
        d_R <- numeric(length(u_j))
        i <- 1
        for (u in u_j) {
          d_R[i] <- sum( (residual == u) & (fevent == event[e]) )
          i <- i + 1
        }
        
        
        Y_R <- numeric(length(u_j))
        i <- 1
        for (u in u_j) {
          # i)  not weight
          Y_R[i] <- Y_R[i] + sum( (u <= residual) )
          
          
          # ii) weight
          weight_subject <- (residual < u) & (u <= L) & (fevent != event[e]) & (fevent != cencode)
          if (!any(weight_subject)) {
            i <- i + 1
            next
          }
          
          # G(inv(u)) / G(X_i)
          if (case_all[j,k] == 0) {  # PH
            if (distribution != "logistic") {
              X_tilde <- inv_res_PH(u, cov[weight_subject,,drop = FALSE], rho, tau, eta, beta, tol)
            } else {
              X_tilde <- inv_res_PH(u, cov[weight_subject,,drop = FALSE], b, c, p, beta, tol)
            }
          } else if (case_all[j,k] == 1) {  # PO
            if (distribution != "logistic") {
              X_tilde <- inv_res_PO(u, cov[weight_subject,,drop = FALSE], rho, tau, eta, beta, tol)
            } else {
              X_tilde <- inv_res_PO(u, cov[weight_subject,,drop = FALSE], b, c, p, beta, tol)
            }
          }
          
          X <- ftime[weight_subject]
          
          q <- G_KM(X_tilde) / G_KM(X)
          
          Y_R[i] <- Y_R[i] + sum(q)
          
          i <- i + 1
        }
        
        delta_H_R <- d_R / Y_R
        
        H_R <- cumsum(delta_H_R)
        
        cox_snell_residual_list[[l]][[e]] <- cbind(residual = u_j, NA_estimated_cumhaz = H_R)
        
      }
    }
  }
  
  
  ## ------------------------------------------------------------------
  ## d1. arguments that only the pages use
  ## ------------------------------------------------------------------
  ## one position per page : a rising baseline leaves the lower right free,
  ## a Cox-Snell curve that drops below y = x leaves the upper left free
  legend.pos <- rep(legend.pos, length.out = 2)
  
  ## time axis of the baseline page
  t_end <- max(ftime)                       # end of follow-up
  if (!is.numeric(xmin) || length(xmin) != 1 || !is.finite(xmin) || xmin < 0)
    stop("`xmin` must be a single non-negative number.", call. = FALSE)
  if (is.null(xmax)) xmax <- t_end
  if (!is.numeric(xmax) || length(xmax) != 1 || !is.finite(xmax) || xmax <= xmin)
    stop("`xmax` must be a single number larger than `xmin`.", call. = FALSE)
  
  n_event <- length(event)
  k_of    <- match(event, mapping)          # internal cause index of each event
  
  ## the link a case fixes for cause k : alpha = 0 is PH, alpha = 1 is PO
  link_of <- function(i_case, k) if (case_all[i_case, k] == 0) "PH" else "PO"
  
  ## ------------------------------------------------------------------
  ## d2. panel grid : the selected events on one page
  ## ------------------------------------------------------------------
  grid_dim <- function(n) {
    nr <- max(1L, as.integer(floor(sqrt(n))))
    c(nr, as.integer(ceiling(n / nr)))
  }
  
  user_split <- !identical(as.integer(par("mfrow")), c(1L, 1L))
  if (is.null(mfrow)) {
    auto_layout <- !user_split
    dims <- if (user_split) as.integer(par("mfrow")) else grid_dim(n_event)
  } else {
    if (!is.numeric(mfrow) || length(mfrow) != 2 || any(!is.finite(mfrow)) ||
        any(mfrow < 1))
      stop("`mfrow` must be a numeric vector of length 2, e.g. c(1, 2).",
           call. = FALSE)
    auto_layout <- TRUE
    dims <- as.integer(mfrow)
  }
  per_page <- dims[1] * dims[2]
  
  ## ------------------------------------------------------------------
  ## d3. device set-up
  ## ------------------------------------------------------------------
  n_quant <- baseline + coxsnell
  n_page  <- n_case * n_quant * ceiling(n_event / per_page)
  
  op <- par(no.readonly = TRUE)
  on.exit(par(op), add = TRUE)
  
  if (is.null(ask)) ask <- (n_page > 1) && dev.interactive()
  if (isTRUE(ask)) {
    oask <- devAskNewPage(TRUE)
    on.exit(devAskNewPage(oask), add = TRUE)
  }
  
  ## ------------------------------------------------------------------
  ## d4. small helpers
  ## ------------------------------------------------------------------
  ## largest finite value, or NA when there is none
  top_of <- function(v) {
    v <- v[is.finite(v)]
    if (length(v) == 0 || max(v) <= 0) NA_real_ else max(v)
  }
  
  ## the semiparametric estimate of one event under one link, drawn as a
  ## step function running from 0 to the end of follow-up. NULL when the
  ## list holds no curve there.
  semi_curve <- function(e, lnk) {
    m <- semi_base_cum_haz_list[[e]][[lnk]]
    if (is.null(m) || NROW(m) == 0) return(NULL)
    tt <- m[, 1]; uu <- m[, 2]
    ok <- is.finite(tt) & is.finite(uu)
    tt <- tt[ok]; uu <- uu[ok]
    if (length(tt) == 0) return(NULL)
    if (tt[1] > 0) { tt <- c(0, tt); uu <- c(0, uu) }
    n <- length(tt)
    if (t_end > tt[n]) { tt <- c(tt, t_end); uu <- c(uu, uu[n]) }
    list(time = tt, u = uu)
  }
  
  base_top <- function(pm, sp) {
    top_of(c(pm[pm[, 1] <= xmax, 2], if (!is.null(sp)) sp$u[sp$time <= xmax]))
  }
  
  
  # cs_top <- function(g) {
  #   if (is.null(g) || NROW(g) == 0) return(NA_real_)    
  #   top_of(c(g[, 1], g[, 2]))  
  # }
  ## the last jump is not drawn (see draw_cs_panel), so it does not set the
  ## axis range either
  cs_top <- function(g) {
    if (is.null(g) || NROW(g) == 0) return(NA_real_)
    if (NROW(g) > 1) g <- g[-NROW(g), , drop = FALSE]
    top_of(c(g[, 1], g[, 2]))
  }
  
  ## A legend is laid out once, at the device size of the moment, and its
  ## position is stored in data coordinates. A screen device redraws at
  ## another size (RStudio on every pane resize, Zoom and Export) by
  ## replaying that layout instead of redoing it, so a legend laid out in a
  ## narrow pane lands left of the plot region and is clipped away. On a
  ## screen device the layout code itself is therefore recorded with
  ## recordGraphics(), so every redraw lays the legend out afresh, shrinking
  ## it when the panel is narrower than its text. Off-screen devices (png,
  ## pdf, and the one gridGraphics draws on for patchwork, cowplot and
  ## ggplotify) are never resized, and gridGraphics cannot read a
  ## recordGraphics() entry, so there the legend is drawn exactly as before.
  draw_legend <- function(pos, labels, cols, ltys, lwds) {
    if (!dev.interactive()) {
      legend(pos, legend = labels, col = cols, lty = ltys, lwd = lwds,
             title = legend.title, bty = "n")
      return(invisible(NULL))
    }
    recordGraphics({
      w   <- legend(pos, legend = labels, col = cols, lty = ltys, lwd = lwds,
                    title = ttl, bty = "n", plot = FALSE)$rect$w
      usr <- par("usr")
      legend(pos, legend = labels, col = cols, lty = ltys, lwd = lwds,
             title = ttl, bty = "n",
             cex = min(1, (usr[2] - usr[1]) / w))
    },
    list(pos = pos, labels = labels, cols = cols, ltys = ltys, lwds = lwds,
         ttl = legend.title),
    getNamespace("graphics"))
    invisible(NULL)
  }
  
  ## the case label, and the quantity beside it, in the outer margin
  page_header <- function(t1, t2) {
    cx <- par("cex")
    if (!is.finite(cx) || cx <= 0) cx <- 1
    mtext(t1, side = 3, outer = TRUE, line = 1.50 / cx,
          font = 2, cex = 1.15 / cx)
    if (!is.null(t2) && nzchar(t2))
      mtext(t2, side = 3, outer = TRUE, line = 0.35 / cx, cex = 0.95 / cx)
  }
  
  ## ------------------------------------------------------------------
  ## d5. panels
  ## ------------------------------------------------------------------
  ## one baseline panel : parametric u_k against the semiparametric one
  draw_base_panel <- function(pm, sp, lnk, main_p, ylim_p) {
    plot(c(xmin, xmax), ylim_p, type = "n",
         xlab = xlab, ylab = ylab, main = main_p, ...)
    
    if (!is.null(sp))
      lines(sp$time, sp$u, type = "s", col = col.semi, lty = lty.semi,
            lwd = lwd)
    lines(pm[, 1], pm[, 2], col = col.par, lty = lty.par, lwd = lwd)
    
    if (isTRUE(legend)) {
      has_sp   <- !is.null(sp)
      semi_lab <- if (lnk == "PH") "Fine-Gray (semiparametric)"
      else "Proportional odds (semiparametric)"
      draw_legend(legend.pos[1],
                  labels = c(paste0(
                    paste0(toupper(substr(distribution, 1, 1)),
                           substr(distribution, 2, nchar(distribution)))," (parametric)"
                  ), if (has_sp) semi_lab),
                  cols = c(col.par, if (has_sp) col.semi),
                  ltys = c(lty.par, if (has_sp) lty.semi),
                  lwds = lwd)
    }
  }
  
  ## one Cox-Snell panel : Nelson-Aalen of the residuals against y = x
  draw_cs_panel <- function(g, main_p, lim_p) {
    if (is.null(g) || NROW(g) == 0 || !any(is.finite(g[, 2]))) {
      plot.new()
      title(main = main_p)
      mtext("no usable events for this cause", side = 3, line = -2,
            cex = 0.8)
      return(invisible(NULL))
    }
    
    plot(lim_p, lim_p, type = "n",
         xlab = xlab.coxsnell, ylab = ylab.coxsnell, main = main_p, ...)
    
    abline(0, 1, col = ref.col, lty = ref.lty)
    
    
    # lines(g[, 1], g[, 2], type = "s", col = col.par, lwd = lwd)
    
    g_plot <- if (NROW(g) > 1) g[-NROW(g), , drop = FALSE] else g
    lines(g_plot[, 1], g_plot[, 2], type = "s", col = col.par, lwd = lwd)
    
    
    if (isTRUE(legend))
      draw_legend(legend.pos[2],
                  labels = c("Nelson-Aalen of residuals", "45-degree line"),
                  cols = c(col.par, ref.col), ltys = c(1, ref.lty),
                  lwds = c(lwd, 1))
  }
  
  ## ------------------------------------------------------------------
  ## d6. one page : every selected event of one case, for one quantity
  ## ------------------------------------------------------------------
  draw_page <- function(ci, page_title, kind) {
    
    ## re-setting mfrow rewinds the panel counter, so the page below
    ## always starts on a fresh device page
    par(mfrow = dims)
    cx <- par("cex")
    if (!is.finite(cx) || cx <= 0) cx <- 1
    par(oma = c(0, 0, 3.4 / cx, 0))
    if (auto_layout) par(mar = c(4.1, 4.1, 2.6, 1.6))
    
    is_base <- identical(kind, "baseline")
    quant   <- if (is_base) ylab else xlab.coxsnell
    lnk     <- vapply(k_of, function(k) link_of(case[ci], k), character(1))
    
    ## what each panel draws
    if (is_base) {
      pm <- para_base_cum_haz_list[[ci]]
      sp <- lapply(seq_len(n_event), function(e) semi_curve(e, lnk[e]))
    } else {
      g <- cox_snell_residual_list[[ci]]
    }
    
    ## axis limits : fixed by the user, shared by the page, or per panel.
    ## On the Cox-Snell page one limit serves both axes, so y = x stays
    ## the diagonal of a square panel.
    fixed <- if (is_base) ylim else lim.coxsnell
    pad   <- if (is_base) 1.05 else 1
    tops  <- vapply(seq_len(n_event), function(e)
      if (is_base) base_top(pm[[e]], sp[[e]]) else cs_top(g[[e]]),
      numeric(1))
    
    if (!is.null(fixed)) {
      lims <- rep(list(fixed), n_event)
    } else if (isTRUE(common.ylim)) {
      v <- suppressWarnings(max(tops, na.rm = TRUE))
      if (!is.finite(v)) v <- 1
      lims <- rep(list(c(0, v * pad)), n_event)
    } else {
      lims <- lapply(tops, function(v) c(0, (if (is.finite(v)) v else 1) * pad))
    }
    
    for (e in seq_len(n_event)) {
      main_p <- paste0("event ", event[e], " (", lnk[e], ")")
      if (is_base) draw_base_panel(pm[[e]], sp[[e]], lnk[e], main_p, lims[[e]])
      else         draw_cs_panel(g[[e]], main_p, lims[[e]])
      
      ## the header belongs to the device page, so it is redrawn
      ## whenever the panels spill over onto another one
      if ((e - 1L) %% per_page == 0L) page_header(page_title, quant)
    }
  }
  
  ## ------------------------------------------------------------------
  ## d7. (case, baseline) then (case, Cox-Snell)
  ## ------------------------------------------------------------------
  for (ci in seq_len(n_case)) {
    ttl <- if (is.null(main)) model[case[ci]]
    else rep(main, length.out = n_case)[ci]
    
    if (baseline) {
      draw_page(ci, ttl, "baseline")
    }
    
    if (coxsnell) {
      draw_page(ci, ttl, "coxsnell")
    }
  }
  
  
  
  out <- list(semi_base_cum_haz_list = semi_base_cum_haz_list,
              para_base_cum_haz_list = para_base_cum_haz_list,
              cox_snell_residual_list = cox_snell_residual_list)
  
}


#' Summarize a Fitted Parametric Competing Risks Regression Model
#'
#' Produces a summary of a fitted \code{"pcrr"} object, including
#' regression coefficients, confidence intervals and baseline parameter
#' estimates, separately for each selected model case.
#'
#' @details
#' A \code{"pcrr"} object holds one fit per model case (see \code{\link{pcrr}}),
#' and \code{summary} reports the tables below for each case selected by
#' \code{case}, which defaults to all of them. Every component of the returned
#' object is therefore a list indexed by the selected cases.
#'
#' Because \eqn{\alpha_k} is held fixed within a case, it is not an estimated
#' parameter of these fits and does not appear in the tables. The assumption
#' tests that chose the cases are reported by \code{\link{print.pcrr}}.
#'
#' Three tables are produced, with a fourth under \code{"gompertz3"}.
#'
#' \emph{Regression coefficients.} For each covariate and cause the estimate
#' \eqn{\hat{\beta}}, its standard error, the Wald statistic
#' \eqn{z = \hat{\beta}/\mathrm{se}(\hat{\beta})} and the two-sided p-value are
#' reported. The Wald test assesses \eqn{H_0\!: \beta = 0}, that is, whether the
#' covariate affects the cumulative incidence of that cause. Standard errors come
#' from the inverse observed information matrix; with \code{variance = FALSE}
#' only the estimates are shown, here and in the other tables.
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
#' \eqn{c_k}, and \eqn{p_k} are reported.
#'
#' These parameters describe the shape of the fitted curve and are not
#' themselves quantities of the fitted model in the sense their names suggest.
#' In particular \eqn{\widehat{p}_k} under \code{"logistic"} is not the
#' long-term event probability unless \eqn{\alpha_k = 0}; use
#' \code{\link{cure.pcrr}} for that, and \code{\link{predict.pcrr}} for the
#' hazard and its turning points.
#'
#' \emph{Baseline shape test.} For \code{distribution = "gompertz3"} an
#' additional Wald test of \eqn{H_0: \eta_k = 0} is reported. Under this null the
#' three-parameter baseline collapses to the two-parameter one, so the test asks
#' whether the extra shape parameter is needed at all: a small p-value supports
#' \code{"gompertz3"}, while a large one indicates that \code{"gompertz2"} is
#' adequate. The test is omitted for \code{"gompertz2"} and \code{"logistic"} fits.
#' 
#' @param object object of class \code{"pcrr"}.
#' @param case integer vector selecting the transformation-model cases to use,
#'  numbered as shown by \code{\link{print.pcrr}}. If \code{NULL} (default),
#'  every available case is used.
#' @param conf.level confidence level for the confidence intervals. The
#'  default is \code{0.95}. It affects the interval table only; the p-values in
#'  the coefficient table always test \eqn{H_0\!: \beta = 0} and do not move
#'  with it, so a 90\% interval excluding 1 alongside a p-value above 0.05 is
#'  expected rather than contradictory.
#' @param digits number of digits to print.
#' @param ... further arguments passed to or from other methods.
#'
#' @return
#' An object of class \code{"summary.pcrr"}. The following components are lists
#' with one element per selected model case, named by the case labels:
#' \describe{
#'   \item{\code{coef}}{regression coefficients, with standard errors, Wald
#'   statistics and p-values.}
#'   \item{\code{conf_int}}{the same coefficients on the exponentiated scale
#'   with confidence limits at \code{conf.level}.}
#'   \item{\code{base_param}}{the distribution-specific baseline parameters:
#'   \code{rho} and \code{tau} for \code{"gompertz2"}, \code{rho}, \code{tau}
#'   and \code{eta} for \code{"gompertz3"}, and \code{b}, \code{c} and
#'   \code{p} for \code{"logistic"}. \code{alpha} is fixed within a case and
#'   is not included.}
#'   \item{\code{shape}}{the \eqn{H_0\!: \eta_k = 0} test under
#'   \code{"gompertz3"}; \code{NULL} for the other baselines and whenever
#'   \code{variance = FALSE}.}
#'   \item{\code{inf}, \code{invinf}}{observed information matrix and its
#'   inverse. A singular matrix is reported with a warning and filled with
#'   \code{NA}. Both are \code{NULL} when the model was fitted with
#'   \code{variance = FALSE}.}
#'   \item{\code{converged}, \code{message}, \code{loglik}, \code{iter}}{
#'   convergence status, optimizer message, maximised log-likelihood and
#'   iteration count of each case. \code{loglik} is \code{NaN} where the
#'   objective terminated on its internal penalty, which is not a fit even
#'   when \code{converged} is \code{TRUE}.}
#' }
#'
#' The remaining components are scalar: \code{call}, \code{n},
#' \code{n_missing}, \code{distribution}, \code{mapping}, \code{digits},
#' \code{case} (the indices selected) and \code{case_model} (the labels of all
#' available cases).
#'
#' @seealso
#' \code{\link{pcrr}},
#' \code{\link{print.summary.pcrr}},
#' \code{\link{plot.pcrr}}
#'
#' @export
summary.pcrr <- function(object, case = NULL, conf.level = 0.95, digits = max(options()$digits - 4, 3), ...){
  
  # if (is.null(object$hess_case_all)) {
  #   message("Please set `variance = TRUE` to perform this operation.")
  #   return(invisible())
  # }
  
  if (!is.numeric(conf.level) || length(conf.level) != 1 || !is.finite(conf.level) ||
      conf.level <= 0 || conf.level >= 1) {
    warning("Invalid confidence level. Using the default value (conf.level = 0.95).")
    conf.level <- 0.95
  }
  
  if (is.null(object$hess_case_all)) {
    variance <- FALSE
  } else {
    variance <- TRUE
  }
  
  mapping <- object$mapping
  
  case_all <- object$case_all
  n_case_all <- nrow(case_all)

  avail_case <- seq_len(n_case_all)
  
  if (is.null(case)) {
    case <- avail_case
  } else {
    if (!is.numeric(case) || length(case) == 0 || any(!is.finite(case)) ||
        any(case != floor(case)) || !all(case %in% avail_case))
      stop("`case` must be one or more of: ",
           paste(avail_case, collapse = ", "),
           "\n\nThe possible model cases are as follows :\n",
           paste(object$case_model, collapse = "\n"), call. = FALSE)
    case <- as.integer(unique(case))
  }
  
  n_case <- length(case)
  
  
  model <- object$case_model
  
  K <- object$k
  P <- object$p
  
  if (object$distribution == "gompertz2"){
    n_param <- K * (3 + P)
    idx_alpha <- (seq_len(K) - 1) * (3 + P) + 1
    idx_rho   <- (seq_len(K) - 1) * (3 + P) + 2
    idx_tau   <- (seq_len(K) - 1) * (3 + P) + 3
    idx_nonbeta <- c(idx_alpha, idx_rho, idx_tau)
  } else if (object$distribution == "gompertz3"){
    n_param <- K * (4 + P)
    idx_alpha <- (seq_len(K) - 1) * (4 + P) + 1
    idx_rho   <- (seq_len(K) - 1) * (4 + P) + 2
    idx_tau   <- (seq_len(K) - 1) * (4 + P) + 3
    idx_eta   <- (seq_len(K) - 1) * (4 + P) + 4
    idx_nonbeta <- c(idx_alpha, idx_rho, idx_tau, idx_eta)
  } else if (object$distribution == "logistic"){
    n_param <- K * (4 + P)
    idx_alpha <- (seq_len(K) - 1) * (4 + P) + 1
    idx_b   <- (seq_len(K) - 1) * (4 + P) + 2
    idx_c   <- (seq_len(K) - 1) * (4 + P) + 3
    idx_p   <- (seq_len(K) - 1) * (4 + P) + 4
    idx_nonbeta <- c(idx_alpha, idx_b, idx_c, idx_p)
  }
  is_alpha <- logical(n_param)
  is_alpha[idx_alpha] <- TRUE
  is_beta  <- rep(TRUE, n_param)
  is_beta[idx_nonbeta] <- FALSE
  is_beta <- is_beta[!is_alpha]
  

  
  
  
  
  if (variance) {
    inf <- vector("list", n_case)
    invinf <- vector("list", n_case)
    
    for (j in 1:n_case) {
      i <- case[j]
      inf[[j]] <- -object$hess_case_all[[i]]
      invinf[[j]] <- tryCatch(solve(inf[[j]]), error = function(e) {
        warning("Hessian matrix is singular or non-invertible. (", model[i], ")")
        matrix(NA_real_, nrow(inf[[j]]), ncol(inf[[j]]),
               dimnames = dimnames(inf[[j]]))
      })
    }
  } else {
    inf <- NULL
    invinf <- NULL
  }

  
  
  
  mle_case_all <- object$mle_case_all
  
  converged <- vector("list", n_case)
  msg <- vector("list", n_case)
  loglik <- vector("list", n_case)
  iter <- vector("list", n_case)
  names(converged) <- model[case]
  names(msg) <- model[case]
  names(loglik) <- model[case]
  names(iter) <- model[case]
  for (j in 1:n_case) {
    i <- case[j]
    converged[[j]] <- ifelse(mle_case_all[[i]]$convergence, FALSE, TRUE) 
    msg[[j]] <- mle_case_all[[i]]$message
    loglik[[j]] <- if(!is.finite(mle_case_all[[i]]$objective) || mle_case_all[[i]]$objective == 1e+100) NaN else -mle_case_all[[i]]$objective
    iter[[j]] <- mle_case_all[[i]]$iterations
  }

  
  coef_tab <- vector("list", n_case)
  ci_tab <- vector("list", n_case)
  base_tab <- vector("list", n_case)

  names(coef_tab) <- model[case]
  names(ci_tab) <- model[case]
  names(base_tab) <- model[case]
  
  if (object$distribution == "gompertz3" && variance) {
    idx_eta <- (seq_len(K) - 1) * (3 + P) + 3
    shape_tab <- vector("list", n_case)
    names(shape_tab) <- model[case]
  } else shape_tab <- NULL
  
  
  if (variance) {
    for (j in 1:n_case) {
      i <- case[j]
      est <- mle_case_all[[i]]$par[!is_alpha]
      
      d <- diag(invinf[[j]])
      names(d) <- rownames(invinf[[j]])
      
      se <- withCallingHandlers(
        sqrt(d), warning = function(w) {
          warning(paste0(model[i], "\n","  sqrt(diag()): ", conditionMessage(w), "\n",
                         "  diag():\n  ", paste(names(d), "=", format(d, digits = 6),collapse = "\n  ")),call. = FALSE)
          invokeRestart("muffleWarning")
        })
      
      zst <- est / se
      pv  <- 2 * (1 - pnorm(abs(zst)))
      
      
      # Regression coefficient table
      coef_tab[[j]] <- cbind(est[is_beta], exp(est[is_beta]), se[is_beta], zst[is_beta], pv[is_beta])
      dimnames(coef_tab[[j]]) <- list(names(est)[is_beta], c("coef", "exp(coef)", "se(coef)", "z value", "Pr(>|z|)"))
      
      # Confidence interval
      a  <- (1 - conf.level) / 2
      a  <- c(a, 1 - a)
      zq <- qnorm(a)
      ci_tab[[j]] <- cbind(exp(est[is_beta]), exp(-est[is_beta]),
                           exp(est[is_beta] + zq[1] * se[is_beta]),
                           exp(est[is_beta] + zq[2] * se[is_beta]))
      dimnames(ci_tab[[j]]) <- list(names(est)[is_beta], c("exp(coef)", "exp(-coef)",
                                                           paste0(format(100 * a, trim = TRUE, digits = 4), "%")))
      
      # Basis parameters
      base_tab[[j]] <- cbind(est[!is_beta], se[!is_beta], zst[!is_beta], pv[!is_beta])
      dimnames(base_tab[[j]]) <- list(names(est)[!is_beta], c("est", "se", "z value", "Pr(>|z|)"))
      
      
      # Baseline shape test : H0 : eta = 0, under which gompertz3 reduces to gompertz2
      if (object$distribution == "gompertz3"){
        z_eta <- est[idx_eta] / se[idx_eta]
        shape_tab[[j]] <- cbind(est = est[idx_eta], se = se[idx_eta], `z value` = z_eta, `Pr(>|z|)` = 2 * (1 - pnorm(abs(z_eta))))
      }
    }
    
  } else {
    for (j in 1:n_case) {
      i <- case[j]
      est <- mle_case_all[[i]]$par[!is_alpha]
      
      
      # Regression coefficient table
      coef_tab[[j]] <- cbind(est[is_beta])
      dimnames(coef_tab[[j]]) <- list(names(est)[is_beta], c("coef"))
      
      # Confidence interval
      ci_tab[[j]] <- cbind(exp(est[is_beta]), exp(-est[is_beta]))
      dimnames(ci_tab[[j]]) <- list(names(est)[is_beta], c("exp(coef)", "exp(-coef)"))
      
      # Basis parameters
      base_tab[[j]] <- cbind(est[!is_beta])
      dimnames(base_tab[[j]]) <- list(names(est)[!is_beta], c("est"))
      
    }
    
    
  }
  

  
  out <- list(call = object$call,
              n = object$n,
              n_missing = object$n_missing,
              distribution = object$distribution,
              mapping = object$mapping,
              digits = digits,
              case = case,
              case_model = object$case_model,
              converged = converged,
              message = msg,
              loglik = loglik,
              iter = iter,
              inf = inf,
              invinf = invinf,
              coef = coef_tab,
              conf_int = ci_tab,
              base_param = base_tab,
              shape = shape_tab
  )
  class(out) <- "summary.pcrr"
  out
}

#' Print a Summary of a Parametric Competing Risks Regression Model
#'
#' Prints a summary object produced by \code{\link{summary.pcrr}}: for each
#' selected model case, the regression coefficients, confidence intervals,
#' baseline parameter estimates and, under \code{"gompertz3"}, the baseline
#' shape test, followed by the convergence status and log-likelihood. A case
#' that did not converge is flagged before its tables. The link tests that
#' selected the cases are printed by \code{\link{print.pcrr}}, not here.
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
  savedig <- options(digits = digits)
  on.exit(options(savedig))
  
  if (is.null(x$invinf)) {
    variance <- FALSE
  } else {
    variance <- TRUE
  }
  
  cat("Parametric Competing Risks Regression\n\n")
  
  if (!is.null(x$call)){ 
    cat("Call:\n")
    dput(x$call)}
 
  cat("\n========================================\n")
  cat("Baseline distribution:", x$distribution, "\n")
  cat("Event codes: ", x$mapping[1], " = failure of interest",
      if (length(x$mapping) > 1)
        paste0("\n             ",paste0(x$mapping[-1], " = competing", collapse = "\n             ")),"\n",sep = "")
  cat("========================================\n")

  
  case <- x$case
  for (j in 1:length(case)){
    i <- case[j]
    cat("\n\n\n========================================\n")
    cat(x$case_model[i])
    cat("\n========================================\n")
    
    if (!x$converged[[j]]){ 
      cat("convergence : ", x$converged[[j]], "\n", "[", x$message[[j]], "]\n",
          "The estimates below are computed at these parameter values ",
          "and should not be interpreted.\n", sep = "")
    }
    
    cat("\nRegression coefficients :\n\n")
    if (variance) {
      printCoefmat(x$coef[[j]], digits = digits, signif.stars = TRUE,
                   has.Pvalue = TRUE, P.values = TRUE, cs.ind = 1:3, tst.ind = 4)
    } else {
      print(x$coef[[j]], digits = digits)
    }
    cat("\n")
    cat("----------------------------------------\n\n")
    print(x$conf_int[[j]]);        cat("\n");
    cat("----------------------------------------\n\n")
    cat("Parameters :\n\n")
    if (variance) {
    printCoefmat(x$base_param[[j]], digits = digits, signif.stars = TRUE,
                 has.Pvalue = TRUE, P.values = TRUE, cs.ind = 1:2, tst.ind = 3)
    } else {
      print(x$base_param[[j]], digits = digits)
    }
    cat("\n")
    
    
    # Baseline shape test
    if (!is.null(x$shape)){
      cat("----------------------------------------\n\n")
      cat("Baseline shape test :\n")
      #cat("[Reduced to two-parameter Gompertz]\n")
      cat("H0 : eta = 0\n\n")
      
      printCoefmat(x$shape[[j]], digits = digits, signif.stars = TRUE, has.Pvalue = TRUE,
                   P.values = TRUE, cs.ind = 1:2, tst.ind = 3)
      
      cat("\n")
    }
    
    cat("----------------------------------------\n\n")
    cat("convergence : ", x$converged[[j]], "  (iteration : ", x$iter[[j]], ")\n", sep = "")
    cat("Log-likelihood =", format(x$loglik[[j]], nsmall = 4), "\n")
  }

  cat("\n\n========================================\n")
  cat("Num. cases =", x$n - x$n_missing)

  if (x$n_missing > 0)
    cat(" (", x$n_missing,
        " cases omitted due to missing values)", sep = "")

  cat("\n")
  
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
#' baseline hazard corresponding to the selected distribution, and
#' \eqn{\widehat{\alpha}_k} is the value at which the link was fixed in the
#' model case being used, that is 0 under PH or 1 under PO.
#'
#' Predictions are obtained separately for each row of \code{cov}, for each
#' event selected by \code{event} and for each model case selected by
#' \code{case}, with one cumulative incidence curve per combination. Both
#' \code{case} and \code{event} default to everything available, so a bare
#' call returns the full set and narrowing it is the deliberate act.
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
#' \eqn{\log\left(\frac{1+p_ke^{-b_kc_k}}{1-p_k}\right)\geq-2b_kc_k.}
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
#' The two coincide only as \eqn{\alpha_k \to 0}, where the derivative
#' condition reduces to \eqn{u_k''(t)=0} and so does not involve the
#' covariates: under PH the covariates scale the hazard and leave its peak
#' where it is. The search is confined to the range of \code{times}, so a
#' turning point beyond the last prediction time is not found, and a value
#' sitting at the edge of the range should be read as a sign of that rather
#' than as a peak.
#'
#' For \code{distribution = "gompertz2"} no turning point is reported.
#' Substituting the two-parameter baseline into the derivative above leaves
#' \eqn{\tau_k e^{\rho_k t}\{\rho_k -
#' \alpha_k \exp(\mathbf{z}^{\top}\boldsymbol{\beta}_k)\tau_k\}}, whose sign
#' does not depend on \eqn{t}: that subdistribution hazard is monotone for
#' every covariate profile, increasing when
#' \eqn{\rho_k > \alpha_k \exp(\mathbf{z}^{\top}\boldsymbol{\beta}_k)\tau_k}
#' and decreasing otherwise.
#'
#' @param object object of class \code{"pcrr"}.
#' @param cov numeric matrix of covariate values. Rows correspond to
#'  covariate profiles and columns correspond to covariates included
#'  in the fitted model. A vector is accepted and is reshaped row-wise when
#'  its length is a multiple of the number of covariates.
#' @param times optional vector of time points at which predictions are
#'  evaluated. If omitted, 200 equally spaced time points over the
#'  observed time range are used.
#' @param case integer vector selecting the transformation-model cases to use,
#'  numbered as shown by \code{\link{print.pcrr}}. If \code{NULL} (default),
#'  every available case is used.
#' @param event vector of event codes, as stored in \code{object$mapping},
#'  identifying the event types for which predictions are computed. If
#'  \code{NULL} (default), every event type in the fitted model is used. A code
#'  that does not appear in \code{object$mapping} is an error.
#' @param ... further arguments passed to or from other methods.
#'
#' @return
#' An object of class \code{"predict.pcrr"}. The first six components are
#' nested lists indexed first by selected model case, named by the case
#' labels, and then by selected event:
#' \describe{
#'   \item{\code{pred}}{matrix with the time points in the first column and one
#'   column of predicted cumulative incidence per covariate profile.}
#'   \item{\code{baseline_hazard}}{matrix with the time points and the baseline
#'   hazard, which carries no covariates.}
#'   \item{\code{subdistribution_hazard}}{matrix with the time points and one
#'   column of predicted subdistribution hazard per covariate profile. Beyond a
#'   finite \code{t_boundary} the entries are \code{NA} rather than zero.}
#'   \item{\code{xmh_base}}{turning point of the baseline hazard under
#'   \code{"gompertz3"} and \code{"logistic"} when one exists at a positive
#'   time, and \code{NA} otherwise. It is always \code{NA} under
#'   \code{"gompertz2"}, whose baseline hazard is monotone.}
#'   \item{\code{xmh_obs}}{vector of profile-specific turning points of the
#'   subdistribution hazard, one per row of \code{cov} and found numerically,
#'   with \code{NA} for a profile that has none. Every entry is \code{NA}
#'   under \code{"gompertz2"}.}
#'   \item{\code{t_boundary}}{named vector of the times at which
#'   \eqn{1+\widehat{\alpha}_k\exp(\mathbf{z}^{\top}\boldsymbol{\widehat{\beta}}_k)
#'   \widehat{u}_k(t)} reaches zero, and \code{Inf} where it does not. Since a
#'   model case fixes \eqn{\alpha_k} at 0 or 1, that quantity is at least 1
#'   throughout and these entries are \code{Inf}; a finite value can arise only
#'   in the unconstrained fit.}
#' }
#'
#' \code{converged} and \code{message} are lists indexed by selected case,
#' giving the convergence status and optimizer message of each case's fit; a
#' case that did not converge is still predicted, and flagged when printed or
#' plotted.
#'
#' The remaining components are scalar: \code{distribution}, \code{eta} (the
#' fitted \eqn{\eta_k}, nested in the same way, for \code{"gompertz3"} only and
#' \code{NULL} otherwise), \code{labels} (the covariate profile labels),
#' \code{case} (the indices selected), \code{case_model} (the labels of all
#' available cases) and \code{event} (the event codes selected).
#'
#' @seealso
#' \code{\link{pcrr}},
#' \code{\link{print.pcrr}},
#' \code{\link{print.predict.pcrr}},
#' \code{\link{plot.predict.pcrr}},
#' \code{\link{cure.pcrr}}
#'
#' @export
predict.pcrr <- function(object, cov, times = NULL, case = NULL, event = NULL, ...){
  P <- object$p
  
  if (is.null(times)) times <- seq(0, object$maxtime, length.out = 200)
  
  if (!is.numeric(times) || any(!is.finite(times)) || any(times < 0))
    stop("`times` must be finite and non-negative.")
  
  if (!is.matrix(cov)) {
    if (is.vector(cov)){
      if (length(cov) %% P == 0) cov <- matrix(cov, ncol = P, byrow = TRUE)
      else stop("cov must have ", P, " column(s).")
    } else cov <- as.matrix(cov)
  }
  if (!is.numeric(cov)) stop("cov must be numeric.")
  if (ncol(cov) != P) stop("cov must have ", P, " column(s).")
  
  
  
  # curve labels for the legend, e.g. "obs 1", "obs 2", ...
  labs <- paste0("obs ", seq_len(nrow(cov)))
  
  case_all <- object$case_all
  n_case_all <- nrow(case_all)
  
  avail_case <- seq_len(n_case_all)
  
  if (is.null(case)) {
    case <- avail_case
  } else {
    if (!is.numeric(case) || length(case) == 0 || any(!is.finite(case)) ||
        any(case != floor(case)) || !all(case %in% avail_case))
      stop("`case` must be one or more of: ",
           paste(avail_case, collapse = ", "),
           "\n\nThe possible model cases are as follows :\n",
           paste(object$case_model, collapse = "\n"), call. = FALSE)
    case <- as.integer(unique(case))
  }
  
  n_case <- length(case)
  
  
  
  mapping <- object$mapping
  if (is.null(event)) {
    event <- mapping
  } else {
    if (!is.numeric(event) || length(event) == 0 || !all(event %in% mapping))
      stop("`event` must be one or more of: ", paste(mapping, collapse = ", "), call. = FALSE)
    event <- unique(event)
  }
  
  
  
  model <- object$case_model
  
  
  
  
  
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
      # u <- rep(NaN, length(t))
      # valid <- 1 + inner > 0
      # u[valid] <- -log1p(inner[valid])
      u <- -log1p(inner)
      u
    }
    du_k <- function(t, b, c, p, tol){
      du <- b * p * (1 + exp(-b * c)) * exp(b * (t - c))/((1 + exp(b * (t - c)))*(1 + p * exp(-b * c)+(1 - p)*exp(b * (t - c))))
      du
    }
  }

  
  tol <- 1e-12
  
  

  
  
  pred_list <- vector("list", n_case)
  #base_cum_haz_list <- vector("list", n_case)
  base_haz_list <- vector("list", n_case)
  sub_haz_list <- vector("list", n_case)
  xmh_base_list <- vector("list", n_case)
  xmh_obs_list <- vector("list", n_case)
  t_boundary_list <- vector("list", n_case)
  eta_list <- vector("list", n_case)
  converged <- vector("list", n_case)
  msg <- vector("list", n_case)
  
  names(pred_list) <- model[case]
  #names(base_cum_haz_list) <- model[case]
  names(base_haz_list) <- model[case]
  names(sub_haz_list) <- model[case]
  names(xmh_base_list) <- model[case]
  names(xmh_obs_list) <- model[case]
  names(t_boundary_list) <- model[case]
  names(eta_list) <- model[case]
  names(converged) <- model[case]
  names(msg) <- model[case]
  
  for (l in 1:n_case) {
    i <- case[l]
    coef <- object$mle_case_all[[i]]$par
    
    pred_list[[l]] <- vector("list", length(event))
    #base_cum_haz_list[[l]] <- vector("list", length(event))
    base_haz_list[[l]] <- vector("list", length(event))
    sub_haz_list[[l]] <- vector("list", length(event))
    xmh_base_list[[l]] <- vector("list", length(event))
    xmh_obs_list[[l]] <- vector("list", length(event))
    t_boundary_list[[l]] <- vector("list", length(event))
    converged[[l]] <- object$mle_case_all[[i]]$convergence == 0
    msg[[l]] <- object$mle_case_all[[i]]$message
    
    names(pred_list[[l]]) <- paste0("event ", event)
    #names(base_cum_haz_list[[l]]) <- paste0("event ", event)
    names(base_haz_list[[l]]) <- paste0("event ", event)
    names(sub_haz_list[[l]]) <- paste0("event ", event)
    names(xmh_base_list[[l]]) <- paste0("event ", event)
    names(xmh_obs_list[[l]]) <- paste0("event ", event)
    names(t_boundary_list[[l]]) <- paste0("event ", event)
    
    if (object$distribution == "gompertz3") {
      eta_list[[l]] <- vector("list", length(event))
      names(eta_list[[l]]) <- paste0("event ", event)
    }
    
    for (e in 1:length(event)) {
      k <- match(event[e], mapping)
      
      if (object$distribution == "gompertz2"){
        alpha <- coef[(k - 1) * (3 + P) + 1]
        rho   <- coef[(k - 1) * (3 + P) + 2]
        tau   <- coef[(k - 1) * (3 + P) + 3]
        eta   <- 0
        beta  <- coef[((k - 1) * (3 + P) + 4):(k * (3 + P))]
        
        u <- u_k(times, rho, tau, eta, tol)
        base_haz <- du_k(times, rho, tau, eta, tol)
      } else if (object$distribution == "gompertz3"){
        alpha <- coef[(k - 1) * (4 + P) + 1]
        rho   <- coef[(k - 1) * (4 + P) + 2]
        tau   <- coef[(k - 1) * (4 + P) + 3]
        eta   <- coef[(k - 1) * (4 + P) + 4]
        beta  <- coef[((k - 1) * (4 + P) + 5):(k * (4 + P))]
        
        u <- u_k(times, rho, tau, eta, tol)
        base_haz <- du_k(times, rho, tau, eta, tol)
      } else if (object$distribution == "logistic"){
        alpha <- coef[(k - 1) * (4 + P) + 1]
        b   <- coef[(k - 1) * (4 + P) + 2]
        c   <- coef[(k - 1) * (4 + P) + 3]
        p   <- coef[(k - 1) * (4 + P) + 4]
        beta  <- coef[((k - 1) * (4 + P) + 5):(k * (4 + P))]
        
        u <- u_k(times, b, c, p, tol)
        base_haz <- du_k(times, b, c, p, tol)
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
        valid <- (!is.nan(u)) & (base > tol)
        cif[valid, j] <- 1.0 - base[valid]^(-1.0 / alpha)
        cif[!valid, j] <- 1.0 # treat undefined region as CIF = 1
        if (any(!valid)){
          undef_t <- times[which(!valid)[1]]
          warning(model[i], "\nCIF cannot be evaluated after time ", undef_t, " for observation ", j, " (event ", event[e], ").", call. = FALSE)
        }
      }
      
      
      
      
      # maximum baseline hazard rate
      xmh_base <- NA_real_
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
        valid <- (!is.nan(u)) & (base > tol)
        sub_haz[valid, j] <- ezb * base_haz[valid] / base[valid]
        sub_haz[!valid, j] <- NA_real_ # treat undefined region as NA for prediction
        if (any(!valid)){
          undef_t <- times[which(!valid)[1]]
          warning(model[i], "\nSubdistribution Hazard cannot be evaluated after time ", undef_t, " for observation ", j, " (event ", event[e], ").", call. = FALSE)
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
                || (rho < 0 && 1 + alpha * ezb * u_k(Inf, rho, tau, eta, tol) < 0)
                || (rho > 0 && eta < 0 && 1 + alpha * ezb * u_k(Inf, rho, tau, eta, tol) < 0)){
              inner <- function(t, alpha, rho, tau, eta, ezb){
                u <- u_k(t, rho, tau, eta, tol)
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
                res <- tryCatch(uniroot(inner, interval = c(0, time_inner), tol = tol, alpha = alpha, rho = rho, tau = tau, eta = eta, ezb = ezb), 
                                error = function(e) NULL)
                if (!is.null(res)) t_boundary[j] <- res$root
              } 
            }
          }
        }
      } else if (object$distribution == "logistic"){
        for (j in seq_len(nrow(cov))){
          ezb <- exp(sum(cov[j, ] * beta))
          if (alpha < 0 && (p > 1 - exp(1/(alpha * ezb)))){
            t_boundary[j] <- (log1p(p*exp(-b * c) - exp(1/(alpha * ezb))) - log(exp(1/(alpha * ezb)) + p - 1)) / b + c
          } # else if (alpha > 0 && p > 1) {
          #  t_boundary[j] <- (log1p(p*exp(-b * c)) - log(p - 1)) / b + c
          #}
        }
      }
      
      
      
      # subdistribution hazard peak
      xmh_obs <- rep(NA_real_, nrow(cov))
      if (object$distribution == "gompertz3"){

        lambda_prime <- function(t, alpha, rho, tau, eta, ezb){
          u <- u_k(t, rho, tau, eta, tol)
          du <- du_k(t, rho, tau, eta, tol)
          rho * (1 + eta * exp(rho * t)) * (1 + alpha * ezb * u) - alpha * ezb * du
        }
        
        for (j in seq_len(nrow(cov))){
          ezb <- exp(sum(cov[j, ] * beta))
          
          lambda_prime_t0 <- lambda_prime(0, alpha, rho, tau, eta, ezb)
          if (abs(lambda_prime_t0) < tol) {xmh_obs[j] <- 0; next}
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
          res <- tryCatch(uniroot(lambda_prime, c(0, t_max), tol = tol, alpha = alpha, rho = rho, tau = tau, eta = eta, ezb = ezb),
                          error = function(e) {warning("Failed to find root: ", conditionMessage(e)); NULL})
          if (is.null(res)) next
          xmh_obs[j] <- res$root
        }
        
        names(xmh_obs) <- paste0("obs ", seq_len(nrow(cov)))
        
      } else if (object$distribution == "logistic"){

        lambda_prime <- function(t, alpha, b, c, p, ezb){
          u <- u_k(t, b, c, p, tol)
          du <- du_k(t, b, c, p, tol)
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
            res <- tryCatch(uniroot(lambda_prime, c(0, t_max), tol = tol, alpha = alpha, b = b, c = c, p = p, ezb = ezb),
                            error = function(e) {warning("Failed to find root: ", conditionMessage(e)); NULL})
            if (is.null(res)) next
            xmh_obs[j] <- res$root
          } else if (abs(exp(b*c) + p *(1 - alpha * ezb) -1) < tol) {
            xmh_obs[j] <- 0
          }
        }
        names(xmh_obs) <- paste0("obs ", seq_len(nrow(cov)))
      }
      
      colnames(sub_haz) <- paste0("obs ", seq_len(nrow(cov)))
      names(t_boundary) <- paste0("obs ", seq_len(nrow(cov)))
      

      pred <- cbind(times, cif)
      colnames(pred) <- c("time", labs)
      base_cum_haz <- cbind(time = times, `baseline cumulative hazard` = u)
      base_haz <- cbind(time = times, `baseline hazard` = base_haz)
      sub_haz <- cbind(time = times, sub_haz)
      
      
      
      pred_list[[l]][[e]] <- pred
      #base_cum_haz_list[[l]][[e]] <- base_cum_haz
      base_haz_list[[l]][[e]] <- base_haz
      sub_haz_list[[l]][[e]] <- sub_haz
      xmh_base_list[[l]][[e]] <- xmh_base
      xmh_obs_list[[l]][[e]] <- xmh_obs
      t_boundary_list[[l]][[e]] <- t_boundary
      if (object$distribution == "gompertz3") eta_list[[l]][[e]] <- eta
      
    }

  }


  out <- list(pred = pred_list,
              baseline_hazard = base_haz_list,
              #baseline_cumulative_hazard = base_cum_haz_list,
              subdistribution_hazard = sub_haz_list,
              xmh_base = xmh_base_list,
              xmh_obs = xmh_obs_list,
              t_boundary = t_boundary_list,
              distribution = object$distribution,
              eta = if(object$distribution=="gompertz3") eta_list else NULL,
              labels = labs,
              case_model = object$case_model,
              case = case,
              event = event,
              converged = converged,
              message = msg
              )
  
  class(out) <- "predict.pcrr"
  
  out
}


#' Print Predicted Cumulative Incidence Functions
#'
#' Prints a summary of predicted cumulative incidence functions produced by
#' \code{\link{predict.pcrr}}, for each selected case and event: the hazard
#' turning points (for \code{"gompertz3"} and \code{"logistic"}) and the
#' predicted values. A case that did not converge is flagged before its
#' output.
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
#' \code{\link{predict.pcrr}},
#' \code{\link{plot.predict.pcrr}}
#'
#' @export
print.predict.pcrr <- function(x, digits = 4, ...){
  
  
  cat("Predicted Cumulative Incidence")
  
  
  for (l in 1:length(x$case)){
    i <- x$case[l]
    cat("\n\n========================================\n")
    cat(x$case_model[i])
    cat("\n========================================\n")
    
    if (!x$converged[[l]]) {
      cat("\nconvergence : ", x$converged[[l]], "\n",
          "[", x$message[[l]], "]\n",
          "The predictions below are computed at these parameter values ",
          "and should not be interpreted.\n", sep = "")
    }
    

    for (e in 1:length(x$event)) {
      
      xmh_base_le   <- x$xmh_base[[l]][[e]]
      xmh_obs_le <- x$xmh_obs[[l]][[e]]
      
      cat("\n[event : ", x$event[e], "]\n", sep = "")
      
      
      if (x$distribution == "gompertz3"){
      
        if (is.na(xmh_base_le)) {
          cat("\nMaximum hazard rate time (baseline) : None\n")
        } else {
          cat("\nMaximum hazard rate time (baseline) :\n")
          cat("  baseline :", format(round(xmh_base_le, digits)), "\n")
        }
        cat("\nStationary point of subdistribution hazard :\n")
        if (any(is.finite(xmh_obs_le))){
          if (x$eta[[l]][[e]] < 0){
            cat("  Type : Maximum (unimodal)\n")
          } else if (x$eta[[l]][[e]] > 0){
            cat("  Type : Minimum (U-shape)\n")
          } else {
            cat("  Type : Gompertz2 case (no stationary point)\n")
          }
          
          for (j in seq_along(xmh_obs_le)){
            cat(sprintf("  %-8s : %s\n",
                        names(xmh_obs_le)[j], ifelse(is.finite(xmh_obs_le[j]), 
                                                              format(round(xmh_obs_le[j], digits)),"None")))}
        } else {
          cat("  (no stationary point found for given covariates)\n")
        }
        

      
      } else if (x$distribution == "logistic"){
        
        if (is.na(xmh_base_le)) {
          cat("\nMaximum hazard rate time (baseline) : None\n")
        } else {
          cat("\nMaximum hazard rate time (baseline) :\n")
          cat("  baseline :", format(round(xmh_base_le, digits)), "\n")
        }
        cat("\nStationary point of subdistribution hazard :\n")
        if (any(is.finite(xmh_obs_le))){
          cat("  Type : Maximum (unimodal)\n")
          
          for (j in seq_along(xmh_obs_le)){
            cat(sprintf("  %-8s : %s\n", names(xmh_obs_le)[j],
                        ifelse(is.finite(xmh_obs_le[j]),
                               format(round(xmh_obs_le[j], digits)), "None")))
          }
        } else {
          cat("  (no stationary point found for given covariates)\n")
        }
      }
      
      cat("\n\n")
      print(x$pred[[l]][[e]], ...)
      
      if (e != length(x$event)){
        cat("\n----------------------------------------\n")
      }
    }
  
  }

  invisible(x)
}

#' Plot Predicted Cumulative Incidence Functions
#'
#' Plots the cumulative incidence functions and subdistribution hazards
#' produced by \code{\link{predict.pcrr}}. For each selected model case a page
#' of cumulative incidence panels is drawn, followed, when the hazards are
#' available, by a page of subdistribution hazard panels. Each page has one
#' panel per selected event and each panel one curve per covariate profile.
#'
#' @details
#' A \code{"predict.pcrr"} object may hold several transformation-model cases,
#' several event types and several covariate profiles, and each has its own
#' place in the display. \emph{Cases} are never overlaid: the selected cases
#' are drawn in turn, each case's cumulative incidence page followed by its
#' subdistribution hazard page (see \code{hazard}). \emph{Events} are the
#' panels of a page, each titled by its event code. \emph{Covariate profiles}
#' are the curves within a panel, distinguished by \code{color} (and, if
#' supplied, \code{lty}), so that a shared color means a shared profile across
#' panels and pages.
#'
#' Each page carries a header giving the case label (or \code{main}) and,
#' beneath it, the quantity shown (\code{ylab} or \code{ylab.hazard}). If the
#' events do not fit the grid (see \code{mfrow}) the page continues on further
#' pages, each with the header repeated. Graphical parameters are restored on
#' exit. A case whose optimization did not converge is still drawn, with a
#' warning that its curves should not be interpreted.
#'
#' With \code{mark.xmh = TRUE} the turning points of the subdistribution hazard
#' are marked on both pages: each profile's turning point by a filled point on
#' its curve and a vertical line in its color, and the baseline turning point
#' by a vertical line in \code{xmh.col}. Turning points that are missing or lie
#' outside the range from \code{xmin} to \code{xmax} are not drawn. The legend
#' gives each profile's turning point and, when it is drawn, the baseline one.
#'
#' When the number of covariate profiles exceeds \code{max.marks} the
#' annotation is simplified automatically: the profile-specific vertical lines
#' are dropped (the points on the curves are kept) and the numeric
#' \eqn{x_{mh}} values are removed from the legend. The baseline line stays,
#' labelled without its value.
#'
#' @param x object of class \code{"predict.pcrr"}.
#' @param case integer vector of model cases to draw, using the case numbers
#'  stored in \code{x$case}. Cases are drawn in the order given, each on its
#'  own pages, and duplicates are dropped. If \code{NULL} (default) every
#'  stored case is drawn.
#' @param event vector of event types to draw, using the codes stored in
#'  \code{x$event}, one panel per event in the order given. If \code{NULL}
#'  (default) every stored event is drawn.
#' @param color line colors for the covariate profiles, recycled to the number
#'  of profiles. Defaults to \code{seq_len(ncurve) + 1}, where \code{ncurve} is
#'  the number of profiles. The profiles are labelled by \code{x$labels}, or
#'  \code{"obs 1"}, \code{"obs 2"}, \dots when these are absent.
#' @param lty line types for the covariate profiles, recycled to the number of
#'  profiles. Defaults to 1 for every profile.
#' @param ylim,ylim.hazard y-axis limits of the cumulative incidence panels and
#'  of the subdistribution hazard panels respectively, applied to every panel
#'  of that page. If \code{NULL} (default) each panel runs from 0 to the
#'  largest value of its curves over all prediction times, or to 1.05 times
#'  that value on the hazard page; see also \code{common.ylim}.
#' @param xmin,xmax range of the x-axis, shared by every panel on both pages.
#'  If \code{xmax} is \code{NULL} (default), the largest prediction time of the
#'  selected cases is used. An \code{xmax} that is not larger than \code{xmin}
#'  is replaced by \code{xmin + 1}.
#' @param xlab label for the x-axis.
#' @param ylab,ylab.hazard y-axis labels of the cumulative incidence panels and
#'  of the subdistribution hazard panels; each also appears in the header of
#'  its page.
#' @param legend logical value indicating whether a legend is drawn on each
#'  panel.
#' @param legend.pos position of the legend on both pages, passed to
#'  \code{\link[graphics]{legend}}. The default is \code{"topleft"}.
#' @param legend.title optional title for the legend.
#' @param lwd line width for the curves.
#' @param main title(s) replacing the case label in the page headers, recycled
#'  over the selected cases so that each element titles one case's pages. If
#'  \code{NULL} (default) the case label is used. Panel titles, which name the
#'  event, are not affected.
#' @param hazard logical value. If \code{TRUE} (default) and \code{x} stores
#'  subdistribution hazards, each case's cumulative incidence page is followed
#'  by a page of subdistribution hazards; otherwise only the cumulative
#'  incidence is drawn.
#' @param mark.xmh logical value. If \code{TRUE} (default), the profile-specific
#'  hazard turning points (\code{x$xmh_obs}) and the baseline turning point
#'  (\code{x$xmh_base}) are marked when they exist; see Details.
#' @param max.marks maximum number of covariate profiles for which the
#'  profile-specific vertical lines and the numeric legend entries are kept.
#'  Above this threshold the annotation is simplified. Default is 4.
#' @param xmh.lty line type for the vertical reference lines at the hazard
#'  turning points. Default is 3 (dotted).
#' @param xmh.col color of the baseline \eqn{x_{mh}} reference line. Default is
#'  \code{"black"}, so that the covariate-free baseline turning point stands out
#'  from the profile-specific lines, which follow the color of their own curve.
#' @param mfrow panel grid of a page as \code{c(rows, columns)}. If \code{NULL}
#'  (default), a grid already set with \code{par(mfrow = )} is used together
#'  with the current margins; failing that, the grid has \code{floor(sqrt(n))}
#'  rows and enough columns for the \code{n} selected events.
#' @param common.ylim logical value. If \code{TRUE}, all panels of a page share
#'  one y-axis range instead of being scaled separately. Ignored where
#'  \code{ylim} or \code{ylim.hazard} fixes the range. Default \code{FALSE}.
#' @param ask logical value. If \code{TRUE}, the user is prompted before each
#'  new page. If \code{NULL} (default), prompting is enabled when more than one
#'  page is drawn on an interactive device (see
#'  \code{\link[grDevices]{dev.interactive}}).
#' @param ... additional graphical parameters passed to
#'  \code{\link[graphics]{plot.default}} when each panel is set up; they do not
#'  reach the curves or the legends.
#'
#' @return
#' Produces the plots and returns the input object invisibly.
#'
#' @seealso
#' \code{\link{pcrr}},
#' \code{\link{predict.pcrr}}
#'
#' @export
plot.predict.pcrr <- function(x, case = NULL, event = NULL,
                              color = NULL, lty = NULL,
                              ylim = NULL, ylim.hazard = NULL,
                              xmin = 0, xmax = NULL,
                              xlab = "Time",
                              ylab = "Cumulative Incidence",
                              ylab.hazard = "Subdistribution Hazard",
                              legend = TRUE, legend.pos = "topleft",
                              legend.title = NULL, lwd = 2, main = NULL,
                              hazard = TRUE, mark.xmh = TRUE, max.marks = 4,
                              xmh.lty = 3, xmh.col = "black",
                              mfrow = NULL, common.ylim = FALSE,
                              ask = NULL, ...) {
  
  ## ------------------------------------------------------------------
  ## 1. case selection : one page pair per case, never overlaid
  ## ------------------------------------------------------------------
  avail_case <- x$case
  if (is.null(case)) {
    case <- avail_case
    if (is.null(avail_case))
      message("Use `case = ` to select cases :\n",
              paste(x$case_model[avail_case], collapse = "\n"))
  } else {
    if (!is.numeric(case) || length(case) == 0 || any(!is.finite(case)) ||
        any(case != floor(case)) || !all(case %in% avail_case))
      stop("`case` must be one or more of: ",
           paste(avail_case, collapse = ", "), "\n",
           paste(x$case_model[avail_case], collapse = "\n"), call. = FALSE)
    case <- as.integer(unique(case))
  }
  lpos   <- match(case, avail_case)      # position inside the stored lists
  n_case <- length(case)
  
  if (!is.null(x$converged)) {
    bad <- which(!vapply(lpos, function(l) isTRUE(x$converged[[l]]), logical(1)))
    if (length(bad) > 0)
      warning("The following model case(s) did not converge. The curves are ",
              "drawn at these parameter values and should not be interpreted :\n",
              paste(x$case_model[case[bad]], collapse = "\n"), call. = FALSE)
  }
  
  ## ------------------------------------------------------------------
  ## 2. event selection : one panel each, inside every page
  ## ------------------------------------------------------------------
  avail_event <- x$event
  if (is.null(event)) {
    event <- avail_event
  } else {
    if (length(event) == 0 || !all(event %in% avail_event))
      stop("`event` must be one or more of: ",
           paste(avail_event, collapse = ", "), call. = FALSE)
    event <- unique(event)
  }
  epos    <- match(event, avail_event)
  n_event <- length(event)
  
  ## ------------------------------------------------------------------
  ## 3. curves, colours and line types
  ## ------------------------------------------------------------------
  labs   <- x$labels
  ncurve <- length(labs)
  if (is.null(labs) || ncurve == 0) {
    ncurve <- ncol(x$pred[[lpos[1]]][[epos[1]]]) - 1
    labs   <- paste("obs", seq_len(ncurve))
  }
  
  if (is.null(color)) color <- seq_len(ncurve) + 1
  color <- rep(color, length.out = ncurve)
  if (is.null(lty))   lty <- 1
  lty   <- rep(lty, length.out = ncurve)
  
  if (!is.numeric(xmin) || length(xmin) != 1 || !is.finite(xmin) || xmin < 0)
    stop("`xmin` must be a single non-negative number.", call. = FALSE)
  if (is.null(xmax))
    xmax <- suppressWarnings(max(unlist(lapply(lpos, function(l)
      lapply(epos, function(e) x$pred[[l]][[e]][, 1]))), na.rm = TRUE))
  if (!is.numeric(xmax) || length(xmax) != 1 || !is.finite(xmax) ||
      xmax <= xmin)
    stop("`xmax` must be a single number larger than `xmin`.", call. = FALSE)
  
  ## simplify the annotation when the panel is crowded
  crowded     <- ncurve > max.marks
  show_lines  <- isTRUE(mark.xmh) && !crowded
  show_values <- isTRUE(mark.xmh) && !crowded
  
  show_haz <- isTRUE(hazard) && !is.null(x$subdistribution_hazard)
  
  ## ------------------------------------------------------------------
  ## 4. panel grid : K events on one page
  ## ------------------------------------------------------------------
  grid_dim <- function(n) {
    nr <- max(1L, as.integer(floor(sqrt(n))))
    c(nr, as.integer(ceiling(n / nr)))
  }
  
  user_split <- !identical(as.integer(par("mfrow")), c(1L, 1L))
  if (is.null(mfrow)) {
    auto_layout <- !user_split
    dims <- if (user_split) as.integer(par("mfrow")) else grid_dim(n_event)
  } else {
    if (!is.numeric(mfrow) || length(mfrow) != 2 || any(!is.finite(mfrow)) ||
        any(mfrow < 1))
      stop("`mfrow` must be a numeric vector of length 2, e.g. c(1, 2).",
           call. = FALSE)
    auto_layout <- TRUE
    dims <- as.integer(mfrow)
  }
  per_page <- dims[1] * dims[2]
  
  ## ------------------------------------------------------------------
  ## 5. device set-up
  ## ------------------------------------------------------------------
  n_quant <- 1L + isTRUE(show_haz)
  n_page  <- n_case * n_quant * ceiling(n_event / per_page)
  
  op <- par(no.readonly = TRUE)
  on.exit(par(op), add = TRUE)
  
  if (is.null(ask)) ask <- (n_page > 1) && dev.interactive()
  if (isTRUE(ask)) {
    oask <- devAskNewPage(TRUE)
    on.exit(devAskNewPage(oask), add = TRUE)
  }
  
  ## ------------------------------------------------------------------
  ## 6. small helpers
  ## ------------------------------------------------------------------
  ## largest finite value of the curve columns, NA-safe
  col_max <- function(m) {
    v <- m[, -1, drop = FALSE]
    if (!any(is.finite(v))) return(NA_real_)
    max(v[is.finite(v)])
  }
  ## height of a curve at xout, NA-safe
  curve_y <- function(tt, yy, xout) {
    ok <- is.finite(tt) & is.finite(yy)
    if (sum(ok) < 2) return(NA_real_)
    approx(tt[ok], yy[ok], xout = xout)$y
  }
  ## is this turning point drawable?
  usable <- function(v) length(v) == 1L && is.finite(v) && v >= xmin && v <= xmax
  
  ## legend of a single panel : the profiles, then the baseline marker
  make_legend <- function(mk, base_v) {
    lg <- labs
    if (show_values && !is.null(mk) && any(is.finite(mk)))
      lg <- ifelse(is.finite(mk), sprintf("%s   x_mh = %.2f", labs, mk), labs)
    out <- list(legend = lg, col = color, lty = lty, lwd = rep(lwd, ncurve))
    
    if (isTRUE(mark.xmh) && usable(base_v)) {
      lab_b <- if (!crowded) sprintf("baseline x_mh = %.2f", base_v)
      else "baseline x_mh"
      out$legend <- c(out$legend, lab_b)
      out$col    <- c(out$col, xmh.col)
      out$lty    <- c(out$lty, xmh.lty)
      out$lwd    <- c(out$lwd, 1.6)
    }
    out
  }
  
  ## the case label, and the quantity beside it, in the outer margin
  page_header <- function(t1, t2) {
    cx <- par("cex")
    if (!is.finite(cx) || cx <= 0) cx <- 1
    mtext(t1, side = 3, outer = TRUE, line = 1.50 / cx,
          font = 2, cex = 1.15 / cx)
    if (!is.null(t2) && nzchar(t2))
      mtext(t2, side = 3, outer = TRUE, line = 0.35 / cx, cex = 0.95 / cx)
  }
  
  ## one panel : one event of one case, one quantity
  draw_panel <- function(mat, mk, base_v, ylim_p, ylab_p, main_p) {
    plot(c(xmin, xmax), ylim_p, type = "n",
         xlab = xlab, ylab = ylab_p, main = main_p, ...)
    
    ## baseline turning point
    if (isTRUE(mark.xmh) && usable(base_v))
      abline(v = base_v, lty = xmh.lty, col = xmh.col, lwd = 1.6)
    
    ## profile-specific vertical reference lines
    if (show_lines && !is.null(mk))
      for (j in seq_len(ncurve))
        if (usable(mk[j])) abline(v = mk[j], lty = xmh.lty, col = color[j])
    
    ## the curves themselves
    for (j in seq_len(ncurve))
      lines(mat[, 1], mat[, j + 1], lty = lty[j], col = color[j], lwd = lwd)
    
    ## turning points marked on the curves (kept even when crowded)
    if (isTRUE(mark.xmh) && !is.null(mk))
      for (j in seq_len(ncurve)) {
        if (!usable(mk[j])) next
        yy <- curve_y(mat[, 1], mat[, j + 1], mk[j])
        if (is.finite(yy))
          points(mk[j], yy, pch = 19, col = color[j], cex = 0.9)
      }
    
    if (isTRUE(legend)) {
      lg <- make_legend(mk, base_v)
      legend(legend.pos, legend = lg$legend, lty = lg$lty, col = lg$col,
             lwd = lg$lwd, title = legend.title, bty = "n")
    }
  }
  
  ## one page : every selected event of one case, for one quantity
  draw_page <- function(l, page_title, kind) {
    
    ## re-setting mfrow rewinds the panel counter, so the page below
    ## always starts on a fresh device page
    par(mfrow = dims)
    cx <- par("cex")
    if (!is.finite(cx) || cx <= 0) cx <- 1
    par(oma = c(0, 0, 3.4 / cx, 0))
    if (auto_layout) par(mar = c(4.1, 4.1, 2.6, 1.6))
    
    is_cif  <- identical(kind, "cif")
    quant   <- if (is_cif) ylab else ylab.hazard
    mats    <- lapply(epos, function(e)
      if (is_cif) x$pred[[l]][[e]] else x$subdistribution_hazard[[l]][[e]])
    
    ## y limits : fixed by the user, shared by the page, or per panel
    fixed <- if (is_cif) ylim else ylim.hazard
    pad   <- if (is_cif) 1 else 1.05
    top   <- function(m) {
      v <- col_max(m)
      if (!is.finite(v) || v <= 0) v <- 1
      c(0, v * pad)
    }
    if (!is.null(fixed)) {
      ylims <- rep(list(fixed), n_event)
    } else if (isTRUE(common.ylim)) {
      v <- suppressWarnings(max(vapply(mats, col_max, numeric(1)), na.rm = TRUE))
      if (!is.finite(v) || v <= 0) v <- 1
      ylims <- rep(list(c(0, v * pad)), n_event)
    } else {
      ylims <- lapply(mats, top)
    }
    
    for (ei in seq_len(n_event)) {
      e  <- epos[ei]
      ev <- event[ei]
      
      mk <- if (isTRUE(mark.xmh)) x$xmh_obs[[l]][[e]] else NULL
      bv <- if (isTRUE(mark.xmh)) x$xmh_base[[l]][[e]] else NA_real_
      bv <- if (is.null(bv) || length(bv) == 0) NA_real_ else as.numeric(bv)[1]
      
      draw_panel(mats[[ei]], mk, bv, ylims[[ei]], quant, paste("event", ev))
      
      ## the header belongs to the device page, so it is redrawn
      ## whenever the panels spill over onto another one
      if ((ei - 1L) %% per_page == 0L) page_header(page_title, quant)
    }
  }
  
  ## ------------------------------------------------------------------
  ## 7. (case, cumulative incidence) then (case, subdistribution hazard)
  ## ------------------------------------------------------------------
  for (ci in seq_len(n_case)) {
    l   <- lpos[ci]
    ttl <- if (is.null(main)) x$case_model[case[ci]]
    else rep(main, length.out = n_case)[ci]
    
    draw_page(l, ttl, "cif")
    if (show_haz) draw_page(l, ttl, "hazard")
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
#' The third state cannot occur in a fit obtained from \code{\link{pcrr}}.
#' \eqn{A_k(t;\mathbf{z})} can decrease only when \eqn{\alpha_k < 0}, whereas
#' every model case fixes \eqn{\alpha_k} at 0 or 1, so
#' \eqn{A_k(t;\mathbf{z}) \geq 1} throughout and \code{t_boundary} is
#' \code{Inf}. It is retained because it is a genuine state of the
#' transformation model for negative \eqn{\alpha_k}.
#'
#'
#' @param object a fitted model object.
#' @param ... further arguments passed to or from methods.
#'
#' @return
#' An object whose class depends on the method; for a \code{"pcrr"} fit see
#' \code{\link{cure.pcrr}} below.
#'
#'
#' @export
cure <- function(object, ...) UseMethod("cure")

#' @rdname cure
#'
#' @param cov numeric matrix of covariate values, one row per profile. A
#'  vector is accepted and is reshaped row-wise when its length is a multiple
#'  of the number of covariates.
#' @param case integer vector selecting the transformation-model cases to use,
#'  numbered as shown by \code{\link{print.pcrr}}. If \code{NULL} (default),
#'  every available case is used.
#' @param event vector of event codes, as stored in \code{object$mapping}, for
#'  which the cure fraction is computed. If \code{NULL} (default), every event
#'  type in the fitted model is used. Reporting every cause is the sensible
#'  default here, since all of them are estimated from the same likelihood.
#'
#' @return
#' An object of class \code{"cure.pcrr"}. The first three components are nested
#' lists indexed first by selected model case, named by the case labels, and
#' then by selected event:
#' \describe{
#'   \item{\code{cure}}{estimated cure fractions, one entry per row of
#'   \code{cov}. Positive for \code{"cure"} profiles, zero for
#'   \code{"no cure (asymptotic)"} profiles and \code{NA} for
#'   \code{"no cure (finite support)"} profiles.}
#'   \item{\code{status}}{the classification of each profile, as described
#'   above.}
#'   \item{\code{t_boundary}}{the time at which \eqn{A_k(t; \mathbf{z})}
#'   reaches zero where such a finite time exists, and \code{Inf} otherwise.
#'   In a fit from \code{\link{pcrr}} it is always \code{Inf}, for the reason
#'   given above.}
#' }
#'
#' The remaining components are \code{case} (the indices selected),
#' \code{case_model} (the labels of all available cases), \code{event} (the
#' event codes selected), and \code{converged} and \code{message}, lists
#' indexed by selected case giving the convergence status and optimizer
#' message of each case's fit.
#'
#' Because the estimate is a model-based extrapolation beyond the observed
#' follow-up, it is a statement about the fitted distribution rather than an
#' observed proportion. Where more than one model case survived, comparing the
#' cure fractions across cases is informative: the cases differ most where the
#' data run out, which is exactly where this quantity is read.
#'
#' @seealso
#' \code{\link{pcrr}},
#' \code{\link{print.cure.pcrr}},
#' \code{\link{predict.pcrr}}
#'
#' @export
cure.pcrr <- function(object, cov, case = NULL, event = NULL, ...){
  P <- object$p
  
  if (!is.matrix(cov)) {
    if (is.vector(cov)){
      if (length(cov) %% P == 0) cov <- matrix(cov, ncol = P, byrow = TRUE)
      else stop("cov must have ", P, " column(s).")
    } else cov <- as.matrix(cov)
  }
  if (!is.numeric(cov)) stop("cov must be numeric.")
  if (ncol(cov) != P) stop("cov must have ", P, " column(s).")
  
  
  case_all <- object$case_all
  n_case_all <- nrow(case_all)
  
  avail_case <- seq_len(n_case_all)
  
  if (is.null(case)) {
    case <- avail_case
  } else {
    if (!is.numeric(case) || length(case) == 0 || any(!is.finite(case)) ||
        any(case != floor(case)) || !all(case %in% avail_case))
      stop("`case` must be one or more of: ",
           paste(avail_case, collapse = ", "),
           "\n\nThe possible model cases are as follows :\n",
           paste(object$case_model, collapse = "\n"), call. = FALSE)
    case <- as.integer(unique(case))
  }
  
  n_case <- length(case)
  

  mapping <- object$mapping
  if (is.null(event)) {
    event <- mapping
  } else {
    if (!is.numeric(event) || length(event) == 0 || !all(event %in% mapping))
      stop("`event` must be one or more of: ", paste(mapping, collapse = ", "), call. = FALSE)
    event <- unique(event)
  }
  
  

  
  model <- object$case_model
  
  
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
  
  
  tol <- 1e-12
  
  pred <- suppressWarnings(predict.pcrr(object, cov, case = case, event = event))
  t_boundary <- pred$t_boundary
  
  cure_list <- vector("list", n_case)
  status_list <- vector("list", n_case)
  
  names(cure_list) <- model[case]
  names(status_list) <- model[case]
  
  for (l in 1:n_case) {
    i <- case[l]
    coef <- object$mle_case_all[[i]]$par
    
    cure_list[[l]] <- vector("list", length(event))
    status_list[[l]] <- vector("list", length(event))
    
    names(cure_list[[l]]) <- paste0("event ", event)
    names(status_list[[l]]) <- paste0("event ", event)
  
  
    for (e in 1:length(event)) {
      k <- match(event[e], mapping)
  
      
      if (object$distribution == "gompertz2"){
        alpha <- coef[(k - 1) * (3 + P) + 1]
        rho   <- coef[(k - 1) * (3 + P) + 2]
        tau   <- coef[(k - 1) * (3 + P) + 3]
        eta   <- 0
        beta  <- coef[((k - 1) * (3 + P) + 4):(k * (3 + P))]
        
        u_Inf <- u_k(Inf, rho, tau, eta, tol)
      } else if (object$distribution == "gompertz3"){
        alpha <- coef[(k - 1) * (4 + P) + 1]
        rho   <- coef[(k - 1) * (4 + P) + 2]
        tau   <- coef[(k - 1) * (4 + P) + 3]
        eta   <- coef[(k - 1) * (4 + P) + 4]
        beta  <- coef[((k - 1) * (4 + P) + 5):(k * (4 + P))]
        
        u_Inf <- u_k(Inf, rho, tau, eta, tol)
      } else if (object$distribution == "logistic"){
        alpha <- coef[(k - 1) * (4 + P) + 1]
        b   <- coef[(k - 1) * (4 + P) + 2]
        c   <- coef[(k - 1) * (4 + P) + 3]
        p   <- coef[(k - 1) * (4 + P) + 4]
        beta  <- coef[((k - 1) * (4 + P) + 5):(k * (4 + P))]
        
        u_Inf <- u_k(Inf, b, c, p, tol)
      }

      
      
      cure <- rep(NA_real_, nrow(cov))
      status <- character(nrow(cov))
      
      for (j in seq_len(nrow(cov))){
        ezb <- exp(sum(cov[j, ] * beta))
        
        if (is.finite(t_boundary[[l]][[e]][j])){
          status[j] <- "no cure (finite support)"
          cure[j] <- NA_real_
        } else{
          if (is.finite(u_Inf)){
            if (abs(alpha) < 1e-8) cure[j] <- exp(-ezb * u_Inf)
            else cure[j] <- (1.0 + alpha * ezb * u_Inf)^(-1.0 / alpha)
            status[j] <- "cure"
          } else {
            cure[j] <- 0
            status[j] <- "no cure (asymptotic)"
          }
        }
      }
      
      
      cure_list[[l]][[e]] <- unname(cure)
      status_list[[l]][[e]] <- unname(status)
      t_boundary[[l]][[e]] <- unname(t_boundary[[l]][[e]])
    }
  }


  result <- list(
    cure = cure_list,
    status = status_list,
    t_boundary = t_boundary,
    case_model = object$case_model,
    case = case,
    event = event,
    converged = pred$converged,
    message = pred$message
  )
  
  class(result) <- "cure.pcrr"
  
  result
}


#' Print Estimated Cure Fractions
#'
#' Prints a summary of estimated cure fractions produced by
#' \code{\link{cure.pcrr}}, for each selected case and event: the estimated
#' cure fraction of each profile, its classification, and the boundary time
#' of the cumulative incidence function when a finite boundary exists. A case
#' that did not converge is flagged before its output.
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
  
  cat("Cure Fraction")
  
  for (l in 1:length(x$case)) {
    i <- x$case[l]
    cat("\n\n========================================\n")
    cat(x$case_model[i])
    cat("\n========================================\n")
    
    if (!x$converged[[l]]) {
      cat("\nconvergence : ", x$converged[[l]], "\n",
          "[", x$message[[l]], "]\n",
          "The cure fractions below are computed at these parameter values ",
          "and should not be interpreted.\n", sep = "")
    }
    
    for (e in 1:length(x$event)) {
      cat("\n[event : ", x$event[e], "]\n\n", sep = "")
      
      cure_le <- x$cure[[l]][[e]]
      status_le <- x$status[[l]][[e]]
      t_boundary_le <- x$t_boundary[[l]][[e]]

      width <- digits + 4
      for (j in seq_along(cure_le)) {
        
        cure_txt <- if (is.na(cure_le[j])) {
          sprintf(paste0("%-", width, "s"), "NA")
        } else {
          sprintf(paste0("%-", width, ".", digits, "f"), cure_le[j])
        }
        
        if (is.finite(t_boundary_le[j])) {
          interval_txt <- sprintf("[0, %s]", format(round(t_boundary_le[j], digits)))
        } else {
          interval_txt <- "[0, Inf)"
        }
        
        cat("  obs ", j, " : ", cure_txt, "<", status_le[j], ">  ", interval_txt, "\n", sep = "")
      }
      
      if (e != length(x$event)){
        cat("\n----------------------------------------\n")
      }
    }
  }
  
  invisible(x)
}

#' Parameter Vector Order
#'
#' This section describes the order of parameters in the parameter vector
#' used by the \code{init} argument of \code{pcrr} and the \code{coef}
#' component of objects of class \code{pcrr}.
#'
#' The parameter vector is generally structured as follows:
#'
#' \code{c(alpha1, parameters of u1(t), beta11, beta12, ...,
#'        alpha2, parameters of u2(t), beta21, beta22, ...)}
#'
#' Here, \code{alpha_k} denotes the parameter of the GOR transformation,
#' and \code{beta_k} denotes the regression coefficients for the covariates.
#' The coefficients are ordered as \code{beta_k1, beta_k2, ..., beta_kP}
#' and follow the column order of the input covariate matrix.
#' 
#' The blocks follow the internal cause order rather than the order in which
#' the codes appear in the data: \code{failcode} occupies the first block and
#' the remaining event codes follow in ascending order. This ordering is
#' recorded in the \code{mapping} component of a fitted object, which is worth
#' inspecting before writing an \code{init} vector by hand. The censoring code
#' has no block.
#' 
#' \subsection{Gompertz2}{
#'
#' For \code{distribution = "gompertz2"}, the parameters of \eqn{u_k(t)}
#' are ordered as \eqn{(\rho_k, \tau_k)}. Therefore, the parameter vector is
#' interpreted as:
#'
#' \code{c(alpha1, rho1, tau1, beta11, beta12, ...,
#'        alpha2, rho2, tau2, beta21, beta22, ...)}
#'
#' The expected length of the parameter vector is \code{K * (3 + P)},
#' where \code{K} is the number of event types, excluding censoring events,
#' and \code{P} is the number of covariates.
#' }
#'
#' \subsection{Gompertz3}{
#'
#' For \code{distribution = "gompertz3"}, the parameters of \eqn{u_k(t)}
#' are ordered as \eqn{(\rho_k, \tau_k, \eta_k)}. Therefore, the parameter
#' vector is interpreted as:
#'
#' \code{c(alpha1, rho1, tau1, eta1, beta11, beta12, ...,
#'        alpha2, rho2, tau2, eta2, beta21, beta22, ...)}
#'
#' The expected length of the parameter vector is \code{K * (4 + P)}.
#' }
#'
#' \subsection{Logistic}{
#'
#' For \code{distribution = "logistic"}, the parameters of \eqn{u_k(t)}
#' are ordered as \eqn{(b_k, c_k, p_k)}. Therefore, the parameter vector is
#' interpreted as:
#'
#' \code{c(alpha1, b1, c1, p1, beta11, beta12, ...,
#'        alpha2, b2, c2, p2, beta21, beta22, ...)}
#'
#' The expected length of the parameter vector is \code{K * (4 + P)}.
#' Here, \code{p} and \code{P} denote different quantities: \code{p}
#' is a parameter of the modified logistic baseline function, whereas
#' \code{P} denotes the number of covariates.
#' }
#'
#' The functions in \pkg{pcmprsk} validate the parameter vector based on
#' its length and interpret its elements according to the ordering
#' described above.
#'
#' @seealso
#' \code{\link{pcrr}}
#'
#' @name pcrr-parameter-order
#' @rdname pcrr-parameter-order
NULL



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
.estimate_mle_gom2 <- function(x, delta, z, theta_init, gtol, maxiter, fixed_alpha){

  K <- ncol(delta)
  P <- ncol(z)
  
  # Estimate
  lower <- rep(-Inf, K * (3 + P))
  upper <- rep(Inf, K * (3 + P))
  for (k in 1:K){
    lower[(k - 1) * (3 + P) + 3] <- 1e-8   # tau > 0
  }
  
  if (!is.null(fixed_alpha)) {
    for (k in 1:K){
      if (fixed_alpha[k] == 0 || fixed_alpha[k] == 1) {
        theta_init[(k - 1) * (3 + P) + 1] <- fixed_alpha[k]
        lower[(k - 1) * (3 + P) + 1] <- fixed_alpha[k]
        upper[(k - 1) * (3 + P) + 1] <- fixed_alpha[k]
      }
    }
  }
  
  control <- list(rel.tol = gtol, iter.max  = maxiter, eval.max = 3 * maxiter)
  mle <- suppressWarnings(nlminb(start = theta_init, objective = .log_lik_gom2,
                                 x = x, delta = delta, z = z, control = control, lower = lower, upper = upper))
  
  mle
}
.score_hessian_gom2 <- function(x, delta, z, theta_mle, variance, fixed_alpha){

  K <- ncol(delta)
  P <- ncol(z)
  
  # symbol name : alpha1, rho1, tau1, beta11, beta12, ..., alpha2, ...
  # symbol name : delta1, delta2, ...
  # symbol name : z1, z2, ...
  
  tol <- 1e-4
  # Make Expression Text
  txt_zb <- character(K)
  txt_u <- character(K)
  for (k in 1:K){
    txt_zb[k] <- "("
    for (p in 1:P){
      txt_zb[k] <- paste0(txt_zb[k], "z", p, " * beta", k, p)
      if (p != P){
        txt_zb[k] <- paste0(txt_zb[k], " + ")
      } else {
        txt_zb[k] <- paste0(txt_zb[k], ")")
      }
    }
    
    if (abs(theta_mle[(k - 1) * (3 + P) + 2] * max(x)) < tol) {
      txt_rho_term_sequence <- paste0("(x + rho", k, " * x^2 / 2 + rho", k, "^2 * x^3 / 6 + rho", k, "^3 * x^4 / 24)")
      txt_u[k] <- paste0("(tau", k, " * ", txt_rho_term_sequence, ")")
    } else {
      txt_u[k] <- paste0("(tau", k, " * expm1(rho", k, " * x) / rho", k, ")")
    }
  }
  
  txt_log_term <- character(K)
  txt_F <- character(K)
  txt_sum_F <- ""
  for (k in 1:K){
    
    if (is.null(fixed_alpha)){
      if (abs(theta_mle[(k - 1) * (3 + P) + 1]) < tol) {
        txt_ezb_u <- paste0("(exp", txt_zb[k], " * ", txt_u[k], ")")
        txt_alpha_term_sequence <- paste0("(1 - alpha", k, " * ", txt_ezb_u, " / 2",
                                          " + alpha", k, "^2 * ", txt_ezb_u, "^2 / 3",
                                          " - alpha", k, "^3 * ", txt_ezb_u, "^3 / 4)")
        txt_log_term[k] <- paste0("- (1 + alpha", k, ") * ", txt_ezb_u, " * ", txt_alpha_term_sequence)
        txt_F[k] <- paste0("1 - exp(-", txt_ezb_u, " * ", txt_alpha_term_sequence, ")")
      } else {
        txt_log_term[k] <- paste0("- (1 / alpha", k, " + 1) * ", "log1p(alpha", k, " * exp", txt_zb[k], " * ", txt_u[k], ")")
        txt_F[k] <- paste0("1 - (1 + alpha", k, " * exp", txt_zb[k], " * ", txt_u[k], ")^(-1 / alpha", k, ")")
      }
      
    } else {
      if (fixed_alpha[k]) { # alpha = 1
        txt_log_term[k] <- paste0("- 2 * ", "log1p(exp", txt_zb[k], " * ", txt_u[k], ")")
        txt_F[k] <- paste0("1 - (1 + exp", txt_zb[k], " * ", txt_u[k], ")^(-1)")
      } else { # alpha = 0
        txt_log_term[k] <- paste0("- exp", txt_zb[k], " * ", txt_u[k])
        txt_F[k] <- paste0("1 - exp(-exp", txt_zb[k], " * ", txt_u[k], ")")
      }
    }
    
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
  
  if (!is.null(fixed_alpha)) {
    idx_alpha <- (seq_len(K) - 1) * (3 + P) + 1
    is_alpha <- logical(n_param)
    is_alpha[idx_alpha] <- TRUE
    
    param_names <- param_names[!is_alpha]
    
    theta_mle <- theta_mle[!is_alpha]
    
    n_param <- K * (2 + P)
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
    m2 <- .estimate_mle_gom2(x, delta, z, t2, 1e-6, 200, NULL)
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
                                     lower = c(-Inf, 1e-8), x = x, delta = delta, k_num = k))
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
.estimate_mle_gom3 <- function(x, delta, z, theta_init, gtol, maxiter, fixed_alpha){

  K <- ncol(delta)
  P <- ncol(z)
  
  lower <- rep(-Inf, K * (4 + P))
  upper <- rep(Inf, K * (4 + P))
  for (k in 1:K){
    lower[(k - 1) * (4 + P) + 3] <- 1e-8   # tau > 0
  }
  
  if (!is.null(fixed_alpha)) {
    for (k in 1:K){
      if (fixed_alpha[k] == 0 || fixed_alpha[k] == 1) {
        theta_init[(k - 1) * (4 + P) + 1] <- fixed_alpha[k]
        lower[(k - 1) * (4 + P) + 1] <- fixed_alpha[k]
        upper[(k - 1) * (4 + P) + 1] <- fixed_alpha[k]
      }
    }
  }
  
  # Coarse control
  control_coarse <- list(rel.tol = 1e-4, iter.max = 50, eval.max = 150)
  # Fine control
  control_fine <- list(rel.tol = gtol, iter.max = maxiter, eval.max = 3 * maxiter)
  
  fit_coarse <- function(st){
    suppressWarnings(nlminb(start = st, objective = .log_lik_gom3,
                            x = x, delta = delta, z = z,
                            lower = lower, upper = upper, control = control_coarse))
  }
  fit_fine <- function(st){
    suppressWarnings(nlminb(start = st, objective = .log_lik_gom3,
                            x = x, delta = delta, z = z,
                            lower = lower, upper = upper, control = control_fine))
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
.score_hessian_gom3 <- function(x, delta, z, theta_mle, variance, fixed_alpha){

  K <- ncol(delta)
  P <- ncol(z)
  
  # symbol name : alpha1, rho1, tau1, eta1, beta11, beta12, ..., alpha2, ...
  # symbol name : delta1, delta2, ...
  # symbol name : z1, z2, ...
  
  tol <- 1e-4
  # Make Expression Text
  txt_zb <- character(K)
  txt_u <- character(K)
  for (k in 1:K){
    txt_zb[k] <- "("
    for (p in 1:P){
      txt_zb[k] <- paste0(txt_zb[k], "z", p, " * beta", k, p)
      if (p != P){
        txt_zb[k] <- paste0(txt_zb[k], " + ")
      } else {
        txt_zb[k] <- paste0(txt_zb[k], ")")
      }
    }
    
    
    if (abs(theta_mle[(k - 1) * (4 + P) + 2] * max(x)) < tol) {
      txt_rho_term_sequence <- paste0("(x",
                                      " + (1 + eta", k, ") * rho", k, " * x^2 / 2",
                                      " + (1 + 3 * eta", k, " + eta", k, "^2) * rho", k, "^2 * x^3 / 6",
                                      " + (1 + 7 * eta", k, " + 6 * eta", k, "^2 + eta", k, "^3) * rho", k, "^3 * x^4 / 24)")
      txt_u[k] <- paste0("(tau", k, " * exp(eta", k, ") * ", txt_rho_term_sequence, ")")
    } else if (abs(theta_mle[(k - 1) * (4 + P) + 4] * expm1(theta_mle[(k - 1) * (4 + P) + 2] * max(x))) < tol) {
      txt_A <- paste0("expm1(rho", k, " * x)")
      txt_eta_term_sequence <- paste0("(", txt_A,
                                      " + eta", k, " * ", txt_A, "^2 / 2",
                                      " + eta", k, "^2 * ", txt_A, "^3 / 6",
                                      " + eta", k, "^3 * ", txt_A, "^4 / 24)")
      txt_u[k] <- paste0("(tau", k, " * exp(eta", k, ") / rho", k, " * ", txt_eta_term_sequence, ")")
    } else {
      txt_u[k] <- paste0("(tau", k, " * exp(eta", k, ") * expm1(eta", k, " * expm1(rho", k, " * x)) / (rho", k, " * eta", k, "))")
    }
    
  }
  txt_log_term <- character(K)
  txt_F <- character(K)
  txt_sum_F <- ""
  for (k in 1:K){
    
    if (is.null(fixed_alpha)){
      if (abs(theta_mle[(k - 1) * (4 + P) + 1]) < tol) {
        txt_ezb_u <- paste0("(exp", txt_zb[k], " * ", txt_u[k], ")")
        # txt_alpha_term_sequence <- paste0("(1 - alpha", k, " * ", txt_ezb_u, " / 2 + alpha", k, "^2 * ", txt_ezb_u, "^2 / 3)")
        txt_alpha_term_sequence <- paste0("(1 - alpha", k, " * ", txt_ezb_u, " / 2",
                                          " + alpha", k, "^2 * ", txt_ezb_u, "^2 / 3",
                                          " - alpha", k, "^3 * ", txt_ezb_u, "^3 / 4)")
        txt_log_term[k] <- paste0("- (1 + alpha", k, ") * ", txt_ezb_u, " * ", txt_alpha_term_sequence)
        txt_F[k] <- paste0("1 - exp(-", txt_ezb_u, " * ", txt_alpha_term_sequence, ")")
      } else {
        txt_log_term[k] <- paste0("- (1 / alpha", k, " + 1) * ", "log1p(alpha", k, " * exp", txt_zb[k], " * ", txt_u[k], ")")
        txt_F[k] <- paste0("1 - (1 + alpha", k, " * exp(", txt_zb[k], ") * ", txt_u[k], ")^(-1 / alpha", k, ")")
      }
    } else {
      if (fixed_alpha[k]) { # alpha = 1
        txt_log_term[k] <- paste0("- 2 * ", "log1p(exp", txt_zb[k], " * ", txt_u[k], ")")
        txt_F[k] <- paste0("1 - (1 + exp", txt_zb[k], " * ", txt_u[k], ")^(-1)")
      } else { # alpha = 0
        txt_log_term[k] <- paste0("- exp", txt_zb[k], " * ", txt_u[k])
        txt_F[k] <- paste0("1 - exp(-exp", txt_zb[k], " * ", txt_u[k], ")")
      }
    }
    
    
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
  
  
  if (!is.null(fixed_alpha)) {
    idx_alpha <- (seq_len(K) - 1) * (4 + P) + 1
    is_alpha <- logical(n_param)
    is_alpha[idx_alpha] <- TRUE
    
    param_names <- param_names[!is_alpha]
    
    theta_mle <- theta_mle[!is_alpha]
    
    n_param <- K * (3 + P)
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
.estimate_mle_logi <- function(x, delta, z, theta_init, gtol, maxiter, fixed_alpha){

  K <- ncol(delta)
  P <- ncol(z)
  
  # Estimate
  lower <- rep(-Inf, K * (4 + P))
  upper <- rep(Inf, K * (4 + P))
  for (k in 1:K){
    lower[(k - 1) * (4 + P) + 2] <- 1e-8   # b > 0
    lower[(k - 1) * (4 + P) + 4] <- 1e-8   # p > 0
    upper[(k - 1) * (4 + P) + 4] <- 1-1e-8 # p < 1
  }
  
  if (!is.null(fixed_alpha)) {
    for (k in 1:K){
      if (fixed_alpha[k] == 0 || fixed_alpha[k] == 1) {
        theta_init[(k - 1) * (4 + P) + 1] <- fixed_alpha[k]
        lower[(k - 1) * (4 + P) + 1] <- fixed_alpha[k]
        upper[(k - 1) * (4 + P) + 1] <- fixed_alpha[k]
      }
    }
  }
  
  control <- list(rel.tol = gtol, iter.max  = maxiter, eval.max = 3 * maxiter)
  mle <- suppressWarnings(nlminb(start = theta_init, objective = .log_lik_logi,
                                 x = x, delta = delta, z = z, control = control, lower = lower, upper = upper))
  
  mle
}
.score_hessian_logi <- function(x, delta, z, theta_mle, variance, fixed_alpha){

  K <- ncol(delta)
  P <- ncol(z)
  
  # symbol name : alpha1, b1, c1, p1, beta11, beta12, ..., alpha2, ...
  # symbol name : delta1, delta2, ...
  # symbol name : z1, z2, ...
  
  tol <- 1e-4
  # Make Expression Text
  txt_zb <- character(K)
  txt_u <- character(K)
  for (k in 1:K){
    txt_zb[k] <- "("
    for (p in 1:P){
      txt_zb[k] <- paste0(txt_zb[k], "z", p, " * beta", k, p)
      if (p != P){
        txt_zb[k] <- paste0(txt_zb[k], " + ")
      } else {
        txt_zb[k] <- paste0(txt_zb[k], ")")
      }
    }
    txt_u[k] <- paste0("(-log1p(-p", k, " + p", k, " * (exp(-b", k, " * c", k, ") + 1) / (1 + exp(b", k ," * (x - c", k, ")))))")
  }
  txt_log_term <- character(K)
  txt_F <- character(K)
  txt_sum_F <- ""
  for (k in 1:K){
    
    if (is.null(fixed_alpha)) {
      if (abs(theta_mle[(k - 1) * (4 + P) + 1]) < tol) {
        txt_ezb_u <- paste0("(exp", txt_zb[k], " * ", txt_u[k], ")")
        # txt_alpha_term_sequence <- paste0("(1 - alpha", k, " * ", txt_ezb_u, " / 2 + alpha", k, "^2 * ", txt_ezb_u, "^2 / 3)")
        txt_alpha_term_sequence <- paste0("(1 - alpha", k, " * ", txt_ezb_u, " / 2",
                                          " + alpha", k, "^2 * ", txt_ezb_u, "^2 / 3",
                                          " - alpha", k, "^3 * ", txt_ezb_u, "^3 / 4)")
        txt_log_term[k] <- paste0("- (1 + alpha", k, ") * ", txt_ezb_u, " * ", txt_alpha_term_sequence)
        txt_F[k] <- paste0("1 - exp(-", txt_ezb_u, " * ", txt_alpha_term_sequence, ")")
      } else {
        txt_log_term[k] <- paste0("- (1 / alpha", k, " + 1) * ", "log1p(alpha", k, " * exp(", txt_zb[k], ") * ", txt_u[k], ")")
        txt_F[k] <- paste0("1 - (1 + alpha", k, " * exp(", txt_zb[k], ") * ", txt_u[k], ")^(-1 / alpha", k, ")")
      }
    } else {
      if (fixed_alpha[k]) { # alpha = 1
        txt_log_term[k] <- paste0("- 2 * ", "log1p(exp", txt_zb[k], " * ", txt_u[k], ")")
        txt_F[k] <- paste0("1 - (1 + exp", txt_zb[k], " * ", txt_u[k], ")^(-1)")
      } else { # alpha = 0
        txt_log_term[k] <- paste0("- exp", txt_zb[k], " * ", txt_u[k])
        txt_F[k] <- paste0("1 - exp(-exp", txt_zb[k], " * ", txt_u[k], ")")
      }
    }
    
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
  
  
  if (!is.null(fixed_alpha)) {
    idx_alpha <- (seq_len(K) - 1) * (4 + P) + 1
    is_alpha <- logical(n_param)
    is_alpha[idx_alpha] <- TRUE
    
    param_names <- param_names[!is_alpha]
    
    theta_mle <- theta_mle[!is_alpha]
    
    n_param <- K * (3 + P)
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
  for (i1 in 1:n_param){
    for (i2 in 1:n_param){
      hessian[i1, i2] <- sum( eval( hessian_exprs[i1, i2][[1]] ) )
    }
  }
  
  list(score = score, hessian = hessian)
  
}
