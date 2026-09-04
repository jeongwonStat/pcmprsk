# pcmprsk: Parametric Competing Risks Regression Models for the Cumulative Incidence Function

The **pcmprsk** package provides tools for fitting parametric regression models for the **cumulative incidence function (CIF)** in the presence of competing risks.

The package implements the generalized odds rate (GOR) transformation model proposed by Jeong and Fine (2007), together with several flexible parametric baseline functions. The supported baseline models include two-parameter and three-parameter Gompertz models and a three-parameter modified logistic model.

The package provides:

- Parametric regression modeling of cumulative incidence functions
- Maximum likelihood estimation
- Statistical inference based on the observed information matrix
- Prediction of cumulative incidence functions for specified covariate values
- Estimation of long-term event probabilities and cure fractions
- Support for monotone and unimodal baseline hazard shapes

---

## Installation

### From CRAN

The current released version can be installed from CRAN using:

```R
install.packages("pcmprsk")
```

### Development version from GitHub

The development version can be installed directly from GitHub using `devtools`:

```R
devtools::install_github("jeongwonStat/pcmprsk")
```

---

## Methodology

### Generalized Odds Rate Regression Model

For competing risks data, the cumulative incidence function for cause \(k\) is modeled directly through the generalized odds rate transformation:

```math
F_k(t;\mathbf Z)
=
1-
\left\{
1+
\alpha_k
\exp(\mathbf Z^\top\boldsymbol\beta_k)
u_k(t)
\right\}^{-1/\alpha_k},
```

where:

- $F_k(t;\mathbf Z)$ is the CIF for cause $k$,
- $\mathbf Z$ is a vector of covariates,
- $\boldsymbol\beta_k$ is the regression coefficient vector,
- $\alpha_k$ is the generalized odds-rate link parameter, and
- $u_k(t)$ is a parametric baseline cumulative function.

The generalized odds-rate framework includes important special cases such as proportional hazards and proportional odds models.

Model parameters are estimated by **full maximum likelihood**.

---

## Supported Baseline Models

`pcmprsk` currently supports three parametric baseline models.

| Model | Argument | Parameters | Main feature |
|---|---|---|---|
| Two-parameter Gompertz | `"gompertz2"` | ($\rho,\tau$) | Monotone hazard shapes |
| Three-parameter Gompertz | `"gompertz3"` | ($\rho,\tau,\eta$) | Can represent unimodal hazard shapes |
| Modified logistic | `"logistic"` | (b,c,p) | Flexible sigmoidal and bounded cumulative functions |

Short aliases are also available:

```text
"gom2"  -> "gompertz2"
"gom3"  -> "gompertz3"
"logi"  -> "logistic"
```

## Basic Usage

The main function is `pcrr()`:

```R
pcrr(ftime,
     fstatus,
     cov,
     distribution = "gompertz2",
     dist = NULL,
     failcode = 1,
     cencode = 0,
     na.action = na.omit,
     gtol = 1e-6,
     maxiter = 300,
     init,
     variance = TRUE)
```

The essential inputs are:

- `ftime`: failure or censoring times
- `fstatus`: event status
- `cov`: covariate matrix
- `failcode`: event code of interest
- `cencode`: censoring code

For example:

```R
fit <- pcrr(ftime = time_vector_data,
       fstatus = event_num_vector_data,
       cov = covariate_matrix_data)
```

When `distribution` is omitted, the **two-parameter Gompertz model** (`"gompertz2"`) is used by default.

---

## Choosing a Baseline Model

### Two-parameter Gompertz

```R
fit_gom2 <- pcrr(ftime = time_vector_data,
                 fstatus = event_num_vector_data,
                 cov = covariate_matrix_data,
                 distribution = "gompertz2")
```

### Three-parameter Gompertz

```R
fit_gom3 <- pcrr(ftime = time_vector_data,
                 fstatus = event_num_vector_data,
                 cov = covariate_matrix_data,
                 distribution = "gompertz3")
```

### Modified Logistic

```R
fit_logi <- pcrr(ftime = time_vector_data,
                 fstatus = event_num_vector_data,
                 cov = covariate_matrix_data,
                 distribution = "logistic")
```

The three models provide different levels of flexibility.

---

## Development

Future development of `pcmprsk` includes:

- More numerically stable maximum likelihood optimization
- Additional parametric baseline hazard functions
- Improved diagnostics for model convergence and identifiability
- Additional model comparison tools
---

### References
Cheng, Y. (2009). Modeling Cumulative Incidences of Dementia and Dementia-Free Death Using a Novel Three-Parameter Logistic Function. The International Journal of Biostatistics, 5(1), Article 29.

Dabrowska, D. M. and Doksum, K. A. (1988). Estimation and testing in a two-sample generalized odds-rate model. Journal of the American Statistical Association, 83(403), 744–749.

Haile, S. R., Jeong, J.-H., Chen, X. and Cheng, Y. (2016). A 3-parameter Gompertz distribution for survival data with competing risks, with an application to breast cancer data. Journal of Applied Statistics, 43(12), 2239–2253.

Jeong, J.-H. and Fine, J. P. (2007). Parametric regression on the cumulative incidence function. Biostatistics, 8(2), 184–196.
