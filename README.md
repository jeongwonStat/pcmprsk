# pcmprsk: Parametric Competing Risks Regression Models for the Cumulative Incidence Function

The `pcmprsk` package provides tools for fitting parametric regression models for the **cumulative incidence function (CIF)** in the presence of competing risks.

The package implements the generalized odds rate (GOR) transformation model proposed by Jeong and Fine (2007), together with several flexible parametric baseline functions. The supported baseline models include two-parameter and three-parameter Gompertz models and a three-parameter modified logistic model.

The package provides:

- Parametric regression modeling of cumulative incidence functions
- Maximum likelihood estimation
- Statistical inference based on the observed information matrix
- Prediction of cumulative incidence functions for specified covariate values
- Estimation of long-term event probabilities and cure fractions
- Support for monotone and unimodal baseline hazard shapes


### Installation

The current version can be installed from source using the package `devtools`
```R
devtools::install_github("jeongwonStat/pcmprsk")
```

It can also be found on CRAN
```R
install.packages("pcmprsk")
```


### Usage Examples

#### - `pcrr` function 
다음은 pcrr에 포함된 전체 인자입니다. 필수 입력 인자는 ftime, fstatus, cov, fstatus에 알맞은 failcode, cencode 입니다.
```R
pcrr(ftime, fstatus, cov,
     distribution = "gompertz2", dist = NULL,
     failcode = 1, cencode = 0,
     na.action = na.omit, gtol = 1e-06, maxiter = 300,
     init, variance = TRUE)
```

따라서 기본적으로 다음과 같이 사용할 수 있습니다. 
```R
pcrr(ftime = time_vector_data,
     fstatus = event_num_vector_data,
     cov = covariate_matrix_data)
```
distribution이 생략된 경우, two-parameter Gompertz 분포 ("gompertz2")가 자동으로 선택됩니다.
곰퍼츠2 위험 분포 설명~~

```R
pcrr(ftime = time_vector_data,
     fstatus = event_num_vector_data,
     cov = covariate_matrix_data,
     distribution = "gompertz3")
```
곰퍼츠3 위험 분포 설명~~

```R
pcrr(ftime = time_vector_data,
     fstatus = event_num_vector_data,
     cov = covariate_matrix_data,
     distribution = "logistic")
```
로지스틱 위험 분포 설명~~

### Developing
- Addition of a more stable MLE optimization methodology
- Add another baseline hazard

### References
Cheng, Y. (2009). Modeling Cumulative Incidences of Dementia and Dementia-Free Death Using a Novel Three-Parameter Logistic Function. The International Journal of Biostatistics, 5(1), Article 29.

Dabrowska, D. M. and Doksum, K. A. (1988). Estimation and testing in a two-sample generalized odds-rate model. Journal of the American Statistical Association, 83(403), 744–749.

Haile, S. R., Jeong, J.-H., Chen, X. and Cheng, Y. (2016). A 3-parameter Gompertz distribution for survival data with competing risks, with an application to breast cancer data. Journal of Applied Statistics, 43(12), 2239–2253.

Jeong, J.-H. and Fine, J. P. (2007). Parametric regression on the cumulative incidence function. Biostatistics, 8(2), 184–196.
