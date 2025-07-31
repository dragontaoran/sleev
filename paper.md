---
title: 'sleev: an R Package for Semiparametric Likelihood Estimation with Errors in Variables'
tags:
  - R
  - measurement error
  - "two-phase sampling"
  - nonparametric statistics
  - "B-splines"
date: "24 April 2025"
output: pdf_document
bibliography: paper.bib
affiliations:
  - name: Department of Biostatistics, Vanderbilt University Medical Center, USA
    index: 1
  - name: Department of Statistical Sciences, Wake Forest University, USA
    index: 2
  - name: Brigham Young University, USA
    index: 3
  - name: Vanderbilt Genetics Institute, Vanderbilt University Medical Center, USA
    index: 4
authors:
  - name: Jiangmei Xiong
    corresponding: yes
    orcid: 0000-0002-7481-765X
    affiliation: 1
  - name: Sarah C. Lotspeich
    affiliation: 2
    orcid: 0000-0001-5380-2427
  - name: Joey B. Sherrill
    affiliation: 3
    orcid: 0009-0002-2741-0475
  - name: Gustavo Amorim
    affiliation: 1
    orcid: 0000-0002-2941-5360
  - name: Bryan E. Shepherd
    affiliation: 1
    orcid: 0000-0002-3758-5992
  - name: Ran Tao
    affiliation: '1,4'
    orcid: 0000-0002-1106-2923
editor_options: 
  markdown: 
    wrap: 72
---

# Summary

Data with measurement error in the outcome, covariates, or both are not
uncommon, particularly with the increased use of routinely collected
data for biomedical research. With error-prone data, often only a
subsample of study data is validated; such settings are known as
two-phase studies. The sieve maximum likelihood estimator (SMLE), which
combines the error-prone data on all records with the validated data on
a subsample, is a highly efficient and robust method to analyze such
data. However, given their complexity, a computationally efficient and
user-friendly tool is needed to obtain the SMLEs. The R package `sleev`
fills this gap by making semiparametric likelihood-based inference using
the SMLEs for error-prone two-phase data in settings with binary and
continuous outcomes. Functions from this package can be used to analyze
data with error-prone binary or continuous responses and error-prone
covariates.

# Statement of Need

Routinely collected data are being used frequently in biomedical
research, such as electronic health records. However, these data tend to
be error-prone, and using these data without correcting for their
error-prone nature could lead to biased estimates and misleading
research findings [@duan2016empirical]. To avoid such invalid study
results, trained experts carefully verify and extract data elements.
However, it is usually only feasible to validate data for a subset of
records or variables. After validation, researchers have error-prone,
pre-validation data for all records (phase one) and error-free validated
data on a subset of records (phase two). Analyses aim to combine the two
types of data to obtain estimates that have low bias and are as robust
and efficient as possible.

There are several packages for [R](https://www.r-project.org/) [@Rlanguage] that address
measurement error, including `augSIMEX` [@augSIMEX], `attenuation` [@attenuation], `decon` [@decon],
`eivtools` [@eivtools], `GLSME` [@GLSME], `mecor` [@mecor],
`meerva` [@meerva], `mmc` [@mmc], `refitME` [@refitME], and `simex` [@SIMEX]. The various R packages reflect many different approaches, such
as regression calibration [@decon], SIMEX (i.e.,
                                           simulation-extrapolation) [@SIMEX], and moment-based corrections
[@mecor], to mention a few. Nearly all of these existing R packages deal
with errors in either the outcome or covariates, but not both, and none
of these packages permits efficient inference that
incorporates both the error-prone phase-one data and the validated
phase-two data.

The sieve maximum likelihood estimator (SMLE) is an estimator that
analyzes two-phase data by combining the error-prone data on all records
with the validated data on a subsample. By leveraging all available
data, the SMLE operates with high efficiency [@tao2021efficient;
                                              @lotspeich2022efficient]. Since it does not make any parametric
assumptions on the error model, the SMLE is also robust. For example,
@tao2021efficient performed a set of simulations highlighting the SMLE's
robustness to different error mechanisms including settings where the
errors had non-zero mean or were multiplicative. Moreover, the SMLE
allows error-prone outcome and error-prone covariates in the same model.
Still, in practice these estimators can be difficult to implement, as
they involve approximating nuisance conditional densities using
B-splines [@schumaker2007spline] and then maximizing the semiparametric
likelihood via a sophisticated EM algorithm [@tao2017efficient]. Here,
we present the R package `sleev`, which makes the SMLE readily
applicable for practitioners in a user-friendly way. `sleev` integrates
and extends primitive R packages, `logreg2ph` and `TwoPhaseReg`,
developed with the original methods papers [@tao2021efficient;
@lotspeich2022efficient]. These two packages lacked proper documentation
and were difficult to use. `logreg2ph` was also computationally slow.

To promote the use of the SMLE, extensive work has been done to create
`sleev`, a computationally efficient and user-friendly R package to
analyze two-phase, error-prone data. Specifically, in `sleev` we rewrote
the core algorithms of `logreg2ph` in C++ to speed up the computation,
and we unified the syntax across functions. To compare the computational
times, we set up simulations with the same code in the [package vignette](https://github.com/dragontaoran/sleev/blob/main/inst/article/sleev_vignette.pdf). The
simulations included phase-one and phase-two sample sizes of 2087 and
835, respectively, and were performed on a 64-bit Linux OS machine with
8G memory. Across 100 simulations, the previous `logreg2ph` took an
average of 289.44 seconds with a standard deviation of 8.83 seconds to perform the analysis, while
the corresponding new function in `sleev` only took an average of 122.32 seconds with a
standard deviation of 8.18 seconds.


# SMLE for Linear Regression

In this section, we briefly introduce the SMLE for linear regression.
Suppose that we want to fit a standard linear regression model for a
continuous outcome $Y$ and covariates $\pmb{X}$:
$Y = \alpha + \boldsymbol{\beta}^\mathrm{T}\pmb{X} + \epsilon$, where
$\epsilon\sim N({0,\sigma}^{2})$. Our goal is to obtain estimates of
${\pmb{\theta} = (\alpha,\boldsymbol{\beta}^\mathrm{T},\sigma^{2})}^\mathrm{T}$.
When we have error-prone data, $Y$ and $\pmb{X}$ are unobserved except
for a subset of validated records. For unvalidated records (the
majority), only the error-prone outcome $Y^{*} = Y + W$ and covariates
$\pmb{X}^{*} = \pmb{X} + \pmb{U}$ are observed in place of $Y$ and
$\pmb{X}$, where $W$ and $\pmb{U}$ are the errors for the outcome and
covariates, respectively. We assume that $W$ and $\pmb{U}$ are
independent of $\epsilon$ . With potential errors in our data, a naive
regression analysis using error-prone variables $Y^{*}$ and
$\pmb{X}^{*}$ could render misleading results [@fuller2009measurement].

We assume that the joint density of the complete data
$\left( Y^{*},\pmb{X}^{*},W,\pmb{U} \right)$ takes the form

$$P\left( Y^{*},\pmb{X}^{*},\ W,\ \pmb{U} \right) = P\left( Y^{*}|\pmb{X}^{*},\ W,\ \pmb{U} \right)P\left( W,\ \pmb{U|}\pmb{X}^{*}\right)P\left( \pmb{X}^{*} \right)$$

$$= P_{\pmb{\theta}}\left( Y|X \right)P\left( W,\ \pmb{U}|\pmb{X}^{*} \right)P\left( \pmb{X}^{*} \right),$$

where $P( \cdot )$ and $P\left( \cdot | \cdot \right)$ denote density
and conditional density functions, respectively. Specifically,
$P_{\pmb{\theta}}\left( Y|\pmb{X} \right)$ then refers to the
conditional density function of the linear regression model of $Y$ given $\pmb{X}$. Denote the validation indicator variable by $V$, with $V = 1$
indicating that a record was validated and $V = 0$ otherwise. For
records with $V = 0$, their measurement errors
$\left( W,\pmb{U} \right)$ are missing, and therefore their
contributions to the log-likelihood can be obtained by integrating out
$W$ and $\pmb{U}$.

Let
$\left( Y_{i}^{*},\pmb{X}_{i}^{*},W_{i},\pmb{U}_{i},V_{i},Y_{i},\pmb{X}_{i} \right)$
for $i = 1,\ldots,n$ denote independent and identically distributed
realizations of $\left( Y^{*},\pmb{X}^{*},W,\pmb{U},V,Y,\pmb{X} \right)$
in a sample of $n$ subjects. Then, the observed-data log-likelihood is
proportional to

$$\sum_{i = 1}^{n}{V_{i}\{\log P_{\pmb{\theta}}\left( Y_{i} \middle| \pmb{X}_{i} \right) + \log{P(W_{i},\pmb{U}_{i}|\pmb{X}_{i}^{*})\}}}$$

$$+ \sum_{i = 1}^{n}{\left( 1 - V_{i} \right)\log\left\{ \int\int P_{\pmb{\theta}}\left( Y_{i}^{*} - w|\pmb{X}_{i}^{*} - \pmb{u} \right)P\left( w,\pmb{u} \middle| \pmb{X}_{i}^{*} \right)dwd\pmb{u} \right\}},\ \ \ \ \ \ \ \ (1)$$

where $P(\pmb{X}^*)$ is left out, because the error-prone covariates are
fully observed and thus $P(\pmb{X}^*)$ can simply be estimated
empirically. We estimate the unknown measurement error model,
$P\left( W_{i},\pmb{U}_{i}|\pmb{X}_{i}^{*} \right)$, using B-spline
sieves. Specifically, we approximate
$P\left( w,\pmb{u}|\pmb{X}_{i}^{*} \right)$ and
$\log P\left( W_{i},\pmb{U}_{i}|\pmb{X}_{i}^{*} \right)$ by
$\sum_{k = 1}^{m}\mathrm{I}\left( w = w_{k},\pmb{u} = \pmb{u}_{k} \right)\sum_{j = 1}^{s_{n}}B_{j}^{q}\left( \pmb{X}_{i}^{*} \right)p_{kj}$
  and
$\sum_{k = 1}^{m}\mathrm{I}\left( W_{i} = w_{k},\pmb{U}_{i} = \pmb{u}_{k} \right)\sum_{j = 1}^{s_{n}}B_{j}^{q}\left( \pmb{X}_{i}^{*} \right)\log p_{kj}$,
respectively. Here,
$\{ \left( w_{1},\pmb{u}_{1} \right),...,$ $\left( w_{m},\pmb{u}_{m} \right) \}$
  are the $m$ distinct observed $\left( W,\pmb{U} \right)$ values from the
validation study, $B_{j}^{q}\left( \pmb{X}_{i}^{*} \right)$ is the $j$th
B-spline basis function of order $q$ evaluated at $\pmb{X}_{i}^{*}$,
$s_{n}$ is the dimension of the B-spline basis, and $p_{kj}$ is the
coefficient associated with $B_{j}^{q}\left( \pmb{X}_{i}^{*} \right)$
  and $\left( w_{k},\pmb{u}_{k} \right)$. The expression (1) is now
approximated by 

$$\sum_{i = 1}^{n}V_{i}\left[\log P_{\pmb{\theta}}\left(Y_{i} \middle| \pmb{X}_{i} \right) +\sum_{k = 1}^{m}\left\{{\mathrm{I}(W_{i}=w_{k},\pmb{U}_{i}= \pmb{u}_{k})\sum_{j=1}^{s_{n}}{B_{j}^{q}(\pmb{X}_{i}^{*})}}\log{p_{kj}}\right\}\right]$$
  
  $$+ \sum_{i = 1}^{n}\log{\left[\sum_{k = 1}^{m}\left\{P_{\pmb{\theta}}\left( Y_{i}^{*} - w_{k}|\pmb{X}_{i}^{*} - \pmb{u}_{k} \right)\sum_{j = 1}^{s_{n}}B_{j}^{q}(\pmb{X}_{i}^{*})\log{p_{kj}}\right\}\right]}. \ \ \ \ \ \ \ \ (2)$$
  
  The maximization of expression (2) is carried out through an EM
algorithm to find the SMLEs $\widehat{\pmb{\theta}}$ and
${\widehat{p}}_{kj}$. The covariance matrix of the SMLE
$\widehat{\pmb{\theta}}$ is obtained through the method of profile
likelihood [@murphy2000profile].

The SMLEs for logistic regression are similar to linear regression and
described in the [package vignette](https://github.com/dragontaoran/sleev/blob/main/inst/article/sleev_vignette.pdf),
and the theoretical properties can be found in @lotspeich2022efficient.

# Functionalities of the `sleev` R Package

The `sleev` package provides a user-friendly way to obtain the SMLEs and
their standard errors. The package can be installed from
[CRAN](https://cran.r-project.org/web/packages/sleev/index.html) or
[GitHub](https://github.com/dragontaoran/sleev). The `sleev` package
includes two main functions: `linear2ph()` and `logistic2ph()`, to fit
linear and logistic regressions, respectively, under two-phase sampling
with an error-prone outcome and covariates. The input arguments are
similar for the two functions and listed in Table 1. In addition to the
arguments for error-prone and error-free outcome and covariates, the
user needs to specify the B-spline matrix
$B_{j}^{q}\left( \pmb{X}_{i}^{*} \right)$ to be used in the estimation
of the error densities.

### Table 1: Main arguments for the `linear2ph()` and `logistic2ph()` functions

| Argument | Description|
|:---------|:--------------------------------------------------|
| y_unval  | Column name of unvalidated outcome in the input dataset. |
| y        | Column name of validated outcome in the input dataset. `NA`s in the input will be counted as individuals not selected in phase two.|
| x_unval  | Column names of unvalidated covariates in the input dataset. |
| x        | Column names of validated covariates in the input dataset. `NA`s in the input will be counted as individuals not selected in phase two. |
| z        | Column names of error-free covariates in the input dataset.                                                                             |
| data     | Dataset generated from `sleev::spline2ph()`.                                                                                                                        |
| hn_scale | Scale of the perturbation constant in the variance estimation via the method of profile likelihood. The default is `1`.                 |
| se       | Standard errors of the parameter estimators will be estimated when set to `TRUE`. The default is `TRUE`.                                  |
| tol      | Convergence criterion. The default is `0.0001`.                                                                                         |
| max_iter | Maximum number of iterations in the EM algorithm. The default is `1000`.                                                                |
| verbose  | Print analysis details when set to `TRUE`. The default is `FALSE`.                                                                      |

# Example: Case study with mock data

For demonstration, the `sleev` package includes a dataset constructed to
mimic data from the Vanderbilt Comprehensive Care Clinic (VCCC) patient
records from @giganti2020accounting. Table 2 describes the variables in
this dataset.

### Table 2: Data dictionary for `mock.vccc`

| Name        | Status      | Type       | Description                                         |
|:------------|:------------|:-----------|:----------------------------------------------------|
| ID          | error-free  |            | Patient ID                                          |
| VL_unval    | error-prone | continuous | Viral load (VL) at antiretroviral therapy (ART)     |
| VL_val      | validated   | continuous | initiation                                          |
| ADE_unval   | error-prone | binary     | Had an AIDS-defining event (ADE) within one         |
| ADE_val     | validated   | binary     | year of ART initiation: 1 - yes, 0 -- no            |
| CD4_unval   | error-prone | continuous | CD4 count at ART initiation                         |
| CD4_val     | validated   | continuous |                                                     |
| Prior_ART   | error-free  | binary     | Experienced ART before enrollment: 1 - yes, 0 - no  |
| Sex         | error-free  | binary     | Sex at birth of patient: 1 - male, 0 - female       |
| Age         | error-free  | continuous | Age of patient                                      |

We now illustrate how to obtain the SMLEs using the `sleev` package with
the `mock.vccc` dataset. Specifically, we show how to fit a linear
regression model in the presence of errors in both the outcome and
covariates using the `linear2ph()` function. Situations with more
covariates and examples with logistic regression are included in the
[package vignette](https://github.com/dragontaoran/sleev/blob/main/inst/article/sleev_vignette.pdf).

This example fits a linear regression model with CD4 count at
antiretroviral therapy (ART) initiation regressed on viral load (VL) at
ART initiation, adjusting for sex at birth. Both CD4 and VL are
error-prone, partially validated variables, whereas sex is error-free.
Because of skewness, we often transform both CD4 and VL. In our
analysis, CD4 was divided by 10 and square root transformed, and VL was
$\log_{10}$ transformed:

```R         
library("sleev")
data("mock.vccc")
mock.vccc$CD4_val_sq10 <- sqrt(mock.vccc$CD4_val / 10)
mock.vccc$CD4_unval_sq10 <- sqrt(mock.vccc$CD4_unval / 10)
mock.vccc$VL_val_l10 <- log10(mock.vccc$VL_val)
mock.vccc$VL_unval_l10 <- log10(mock.vccc$VL_unval)
```

To obtain the SMLEs, we first need to set up the B-spline basis for the
error-prone covariate `VL_unval_l10` (the transformed VL variable from
phase one) and `Sex`. The `spline2ph()` function in the `sleev` package
can set up the B-spline basis, and combine it with the input data for
the final analysis. Here, we use a cubic B-spline basis with the
`degree = 3` argument. The size of the basis $s_{n}$ is set to be 20,
specified through the `size = 20` argument. More details regarding order
and size selection, as well as run time comparison of B-spline basis,
are discussed in the
[package vignette](https://github.com/dragontaoran/sleev/blob/main/inst/article/sleev_vignette.pdf).
To allow possible heterogeneity in error distribution between males and
females, we can set up B-spline basis separately and proportionally for
the two `Sex` groups by specifying argument `group = "Sex"`. The
described B-spline basis is constructed as follows.

```R         
sn <- 20
data.linear <- spline2ph(x = "VL_unval_l10", data = mock.vccc, size = sn,
                         degree = 3, group = "Sex")
```

Alternatively, if the investigator has prior knowledge that the errors
in `VL_unval_l10` are likely to be homogeneous, one may fit a simpler
model by not stratifying the B-spline basis by `Sex`.

Having constructed the B-spline basis, the SMLEs can be obtained by
running the `linear2ph()` function on `data.linear`, as shown in the
code below. Again, the inputs are explained in Table 1. The fitted SMLEs
are stored in a list object of class `linear2ph`. Here, we assign the
fitted SMLEs to the variable name `res_linear`. The list of class
`linear2ph` contains five components: `coefficient`, `covariance`,
`sigma`, `converge`, and `converge_cov`.

```R         
res_linear <- linear2ph(y_unval = "CD4_unval_sq10", y = "CD4_val_sq10",
                        x_unval = "VL_unval_l10", x = "VL_val_l10", 
                        z = "Sex", data = data.linear,
                        hn_scale = 1, se = TRUE, tol = 1e-04, 
                        max_iter = 1000, verbose = FALSE)
```

We should first check if the EM algorithms for estimating the regression
coefficients and their covariance matrix converged by using the
`print()` for class `linear2ph` directly.

```R         
> res_linear

Call:
linear2ph(y_unval = "CD4_unval_sq10", y = "CD4_val_sq10",
          x_unval = "VL_unval_l10", x = "VL_val_l10", z = "Sex", 
          data = data.linear, hn_scale = 1, se = TRUE, 
          tol = 1e-04, max_iter = 1000, verbose = FALSE)

The parameter estimation has converged.

Coefficients:
 Intercept VL_val_l10        Sex 
 4.8209166 -0.1413168  0.2727984 
```

The `summary()` function for the object of class `linear2ph` returns the
estimated coefficients, their standard errors, test statistics, and
$p$-values as follows:

```R         
> summary(res_linear)

Call:
linear2ph(y_unval = "CD4_unval_sq10", y = "CD4_val_sq10",
          x_unval = "VL_unval_l10", x = "VL_val_l10", z = "Sex", 
          data = data.linear, hn_scale = 1, se = TRUE, 
          tol = 1e-04, max_iter = 1000, verbose = FALSE)

Coefficients:
             Estimate         SE Statistic      p-value
Intercept   4.8209166 0.15865204 30.386729 0.0000000000
VL_val_l10 -0.1413168 0.03983406 -3.547636 0.0003887047
Sex         0.2727984 0.10888178  2.505455 0.0122294098
```

# Acknowledgement

This research was supported by the National Institute of Health grants
R01AI131771, R01HL094786, and P30AI110527 and the 2022 Biostatistics
Faculty Development Award from the Department of Biostatistics at
Vanderbilt University Medical Center. This work leveraged the resources
provided by the Vanderbilt Advanced Computing Center for Research and
Education (ACCRE), a collaboratory operated by and for Vanderbilt
faculty.

# References
