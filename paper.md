---
title: 'sleev: An R Package for Semiparametric Likelihood Estimation with Errors in
  Variables'
tags:
  - R
  - measurement error
  - "two-phase sampling"
  - nonparametric statistics
  - "B-splines"
date: "16 September 2024"
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
subsample of study data is validated; such setting are known as
two-phase studies. The sieve maximum likelihood estimator (SMLE), which
combines the error-prone data on all records with the validated data on
a subsample, is a highly efficient and robust estimator to analyze such
data. However, given their complexity, a computationally efficient and
user-friendly tool is needed to obtain SMLEs. R package `sleev` fills
this gap by making semiparametric likelihood-based inference using SMLEs
for error-prone two-phase data in settings with binary and continuous
outcomes. Functions from this package can be used to analyze data with
error-prone responses and covariates.

# Statement of Need

Routinely collected data are being used frequently in biomedical
research, such as electronic health records. However, these data tend to
be error-prone, and careless use of these data could lead to biased
estimates and misleading research findings [@duan2016empirical]. To
avoid invalid study results, trained experts carefully verify and
extract data elements. However, it is usually only feasible to validate
data for a subset of records or variables. After validation, researchers
have error-prone pre-validation data for all records (phase one data)
and error-free validated data on a subset of records (phase two data).
Analyses aim to combine the two types of data to obtain estimates that
have low bias and are as efficient as possible.

There are several [R](https://www.r-project.org/) packages that address
measurement error[@Rlanguage]. A few examples are `augSIMEX`,
`attentuation`, `decon`, `eivtools`, `GLSME`, `mecor`, `meerva`, `mmc`,
`refitME`, and `simex`. The various R packages reflect many of the
different approaches for addressing measurement error such as regression
calibration [@wang2011deconvolution], SIMEX (i.e.,
simulation-extrapolation) [@lederer2006short], and moment-based
corrections [@nab2021mecor], to mention a few. Nearly all of these
existing R packages deal with errors in either the outcome or
covariates, but not both, and none of these packages permits
semiparametric likelihood-based inference that incorporates both the
error-prone phase-1 data and the validated phase-2 data.

SMLE is an estimator that analyzes two-phase data by combining the
error-prone data on all records with the validated data on a subsample.
SMLE is highly efficient because it uses all available data from both
phases [@tao2021efficient; @lotspeich2022efficient]. SMLE is also
robust, since it does not make any parametric assumptions on the error
model. For example, @tao2021efficient performed a set of simulations
highlighting the SMLE's robustness to different error mechanisms
including settings where the errors had non-zero mean and were
multiplicative. Moreover, SMLE allows error-prone outcome and
error-prone covariates in the same model. Still, in practice these
estimators can be difficult to implement, as they involve approximating
nuisance conditional densities using B-splines [@schumaker2007spline]
and then maximizing the semiparametric likelihood via a sophisticated EM
algorithm [@tao2017efficient]. Here we present R package `sleev` which
makes this method readily applicable for practitioners in a
user-friendly way. `sleev` integrates and extends R packages `logreg2ph`
and `TwoPhaseReg`, two primitive R packages developed with the original
methods papers [@tao2021efficient; @lotspeich2022efficient]. These two
packages lacked proper documentation and were difficult to use.
`logreg2ph` was computationally slow.

To promote the use of SMLE in two-phase data, extensive work has been
done to create `sleev`, a computationally efficient and user-friendly R
package to analyze two-phase, error-prone data. Specifically, in `sleev`
we rewrote the core algorithms of `logreg2ph` in C++ to speed up the
computation, and we unified the syntax across functions. To compare the
computational times, we set up simulations with the same code in the
vignette. The simulations included phase-1 and phase-2 sampe sizes of
2087 and 835, respectively, and we performed on a *ACCRE information
that I will include once the code finishes running*. Across 100
simulations, the previous `logreg2ph` took on averaage xx seconds with a
standard deviation of yy seconds, while the current `logisitc2ph` in
`sleev` took xx seconds with a standard deviation of yy seconds.

# SMLE for linear regression

In this section, we briefly introduce the SMLEs for linear regression.
The SMLEs for logistic regression are similar to linear regression and
described in the [package
vignette](https://github.com/dragontaoran/sleev/blob/main/vignettes/sleev_vignette.pdf).
Suppose that we want to fit a standard linear regression model for a
continuous outcome $Y$ and covariates $\mathbf{X}$:
$Y = \alpha + \boldsymbol{\beta}^{T}\mathbf{X} + \epsilon$, where
$\epsilon\sim N({0,\sigma}^{2})$. Our goal is to obtain estimates of
${\mathbf{\theta} = (\alpha,\boldsymbol{\beta}^{T},\sigma^{2})}^{T}$.
When we have error-prone data, $Y$ and $\mathbf{X}$ are unobserved
except for a subset of validated records. For the majority of
unvalidated records, only the error-prone outcome $Y^{*} = Y + W$ and
covariates $\mathbf{X}^{\mathbf{*}} = \mathbf{X} + \mathbf{U}$ are
observed in place of $Y$ and $\mathbf{X}$, where $W$ and $\mathbf{U}$
are the errors for the outcome and covariates, respectively. We assume
that $W$ and $\mathbf{U}$ are independent of $\epsilon$ . With potential
errors in our data, a naive regression analysis using error-prone
variables $Y^{*}$ and $\mathbf{X}^{\mathbf{*}}$ could render misleading
results.

We assume that the joint density of the complete data
$\left( Y^{*},\mathbf{X}^{*},W,\mathbf{U} \right)$ takes the form

$$P\left( Y^{*},\mathbf{X}^{*},\ W,\ \mathbf{U} \right) = P\left( Y^{*}|\mathbf{X}^{*},\ W,\ \mathbf{U} \right)P\left( W,\ \mathbf{U|}\mathbf{X}^{*}\  \right)P\left( \mathbf{X}^{*} \right)$$

$$= P_{\mathbf{\theta}}\left( Y|X \right)P\left( W,\ \mathbf{U|}\mathbf{X}^{*} \right)P\left( \mathbf{X}^{*} \right)$$

where $P( \cdot )$ and $P\left( \cdot | \cdot \right)$ denote density
and conditional density functions, respectively.
$P_{\mathbf{\theta}}\left( Y|\mathbf{X} \right)$ then refers to the
conditional density function of
$Y = \alpha + \boldsymbol{\beta}^{T}\mathbf{X} + \epsilon$. Denote the
validation indicator variable by $V$, with $V = 1$ indicating that a
record was validated and $V = 0$ otherwise. For records with $V = 0$,
their measurement errors $\left( W,\mathbf{U} \right)$ are missing, and
therefore their contributions to the log-likelihood can be obtained by
integrating out $W$ and $\mathbf{U}$. Let
$\left( Y_{i}^{*},\mathbf{X}_{i}^{\mathbf{*}},W_{i},\mathbf{U}_{i},V_{i},Y_{i},\mathbf{X}_{i} \right)$
for $i = 1,\ldots,n$ denote independent and identically distributed
realizations of
$\left( Y^{*},\mathbf{X}^{\mathbf{*}},W,\mathbf{U},V,Y,\mathbf{X} \right)$
in a sample of $n$ subjects. Then, the observed-data log-likelihood
takes the form

$$\sum_{i = 1}^{n}{V_{i}\{\log P_{\mathbf{\theta}}\left( Y_{i} \middle| \mathbf{X}_{i} \right) + \log{P(W_{i},\mathbf{U}_{i}|\mathbf{X}_{i}^{*})\}}}$$

$$+ \sum_{i = 1}^{n}{\left( 1 - V_{i} \right)\log\left\{ \int\int P_{\mathbf{\theta}}\left( Y_{i}^{*} - w|\mathbf{X}_{i}^{*} - \mathbf{u} \right)P\left( w,\mathbf{u} \middle| \mathbf{X}_{i}^{\mathbf{*}} \right)dwd\mathbf{u} \right\}}\ \ \ \ \ \ \ \ (1)$$

where $P(\pmb{X^*})$ is left out, because the error-prone covariates are
fully observed and thus $P(\pmb{X^*})$ can simply be estimated
empirically. We estimate the unknown measurement error model,
$P\left( W_{i},\mathbf{U}_{i}|\mathbf{X}_{i}^{\mathbf{*}} \right)$ using
B-spline sieves. Specifically, we approximate
$P\left( w,\mathbf{u}|\mathbf{X}_{i}^{\mathbf{*}} \right)$ and
$\log P\left( W_{i},\mathbf{U}_{i}|\mathbf{X}_{i}^{\mathbf{*}} \right)$
by
$\sum_{k = 1}^{m}I\left( w = w_{k},\mathbf{u} = \mathbf{u}_{k} \right)\sum_{j = 1}^{s_{n}}B_{j}^{q}\left( \mathbf{X}_{i}^{\mathbf{*}} \right)p_{kj}$
and
$\sum_{k = 1}^{m}I\left( W_{i} = w_{k},\mathbf{U}_{i} = \mathbf{u}_{k} \right)\sum_{j = 1}^{s_{n}}B_{j}^{q}\left( \mathbf{X}_{i}^{\mathbf{*}} \right)\log p_{kj}$,
respectively.
$\left\{ \left( w_{1},\mathbf{u}_{1} \right),...,\left( w_{m},\mathbf{u}_{m} \right) \right\}$
are the $m$ distinct observed $\left( W,\mathbf{U} \right)$ values,
$B_{j}^{q}\left( \mathbf{X}_{i}^{\mathbf{*}} \right)$ is the $j$th
B-spline basis function of order $q$ evaluated at
$\mathbf{X}_{i}^{\mathbf{*}}$, $s_{n}$ is the dimension of the B-spline
basis, and $p_{kj}$ is the coefficient associated with
$B_{j}^{q}\left( \mathbf{X}_{i}^{\mathbf{*}} \right)$ and
$\left( w_{k},\mathbf{u}_{k} \right)$. The log-likelihood in expression
(1) is now approximated by

$$\sum_{i = 1}^{n}V_{i}\left\{\log P_{\mathbf{\theta}}\left(Y_{i} \middle| \mathbf{X}_{i} \right) +\sum_{k = 1}^{m}\big\{{I(W_{i}=w_{k},\mathbf{U}_{i}= \mathbf{u}_{k})\sum_{j=1}^{s_{n}}{B_{j}^{q}(\mathbf{X}_{i}^{*})}}\log{p_{kj}}\big\}\right\}$$

$$+ \sum_{i = 1}^{n}{\left( 1 - V_{i} \right)\log\left\{ \sum_{k = 1}^{m}\{{P_{\mathbf{\theta}}\left( Y_{i}^{*} - w_{k}|\mathbf{X}_{i}^{*} - \mathbf{u}_{k} \right)}\sum_{j = 1}^{s_{n}}{B_{j}^{q}(\mathbf{X}_{i}^{*})}\log{p_{kj}\}} \right\}}.\ \ \ \ \ \ \ (2)$$

The maximization of expression (2) is carried out through an EM
algorithm to find the SMLEs $\widehat{\mathbf{\theta}}$ and
${\widehat{p}}_{kj}$. The covariance matrix of the SMLE
$\widehat{\mathbf{\theta}}$ is obtained through the method of profile
likelihood [@murphy2000profile]. Full details on the SMLE method for
logistic regression with error-prone data, including its theoretical
properties, can be found in [@lotspeich2022efficient].

# Functionalities of the `sleev` R package

The `sleev` package provides a user-friendly way to obtain the SMLEs and
their standard errors. The package can be installed from
[CRAN](https://cran.r-project.org/web/packages/sleev/index.html) or
[GitHub](https://github.com/dragontaoran/sleev). The `sleev` package
includes two main functions: linear2ph() and logistic2ph(), to fit
linear and logistic regressions, respectively, under two-phase sampling
with an error-prone outcome and covariates. The input arguments are
similar for the two functions and listed in Table 1. From Table 1, we
see that in addition to the arguments for error-prone and error-free
outcome and covariates, the user needs to specify the B-spline matrix
$B_{j}^{q}\left( \mathbf{X}_{i}^{*} \right)$ to be used in the
estimation of the error densities.

Table 1: Main arguments for `linear2ph()` and `logistic2ph()`

| Argument | Description                                                                                                                           |
|:---------|:-------------------------------------------------------------|
| y_unval  | Column name of unvalidated outcome in the input dataset.                                                                              |
| y        | Column name of validated outcome in the input dataset. NAs in the input will be counted as individuals not selected in phase two.     |
| x_unval  | Column names of unvalidated covariates in the input dataset.                                                                          |
| x        | Column names of validated covariates in the input dataset. NAs in the input will be counted as individuals not selected in phase two. |
| z        | Column names of error-free covariates in the input dataset.                                                                           |
| b_spline | Column names of the B-spline basis in the input dataset.                                                                              |
| data     | Dataset                                                                                                                               |
| hn_scale | Scale of the perturbation constant in the variance estimation via the method of profile likelihood. The default is 1.                 |
| se       | Standard errors of the parameter estimators will be estimated when set to TRUE. The default is TRUE.                                  |
| tol      | Convergence criterion. The default is 0.0001.                                                                                         |
| max_iter | Maximum number of iterations in the EM algorithm. The default is 1000.                                                                |
| verbose  | Print analysis details when set to TRUE. The default is FALSE.                                                                        |

The `sleev` package includes a dataset constructed to mimic data from
the Vanderbilt Comprehensive Care Clinic (VCCC) patient records from
[@giganti2020accounting]. Table 2 lists the variables in this dataset.

Table 2: Data dictionary for `mock.vccc`

| Name      | Status      | Type       | Description                                                                          |
|:----------|:------------|:-----------|:-------------------------------------|
| ID        | error-free  |            | Patient ID                                                                           |
| VL_unval  | error-prone | continuous | Viral load (VL) at antiretroviral therapy                                            |
| VL_val    | validated   |            | (ART) initiation                                                                     |
| ADE_unval | error-prone | binary     | Had an AIDS-defining event (ADE) within                                              |
| ADE_val   | validated   |            | one year of ART initiation: 1 - yes, 0 -- no                                         |
| CD4_unval | error-prone | continuous | CD4 count at ART initiation                                                          |
| CD4_val   | validated   |            |                                                                                      |
| prior_ART | error-free  | binary     | Experienced ART before enrollment: 1 - yes, 0 - no                                   |
| Sex       | error-free  | binary     | Sex at birth of patient: 1 - male, 0 - female                                        |
| Age       | error-free  | continuous | Age of patient                                                                       |

# Example: Case study with mock data

We now illustrate how to obtain SMLEs using the `sleev` package with
dataset `mock.vccc`. Specifically, we show how to fit a linear
regression model in the presence of errors in both the outcome and
covariates using the linear2ph() function. Situations with more
covariates and examples with logistic regression are included in the
[package
vignette](https://github.com/dragontaoran/sleev/blob/main/vignettes/sleev_vignette.pdf).

This example fits a linear regression model with CD4 count at
antiretroviral therapy (ART) initiation regressed on viral load (VL) at
ART initiation, adjusting for sex at birth. Both CD4 and VL are
error-prone, partially validated variables, whereas sex is error-free.

Because of skewness, we often transform both CD4 and VL. In our
analysis, CD4 was divided by 10 and square root transformed and VL was
$\log_{10}$ transformed:

```         
library("sleev")
data("mock.vccc")
mock.vccc$CD4_val_sq10 <- sqrt(mock.vccc$CD4_val / 10)
mock.vccc$CD4_unval_sq10 <- sqrt(mock.vccc$CD4_unval / 10)
mock.vccc$VL_val_l10 <- log10(mock.vccc$VL_val)
mock.vccc$VL_unval_l10 <- log10(mock.vccc$VL_unval)
```

To obtain the SMLEs, we first need to set up the B-spline basis for the
covariates `VL_unval_l10` (the transformed error-prone VL variable from
phase one). The `spline2ph()` function in `sleev` packages can set up
the B-spline basis, and combine it with the data input for the final
analysis. Here, we use a cubic B-spline basis with the `degree=3`
argument. The size of the basis $s_{n}$ is set to be 20, specified
through the `size = 20` argument. More details regarding order and size
selection of B-spline basis are discussed in vignette. The B-spline
basis is set up separately for the two `Sex` groups, and the size of the
B-spline basis is assigned in proportion to the relative size of the two
`Sex` groups. This is specified by argument `group`. This allows the
errors in `VL_unval_l10` to be heterogeneous between males and females.
The described B-spline basis is constructed as follows.

```         
sn <- 20
b_spline_names <- paste0("bs", 1:sn)
data.linear <- spline2ph(x = "VL_unval_l10", data = mock.vccc, size = sn,
                         degree = 3,  bs_names = b_spline_names, group = "Sex")
```

Alternatively, if the investigator has prior knowledge that the errors
in `VL_unval_l10` are likely to be homogeneous, one may fit a simpler
model by not stratifying the B-spline basis by `Sex`.

The SMLEs can be obtained by running function `linear2ph()`, as shown in
the code below. Again, a list explaining the inputs are shown in
Table 1. The fitted SMLEs are stored in a list object of class
`linear2ph`. Here, we assign the fitted SMLEs to the variable name
`res_linear`. The list of class `linear2ph` contains five components:
`coefficient`, `covariance`, `sigma`, `converge`, and `converge_cov`.

```         
start.time <- Sys.time()
res_linear <- linear2ph(y_unval = "CD4_unval_sq10", y = "CD4_val_sq10",
                        x_unval = "VL_unval_l10", x = "VL_val_l10", z = "Sex",
                        b_spline = b_spline_names, data = data.linear,
                        hn_scale = 1, se = TRUE, tol = 1e-04, 
                        max_iter = 1000, verbose = FALSE)
paste0("Run time: ", round(difftime(Sys.time(), start.time, 
                                    units = "secs"), 3), " sec")

[1] "Run time: 1.78 sec"
```

We should first check if the EM algorithms for estimating the regression
coefficients and their covariance matrix converged by the `print()` for
class `linear2ph` directly:

```         
> res_linear
This model has converged.
Coefficients:
 Intercept VL_val_l10        Sex 
     4.821     -0.141      0.273 
```

The `summary()` function for list of class `linear2ph()` returns the
estimated coefficients, their standard errors, test statistics, and
p-values as well as their covariance as follows:

```         
> summary(res_linear)
Model Summary:
Coefficients:
           Estimate     SE Statistic  p-value
Intercept     4.821 0.1587     30.39 0.000000
VL_val_l10   -0.141 0.0398     -3.55 0.000389
Sex           0.273 0.1089      2.51 0.012229
```

# Acknowledgement

This research was supported by the National Institute of Health grants
R01AI131771, R01HL094786, and P30AI110527 and the 2022 Biostatistics
Faculty Development Award from the Department of Biostatistics at
Vanderbilt University Medical Center.
