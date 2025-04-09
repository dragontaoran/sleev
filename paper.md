---
title: 'sleev: An R Package for Semiparametric Likelihood Estimation with Errors in Variables'
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

There are several [R](https://www.r-project.org/) packages that address
measurement error [@Rlanguage]. A few examples are
`augSIMEX`[@augSIMEX], `attenuation`[@attenuation], `decon`[@decon],
`eivtools`[@eivtools], `GLSME`[@GLSME], `mecor`[@mecor],
`meerva`[@meerva], `mmc`[@mmc], `refitME`[@refitME], and `simex`
[@SIMEX]. The various R packages reflect many different approaches, such
as regression calibration [@decon], SIMEX (i.e.,
                                           simulation-extrapolation) [@SIMEX], and moment-based corrections
[@mecor], to mention a few. Nearly all of these existing R packages deal
with errors in either the outcome or covariates, but not both, and none
of these packages permits semiparametric likelihood-based inference that
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
errors had non-zero mean and were multiplicative. Moreover, the SMLE
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
times, we set up simulations with the same code in the vignette. The
simulations included phase-one and phase-two sample sizes of 2087 and
835, respectively, and were performed on a 64-bit Linux OS machine with
8G memory. Across 100 simulations, the previous `logreg2ph` took on
average 173.19 seconds with a standard deviation of 4.40 seconds, while
the corresponding new function in `sleev` took 109.92 seconds with a
standard deviation of 7.91 seconds.


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

$$P\left( Y^{*},\pmb{X}^{*},\ W,\ \pmb{U} \right) = P\left( Y^{*}|\pmb{X}^{*},\ W,\ \pmb{U} \right)P\left( W,\ \pmb{U|}\pmb{X}^{*}\  \right)P\left( \pmb{X}^{*} \right)$$

$$= P_{\pmb{\theta}}\left( Y|X \right)P\left( W,\ \pmb{U}|\pmb{X}^{*} \right)P\left( \pmb{X}^{*} \right),$$

where $P( \cdot )$ and $P\left( \cdot | \cdot \right)$ denote density
and conditional density functions, respectively. Specifically,
$P_{\pmb{\theta}}\left( Y|\pmb{X} \right)$ then refers to the
conditional density function of the linear regression model of $Y$ and
$\pmb{X}$. Denote the validation indicator variable by $V$, with $V = 1$
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
$P\left( W_{i},\pmb{U}_{i}|\pmb{X}_{i}^{*} \right)$ using B-spline
sieves. Specifically, we approximate
$P\left( w,\pmb{u}|\pmb{X}_{i}^{*} \right)$ and
$\log P\left( W_{i},\pmb{U}_{i}|\pmb{X}_{i}^{*} \right)$ by
$\sum_{k = 1}^{m}\mathrm{I}\left( w = w_{k},\pmb{u} = \pmb{u}_{k} \right)\sum_{j = 1}^{s_{n}}B_{j}^{q}\left( \pmb{X}_{i}^{*} \right)p_{kj}$
  and
$\sum_{k = 1}^{m}\mathrm{I}\left( W_{i} = w_{k},\pmb{U}_{i} = \pmb{u}_{k} \right)\sum_{j = 1}^{s_{n}}B_{j}^{q}\left( \pmb{X}_{i}^{*} \right)\log p_{kj}$,
respectively. Here,
$\left\{ \left( w_{1},\pmb{u}_{1} \right),...,$ $\left( w_{m},\pmb{u}_{m} \right) \right\}$
  are the $m$ distinct observed $\left( W,\pmb{U} \right)$ values from the
validation study, $B_{j}^{q}\left( \pmb{X}_{i}^{*} \right)$ is the $j$th
B-spline basis function of order $q$ evaluated at $\pmb{X}_{i}^{*}$,
$s_{n}$ is the dimension of the B-spline basis, and $p_{kj}$ is the
coefficient associated with $B_{j}^{q}\left( \pmb{X}_{i}^{*} \right)$
  and $\left( w_{k},\pmb{u}_{k} \right)$. The expression (1) is now
approximated by
