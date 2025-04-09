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

