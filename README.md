<p style="display:inline-block;">
  <img src="hex.png" width="200">
    <h1>SLEEV: Semiparametric Likelihood Estimation with Errors in Variables</h1>
    </p>

License: file `/inst/LICENSE_sleev.txt`

The complete R package `logreg2ph` and code for the simulation settings included in [this paper](https://doi.org/10.1111/biom.13512). 
  
Part of the R package `TwoPhaseReg` as described in [this paper](https://doi.org/10.1002/sim.8876).

If you wish to make suggestions, report a problem, or need more support with package use, please email Dr. Ran Tao (r DOT tao AT vumc DOT org) or Dr. Sarah Lotspeich (lotspes AT wfu DOT edu).
  
Below is a simple example for running a linear regression with `sleev`. The complete vignette detailing all uses of example in `sleev` can be found as an [article](https://github.com/dragontaoran/sleev/blob/main/inst/article/sleev_vignette.pdf).
  
### Install
To install the package for GitHub, run the following in your `R` console: 
    
```{r}
  devtools::install_github("dragontaoran/sleev")
```
or install from CRAN:
```
install.packages("sleev")
```
  
Load package:
    
```{r}
  library(sleev)
```
  
### `mock.vccc` data
  
This data is a simulated two-phase sampling dataset based on VCCC data. To load this dataset, run:
    
```{r}
  data(mock.vccc)
```
  
### Preprocessing
  Because of skewness, we often transform both CD4 and VL. In our analysis, CD4 was divided by 10 and square-root transformed and VL was log_10 transformed:

```
  mock.vccc$CD4_val_sq10 <- sqrt(mock.vccc$CD4_val/10)
  mock.vccc$CD4_unval_sq10 <- sqrt(mock.vccc$CD4_unval/10)
  mock.vccc$VL_val_l10 <- log10(mock.vccc$VL_val)
  mock.vccc$VL_unval_l10 <- log10(mock.vccc$VL_unval)
```
Load `splines`	package
  
```
  library(splines)
```
  
### Construct B-spline Basis
  
```
sn=20
data.linear <- spline2ph(x = "VL_unval_l10", data = mock.vccc, size = sn,
                         degree = 3,  group = "Sex")
```
  
### Model fitting 
  
The SMLEs can be obtained by running
```
res_linear <- linear2ph(y_unval = "CD4_unval_sq10", y = "CD4_val_sq10", 
                        x_unval = "VL_unval_l10",x = "VL_val_l10",z = "Sex", 
                        b_spline = paste0("bs", 1:sn), data = data.linear,
                        hn_scale = 1, se = TRUE, tol = 1e-04, 
                        max_iter = 1000, verbose = FALSE)
```
Check convergence:
    
```
c(res_linear$converge, res_linear$converge_cov)
```
  
