
<!-- README.md is generated from README.Rmd. Please edit that file -->

# riskscores

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/riskscores)](https://CRAN.R-project.org/package=riskscores)

<!-- badges: end -->

Risk scores are sparse linear models that map an integer linear
combination of covariates to the probability of an outcome occurring.
Unlike regression models, risk score models consist of integer
coefficients for often dichotomous variables. This allows risk score
predictions to be easily computed by adding or subtracting a few small
numbers.

Risk scores developed heuristically by altering logistic regression
models have decreased performance, as there is a fundamental trade-off
between the model’s simplicity and its predictive accuracy. In contrast,
this package presents an optimization approach to learning risk scores,
where the constraints for sparsity and integer coefficients are
integrated into the model-fitting process, rather than implemented
afterward.

## Installation

You can install the development version of riskscores from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("hjeglinton/riskscores", build_vignettes = TRUE)
```

## Example

We’ll fit a risk score model to predict breast cancer from biopsy data.
More details can be found in the package’s vignette.

``` r
library(riskscores)

# Prepare data
y <- breastcancer[,1]
X <- as.matrix(breastcancer[,-1])

# Fit risk score model 
mod <- risk_mod(X, y, lambda0 = 0.058)
```

The integer risk score model can be viewed by calling `mod$model_card`.
An individual’s risk score can be calculated by multiplying each
covariate response by its respective number of points and then adding
all points together. In our example below, a patient with a
ClumpThickness value of 1, a BareNuclei value of 5, and a BlandChromatin
value of 10 would receive a score of $10(1) + 7(5) + 8(10) = 125$.

|                | Points |
|:---------------|-------:|
| ClumpThickness |     10 |
| BareNuclei     |      7 |
| BlandChromatin |      8 |

Each score can then be mapped to a risk probability. The `mod$score_map`
dataframe maps an integer range of scores to their associated risk. We
can see that a patient who received a score of 125 would have a 77.9%
risk of their tissue sample being malignant.

| Score |   Risk |
|------:|-------:|
|    25 | 0.0006 |
|    50 | 0.0054 |
|    75 | 0.0446 |
|   100 | 0.2886 |
|   125 | 0.7788 |
|   150 | 0.9683 |
|   175 | 0.9962 |
|   200 | 0.9996 |
