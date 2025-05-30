
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
mod <- risk_mod(X, y, lambda = 0.0392)
```

The integer risk score model can be viewed by calling `mod$model_card`.
An individual’s risk score can be calculated by multiplying each
covariate response by its respective number of points and then adding
all points together. In our example below, a patient with a
ClumpThickness value of 5, a BareNuclei value of 1, and a BlandChromatin
value of 3 would receive a score of $9(5) + 7(1) + 8(3) = 76$.

|                | Points |
|:---------------|-------:|
| ClumpThickness |      9 |
| BareNuclei     |      7 |
| BlandChromatin |      8 |

Each score can then be mapped to a risk probability. The `mod$score_map`
dataframe maps an integer range of scores to their associated risk. We
can see that a patient who received a score of 120 would have a 78.86%
risk of their tissue sample being malignant.

| Score |   Risk |
|------:|-------:|
|    30 | 0.0012 |
|    60 | 0.0176 |
|    90 | 0.2052 |
|   120 | 0.7886 |
|   150 | 0.9818 |
|   180 | 0.9987 |
|   210 | 0.9999 |
|   240 | 1.0000 |
