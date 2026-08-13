# HTEPower

`HTEPower` provides power calculations for machine-learning-based
heterogeneous treatment effect estimation in social science experiments.

## Installation

```r
# install.packages("pak")
pak::pak("andyjiayuwang/HTEPower")
```

## Simplified power formula

The initial function solves the equality boundary of

```math
n \geq
\frac{\sigma^2\left(z_{1-\alpha/2}+z_{\kappa}\right)^2}
{p(1-p)\operatorname{Var}\{s_0(X)\}}.
```

Pass exactly one of `n`, `alpha`, `kappa`, `p`, `sigma2`, or `var_s0` as
`NULL`. The function solves for that quantity using the other five inputs.
When solving for `n`, it rounds the result up to the nearest integer.

```r
library(HTEPower)

hte_power_simple(
  n = NULL,
  alpha = 0.05,
  kappa = 0.8,
  p = 0.5,
  sigma2 = 10,
  var_s0 = 0.25
)
#> [1] 1256
```

Solving for `p` generally returns two symmetric treatment probabilities
because `p * (1 - p)` is unchanged when `p` is replaced by `1 - p`.

```r
hte_power_simple(
  n = 1256,
  alpha = 0.05,
  kappa = 0.8,
  p = NULL,
  sigma2 = 10,
  var_s0 = 0.25
)
```

For more details, see our paper, "When Does Machine Learning Reliably Detect
HTE? Results from Large-Scale Experiments."

Authors and contributors: Jiawei Fu, Donald P. Green, and Andy J. Wang.
