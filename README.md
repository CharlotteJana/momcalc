Package momcalc
================

<!-- README.md is generated from README.Rmd. Please edit that file -->

# Introduction

Package `momcalc` includes different functions to *calculate moments* of
some distributions. It is possible to calculate the moments of a normal,
lognormal or gamma distribution symbolically. These distributions may be
multivariate. The distribution or moments of the [BEGG
distribution](https://www.researchgate.net/publication/280136422_A_Bimodal_Extension_of_the_Generalized_Gamma_Distribution)
can be calculated numerically. The package also provides a test
concerning the modality of a distribution. It tests, if a one
dimensional distribution with compact support can be unimodal and is
based on the moments of that distribution.

You can install `momcalc` from github with:

``` r
# install.packages("devtools")
devtools::install_github("CharlotteJana/momcalc")
```

# Symbolical calculation of moments

# The BEGG distribution

The *bimodal extension of the generalized Gamma-Distribution* (BEGG) was
first introduced in

  - Bulut, Y. M., & Arslan, O. (2015). [A bimodal extension of the
    generalized gamma
    distribution](https://www.researchgate.net/publication/280136422_A_Bimodal_Extension_of_the_Generalized_Gamma_Distribution).
    *Revista Colombiana de Estadística*, 38(2), 371-384.

It is a scale mixture of the generalized gamma distribution that is
almost always bimodal. The two modes can have different shapes,
depending on the parameters \(α\), \(β\), \(δ_0\), \(δ_1\), \(η\),
\(ε\), \(μ\) and \(σ\).

The density function can be calculated with `dBEGG` and is given by \[ 
f(x) = \begin{cases} 
      \frac{αβ}{2η^\frac{δ₁+1}{α}(1+ε)^{δ₁}Γ(\frac{δ₁+1}{αβ})}~(-x)^{δ₁}~e^{-\frac{(-x)^{αβ}}{η^β(1+ε)^{αβ}}}~~~~ &if~x<0\\
      \frac{αβ}{2η^\frac{δ₀+1}{α}(1-ε)^{δ₀}Γ(\frac{δ₀+1}{αβ})}~x^{δ₀}~e^{-\frac{x^{αβ}}{η^β(1-ε)^{αβ}}}~~~~ &if~x≥0\\
\end{cases} 
\]

The k-th raw moment is given by \[
E(X^k) = \frac{(-1)^kη^{k/α}(1+ε)^{k+1}}{2}~\frac{Γ(\frac{δ₁+k+1}{αβ})}{Γ(\frac{δ₁+1}{αβ})}~+~\frac{η^{k/α}(1-ε)^{k+1}}{2}~\frac{Γ(\frac{δ₀+k+1}{αβ})}{Γ(\frac{δ₀+1}{αβ})}
\] and can be calculated with `mBEGG`.

![](man/figures/README-unnamed-chunk-2-1.png)<!-- -->

# Test if a distribution is unimodal

Function `is.unimoal` checks if an (unknown) distribution with compact
support can be unimodal. It uses several inequalities that were
introduced in

  - Teuscher, F., & Guiard, V. (1995). [Sharp inequalities between
    skewness and kurtosis for unimodal
    distributions](https://www.sciencedirect.com/science/article/pii/016771529400074I).
    *Statistics & probability letters*, 22(3), 257-260.
  - Johnson, N. L., & Rogers, C. A. (1951). [The moment problem for
    unimodal
    distributions](https://www.jstor.org/stable/2236630?seq=1#page_scan_tab_contents).
    *The Annals of Mathematical Statistics*, 433-439.
  - Simpson, J. A., & Welch, B. L. (1960). [Table of the bounds of the
    probability integral when the first four moments are
    given](https://www.jstor.org/stable/2333310?seq=1#page_scan_tab_contents).
    *Biometrika*, 47(3/4), 399-410.

Depending on the inequality, moments up to order 2 or 4 are required. A
distribution that satisfies all inequalities that contain only moments
up to order 2 is called **2-b-unimodal**. A distribution that satisfies
all inequalities that contain only moments up to order 4 is called
**4-b-unimodal**. It is possible that a multimodal distribution
satisfies all inequalities and is therefore 2- and even 4-bimodal. But
if at least one of the inequalities is not satisfied, the distribution
cannot be unimodal. In this case, the test returns `not unimodal` as
result.

Here are some examples using the BEGG distribution:

``` r
# example 1 (bimodal)
example1 <- mBEGG(order = 1:4, alpha = 2, beta = 2, delta0 = 1, delta1 = 4, eta = 1, eps = 0)
is.unimodal(-2, 2, example1)
#> [1] "not unimodal"

# example 2 (bimodal)
example2 <- mBEGG(order = 1:4, alpha = 2, beta = 1, delta0 = 0, delta1 = 2, eta = 1, eps = -0.5)
is.unimodal(-2, 3, example2)
#> [1] "4-b-unimodal"

# example 3 (bimodal)
example3 <- mBEGG(order = 1:4, alpha = 3, beta = 2, delta0 = 4, delta1 = 2, eta = 2, eps = 0.3)
is.unimodal(-2.5, 1.5, example3[1:2]) # test with moments of order 1 and 2
#> [1] "2-b-unimodal"
is.unimodal(-2.5, 1.5, example3) # test with moments of order 1 - 4
#> [1] "not unimodal"

# example 4 (unimodal)
example4 <- mBEGG(order = 1:4, alpha = 2, beta = 1, delta0 = 0, delta1 = 0, eta = 1, eps = 0.7)
is.unimodal(-4, 2, example4)
#> [1] "4-b-unimodal"
```

# License

Package `momcalc` is free open source software licensed under the [GNU
Public License](https://www.gnu.org/licenses/#GPL) (GPL 2.0 or above).
The software is provided as is and comes WITHOUT WARRANTY.
