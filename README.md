
# TINC <img src='man/figures/logo.png' align="right" height="139" />

<!-- badges: start -->

[![Lifecycle:
maturing](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing)
[![pkgdown](https://github.com/caravagnalab/TINC/actions/workflows/pkgdown.yaml/badge.svg)](https://github.com/caravagnalab/TINC/actions/workflows/pkgdown.yaml)
[![check-master](https://github.com/caravagnalab/TINC/actions/workflows/check-master.yml/badge.svg)](https://github.com/caravagnalab/TINC/actions/workflows/check-master.yml)
<!-- badges: end -->

`TINC` is a package to determine the contamination of tumour DNA in a
matched normal sample. The approach uses evolutionary theory applied to
read counts data from whole-genome sequencing assays.

Precisely, the package provides methods to determine, for every matched
pair of normal and tumour sample biopsies, the proportion of cancer
cells, or tumour read fractions, contaminating the normal sample (Tumour
in Normal, TIN). Similarly, it determines the proportion of cancer cells
in the tumour sample (Tumour in Tumour, TIT), also called tumour purity.

If the TIN contamination determined by TINC is high, there are good
chances that a somatic caller will be affected by false negative calls
(i.e., non-called mutations that have read support in the normal). In
these cases TINC can be used to flag samples better suited for a
tumour-only somatic calling pipeline.

#### Help and support

## [![](https://img.shields.io/badge/GitHub%20Pages-https://caravagnalab.github.io/TINC/-yellow.svg)](https://caravagnalab.github.io/TINC)

### Installation

You can install the released version of `TINC` from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("caravagnalab/TINC")
```

------------------------------------------------------------------------

#### Copyright and contacts

Giulio Caravagna. Cancer Data Science (CDS) Laboratory.

[![](https://img.shields.io/badge/Email-gcaravagn@gmail.com-steelblue.svg)](mailto:gcaravagn@gmail.com)
[![](https://img.shields.io/badge/CDS%20Lab%20Github-caravagnalab-seagreen.svg)](https://github.com/caravagnalab)
[![](https://img.shields.io/badge/CDS%20Lab%20webpage-https://www.caravagnalab.org/-red.svg)](https://www.caravagnalab.org/)
