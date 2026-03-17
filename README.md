# multChernoff

This package computes a finite-sample tail bound of the likelihood ratio test (LRT) under multinomial sampling. The tail bounds can be used to obtain conservative p-values and critical values. This is useful for inference when the sample size is comparable to or even smaller than the alphabet size, where the standard chi-square asymptotic ([Wilks' theorem](https://en.wikipedia.org/wiki/Wilks%27_theorem?oldformat=true)) may not hold. 

## Installation

You can install the package from CRAN with  
``` r
> install.packages("multChernoff")
```
or from GitHub with
``` r
> devtools::install_github("richardkwo/multChernoff")
```

## Usage

Please refer to the vignette.
``` r
> vignette("multChernoff")
```

The package can be used with the finite-sample critical value `criticalValue`
to construct a convex confidence region on the underlying probability vector.

## Reference

The method is based on the following work:

F. Richard Guo and Thomas S. Richardson, "[Chernoff-Type Concentration of Empirical Probabilities in Relative Entropy](https://dx.doi.org/10.1109/TIT.2020.3034539)," in IEEE Transactions on Information Theory, vol. 67, no. 1, pp. 549-558, Jan. 2021.
