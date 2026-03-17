
# multChernoff

This package computes finite-sample tail bounds of the likelihood ratio test (LRT) under multinomial sampling. The tail bounds can be used to obtain conservative p-values and critical values. 

This is useful for inference when the sample size is comparable to or even smaller than the alphabet size, where the standard chi-square asymptotic ([Wilks' theorem](https://en.wikipedia.org/wiki/Wilks%27_theorem?oldformat=true)) may not hold. 

## Installation

You can install the released version of `multChernoff` from GitHub with:

``` r
devtools::install_github("richardkwo/multChernoff")
```

## Applications

The package can be used with the finite-sample critical value `criticalValue`
to construct convex confidence regions on the underlying probability vector.
Some larger workflows use the package [CVXR](https://cvxr.rbind.io/) for
writing and solving convex programs. `CVXR` can be installed with:

```R
install.packages("CVXR")
```

## Citation

If you find this useful, consider citing the following article. 

F. Richard Guo and Thomas S. Richardson, "[Chernoff-Type Concentration of Empirical Probabilities in Relative Entropy](https://dx.doi.org/10.1109/TIT.2020.3034539)," in IEEE Transactions on Information Theory, vol. 67, no. 1, pp. 549-558, Jan. 2021.

An open-access preprint is available [here](https://arxiv.org/abs/2003.08614).
