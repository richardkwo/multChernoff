
# multChernoff

This package computes finite-sample tail bounds of the likelihood ratio test (LRT) under multinomial sampling. The tail bounds can be used to obtain conservative p-values and critical values. 

This is useful for inference when the sample size is comparable to or even smaller than the alphabet size, where the standard chi-square asymptotic ([Wilks' theorem](https://en.wikipedia.org/wiki/Wilks%27_theorem?oldformat=true)) may not hold. 

## Installation

You can install the released version of `multChernoff` from GitHub with:

``` r
devtools::install_github("richardkwo/multChernoff")
```

## Examples

Through the following examples, we will use the finite-sample critical value `criticalValue` to construct a convex confidence region on the underlying probability vector. We will use the package [CVXR](https://cvxr.rbind.io/) for writing and solving convex programs. `CVXR` can be installed with:

```R
install.packages("CVXR")
```

### Example 1: What percentage of butterflies were unseen by Corbet?

Naturalist Alexander Steven Corbet spent two years trapping butterflies in Peninsular Malaysia. Here is the data collected by Corbet.

| Frequency | 1    | 2    | 3    | 4    | 5    | 6    | 7    | 8    | 9    | 10   | 11   | 12   | 13   | 14   | 15   |
| --------- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- |
| Species   | 118  | 74   | 44   | 24   | 29   | 22   | 20   | 19   | 20   | 15   | 12   | 14   | 6    | 12   | 6    |

Here we ask the question: *what percentage of butterflies in Malaya belonged to the species that Corbet had not seen?* That is, we want to estimate the proportion of butterflies from all the unseen species. Clearly, the MLE is zero based on the sample. Instead, we ask for an upper bound with 95% confidence.

Let k = 435 + 1, where 435 is the number of species observed by Corbet. The corresponding observed distribution is the distribution from the table above with a zero entry at the end. In the following program, we maximize the last entry of the underlying distribution subject to a tail bound on the LRT, which equals the Kullback-Leibler divergence KL(observed distribution || true distribution) multiplied by twice the sample size.

See the paper at the bottom of this page for details. 

``` R
library(multChernoff)
library(CVXR)
library(plyr)

# Corbet butterfly data
# https://en.wikipedia.org/wiki/Unseen_species_problem
corbet.butterfly <- data.frame(j=1:15, n_j=c(118,74,44,24,29,22,20,19,20,15,12,14,6,12,6))
n.butterfly <- Reduce(c, alply(corbet.butterfly, 1, function(.df) rep(.df$j, .df$n_j)))

alpha <- 0.05
n.observed <- c(n.butterfly, 0)    # the last one is the unseen
n <- sum(n.observed)
k <- length(n.observed)
p.observed <- n.observed / n

# critical value
t.alpha <- criticalValue(k, n, p=alpha, verbose = TRUE)
cat(sprintf("critical value = %f\n", t.alpha))
```

We get printout `critical value = 962.402970`. This value is then used in the following convex program.

```R
# convex program
p <- Variable(k)
obj <- p[k]
constr <- list(p>=0, 
               sum(p) == 1, 
               2 * n * sum(p.observed * (log(p.observed) - log(p))) <= t.alpha)
prob <- Problem(Maximize(obj), constr)
result <- solve(prob, verbose=FALSE) # you should try other solvers

# result
print(result$status)
unseen <- result$value
p.maximizer <- c(result$getValue(p))
cat(sprintf("unseen <= %f\n", unseen))
```

With [MOSEK](https://cran.r-project.org/web/packages/Rmosek/) solver, the percentage of unseen species is at most 21.1%.

### Example 2: Analysis of a binary instrumental variable model

Consider the following data from the Lipid Research Coronary Primary Prevention Trial.

|                       | X=0, Y=0 | X=0, Y=1 | X=1, Y=0 | X=1, Y=1 |
| --------------------- | -------- | -------- | -------- | -------- |
| Z=0 (placebo group)   | 158      | 14       | 0        | 0        |
| Z=1 (treatment group) | 52       | 12       | 23       | 78       |

Here X is the actual treatment received and Y is the outcome. The binary instrumental variable (IV) model imposes constraints on the underlying probability P(X, Y| Z=0) and P(X, Y|Z=1). These can be expressed in terms of the [generalized instrumental inequalities](https://doi.org/10.1093/biomet/asaa003). At the same time, these probabilities lie in a convex confidence region characterized by the LRT critical value. In particular, this is an LRT on two independent multinomial samples (Z=0 arm and Z=1 arm). Critical value of this kind can be computed by `criticalValue(c(k1, k2), c(n1, n2), p=0.05)`. The true P(X, Y| Z=0) and P(X, Y|Z=1) are identified by intersecting the aforementioned two regions. If the intersection is empty, the IV model is rejected.

Further, this intersected region implies upper and lower bounds on the average causal effect (ACE) through the Balke-Pearl bounds. The following script executes such analysis. 

```R
library(CVXR)
library(multChernoff)

# each row is a level of Z
data <- matrix(c(158,14,0,0,52,12,23,78), nrow=2, byrow = TRUE)
n.z <- apply(data, 1, sum)
p.empirical <- data / n.z

# get critical value ----
alpha <- 0.05
t.alpha <- criticalValue(rep(4, 2), n.z, p=alpha, verbose = TRUE) 
cat(sprintf("critical value = %f\n", t.alpha))

# make objective -----
p <- Variable(2, 4)
g <- function(x, y) {
  a <- rep(1, 4)
  b <- rep(0, 4)
  d <- rep(0, 4)
  if (x==0 && y==0) {
    a[2] <- 0
    b[1] <- b[3] <- 1
    d[1] <- d[4] <- 1
  } else if (x==0 && y==1) {
    a[1] <- 0
    b[2] <- b[3] <- 1
    d[2] <- d[4] <- 1
  } else if (x==1 && y==0) {
    a[4] <- 0
    b[1] <- b[3] <- 1
    d[2] <- d[3] <- 1
  } else {
    a[3] <- 0
    b[1] <- b[4] <- 1
    d[2] <- d[4] <- 1
  }
  A <- min_entries(p %*% a)
  H <- p %*% replicate(2, b)
  I <- p %*% replicate(2, d)
  B <- min_entries(H + t(I))
  return(min_elemwise(A, B))
}

g.empirical <- function(x, y, p) {
  a <- rep(1, 4)
  b <- rep(0, 4)
  d <- rep(0, 4)
  if (x==0 && y==0) {
    a[2] <- 0
    b[1] <- b[3] <- 1
    d[1] <- d[4] <- 1
  } else if (x==0 && y==1) {
    a[1] <- 0
    b[2] <- b[3] <- 1
    d[2] <- d[4] <- 1
  } else if (x==1 && y==0) {
    a[4] <- 0
    b[1] <- b[3] <- 1
    d[2] <- d[3] <- 1
  } else {
    a[3] <- 0
    b[1] <- b[4] <- 1
    d[2] <- d[4] <- 1
  }
  A <- min(p %*% a)
  H <- p %*% replicate(2, b)
  I <- p %*% replicate(2, d)
  B <- min(H + t(I))
  return(min(A, B))
}

# constraints ------
KL <- sum(sum_entries(p.empirical * (log(p.empirical) - log(p)), axis=1) * 2 * n.z)
IV.fun.1 <- function(x) {
  if (x==0) {
    max_entries(p[,1]) + max_entries(p[,2])
  } else {
    max_entries(p[,3]) + max_entries(p[,4])
  }
}
IV.fun.2 <- function(x, y) {
  if (x==0 && y==0) {
    a <- max_entries(p[,1])
    b <- min_entries(p[,1] + p[,3])
    d <- min_entries(p[,1] + p[,4])
  } else if (x==0 && y==1) {
    a <- max_entries(p[,2])
    b <- min_entries(p[,2] + p[,4])
    d <- min_entries(p[,2] + p[,3])
  } else if (x==1 && y==0) {
    a <- max_entries(p[,3])
    b <- min_entries(p[,1] + p[,3])
    d <- min_entries(p[,2] + p[,3])
  } else {
    a <- max_entries(p[,4])
    b <- min_entries(p[,2] + p[,4])
    d <- min_entries(p[,1] + p[,4])
  }
  return(a - b - d)
}
IV.fun.3 <- function() {
  min_entries(p[,2] + p[,4]) + min_entries(p[,2] + p[,3]) + 
    min_entries(p[,1] + p[,3]) + min_entries(p[,1] + p[,4])
}

constr <- list(p >= 0, 
               sum_entries(p, axis=1) == 1, 
               KL <= t.alpha, 
               IV.fun.1(0) <= 1, IV.fun.1(1) <= 1, 
               IV.fun.2(0,0) <= 0, IV.fun.2(0,1) <= 0, IV.fun.2(1,0) <= 0, IV.fun.2(1,1) <= 0,
               IV.fun.3() >= 1)

# ACE lower bound ----
obj <- 1 - g(1,0) - g(0,1)
prob <- Problem(Minimize(obj), constr) 
result <- solve(prob)

print(result$status)
ACE.lb <- result$value
p.lb <- result$getValue(p)

# ACE upper bound ----
obj <- g(0,0) + g(1,1) - 1
prob <- Problem(Maximize(obj), constr) 
result <- solve(prob, solver="ECOS")

print(result$status)
ACE.ub <- result$value
p.ub <- result$getValue(p)

cat(sprintf("\n%.3f <= ACE <= %.3f (naive plugin gives [%.3f, %.3f])\n", 
            ACE.lb, ACE.ub, 
            1 - g.empirical(1,0,p.empirical) - g.empirical(0,1,p.empirical),
            g.empirical(0,0,p.empirical) + g.empirical(1,1,p.empirical) - 1))

# p-value ------
prob.pval <- Problem(Minimize(KL), list(
  p >= 0, 
  sum_entries(p, axis=1) == 1, 
  IV.fun.1(0) <= 1, IV.fun.1(1) <= 1, 
  IV.fun.2(0,0) <= 0, IV.fun.2(0,1) <= 0, IV.fun.2(1,0) <= 0, IV.fun.2(1,1) <= 0,
  IV.fun.3() >= 1
))
result.pval <- solve(prob.pval)
print(result.pval$status)
KL.ml <- max(result.pval$value, 0)
p.mle <- result.pval$getValue(p)
p.value <- tailProbBound(KL.ml, rep(4, 2), n.z, verbose = TRUE)

cat(sprintf("\np-value (conservative) = %g\n", p.value))
```

Running the script, we get

```R
0.151 <= ACE <= 0.907 (naive plugin gives [0.391, 0.779])
p-value (conservative) = 1
```

## Citation

If you find this useful, consider citing the following article. 

F. Richard Guo and Thomas S. Richardson, "[Chernoff-Type Concentration of Empirical Probabilities in Relative Entropy](https://dx.doi.org/10.1109/TIT.2020.3034539)," in IEEE Transactions on Information Theory, vol. 67, no. 1, pp. 549-558, Jan. 2021.

An open-access preprint is available [here](https://arxiv.org/abs/2003.08614).

