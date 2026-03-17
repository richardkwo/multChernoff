#' An upper bound on the moment generating function of LRT
#'
#' The LRT is the log-likelihood ratio test statistic, which can be written as
#' \deqn{LRT = n \, KL(\hat{p} \| p),}
#' namely the Kullback-Leibler divergence from the empirical probabilities to
#' the true probabilities multiplied by the sample size. \eqn{G(\lambda; k,n)} is
#' a polynomial in \eqn{\lambda} such that
#' \deqn{MGF(\lambda; LRT) := E_{p}[\exp(\lambda \times LRT)] \leq G(\lambda; k,n)}
#' holds for every \eqn{\lambda \in [0,1]} and every \eqn{p}.
#'
#' @param k number of categories
#' @param n sample size
#' @param lambda number between 0 and 1
#'
#' @return A numeric upper bound on the MGF of LRT evaluated at \code{lambda}.
#' @seealso \code{\link{tailProbBound}}, \code{\link{criticalValue}}
#' @export
#'
#' @examples
#' mgfBound(k = 5, n = 20, lambda = 0.5)
#' mgfBound(k = 5, n = 20, lambda = 0)  # always 1 at lambda = 0
mgfBound <- function(k, n, lambda) {
  .S <- sapply(0:n, function(m) {
    lgamma(m+1) + lchoose(n, m) + lchoose(m + k - 2, m) + m *  (log(lambda) - log(n))
  })
  if (length(lambda) == 1) {
    if (lambda == 0) {
      return(1)
    }
    maxS <- max(.S)
    if (is.infinite(maxS)) {
      return(Inf)
    }
    result <- exp(maxS + log(sum(exp(.S - maxS))))
  } else {
    maxS <- apply(.S, 1, max)
    result <- exp(maxS + log(rowSums(exp(.S - maxS))))
  }
  result[lambda == 0] <- 1
  return(result)
}

#' Tail bound on P(2 LRT > x).
#'
#' The LRT is the log-likelihood ratio test statistic, which can be written as
#' \deqn{LRT = n \, KL(\hat{p} \| p).}
#' By the Wilks' theorem, for a fixed k-dimensional
#' probability vector, it holds that
#' \deqn{2 \times LRT \to_{d} \chi^2_{k-1}.}
#' This function returns a finite-sample counterpart to
#' \code{pchisq(x, k-1, lower.tail=FALSE)}.
#' The LRT is also extended to multiple independent multinomial trials.
#' For example, for a \eqn{(k_1, n_1)}-trial and a \eqn{(k_2, n_2)}-trial,
#' we have
#' \deqn{LRT = n_1 \, KL(\hat{p}_1 \| p_1) + n_2 \, KL(\hat{p}_2 \| p_2).}
#'
#' @param x the value of \eqn{2 \times LRT}
#' @param k number of categories (a vector for independent multinomial draws)
#' @param n sample size (a vector for independent multinomial draws)
#' @param verbose draw the minimizer if \code{TRUE}
#' @note For independent multinomial samples, k and n must be of the same length.
#'
#' @seealso \code{\link{criticalValue}}, \code{\link{mgfBound}}
#'
#' @return An upper bound on P(2 LRT > x), which can be used as a conservative p-value.
#' @export
#'
#' @examples
#' tailProbBound(20, 7, 50)
#' pchisq(20, 6, lower.tail=FALSE) # compare with the standard chi-square asymptotic
#' # two independent multinomial trials (k=3, n=4) and (k=12, n=20)
#' tailProbBound(12, c(3, 4), c(12, 20))
#' pchisq(12, 5, lower.tail=FALSE) # compare with the standard chi-square asymptotic
tailProbBound <- function(x, k, n, verbose=FALSE) {
  stopifnot("k and n must be of the same length" = length(k) == length(n))
  if (length(k)==1) {
    f <- function(lambda) {
      lambda <- pmin(pmax(lambda, 0), 1)
      Gval <- mgfBound(k, n, lambda)
      if (!is.finite(Gval) || Gval <= 0) {
        return(1e10)
      }
      - lambda * x / 2 + log(Gval)
    }
  } else {
    f <- function(lambda) {
      lambda <- pmin(pmax(lambda, 0), 1)
      result <- - lambda * x / 2
      for (i in 1:length(k)) {
        Gval <- mgfBound(k[i], n[i], lambda)
        if (!is.finite(Gval) || Gval <= 0) {
          return(1e10)
        }
        result <- result + log(Gval)
      }
      return(result)
    }
  }
  if (length(n) == 1) {
    lambda.init <- max(min(1 - (k-1) / (x / 2) + k / (k-1) * (x / 2 - k + 1) / n, 1), 0)
    lambda.init <- c(lambda.init, seq(0, 1, length.out = 3))
  } else {
    lambda.init <- seq(0, 1, length.out = 4)
  }
  sols <- plyr::laply(lambda.init, function(.x) {
    stats::optim(.x, f, method="L-BFGS-B", lower = 1e-8, upper=0.99999)
  })
  m <- which.min(sols[,2])
  sol.par <- unlist(sols[,1])[m]
  sol.val <- unlist(sols[,2])[m]
  if (verbose) {
    lambda.vec <- seq(0, 1, length.out = 100)
    plot_vals <- sapply(lambda.vec, function(lam) {
      Gval <- mgfBound(k, n, lam)
      if (length(Gval) > 1) {
        Gval <- Gval[1]
      }
      if (!is.finite(Gval) || Gval <= 0) {
        return(1e10)
      }
      -(lam * x / 2) + log(Gval)
    })
    plot(lambda.vec, plot_vals, type="l", xlab="lambda", ylab="bound")
    graphics::abline(v=sol.par, col="red")
  }
  return(min(exp(sol.val), 1))
}

#' Critical value x such that \eqn{P(2 LRT > x) \le p}
#'
#' The LRT is the log-likelihood ratio test statistic, which can be written as
#' \deqn{LRT = n \, KL(\hat{p} \| p).}
#' By the Wilks' theorem, for a fixed k-dimensional
#' probability vector, it holds that
#' \deqn{2 \times LRT \to_{d} \chi^2_{k-1}.}
#' This function returns a finite-sample counterpart to
#' \code{qchisq(p, k-1, lower.tail=FALSE)}.
#' The LRT is also extended to multiple independent multinomial trials.
#' For example, for a \eqn{(k_1, n_1)}-trial and a \eqn{(k_2, n_2)}-trial,
#' we have
#' \deqn{LRT = n_1 \, KL(\hat{p}_1 \| p_1) + n_2 \, KL(\hat{p}_2 \| p_2).}
#'
#' @param k number of categories (a vector for independent multinomial draws)
#' @param n sample size (a vector for independent multinomial draws)
#' @param p significance level (e.g., 0.05)
#' @param verbose draw the minimizer if \code{TRUE}
#' @note For independent multinomial samples, k and n must be of the same length.
#'
#' @return A finite-sample critical value \eqn{x} such that the bound on
#'   \eqn{P(2 \times LRT > x)} is at most \code{p}.
#' @seealso \code{\link{tailProbBound}}, \code{\link{mgfBound}}
#' @export
#'
#' @examples
#' n <- 1:40
#' crit <- sapply(n, function(.n) criticalValue(20, .n, p=0.01))
#' plot(n, crit)
#' # chi-squared asymptotic by Wilks' theorem
#' abline(h=qchisq(0.01, df=20-1, lower.tail = FALSE))
#' criticalValue(10, 40, p=0.05)
#' # two independent multinomial trials (k=3, n=4) and (k=12, n=20)
#' criticalValue(c(3, 4), c(12, 20), p=0.05)
criticalValue <- function(k, n, p=0.05, verbose=FALSE) {
  sol <- stats::uniroot(function(x) tailProbBound(x, k, n) - p, c(0, 1e3))
  if (verbose) {
    tailProbBound(sol$root, k, n, verbose = TRUE)
  }
  return(sol$root)
}
