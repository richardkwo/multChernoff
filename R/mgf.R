#' G(lambda; k,n) upper bound of the moment generating function on (LRT / 2)
#'
#' @param k number of categories
#' @param n sample size
#' @param lambda number between 0 and 1
#'
#' @export
mgfBound <- function(k, n, lambda) {
  # compute on the log scale
  .S <- sapply(0:n, function(m) {
    lgamma(m+1) + lchoose(n, m) + lchoose(m + k - 2, m) + m *  (log(lambda) - log(n))
  })
  if (length(lambda) == 1) {
    result <- sum(exp(.S))
  } else {
    result <- apply(exp(.S), 1, sum)
  }
  result[lambda == 0] <- 1
  return(result)
}

#' Tail bound on P(LRT > x)
#'
#' @param x the value of the LRT
#' @param k number of categories (a vector for independent multinomial draws)
#' @param n sample size (a vector for independent multinomial draws)
#' @param verbose draw the minimizer if \code{TRUE}
#' @note For independent multinomial samples, k and n must be of the same length.
#'
#' @return an upper bound on P(LRT > x), which can be used as a conservative p-value
#' @export
#'
#' @examples
#' tailProbBound(20, 7, 50)
#' pchisq(20, 6, lower.tail=FALSE) # compare with the standard chi-square asymptotic
#' tailProbBound(30, c(5,10), c(20,50))
#' pchisq(20, 13, lower.tail=FALSE) # compare with the standard chi-square asymptotic
tailProbBound <- function(x, k, n, verbose=FALSE) {
  stopifnot("k and n must be of the same length" = length(k) == length(n))
  if (length(k)==1) {
    f <- function(lambda) {
      lambda <- pmin(pmax(lambda, 0), 1)
      - lambda * x / 2 + log(mgfBound(k, n, lambda))
    }
  } else {
    f <- function(lambda) {
      lambda <- pmin(pmax(lambda, 0), 1)
      result <- - lambda * x / 2
      for (i in 1:length(k)) {
        result <- result + log(mgfBound(k[i], n[i], lambda))
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
    stats::optim(.x, f, method="L-BFGS-B", lower = 0, upper=1)
  })
  m <- which.min(sols[,2])
  sol.par <- unlist(sols[,1])[m]
  sol.val <- unlist(sols[,2])[m]
  if (verbose) {
    lambda.vec <- seq(0, 1, length.out = 100)
    plot(lambda.vec, f(lambda.vec), type="l", xlab="lambda", ylab="bound")
    graphics::abline(v=sol.par, col="red")
  }
  return(exp(sol.val))
}

#' Critical value x such that \eqn{P(LRT > x) \le p}
#'
#' @param k number of categories (a vector for independent multinomial draws)
#' @param n sample size (a vector for independent multinomial draws)
#' @param p significance level (e.g., 0.05)
#' @param verbose draw the minimizer if \code{TRUE}
#'
#' @seealso tailProbBound
#' @export
#'
#' @examples
#' n <- 1:100
#' crit <- sapply(n, function(.n) criticalValue(20, .n, p=0.01))
#' plot(n, crit)
#' abline(h=qchisq(0.01, df=20-1, lower.tail = FALSE)) # standard chi-square asymptotic
#' criticalValue(20, 100, p=0.01, verbose=TRUE)
criticalValue <- function(k, n, p=0.05, verbose=FALSE) {
  sol <- stats::uniroot(function(x) tailProbBound(x, k, n) - p, c(0, 1e3))
  if (verbose) {
    tailProbBound(sol$root, k, n, verbose = TRUE)
  }
  return(sol$root)
}
