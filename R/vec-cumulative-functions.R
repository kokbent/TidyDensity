#' Cumulative Variance
#'
#' @family Vector Function
#'
#' @author Steven P. Sanderson II, MPH
#'
#' @details
#' A function to return the cumulative variance of a vector.
#' `exp(cummean(log(.x)))`
#'
#' @description
#' A function to return the cumulative variance of a vector.
#'
#' @param .x A numeric vector
#'
#' @examples
#' x <- mtcars$mpg
#'
#' cvar(x)
#'
#' @return
#' A numeric vector. Note: The first entry will always be
#' NaN.
#'
#' @export
#'

cvar <- function(.x) {
  n <- length(.x)
  csumx <- base::cumsum(.x)
  cmeanx <- cmean(.x)

  p1 <- base::cumsum(.x^2)
  p2 <- -2 * cmeanx * csumx
  p3 <- (1:n) * cmeanx^2
  (p1 + p2 + p3) / ((1:n) - 1)
}

#' Cumulative Skewness
#'
#' @family Vector Function
#'
#' @author Steven P. Sanderson II, MPH
#'
#' @details
#' A function to return the cumulative skewness of a vector.
#'
#' @description
#' A function to return the cumulative skewness of a vector.
#'
#' @param .x A numeric vector
#'
#' @examples
#' x <- mtcars$mpg
#'
#' cskewness(x)
#'
#' @return
#' A numeric vector
#'
#' @export
#'

cskewness <- function(.x) {
  n <- length(.x)
  csumx <- base::cumsum(.x)
  csumx2 <- base::cumsum(.x^2)
  csumx3 <- base::cumsum(.x^3)
  cmeanx <- cmean(.x)

  order3 <- csumx3 - (1:n) * cmeanx^3
  order2 <- -3 * cmeanx * csumx2 + 3 * cmeanx^2 * csumx
  num <- order3 + order2

  order2 <- csumx2 + (1:n) * cmeanx^2
  order1 <- -2 * cmeanx * csumx
  den <- order2 + order1

  sqrt(1:n) * num / den^(3/2)
}

#' Cumulative Kurtosis
#'
#' @family Vector Function
#'
#' @author Steven P. Sanderson II, MPH
#'
#' @details
#' A function to return the cumulative kurtosis of a vector.
#'
#' @description
#' A function to return the cumulative kurtosis of a vector.
#'
#' @param .x A numeric vector
#'
#' @examples
#' x <- mtcars$mpg
#'
#' ckurtosis(x)
#'
#' @return
#' A numeric vector
#'
#' @export
#'

ckurtosis <- function(.x) {
  n <- length(.x)
  csumx <- base::cumsum(.x)
  csumx2 <- base::cumsum(.x^2)
  csumx3 <- base::cumsum(.x^3)
  csumx4 <- base::cumsum(.x^4)
  cmeanx <- cmean(.x)

  order4 <- csumx4 + (1:n) * cmeanx^4
  order3 <- -4 * csumx3 * cmean(.x) - 4 * csumx * cmean(.x)^3
  order2 <- 6 * csumx2 * cmean(.x)^2
  num <- order4 + order3 + order2

  order2 <- csumx2 + (1:n) * cmeanx^2
  order1 <- -2 * cmeanx * csumx
  den <- (order2 + order1)^2

  (1:n) * num / den
}

#' Cumulative Standard Deviation
#'
#' @family Vector Function
#'
#' @author Steven P. Sanderson II, MPH
#'
#' @details
#' A function to return the cumulative standard deviation of a vector.
#'
#' @description
#' A function to return the cumulative standard deviation of a vector.
#'
#' @param .x A numeric vector
#'
#' @examples
#' x <- mtcars$mpg
#'
#' csd(x)
#'
#' @return
#' A numeric vector. Note: The first entry will always be
#' NaN.
#'
#' @export
#'

csd <- function(.x) {
  sqrt(cvar(.x))
}

#' Cumulative Median
#'
#' @family Vector Function
#'
#' @author Steven P. Sanderson II, MPH
#'
#' @details
#' A function to return the cumulative median of a vector.
#'
#' @description
#' A function to return the cumulative median of a vector.
#'
#' @param .x A numeric vector
#'
#' @examples
#' x <- mtcars$mpg
#'
#' cmedian(x)
#'
#' @return
#' A numeric vector
#'
#' @export
#'

cmedian <- function(.x) {
  sapply(seq_along(.x), function(k, z) stats::median(z[1:k]), z = .x)
}

#' Cumulative Geometric Mean
#'
#' @family Vector Function
#'
#' @author Steven P. Sanderson II, MPH
#'
#' @details
#' A function to return the cumulative geometric mean of a vector.
#' `exp(cummean(log(.x)))`
#'
#' @description
#' A function to return the cumulative geometric mean of a vector.
#'
#' @param .x A numeric vector
#'
#' @examples
#' x <- mtcars$mpg
#'
#' cgmean(x)
#'
#' @return
#' A numeric vector
#'
#' @export
#'

cgmean <- function(.x) {
  exp(cmean(log(.x)))
}

#' Cumulative Harmonic Mean
#'
#' @family Vector Function
#'
#' @author Steven P. Sanderson II, MPH
#'
#' @details
#' A function to return the cumulative harmonic mean of a vector.
#' `1 / (cumsum(1 / .x))`
#'
#' @description
#' A function to return the cumulative harmonic mean of a vector.
#'
#' @param .x A numeric vector
#'
#' @examples
#' x <- mtcars$mpg
#'
#' chmean(x)
#'
#' @return
#' A numeric vector
#'
#' @export
#'

chmean <- function(.x) {
  1 / (cumsum(1 / .x))
}

#' Cumulative Mean
#'
#' @family Vector Function
#'
#' @author Steven P. Sanderson II, MPH
#'
#' @details
#' A function to return the cumulative mean of a vector. It uses [dplyr::cummean()]
#' as the basis of the function.
#'
#' @description
#' A function to return the cumulative mean of a vector.
#'
#' @param .x A numeric vector
#'
#' @examples
#' x <- mtcars$mpg
#'
#' cmean(x)
#'
#' @return
#' A numeric vector
#'
#' @export
#'

cmean <- function(.x) {
  dplyr::cummean(.x)
}
