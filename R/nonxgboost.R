#' Extract a wide form of the adverse event, conmeds, or other dataset
#' @export
#' @param data The adverse event dataset
#' @param dm The demography dataset
#' @param subject The unique subject identifier. Defaults to \code{id="usubjid"}
#' @param term The name of the column identifying AE terms. Defaults to \code{term="aeterm"}
#' @param arm The name of the treatment group column. Defaults to \code{trt="arm"}
#' @param arm.col Whether to include the treatment arm as a column in the output. Defaults to \code{arm.col=FALSE}
#' @param subject.col Whether to include the subject identifier as a column in the output (the rownames get the subject identifiers). Defaults to \code{id.col=TRUE}
#' @param drop Drop columns in which fewer than \code{drop} of cases are 1. Defaults to \code{drop=0}
#' @export
wideData <- function(data, dm, subject="subject", term="PT", arm="arm",
                     arm.col=FALSE, subject.col=TRUE, drop=0){
  if (class(data) != "data.frame" | class(dm) != "data.frame")
    stop("data and dm should be data.frames")

  if(arm.col & !is.element(arm, names(dm)))
    stop("dm doesn't contain arm")

  uae <- unique(data[, term])

  if (length(uae) != length(unique(casefold(uae))))
    stop("term contains mixed case versions of identical entries")

  data[, term] <- casefold(trimws(as.character(data[, term])))
  data[, subject] <- as.character(data[, subject])
  dm[, subject] <- as.character(dm[, subject])

  if (any(data[, term] == ""))
    stop("term includes empty strings")

  uae <- unique(data[, term])

  fun <- function(i, data, subs, demog, term){
    sb <- demog[, subs]
    data <- data[data[, term] == i,  ]
    m <- match(data[, subs], sb)
    res <- rep(0, length(sb))
    res[m] <- 1
    res
  }

  aes <- lapply(uae, fun, data = data, subs = subject, demog = dm, term = term)
  wae <- as.data.frame(do.call("cbind", aes))
  dimnames(wae) <- list(dm[, subject], uae)

  # Drop low-frequency events
  mns <- apply(wae, 2, sum) <= drop
  wae <- wae[, !mns]

  # Add arm and id if requested
  if (arm.col) wae[, arm] <- factor(dm[, arm])
  if (subject.col) wae[, subject] <- rownames(wae)

  invisible(wae)
}

#' Compute relative risk with approximate confidence intervals
#'
#' @param x The number of success in the first group
#' @param nx The number of observations in the first group
#' @param y The number of success in the second group
#' @param ny The number of observations in the second group
#' @param alpha Determines the confidence level. Defaults to \code{alpha=c(0.5, 0.1)}
#'   and both 90\% and 50\% intervals are produced.
#' @param method Character string giving the method. Allowable options are
#'   \code{method="MLE"} and \code{method="Jeffreys"}.(the default). It isn't case-sensitive
#' @details If \code{method="Jeffreys"} (the default), the function adds 1/2 to
#'   the numerator and 1 to the denominator when estimating the
#'   proporions of successes in each group, thus avoiding infinite values. Then,
#'   estimate the standard error of the log relative risks, compute the confidence
#'   intervals, and transform back to the original scale.
#' @export
rrci <- function(x, nx, y, ny, alpha=.050, method="Jeffreys"){
  method <- casefold(method)

  z <- qnorm(1 - alpha/2)

  if (method == "jeffreys"){
    px <- (x + .50) / (nx + 1)
    py <- (y + .50) / (ny + 1)
    sel <- sqrt((1/(x + .50) + 1/(y + .50)) - (1/(nx + 1) + 1/(ny + 1)))
  } else if (method == "mle"){
    px <- x / nx
    py <- y / ny
    sel <- sqrt((1/x + 1/y) - (1/nx + 1/ny))
  } else {
    stop("method can be either 'Jeffreys' or 'MLE'")
  }

  lrr <- log(px) - log(py)

  data.frame(RR=lrr, Lower=lrr - z * sel, Upper = lrr + z * sel) %>%
    mutate(event = rownames(.), across(RR:Upper, exp), px = px, py = py, x = x, y = y)
}

#' Get similar strings in a vector of strings
#' @param x The vector of strings to be searched.
#' @param terms Vetor of strings to be approximately matched.
#' @param maxdist Passed to \code{agrepl} as \code{max.distance}.
#' @export
getSynonyms <- function(x, terms, maxdist = .1){
  u <- sort(unique(x))

  theTerms <- c()

  for (term in terms){
    ii <- agrepl(term, x, max.distance = maxdist)
    theTerms <- c(theTerms, x[ii])
  }

  unique(theTerms)
}
