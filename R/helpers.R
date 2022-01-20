#' Create a model frame when there are many predictors.
#'
#' @param fo A formula.
#' @param data A data frame.
#' @param n.cores The number of cores to use in parallel processing.
#' @details The function breaks the formula's right had side into it's
#'   component part, and uses each to construct a model matrix without an
#'   intercept. So the returned object has no intercept, but \code{xgboost}
#'   ought to be able to cope.
#' @export
big.model.matrix <- function(fo, data, n.cores){
  nms <- make.names(names(data))
  if (any(names(data) != nms)){
    stop("Do names(data) <- make.names(names(data))")
  }

  cf <- as.character(fo)
  y <- cf[2]
  lp <- cf[3]
  Xnms <- strsplit(lp, " \\+ ")[[1]]

  if (Xnms == "."){
    Xnms <- names(data)[names(data) != y]
  }

  mm.na <- function(fo, d, X){
    o <- model.matrix(fo, d)
    nas <- is.na(d[, X])
    o2 <- matrix(NA, ncol = ncol(o), nrow = nrow(d))
    colnames(o2) <- colnames(o)
    o2[!nas, ] <- o
    as.data.frame(o2)
  }

  if (n.cores == 1){
    warning("Untested with n.cores = 1")
    X <- lapply(Xnms, function(X){
      f <- formula(paste("~ -1 +", X))
      as.data.frame(mm.na(f, data, X))
    }) %>% bind_cols()
  } else {
    cl <- parallel::makeCluster(n.cores)
    on.exit(parallel::stopCluster(cl))

    X <- parallel::parLapply(cl, Xnms, function(X){
      f <- formula(paste("~ -1 +", X))
      as.data.frame(mm.na(f, data, X))
    }) %>% bind_cols()

    X
  }
}



#' Plot the training and testing set error over iterations
#' @param object An object of class 'xgbm'.
#' @param plot.it Whether to plot, defaulting to \code{plot.it = TRUE}.
#'   If \code{plot.it = FALSE}, the best number of trees, as determined by
#'   cross validation, is returned.
#' @export
xgbm.perf <- function(object, plot.it = TRUE){
  if (!plot.it){
    return(bestCViter(object))
  } else {
    d <- as.data.frame(object$cv$evaluation_log)
    bt <- bestCViter(object)

    test <- getErrorName(d)
    train <- getErrorName(d, "train")
    d[, "testerror"] <- d[, test]
    d[, "trainerror"] <- d[, train]
    d <- d[, c("iter", "testerror", "trainerror")]
    names(d) <- c("iter", "CV error", "Training error")
    d <- gather(d, what, value, -iter)

    ggplot(d, aes(iter, value, color = what)) +
      geom_line() +
      xlab("Iteration") +
      ylab("Error") +
      theme(legend.title = element_blank())
  }
}

getErrorName <- function(el, what = "test_"){
  nms <- names(el)[substring(names(el), 1, 5) == what]
  nms <- nms[substring(nms, nchar(nms) - 3) == "mean"]
  nms
}



#' Get the best interation in terms of CV error (log-likelihood)
#' @param object An object of class "xgbm"
#' @export
bestCViter <- function(object){
  d <- as.data.frame(object$cv$evaluation_log)

  iter <- d$iter
  what <- getErrorName(d)

  d <- d[, what]

  minCV <- iter[d == min(d)]
  min(minCV)
}

#' Deal with fallout from R 4.0 putting stringsAsFactors = FALSE
#' @param data A data frame.
#' @return A data frame in which all character columns have been turned into factors.
#' @export
stringsAsFactors <- function (data){
  data[sapply(data, is.character)] <- lapply(data[sapply(data,
                                                         is.character)], as.factor)
  data
}
