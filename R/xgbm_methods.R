#' @method print xgbm
#' @export
print.xgbm <- function(x, ...){
  bcv <- bestCViter(x)
  cat("An xgbm object fit by xgboost::xgb.train and xgboost::xgb.cv\n")
  cat("The best cross-validation iteration was", bcv, "\n\nCall:\n")
  print(x$call)
}

#' @method summary xgbm
#' @export
summary.xgbm <- function(object, ...){
  # This is a simple wrapper to \code{relative.influence}, mostly to
  # mimic the gbm package and make existing code easier to reuse.
  res <- as.data.frame(relative.influence(object))
  names(res) <- c("var", "rel.inf")

  if (length(unique(object$y)) == 2){
    p <- as.numeric(predict(object, type = "response") > .5)
    p <- paste("Predicted:", p)
    y <- paste("Observed:", object$y)
    cm <- as.data.frame.matrix(table(p, y))
  }

  res <- list(RI = res, ConfusionMatrix = cm)
  class(res) <- "summary.xgbm"
  res
}

#' @export
print.summary.xgbm <- function(x, ...){
  cat("Relative influence (first 6):\n")
  print(head(x[[1]]))
  cat("...\n\nConfusion matrix (training data, best CV error):\n")
  print(x[[2]])

  invisible()
}

#' Get predictions from an xgbm object
#' @param x A object of class 'xgbm'.
#' @param newdata A data frame. If not specified, the original training data is used.
#' @param n.trees The number of trees to use. If not specified, the number that
#'   optimizes cross-validation error is used.
#' @param type String specifying either "link" or "response". Defaults to
#'   \code{type = "link"} for compatibility with \code{gbm}, but will seldom
#'   be desired.
#' @method predict xgbm
#' @export
predict.xgbm <- function(x, newdata = NULL, n.trees = NULL, type = "link", n.cores = NULL){
  if (is.null(n.trees)){
    n.trees <- bestCViter(x)
  }

  if (type == "link"){
    outputmargin <- TRUE
  } else if (type == "response") {
    outputmargin <- FALSE
  } else {
    stop("type should be either 'link' or 'response'")
  }

  if (is.null(newdata)){
    newdata <- x$X
  } else {
    newdata <- newdata[, x$xnames]
    newdata <- stringsAsFactors(newdata)

    classData <- sapply(x$data[, x$xnames], class)
    classNewdata <- sapply(newdata, class)

    if (any(classData != classNewdata)){
      wrn <- "columns in newdata don't all have the same classes as columns in the training data."
      wrn <- paste(wrn, "Does newdata have a column of all NAs that R has decided is logical?")
      wrn <- paste(wrn, "Attempting coercion, but this is risky!")
      warning(wrn)

      newdata <- lapply(1:ncol(newdata), function(X){
        if (class(newdata[[X]]) != class(x$data[[X]])){
          wh <- try(class(newdata[[X]]) <- class(x$data[[X]]))
          if (inherits(wh, "try-error")){
            stop("Coercion failed. Check the column types of newdata: they need to match those in the training data.")
          }
        }

        newdata[[X]]
      })
    }

    fo <- formula(x$formula)[-2]

    # Set any levels unseen in the training data to NA
    for(fac in names(x$xlevels)){
      tmp <- as.character(newdata[[fac]])
      tmp <- ifelse(tmp %in% x$xlevels[[fac]], tmp, NA)
      newdata[[fac]] <- factor(tmp)
    }

    fo <- as.character(formula(x))
    fo <- as.formula(paste("~", fo[3]))

    ##mf <- model.frame(fo, newdata, na.action = na.pass, xlev = x$xlevels)
    ##newdata <- model.matrix(fo, mf)

    if (is.null(n.cores)){
      n.cores <- parallel::detectCores() - 1
    }
    newdata <- big.model.matrix(formula(x), newdata, n.cores)
    newdata <- newdata[, colnames(newdata) %in% x$model$feature_names]
  }

  mc <- colnames(x$X)[!(colnames(x$X) %in% colnames(newdata))]
  if (length(mc) != 0){
    stop(paste0("Missing column(s) in newdata: ", paste(mc, collapse = ", "), ". (All NAs?)"))
  }

  res <- predict(x$model, newdata = as.matrix(newdata), ntreelimit = n.trees,
                 outputmargin = outputmargin)

  if (x$call$distribution == "multinomial"){
    res <- matrix(res, byrow = TRUE, ncol = length(unique(x$y)))
    yname <- as.character(x$call$formula[2])
    colnames(res) <- levels(x$data[, yname])
  }

  res
}

#' Fitted values from an xgbm object
#' @param object An object of class 'xgbm'.
#' @param type String giving the type of fitted values. Defaults to 'cv' and
#'   the fitted values are the test set cross-validation predictions. The only
#'   other option is 'model' in which case the fitted model, at the best
#'   cross-validation iteration, is run on the data used to train it.
#' @method fitted xgbm
#' @export
fitted.xgbm <- function(object, type = "cv", ...){
  if (type == "cv"){
    object$cv$pred
  } else if (type == "model") {
    nd <- object$X
    invisible(predict(object$m, newdata=nd, ntreelimit = bestCViter(object)))
  } else {
    stop("type should be 'cv' or 'model'")
  }
}

#' @method residuals xgbm
#' @export
residuals.xgbm <- function(object, ...){
  object$y - fitted(object)
}


if (FALSE){
  #' @method ggplot xgbm
  #' @export
  ggplot.xgbm <- function(data, mapping=NULL, ...){
    pd <- data$cv$evaluation_log[, c(1, 2, 4)]
    names(pd) <- c("iter", "train", "test")

    gather(Set, value, -iter)

    minCV <- bestCViter(data)
    minCV <- data.frame(v = min(minCV))

    ggplot(pd, aes(iter, value, color=Set, group=Set)) +
      geom_line() +
      theme(legend.title = element_blank()) +
      xlab("Iteration") +
      ylab("Objective") +
      geom_vline(data=minCV, aes(xintercept = v), color="blue", linetype=2)
  }
}


#' Partial dependence plots for xgbm
#' @param data An object of class 'xgbm'.
#' @param mapping Not used.
#' @param i.var A number giving the predictor variable to plot (the function
#'   assumes there is no intercept), or a string giving the variable name.
#' @param return.grid Whether to return the grid of data points and predicted
#'   values. Defaults to \code{return.grid = FALSE} and the plot is returned.
#' @param type String with value 'auto' (the default), 'regression' or 'classification'.
#'   Passed through to \code{pdp::partial}.
#' @method ggplot xgbm
#' @export
ggplot.xgbm <- function(data, mapping = NULL, i.var = 1, return.grid = FALSE, type = "auto", color = "blue", ...){
  if (is.numeric(i.var)){
    i.var <- data$xnames[i.var]
  }

  # Deal with factors
  fac <- is.factor(data$data[[i.var]])
  if (fac){
    mm <- as.data.frame(unique(model.matrix(formula(paste("~", i.var)), data = data$data)))
    mm <- na.omit(mm)
    what <- colnames(mm)
    #cats <- what
    grid.resolution <- 2
    lvls <- factor(levels(data$data[[i.var]]), levels = levels(data$data[[i.var]]))
  } else {
    what <- i.var
    #cats <- NULL
    grid.resolution <- NULL
  }

  if (fac){
    pp <- pdp::partial(data$model, what, train = data$X, pred.grid = mm, #cats = cats,
                       grid.resolution = grid.resolution, type = type)
  } else {
    pp <- pdp::partial(data$model, what, train = data$X,
                       grid.resolution = grid.resolution, type = type)
  }

  if (!return.grid){
    if (!fac & ncol(pp) == 2){
      # So-called tidy evaluation is extremely untidy
      ggplot(pp, aes(!!ensym(what), yhat)) +
        geom_line(color = color)
    } else if (fac){
      pp[, i.var] <- lvls
      ggplot(pp, aes(!!ensym(i.var), yhat)) +
        geom_point(color = color)
    } else {
      stop("ggplot.xgbm is currently only implemented for a single predictor - try using pdp::partial directly")
    }
  } else { # return.grid = TRUE
    if (fac){
      pp <- data.frame(pp$yhat, lvls)
      names(pp) <- c("yhat", i.var)
    }
    pp
  }
}


#' @method plot xgbm
#' @export
plot.xgbm <- ggplot.xgbm
