#' Wrapper to xgb.cv and xgb.train
#' @param formula,data Formula and data.frame from which to create the model.matrix and response.
#' @param n.trees Number of trees to build.
#' @param interaction.depth,learning.rate Interaction depth (6) and learning rate (0.1).
#' @param bag.fraction Proportion of data in each tree (0.5). If \code{bag.fraction = 1}
#'   and, the whole data set is used as-is.
#' @param col.fraction Proportion of columns of features to use. Defaults to
#'   \code{col.fraction = 1} and the only reason you'd want to reduce it is due
#'   to memory or processing time issues. This is passed as \code{colsample_bynode}
#'   for no reason other than that's what random forest does. I've no idea if
#'   one method is better than another.
#' @param cv.folds Number of cross-validation folds (10).
#' @param cv.class.stratify Whether to stratify cross-validation by the response
#'   values (FALSE).
#' @param n.cores Number of cores to use (1).
#' @param n.minobsinnode Minimum number of observations allowed in a tree node (3).
#' @param leaf.penalty Penalty factor for the total number of leaves in trees (0).
#' @param weight.penalty.L1,weight.penalty.L2 L1 and L2 penalties for leaf weights
#'   (0 for L1 and 1 for L2).
#' @param early.stopping.trees Passed through as \code{early_stopping_rounds} and
#'   defaults to 100.
#' @param distribution The only values allowed are, "gaussian", "huber" (which
#'   uses pseudo-huber loss), "binomial",
#'   "multinomial", "poisson", "quantile" or "coxph". Others should be added as the (my)
#'   need arises. There is no default because experience suggests that leads too
#'   easily to mistakes.
#' @param quant Quantile to be modelled when \code{distribution = "quantile"}. Defaults
#'   to \code{quant = 0.5}.
#' @param event Only used when \code{distribution = "coxph"}. Indicates if the
#'   observational unit has an event (1) or is censored (0).
#' @param verbose Control printing (100). Use \code{verbose = 0} for no printing.
#' @param fail.if.not.converged Defaults to TRUE
#'
#' @details The function takes on the job of turning the data and formula into
#'   the favoured stuff of \code{xgboost} and applies sensible metrics given
#'   the distribution: that is, it does maximum likelihood when a likelihood
#'   function is available (everything but quantile regression).
#'
#'   The response (and any other) variable should be transformed prior to
#'   using \code{xgbm}, if necessary. An apparent bug, somewhere or other,
#'   means that if the transformation is done via the formula, relative
#'   influence goes wrong.
#'
#'   For \code{distribution = "coxph"}, you need to use the \code{event}
#'   argument.
#'
#'   For \code{distribution = "poisson"}, if you need an offset, divide the
#'   response by the exposure and pass exposure in using the weight argument.
#'
#'   Some of the code is quite inefficient and probably annoying. One of the
#'   reasons is that it's best to fail quickly rather than wait for a lot of
#'   processing to be done and then fail.
#' @export
xgbm <- function(formula, data, n.trees=100, interaction.depth=6, learning.rate=.1,
                 weight = NULL,
                 bag.fraction=.5, col.fraction = 1,
                 cv.folds=10, cv.class.stratify=FALSE, n.cores = NULL,
                 n.minobsinnode = 3, leaf.penalty = 0,
                 weight.penalty.L1 = 0, weight.penalty.L2 = 0,
                 early.stopping.trees = 100,
                 distribution = NULL,
                 quant = .5, event = NULL,
                 verbose  =100, fail.if.not.converged = TRUE){

  theCall <- match.call()

  if (is.null(distribution)){
    stop("distribution not specified and I daren't guess")
  }

  nms <- make.names(names(data))
  if (any(names(data) != nms)){
    stop("Do names(data) <- make.names(names(data))")
  }

  formula <- as.formula(formula)

  if (length(as.character(formula[[2]])) != 1){
    stop("You need to do transformations prior to calling xgbm, and the response appears to use a transformation")
  }

  distribution <- tolower(distribution)

  #if (!is.null(metrics) & !is.function(metrics) & length(metrics) != 1){
  #  stop("metrics should be a string (not vector of strings) or a function")
  #}

  if (is.null(n.cores)){
    n.cores <- parallel::detectCores() - 1
  }

  if (distribution == "coxph"){
    event <- deparse(substitute(event))
    if(event == "NULL"){
      stop("For coxph, you need to specify which column indicates an event (unquoted).")
    }

    if (!(event %in% names(data))){
      stop("Specified event column is not in the data.")
    }

    trueevent <- data[, event]
    data <- data[, names(data) != event]
  }

  if (!is.logical(verbose) & verbose > 0 & verbose < Inf){
    print_every_n <- verbose
    verbose <- TRUE
  } else {
    print_every_n <- 1
    verbose <- FALSE
  }

  # Get column classes
  ytxt <- as.character(formula)[2]
  cc <- sapply(data[, names(data) != ytxt], class)

  if ("character" %in% cc){
    if (verbose) message("converting character columns to factors")

    data[sapply(data, is.character)] <- lapply(data[sapply(data, is.character)], as.factor)
  }

  data <- droplevels(data) # otherwise xgb.importance can fail

  X <- big.model.matrix(formula, data, n.cores = n.cores)

  lp <- as.character(formula)[3]
  if (lp == "."){
    Xnames <- names(data)[names(data) != ytxt]
  } else {
    Xnames <- strsplit(as.character(formula)[3], " \\+ ")[[1]]
  }

  xdata <- data[, Xnames]
  data <- data[, c(ytxt, Xnames)]

  xlevels <- lapply(xdata, function(X){
    if (is.factor(X) | is.ordered(X)){
      levels(X)
    } else {
      NULL
    }
  })
  xlevels[sapply(xlevels, is.null)] <- NULL

  y <- data[, ytxt]

  if (distribution == "multinomial"){
    y <- as.numeric(y) - 1
  } else if (distribution %in% c("binomial", "bernoulli")){
    if (inherits(y, "character")){
      y <- as.numeric(as.factor(y)) - 1
    } else if (inherits(y, "factor")){
      y <- as.numeric(y) - 1
    } else {
      y <- as.numeric(y)
    }
  } else if (distribution == "coxph"){
    y <- ifelse(trueevent, y, -y)
  }

  if (any(is.na(y))){
    stop("missing values in the response are not allowed")
  }

  xd <- xgb.DMatrix(data = as.matrix(X), info = list(label = y))

  if (class(distribution) == "character"){
    if (length(unique(y)) == 2 | distribution %in% c("binomial", "bernoulli")){
      distribution <- "reg:logistic"
      metrics <- "logloss"
    } else if (distribution == "gaussian"){
      distribution <- "reg:squarederror"
      metrics <- "rmse"
    } else if (distribution == "huber") {
      distribution <- "reg:pseudohubererror"
      metrics <- "mphe"
    } else if (distribution == "quantile"){
      distribution <- function(preds, dtrain){
        qntregobj(preds, dtrain, quant = quant)
      }
      metrics <- list(function(preds, dtrain) qnterr(preds, dtrain, quant=quant))
    } else if (distribution == "multinomial"){
      distribution <- "multi:softprob"
      metrics <- "mlogloss"
      num_class <- length(unique(y))
    } else if (distribution == "poisson"){
      distribution <- "count:poisson"
      metrics <- "poisson-nloglik"
    } else if (distribution == "coxph"){
      event <- deparse(substitute(event))
      if(event == "NULL"){
        stop("For coxph, you need to specify which column indicates an event (unquoted).")
      }
      distribution <- "survival:cox"
      metrics <- "cox-nloglik"
    } else {
      stop("you're going to have to implement that distribution before you can use it")
    }
  }

  params <- list(eta = learning.rate, max_depth = interaction.depth,
                 nthread = n.cores, subsample = bag.fraction,
                 colsample_bynode = col.fraction,
                 min_child_weight = n.minobsinnode,
                 gamma = leaf.penalty, alpha = weight.penalty.L1,
                 lambda = weight.penalty.L2,
                 objective = distribution,
                 eval_metric = metrics)
  if (distribution == "multi:softprob"){
    params <- c(params, num_class = num_class)
  }

  xm <- xgb.cv(xd, params=params, nfold=cv.folds, nrounds = n.trees,
               early_stopping_rounds = early.stopping.trees,
               weight = weight,
               print_every_n = print_every_n,
               prediction = TRUE,
               verbose = verbose)

  nms <- names(xm$evaluation_log)

  nms <- gsub("rmse", "error", nms)
  nms <- gsub("logloss", "error", nms)

  names(xm$evaluation_log) <- nms

  if (bestCViter(list(cv = xm)) == n.trees){
    msg <- "Not enough trees for CV convergence"
    if (fail.if.not.converged){
      stop(msg)
    } else {
      warning(msg)
    }
  }

  m <- xgb.train(xd, params=params, nrounds=n.trees, print_every_n = print_every_n,
                 verbose = verbose)

  res <- list(call = theCall, model = m, cv = xm, cv.folds = cv.folds,
              xnames = Xnames, X = X, y = y, data = data,
              xlevels=xlevels, n.trees = n.trees)
  class(res) <- "xgbm"
  res
}

