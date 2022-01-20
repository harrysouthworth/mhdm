#' Gradient Boosting wrappers for xgboost
#'
#' Wrappers for xgboost functions, in the hope of making that package more
#' friendly to statisticians, especially those used to using the gbm package.
#' Also includes helper functions for data mining of medical history or
#' adverse event data.
#'
#' @name mhdm-package
#' @aliases mhdm-package
#' @docType package
#'
#' @importFrom utils head tail
#' @importFrom stats model.frame model.matrix model.response na.omit na.pass
#'   predict setNames
#' @importFrom xgboost xgb.DMatrix xgb.cv xgb.train
NULL
