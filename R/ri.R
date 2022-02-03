#' Relative influence of predictors in an xgbm model
#' @param x An object of class 'xgbm'.
#' @details Relative influence is referred to as relative importance in the
#'   xgboost documentation, but I prefer 'influence' because it keeps the idea
#'   distinct from Breiman's variable importance. Also, it should be computed
#'   as the sum over all splits on a variable (according to the xgboost
#'   documentation), but isn't when there are factors in the model. So here
#'   we add the values associated with factors back together. Finally,
#'   the documentation for \code{xgb.importance} (for which this is a
#'   wrapper) states that you need to give it an integer vector of tree numbers
#'   indexed at 0. This function takes care of it all by assuming you
#'   used cross-validation.
#'
#' @export
relative.influence <- function(x){
  if (!inherits(x, "xgbm")){
    stop("x should have class 'xgbm'")
  }
  cv.trees <- bestCViter(x) - 1

  o <- as.data.frame(xgb.importance(model = x$model, trees = 0:cv.trees))

  # Need to lump together RI from levels of any factors here
  if (length(x$xlevels) > 0){
    for (i in 1:length(x$xlevels)){
      fns <- paste0(names(x$xlevels)[i], x$xlevels[[i]])

      if (any(fns %in% o$Feature)){

        so <- data.frame(Feature = names(x$xlevels)[i],
                         Gain = sum(o$Gain[o$Feature %in% fns]),
                         Cover = sum(o$Cover[o$Feature %in% fns]),
                         Frequency = sum(o$Frequency[o$Feature %in% fns]),
                         stringsAsFactors = FALSE) %>% na.omit()

        if (nrow(so) > 0) o <- bind_rows(o, so) %>% filter(!(Feature %in% fns))
      }
    }
  }

  # Remove bloody backticks
  o$Feature <- ifelse(substring(o$Feature, 1, 1) == "`",
                      substring(o$Feature, 2, nchar(o$Feature) - 1),
                      o$Feature)

  o <- arrange(o, desc(Gain))
  o$Gain <- o$Gain / max(o$Gain)

  class(o) <- "relative.influence"
  o
}

#' @method print relative.influence
#' @export
print.relative.influence <- function(x, n=6, ...){
  x <- unclass(x)
  d <- as.data.frame(x)

  if (n > nrow(d)) n <- nrow(d)

  print(d[1:n, ])
  if (nrow(d) > n){
    cat(paste("...plus", nrow(d) - n, "other rows."))
  }
  invisible()
}

#' @method head relative.influence
#' @export
head.relative.influence <- function(x, ...){
  x <- as.data.frame(unclass(x))
  head(x, ...)
}

#' @method as.data.frame relative.influence
#' @export
as.data.frame.relative.influence <- function(x, row.names=NULL, optional=FALSE, ...){
  as.data.frame(unclass(x), stringsAsFactors = FALSE)
}

#' @method ggplot relative.influence
#' @export
ggplot.relative.influence <- function(data, mapping=NULL, ..., n=30){
  data <- as.data.frame(unclass(data), stringsAsFactors=FALSE)

  if (nrow(data) > n) data <- data[1:n, ]

  data$Feature <- ifelse(substring(data$Feature, 1, 1) == "`",
                         substring(data$Feature, 2),
                         data$Feature)

  data$Feature <- ifelse(substring(data$Feature, nchar(data$Feature), nchar(data$Feature)) == "`",
                         substring(data$Feature, 1, nchar(data$Feature) - 1),
                         data$Feature)

  data$Feature <- factor(data$Feature, levels = rev(data$Feature))

  ggplot(data, aes(Gain, Feature)) +
    geom_point(color="blue") +
    xlab("Relative influence") +
    ylab("")
}

#' @method plot relative.influence
#' @export
plot.relative.influence <- ggplot.relative.influence
