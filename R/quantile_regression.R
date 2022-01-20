qntregobj <- function(preds, dtrain, quant){
  y <- getinfo(dtrain, "label")
  r <- preds - y

  grad <- ifelse(r > 0, 1 - quant, -quant)
  hess <- rep(1, length(grad))

  list(grad = grad, hess = hess)
}

qnterr <- function(preds, dtrain, quant){
  y <- getinfo(dtrain, "label")
  r <- preds - y

  err <- mean(ifelse(r < 0, -quant * r, (1 - quant) * r))

  list(metric="error", value=err)
}


qnt.5 <- function(preds, dtrain){
  qntregobj(preds, dtrain, quant=.5)
}

