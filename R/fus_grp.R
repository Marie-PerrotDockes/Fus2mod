group_fl2 <- function(response, regressors, group, modal, seq.a = NULL,
                      lambda = NULL, mina = 0.1, nfold = 5, nrep= 10,
                      nb.cores = 3, plot = TRUE){
  p <- ncol(regressors)
  X  <- model.matrix(~modal + modal:regressors - 1)
  if(is.null(colnames(regressors))) colnames(regressors) <- paste0('regressors', 1:p)
  X2 <- cbind(X,regressors)
  group <- c(group, group)
  b <- (3 * p - a * p + 2) / (2 * p )
  mod <- gglasso(x = X2, group= group, loss = 'ls', y = response, standardize = FALSE,
                pf = c(0, 0, rep(b, (ncol(X) - 2)), rep(a, ncol(regressors))),
                intercept = F, lambda = lambda )
  return(mod)

}
