#' Description of the function
#'
#' @param  response a vector response variable
#' @param  regressors a quantitative matrix of regressor
#' @param  group a vector with two levels. (The group of the ANCOVA)
#' @param  a the parameters that indicate how much the coefficients will be fused
#' @param  lambda if the user wants to use it owns values of lambdas
#' @return The coefficients of the fused lasso ANCOVA for the different value of lambda
#' @examples
#' B <- c(1, -1, 1.5, 1.5, rep(0, 6), 2, 0, 2, 0)
#'group <- c(rep('M1', 10), rep('M2', 10))
#'regressors <- matrix(rnorm(6*20), ncol = 6)
#'X  <- model.matrix(~group + group:regressors - 1)
#'y <- X%*%B + rnorm(20)
#'y <- scale(y)
#'mod1 <- bic.fl2(y, regressors, group)
#'mod <- stab.fl2_fixa(y, regressors, group, a = mod1[[3]])
#'coef(mod, s='lambda.min')
#' @export
stab.fl2_fixa <- function(response, regressors, group, a,
                   lambda = NULL,  nrep= 1000,
                   nb.cores = 3, plot = TRUE){
  p <- ncol(regressors)
  X  <- model.matrix(~group + group:regressors - 1)
  if(is.null(colnames(regressors))) colnames(regressors) <- paste0('regressors', 1:p)
  X2 <- cbind(X,regressors)
  b <- (3 * p - a * p + 2) / (2 * p )
  q <- max(5, min(2 * round(nrow(X)/log(ncol(X))/10)*10, ncol(X)))
  mod <- stabsel(X2, y, q=q, PFER=1, B = nrep,
                 fitfun = glmnet.lasso,
                 args.fitfun = list( penalty.factor=c(0, 0, rep(b, (ncol(X) - 2)), rep(a, ncol(regressors))),
                                     standardize=FALSE), mc.cores = nb.cores, mc.preschedule = TRUE)

  return(mod)
}
