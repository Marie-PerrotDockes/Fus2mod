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
#'mod <- bic.fl2(y, regressors, group)
#'print(mod[[2]])
#' @export
bic.fl2 <- function(response, regressors, group, seq.a = NULL,
                    lambda = NULL, mina = 0.1,
                     plot = TRUE){
  p <- ncol(regressors)
  X  <- model.matrix(~group + group:regressors - 1)
  if(is.null(colnames(regressors))) colnames(regressors) <- paste0('regressors', 1:p)
  X2 <- cbind(X,regressors)
  if (is.null(seq.a)) seq.a <- seq(mina, ((3 * p + 2) / (2 * p)) - 0.001,len=20)
  alla <- lapply(seq.a, function(a) {
    b <- (3 * p - a * p + 2) / (2 * p )
    mod <- glmnet(x = X2, y = response, standardize = FALSE,
                  penalty.factor = c(0, 0, rep(b, (ncol(X) - 2)), rep(a, ncol(regressors))),
                  intercept = F, lambda = lambda )
    Y_hat  <- X2 %*% mod$beta
    errors_sq <- apply((Y_hat - response) ^ 2, 2, sum)
    ddl <- apply(mod$beta, 2, function(x){length(unique(as.numeric(x)))})
    BIC <- errors_sq + log(nrow(regressors)) * ddl
    mini <- which.min(BIC)
    return(list(cbind.data.frame(BIC = BIC[mini], E_sq = errors_sq[mini],
                                 DDL = ddl[mini],a = a), mod$beta[,mini]))

  })
  min.bic <- do.call(rbind, lapply(1:length(seq.a), function(i){alla[[i]][[1]]}))
  if(plot){
    print(ggplot(min.bic, aes(x=a, y=BIC)) + geom_line()+ geom_point() + theme_bw())
  }
  K <- nlevels(as.factor(group))
  Coef <- alla[[which.min(min.bic$BIC)]][[2]]
  Cr <- c(0,0,Coef[rep(((K * (p + 1)) +1) : ((K * (p + 1)) + p), each = 2)])
  Cf <- Coef[1:(K * (p + 1))]+ Cr

  return(list(min.bic = min.bic, Coef = Cf, a = seq.a[which.min(min.bic$BIC)]))
}
