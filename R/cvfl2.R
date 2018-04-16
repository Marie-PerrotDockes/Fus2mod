#' Description of the function
#'
#' @param  response a vector response variable
#' @param  regressors a quantitative matrix of regressor
#' @param  group a vector with two levels. (The group of the ANCOVA)
#' @param  seq.a the parameters that indicate how much the coefficients will be fused. (The more a is small the more the coefficient will be supposed to fused). If NULL the sequence to test will be calculted automaticaly
#' @param  lambda if the user wants to use it owns values of lambdas
#' @param  mina if seq.a is NULL the minimum value of a that is use in seq.a ( a must be higher than 0.)
#' @param  nrep the number of times we perform each CV for the different value of a in seq.a
#' @return The ddl and the error of the lambda that minimise the CV error for differents value of a (we take the mean over nrep  replicats)
#' @examples
#' B <- c(1, -1, 1.5, 1.5, rep(0, 6), 2, 0, 2, 0)
#'group <- c(rep('M1', 10), rep('M2', 10))
#'regressors <- matrix(rnorm(6*20), ncol = 6)
#'X  <- model.matrix(~group + group:regressors - 1)
#'y <- X%*%B + rnorm(20)
#'y <- scale(y)
#'mod <- cv.fl2(y, regressors, group)
#' @export
cv.fl2 <- function(response, regressors, group, seq.a = NULL,
                   lambda = NULL, mina = 0.1, nfold = 5, nrep= 10,
                   nb.cores = 3, plot = TRUE){
  p <- ncol(regressors)
  X  <- model.matrix(~group + group:regressors - 1)
  if(is.null(colnames(regressors))) colnames(regressors) <- paste0('regressors', 1:p)
  X2 <- cbind(X,regressors)
  if (is.null(seq.a)) seq.a <- seq(mina, ((3 * p + 2) / (2 * p)) - 0.001,len=20)
  min.cv <- do.call(rbind, mclapply(1:nrep, function(i) {
    foldid <- sample(rep(seq(nfold), length = nrow(X2)))
    do.call(rbind, lapply(seq.a, function(a) {
      b <- (3 * p - a * p + 2) / (2 * p )
      mod <- cv.glmnet(x = X2, y = response, foldid = foldid, standardize = FALSE,
                       penalty.factor = c(0, 0, rep(b, (ncol(X) - 2)), rep(a, ncol(regressors))),
                       intercept = F, lambda = lambda )
      cvm  <- min(mod$cvm)
      Coef <- coef(mod, s ='lambda.min')
      ddl <- length(unique(as.numeric(Coef)))
      return(cbind(cvm, ddl, a))

    }))}, mc.cores=nb.cores)
  )
  dplot <- min.cv %>% as.data.frame()  %>%   gather(key = 'Criterion', value = 'value', -a) %>%
    group_by(a, Criterion) %>% dplyr::summarize( mean = mean(value), sd  = sd(value) / sqrt(nrep))

  if(plot){

    print(ggplot(dplot, aes(x=a, y=mean, ymin=mean-sd, ymax=mean+sd)) + geom_smooth(stat="identity") + theme_bw() + facet_wrap(~Criterion, scales='free'))
  }


  return(dplot)
}
