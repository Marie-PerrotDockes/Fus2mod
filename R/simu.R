# require(parallel)
# n               <- 400
# p               <- 20
# t               <- rep(c("AC","AD","BC","BD"))
# loc             <- paste0("loc", 1:20)
# Loc_T           <- matrix(rnorm(n * p * 4), ncol = p * 4)
# colnames(Loc_T) <- paste(rep(loc, each = 4), rep (t, 20), sep=":")
# C               <- as.factor(sample(rep(c("I","NI") , round(n / 2)), n))
# X               <- model.matrix(~C + Loc_T:C -1)
# X               <- cbind(X, Loc_T)
# b <- rep(0, (ncol(X) ))
# names(b) <- colnames(X)
# G1 <- c(1, 2)
# G2 <- c(3, 4)
# G3 <- c(5, 6)
# sapply(c(G1, G2, G3), function(i){
#   val <- sample(1:4, 4)
#   if(i %in% G3)  b[grepl(paste0("loc", i,":"), names(b)) & grepl("CI", names(b))]  <<- val
#   if(i %in% G2)  b[grepl(paste0("loc", i,":"), names(b)) & grepl("CNI", names(b))] <<- val
#   if(i %in% G1){
#     b[grepl(paste0("loc", i,":"), names(b)) & !(grepl("CI", names(b)) | grepl("CNI", names(b)))]  <<- val
#   }
# })
# b[1:2] <- c(1, 2)
# b
#
#
#
# simu <-do.call(rbind, lapply(1:20, function(x){
# y <- X %*% b + rnorm(n)
# Loc_T           <- matrix(rnorm(n * p * 4), ncol = p * 4)
# colnames(Loc_T) <- paste(rep(loc, each = 4), rep (t, 20), sep=":")
# mods <- fus_grp(response = y, regressors = Loc_T, grp =NULL, group_by = paste0(loc,":"), uniq_modal = NULL, modal = C, seq.a = NULL,
#                 lambda = NULL, mina = 0.1, nfold = 5, nrep= 10,
#                 nb.cores = 3, plot = TRUE, a = 1, nlambda = 100)
#
# Mod <- mods[[1]]
#
# Mod$ROC <- lapply(mods[[1]]$Coefficients, function(b_hat){ A <- t( apply(b_hat, 2, TP_FP, b = b))
# colnames(A) <-c("TPR", "FPR")
# as.data.frame(A)})
#
#
# Roc <-select(Mod,value,ROC) %>% unnest()
#
# return(Roc)}))
#
#
# head(simu)
# p <- ggplot(data = simu, aes( x= FPR,y = TPR, color = as.factor(value) ))+geom_smooth()
# p
