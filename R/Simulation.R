# rm(list=ls())
# source("R/fus_grp.R")
# n               <- 40
# p               <- 20
# t               <- rep(c("AC","AD","BC","BD"))
#
# require(parallel)
# require(tidyverse)
# require(gglasso)
# require(MASS)
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
# TP_FP <- function(x, b){
#   # x <- x[seq(3, length(x), 4)]
#   # b <- b[seq(3, length(b), 4)]
#   TP  <- sum( x != 0 & b != 0)
#   FP  <- sum( x != 0 & b == 0)
#   TN  <- sum( x == 0 & b == 0)
#   FN  <- sum( x == 0 & b != 0)
#   TPR <- TP / (TP + FN)
#   FPR <- FP / (FP + TN)
#   return(c(TPR, FPR))
# }
#
#
# Simus <- lapply(1:100, function(lalal){
#   loc             <- paste0("loc", 1:20)
#   Loc_T           <- matrix(rnorm(n * p * 4), ncol = p * 4)
#   colnames(Loc_T) <- paste(rep(loc, each = 4), rep (t, 20), sep=":")
#   C               <- as.factor(sample(rep(c("I","NI") , round(n / 2)), n))
#   X               <- model.matrix(~C + Loc_T:C -1)
#   X               <- cbind(X, Loc_T)
#   y               <- X %*% b + rnorm(n)
#   do.call(rbind, mclapply(seq(1, 1.5, len =10), function(a){
#   gfl2 <- group_fl2_adaptive(response = y, regressors = Loc_T, group = NULL, group_by = paste0(loc,":") , modal = C, seq.a = NULL,
#                     lambda = NULL, mina = 0.1, nfold = 5, nrep= 10,
#                     nb.cores = 3, plot = TRUE, a =a)
#   Error <- log(apply(gfl2$Est, 2, function(x){sum((x - y)^2)})/n)
#   Norm_Bet <- apply(gfl2$mod$beta, 2, function(x){
#     tapply(x, gfl2$grp_ord,function(x){sqrt(sum(x^2))})})
#   # Df <- apply(gfl2$mod$beta, 2, function(x){length(unique(x)) + 1}) / 4
#   dj <- 4
#   Df <- apply(Norm_Bet, 2, function(x){sum(x!=0) + (sum(x / gfl2$bOLSgr)*(dj-1)) })
#     #http://hansheng.gsm.pku.edu.cn/pdf/2008/agLasso.pdf
#   BIC <- Error + (Df*log(n) /n)
#   b_hat <- gfl2$mod$beta
#   b <- b[gfl2$ord]
#   Roc <- t( apply(b_hat, 2, TP_FP, b=b) )
#   colnames(Roc) <- c("TPR", "FPR")
#   Roc <- cbind.data.frame(Roc, a = a, BIC =BIC, nb = 1:length(BIC))
#
#   return(Roc)
#   }))
#
# })
# ROC <- do.call(rbind, Simus)
# ROC <- Reduce ( "+",Simus) /100
#  head(ROC)
#
#  p <- ggplot (data = ROC, aes(x= FPR, y =TPR, color = as.character(a))) + geom_smooth()+
#                 theme_bw()
#  p
#
#  p <- ggplot (data = ROC, aes(x= nb, y =BIC, color = as.character(a))) + geom_smooth()+
#    theme_bw()
#  p
