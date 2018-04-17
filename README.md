Install and load the package
----------------------------

``` r
devtools::install_github("Marie-PerrotDockes/Fus2mod")
```

    ## Skipping install of 'Fus2mod' from a github remote, the SHA1 (795fe126) has not changed since last install.
    ##   Use `force = TRUE` to force installation

``` r
require(Fus2mod)
```

    ## Loading required package: Fus2mod

    ## Loading required package: Matrix

    ## Loading required package: glmnet

    ## Loading required package: foreach

    ## Loaded glmnet 2.0-16

    ## Loading required package: parallel

    ## Loading required package: tidyverse

    ## ── Attaching packages ────────────────────────────────────────────────── tidyverse 1.2.1 ──

    ## ✔ ggplot2 2.2.1     ✔ purrr   0.2.4
    ## ✔ tibble  1.4.2     ✔ dplyr   0.7.4
    ## ✔ tidyr   0.8.0     ✔ stringr 1.3.0
    ## ✔ readr   1.1.1     ✔ forcats 0.3.0

    ## ── Conflicts ───────────────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ purrr::accumulate() masks foreach::accumulate()
    ## ✖ tidyr::expand()     masks Matrix::expand()
    ## ✖ dplyr::filter()     masks stats::filter()
    ## ✖ dplyr::lag()        masks stats::lag()
    ## ✖ purrr::when()       masks foreach::when()

    ## Loading required package: stabs

Dataset Simulation
------------------

We start by generate a toy exemple, where we have two modalities H and I, 10 regressors and 30 individuals. For now we only have one response. We generate our dataset such that the 3 first regressors have no effect on the response , then the regressors 4, 5 and 6 have an effect on the response onl y on the sample that comes from the modality I and the regressors 7 and 8 have an effect on the response for the modality N but not for the modlity I.

``` r
n <- 40
p <- 10
K <- 2
q <- 1


# B <- Simul_B(p, q, s, f, k)
B <- c(0, 0,
       0, 0,
       0, 0,
       0, 0,
       0, 3,
       0, 2.5,
       0, 2,
       2.2, 0,
       1.8, 0,
       2, 2,
       2, 2
       )

group <- c(rep("H", n / 2), rep("I", n / 2))
K <- nlevels(as.factor(group))
regressors     <- matrix(rnorm(n * p), ncol = (length(B) / K - 1))

X             <- model.matrix(~group + group:regressors - 1)
y <- as.matrix(X %*% B + matrix(rnorm(n * q ), ncol = q))

y <- scale(y)
```

Model Selection
---------------

Here we foccus on the model *Y* = *X**B* + *E*, where *Y*, *B* and *E* are vector and *X* is a one-way ANCOVA design matrix. We have two objectifs : we want to find which collumns of *X* can explain the response *y* in the first hand. In the second hand we want to observe if the differents modalities influence the values of the coefficient of the regressors on the response, in other terms if the collumn group1:regressorsi will have the same coefficient than the collum group2:regressorsi. To do that we propose to apply a lasso criterion to the model : *Y* = *X*<sub>2</sub>*B* + *E* where *X*2 is the concatenation of *X* and the matrix of the regressors. To be abble to arrange the level of fusion we propose to put weight on the penalties a weight *b* if the coefficients is the coefficient of a couple regressors modalities and a weight *a* if the coefficients is the coefficients of a whole regressors. Like 2*p**b* + *p**a* must be equal to 3*p* + 2 fixing *a* will give us a value of *b*. Thus, the more *a* is small the more the model will encourage the coefficient of two different modalities for the same regressors to be the same. For more detail about this model we confer the reader to the file "2\_modalities.pdf".

We first propose a Cross-validation step. For different values of *a* we will apply a 5-fold CV on our data and keep the minimal error (cvm)and the degree of freedom (ddl) that match this minimal error. We will do that *n**r**e**p* = 10 times for each value of *a* the mean on this 10 replicats are display bellow.

``` r
 CV <- cv.fl2(response = y, regressors = regressors, group = group, mina = 0.1, nfold = 5, nrep= 10,
                   nb.cores = 3, plot = TRUE)
```

![](Vignettes_files/figure-markdown_github/unnamed-chunk-3-1.png)

By watching this plot we propose to take the *a* that minimise the degree of freedom because it is the lower (or really close to ) the lower error of prediction.

``` r
ddl <- CV[CV$Criterion =="ddl", ]
a <- ddl[which.min(ddl$mean), "a"]
a
```

    ## # A tibble: 1 x 1
    ## # Groups:   a [1]
    ##       a
    ##   <dbl>
    ## 1  1.05

If we want we can add a stability selection step to keep the more stable variable. In order to avoid that sometimes the model select as non null value the coefficient of the whole regressors and sometimes it keep it only for one modality we perform the fusion before the stability selection and then apply the stability selection step.

``` r
source('~/Documents/Multivar_selec/Multivar_selec/FusedLasso/Fus2mod/R/stab.fl2_fixa.R', echo=TRUE)
```

    ## 
    ## > stab.fl2_fixa <- function(response, regressors, group, 
    ## +     a, lambda = NULL, nrep = 1000, nb.cores = 3, plot = TRUE) {
    ## +     p <- ncol(regressors .... [TRUNCATED]

``` r
stab <- stab.fl2_fixa(response = y, regressors, group, a,
                   lambda = NULL,  nrep= 1000,
                   nb.cores = 3, plot = TRUE)
stab
```

    ##  Stability Selection with unimodality assumption
    ## 
    ## Selected variables:
    ##             groupH             groupI groupI:regressors4 
    ##                  1                  2                  7 
    ##        regressors9 
    ##                 18 
    ## 
    ## Selection probabilities:
    ## groupI:regressors8 groupI:regressors7 groupI:regressors3 
    ##             0.2285             0.2330             0.2390 
    ## groupI:regressors2 groupH:regressors3 groupH:regressors2 
    ##             0.2565             0.2900             0.3035 
    ## groupH:regressors5        regressors1 groupH:regressors6 
    ##             0.3170             0.5660             0.5960 
    ##        regressors4 groupH:regressors8 groupH:regressors7 
    ##             0.6375             0.8460             0.9260 
    ## groupI:regressors6       regressors10 groupI:regressors5 
    ##             0.9275             0.9560             0.9575 
    ## groupI:regressors4        regressors9             groupH 
    ##             0.9870             0.9955             1.0000 
    ##             groupI 
    ##             1.0000 
    ## 
    ## ---
    ## Cutoff: 0.9765; q: 14; PFER (*):  0.989 
    ##    (*) or expected number of low selection probability variables
    ## PFER (specified upper bound):  1 
    ## PFER corresponds to signif. level 0.0521 (without multiplicity adjustment)

``` r
par(mf.row=c(1,2))
```

    ## Warning in par(mf.row = c(1, 2)): "mf.row" n'est pas un paramètre graphique

``` r
plot(stab, type="maxsel", main="Maximum selection frequencies")
```

![](Vignettes_files/figure-markdown_github/unnamed-chunk-6-1.png)

``` r
plot(stab, type="path", main="Stability paths")
```

![](Vignettes_files/figure-markdown_github/unnamed-chunk-6-2.png)

``` r
sel <-names(stabsel(stab, cutoff=0.8)$selected)[-c(1,2)]
sel 
```

    ## [1] "groupI:regressors4" "groupI:regressors5" "groupI:regressors6"
    ## [4] "groupH:regressors7" "groupH:regressors8" "regressors9"       
    ## [7] "regressors10"

``` r
names(B) <- colnames(X)

B_long <- c(0, 0,
       0, 0,
       0, 0,
       0, 0,
       0, 1.4,
       0, 1.5,
       0, 2,
       1.2, 0,
       1.8, 0,
       0, 0,
       0, 0, 
       rep(0 ,8),1,1
       )
colnames(regressors) <- paste0('regressors', 1:p)
names(B_long) <- c(colnames(X), colnames(regressors))
names(B_long[B_long != 0])
```

    ## [1] "groupI:regressors4" "groupI:regressors5" "groupI:regressors6"
    ## [4] "groupH:regressors7" "groupH:regressors8" "regressors9"       
    ## [7] "regressors10"

If we take a threshold of 0.8, we have a True Positive Rate equal to 1 and a False Positive Rate equal to 0.
