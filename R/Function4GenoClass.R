#Author : Charlotte Brault


####################  plotResponseVsGclasses    ###############################
#Format : 
# - clgeno = genotypic class coded in updated Joinmap format
#            individuals in columns / markers in rows
#
# - response = phenotypes in column pheno
#          trait in column trait
#          scenario in column scenario
#          year in column year
#
# - marker = marker in clgeno to select
# ylab : y label to plot : whether BLUP or phenotypes is plotted
# plot.type : whether boxplot or beanplot
# NA.value : how are coded missing values in clgeno data frame
# cat.response = name of the blup / pheno column
#'
#' @param clgeno data frame, individuals are in column (colnames = geno), markers are in row (rownames = marker name)
#' @param marker a specific rowname from clgeno
#' @param NA.value how are coded NA in the clgeno data frame ? Defaut "--"
#' @param response data frame with 
#'  - the response : colname = cat.response
#'  - the  trait : colname = trait
#'  - the scenario : colname = scenario
#'  - yhe year : colname = year
#' @param cat.response category of response (use as ylab in the plot)
#' @param trait which trait is considered in response 
#' @param scenario which scenario is considered in response
#' @param year which year is considered in response
#' @param plot.type whether "boxplot" or "beanplot" to plot the results
#' @param ylim extrem limits of the y axis, for example c(0,0.8)
#' @param logy should the y axis be transformed in logarithm ?
#' @param verbose = 1 by default
#'
#' @return plot of the response (for a trait -scenario - year) according the genotypic classe 
#' at the marker
#' 
#' 
plotResponseVsGclasses <- function(clgeno, marker, NA.value = "--", response, cat.response = "blup", trait, scenario,
                                   year, plot.type="boxplot", ylab="blup", ylim=NULL, logy="", verbose=1 ){
  
  stopifnot(cat.response %in% colnames(response))
  stopifnot("ind" %in% colnames(response))
  stopifnot("year" %in% colnames(response))
  stopifnot("trait" %in% colnames(response))
  stopifnot("scenario" %in% colnames(response))
  stopifnot(!all(!rownames(clgeno) == marker)) #check if the marker exists
  
  # Verifiy if differences in individuals 
  setd1 = c(setdiff(colnames(clgeno),response$ind))
  setd2 = c(setdiff(response$ind,colnames(clgeno)))
  
  clgeno2 = clgeno[, !(names(clgeno) %in% setd1)]
  response2 = response[ ! (response$ind %in% setd2),]
  
  
  #Select the marker and the response to plot
  clgeno2 = as.matrix(clgeno2[rownames(clgeno2) == marker,])
  response2 = response2[response2$trait == trait & response2$year == year & 
                          response2$scenario == scenario,]
  
  stopifnot(! all(is.na(response2)))
  
  clgeno2[clgeno2 == NA.value] <- NA
  
  if(plot.type == "boxplot"){
    boxplot(response2[[cat.response]] ~ t(clgeno2), ylim=ylim, xlab = "genotypic class", ylab = cat.response, 
            main = paste0(trait,"-",scenario,"-",year,", ",marker), las=1, varwidth=TRUE)
  }
  if(plot.type == "beanplot"){
    beanplot(response2[[cat.response]] ~ t(clgeno2), ylim=ylim, xlab = "genotypic class", ylab = cat.response, 
             main = paste0(trait,"-",scenario,"-",year,", ",marker), las=1, log=logy,bw="nrd0", 
             what=c(1,1,1,0), col=c("grey",0,0,"red"), border=NA)
  }
  
}


################################ MeanResponseOnGenoclass ##############################


#' For each allele i, compute its additive effect (alpha_i), where i is in {a,b,c,d} if segregation is abxcd, etc.
#' Moreover, for each parent (par1 of genotype ab and par2 of genotype cd), 
#' compute half of its "average effect of allelic substitution" (alpha).
#' Ref: Lynch & Walsh (1998) and Doligez et al (2013)
#' If <abxcd>:
#' Af = half.alpha.par1 = (alphaA - alphaB) / 2 = ((muAD + muAC) - (muBD + muBC)) / 4
#' Am = half.alpha.par2 = (alphaC - alphaD) / 2 = ((muAC + muBC) - (muAD + muBD)) / 4
#' Indeed, muAD = alphaA + alphaD
#' If <efxeg>:
#' Af = (alphaE - alphaF) / 2 = ((muEG + muEE) - (muFG + muFE)) / 4
#' Am = (alphaE - alphaG) / 2 = ((muEE + muFE) - (muEG + muFG)) / 4
#' If <hkxhk>:
#' Af = Am = (alphaH - alphaK) / 2 = (muHH - muKK) / 4
#' If <lmxll>:
#' Af = (alphaL _ alphaM) / 2 = (muLL - muLM) / 2
#' Am = 0, no effect of allelic substitution for homozygous parents
#' If <nnxnp>:
#' Af = 0
#' Am = (alphaN - alphaP) / 2 = (muNN - muNP) / 2
#' 
#' @param clgeno data frame, individuals are in column (colnames = geno), markers are in row (rownames = marker name)
#' @param marker a specific rowname from clgeno
#' @param response a data frame containing the response and the individuals (it must be in a column)
#' @param cat.response category of response (use as ylab in the plot)
#' @param NA.value  how are coded NA in the clgeno data frame ? Defaut "--" 
#'
#' @return data frame with a column Am for the male additive effect and Af for the female additive effect
#' 
#' @seealso AD2AlphaDelta, AlphasDeltas2AD

MeanResponseOnGenoclass = function(clgeno, marker, response, cat.response, NA.value = "--"){
  stopifnot("ind" %in% colnames(response))
  stopifnot(!all(!rownames(clgeno)==marker)) #check if the marker exists
  
  # Verifiy if differences in individuals 
  setd1 = c(setdiff(colnames(clgeno),response$ind))
  setd2 = c(setdiff(response$ind,colnames(clgeno)))
  
  clgeno2 = clgeno[, !(names(clgeno) %in% setd1)]
  response = response[ ! (response$ind %in% setd2),]
  
  #Select the marker and the response to plot
  clgeno2 = t(as.matrix(clgeno2[rownames(clgeno2) == marker,]))

  clgeno2[clgeno2 == NA.value] <- NA
  response$clgeno <- NA
  
  for (i in 1: nrow(response)){
    idx = which(rownames(clgeno2) == response$ind[i])
    response$clgeno[i] <- clgeno2[idx]
  }
  
  mgeno = as.matrix(tapply(response[[cat.response]], INDEX = response$clgeno, 
                                FUN = mean, na.rm = T))

  sdclgeno = as.matrix(tapply(response[[cat.response]],INDEX = response$clgeno,
                              FUN = sd, na.rm = T))
  
  if(clgeno[rownames(clgeno) == marker,"seg"] == "<abxcd>"){ # 4 alleles 
  # Additive effect parent 1&2
    #' Af = half.alpha.par1 = ((muAD + muAC) - (muBD + muBC)) / 4
    #' Am = half.alpha.par2 = ((muAC + muBC) - (muAD + muBD)) / 4
    A1 = ((mgeno[rownames(mgeno) == "ad",] + mgeno[rownames(mgeno) == "ac",]) -
            (mgeno[rownames(mgeno) == "bd",] + mgeno[rownames(mgeno) == "bc",])) /4
    A2 = ((mgeno[rownames(mgeno) == "ac",] + mgeno[rownames(mgeno) == "bc",]) - 
            (mgeno[rownames(mgeno) == "ad",] + mgeno[rownames(mgeno) == "bd",])) / 4
  # Dominance effect
    # D = ((muAC + muBD) - (muBC + muAD)) / 4
    D = ((mgeno[rownames(mgeno) =="ac",] + mgeno[rownames(mgeno) =="bd",]) - 
           (mgeno[rownames(mgeno) =="bc",] + mgeno[rownames(mgeno) =="ad",]))/4
    
  } else if(clgeno[rownames(clgeno) == marker,"seg"] == "<efxeg>"){ # 3 alleles 
    # Additive effect parent 1&2
    #' Af = ((muEG + muEE) - (muFG + muFE)) / 4
    #' Am = ((muEE + muFE) - (muEG + muFG)) / 4
    A1 = ((mgeno[rownames(mgeno) == "eg",] + mgeno[rownames(mgeno) == "ee",]) - 
      (mgeno[rownames(mgeno) == "fg",] + mgeno[rownames(mgeno) == "ef",])) / 4
    A2 = ((mgeno[rownames(mgeno) == "ee",] + mgeno[rownames(mgeno) == "ef",]) - 
            (mgeno[rownames(mgeno) == "eg",] + mgeno[rownames(mgeno) == "fg",])) / 4
    # Dominance effect 
      # D = ((muEE + muFG) - (muEF + muEG)) / 4
    D = ((mgeno[rownames(mgeno) == "ee",] + mgeno[rownames(mgeno) == "fg",]) - 
           (mgeno[rownames(mgeno) == "ef",] + mgeno[rownames(mgeno) == "eg",])) / 4


  } else if(clgeno[rownames(clgeno) == marker,"seg"] == "<hkxhk>"){ # 2 alleles, both parents homozygous
    # Additive effect parent 1&2
    #' Af = Am = (muHH - muKK) / 4
    A1 = (mgeno[rownames(mgeno) =="hh",]  - mgeno[rownames(mgeno) =="kk",]) / 4
    A2 = A1
    # Dominance effect : D = ((muHH - muKK) - 2*muHK)/4
    D = ((mgeno[rownames(mgeno) =="hh",] + mgeno[rownames(mgeno) =="kk",]) - 
           2*mgeno[rownames(mgeno) =="hk",]) / 4

    
  } else if(clgeno[rownames(clgeno) == marker,"seg"] == "<lmxll>"){ # 2 alleles, 2nd parent homozygous
    # Additive effect prent 1&2
    #' Af = (muLL - muLM) / 2
    #' Am = 0
    A1 = (mgeno[rownames(mgeno) =="ll",] - mgeno[rownames(mgeno) =="lm",]) / 2
    A2 = 0
    # Dominance effect
    D = 0
    
  } else if (clgeno[rownames(clgeno) == marker,"seg"] == "<nnxnp>"){ # 2 alles, 1st parent homozygous
    # Additive effect parent 1&2
    #' Af = 0
    #' Am = (muNN - muNP) / 2
    A1 = 0 
    A2 = (mgeno[rownames(mgeno) =="nn",] - mgeno[rownames(mgeno) =="np",]) / 2
    #Dominance effect
    D = 0
  }
  A1 = as.numeric(A1)
  A2 = as.numeric(A2)
  D = as.numeric(D)
  
  out = list(seg = clgeno[rownames(clgeno) == marker,"seg"],
             mean.on.genoclass = mgeno,
             sd.on.genoclass = sdclgeno,
             half.alpha.par1 = A1, half.alpha.par2 = A2, Dominance = D )
  return(out)
}

#' 
#' @param seg segregation type
#' @param A1 half of the "average effect of allelic substitution" of parent 1
#' @param A2 half of the "average effect of allelic substitution" of parent 2
#' @param D dominance effect
#' @param constraints.deltas constraints on the deltas
#'
#' @return list of two named vectors
#' @seealso AD2AlphaDelta

AD2AlphaDelta = function(seg="<abxcd>", A1=0.05, A2=0.02, D=0, constraints.deltas=c(0.5, 0.5)){
  
  out = list(alphas=c(),
             deltas=c())

  if(seg == "<abxcd>"){
    alphaA = 0  
    alphaC = 0
    alphaB = -2 * A1 
    alphaD = -2 * A2
    out$alphas = c(A=alphaA, B=alphaB, C=alphaC, D=alphaD)
    
    deltaAC = constraints.deltas[1] * D
    deltaAD = 0
    deltaBC = 0
    deltaBD = constraints.deltas[2] * D
    out$deltas = c(AC=deltaAC, AD=deltaAD, BC=deltaBC, BD=deltaBD)
  } else if(seg == "<efxeg>"){
    alphaE = 0
    alphaF = -2*A1
    alphaG = -2*A2
    out$alphas = c(E=alphaE, F=alphaF, G=alphaG)
    
    deltaEF = 0
    deltaEG = 0
    deltaFG = D
    out$deltas = c(EF=deltaEF, EG=deltaEG, FG=deltaFG)
    
  } else if(seg == "<hkxhk>"){
    alphaH = 0
    alphaK = -2*A1
    out$alphas = c(H=alphaH, K=alphaK)
    
    deltaHK = 0
    out$deltas = c(HK=deltaHK)
    
  }  else if(seg == "<lmxll>"){
    alphaL = 0
    alphaM = -2*A1
    out$alphas = c(N=alphaN, P=alphaP)
    
    deltaNP = 0
    out$deltas = c(NP=deltaNP)

  }else if(seg == "<nnxnp>"){
    alphaN = 0
    alphaP = -2*A2
    out$alphas = c(N=alphaN, P=alphaP)
    
    deltaNP = 0
    out$deltas = c(NP=deltaNP)
  }

  return(out)
}

#' Convert  the alpha's and delta's into A1, A2 and D.
#' See Xie & Xu (1999) and Doligez et al (2013).
#' By definition:
#' A1 = (alphaA - alphaB) / 2
#' A2 = (alphaC - alphaD) / 2
#' We fix alphaA = alphaC = 0 as constraints.
#' Thus: 
#' alphaB = - 2 * A1 
#' alphaD = - 2 * A2
#' Also by definition:
#' D = deltaAC - deltaAD - deltaBC + deltaBD
#' @param seg segregation type
#' @param alphas named vector
#' @param deltas named vector
#'
#' @return a vector of additive effect male and female and dominance effect

AlphasDeltas2AD = function(seg, alphas, deltas){
  out = c()

  if(seg == "<abxcd>"){
    A1 = (alphas["A"] - alphas["B"]) / 2
    A2 = (alphas["C"] - alphas["D"]) / 2
    D = deltas["AC"] - deltas["AD"] - deltas["BC"] + deltas["BD"]
    
  }else if(seg =="<efxeg>"){
    A1 = (alphas["E"] - alphas["F"]) / 2
    A2 = (alphas["E"] - alphas["G"]) / 2
    D = deltas["FG"]
    
  }else if(seg =="<hkxhk"){
    A1 = (alphas["H"] - alphas["K"]) / 2
    A2 = (alphas["H"] - alphas["K"]) / 2
    D = 0
    
  }else if(seg == "<lmxll>"){
    A1 = (alphas["L"] - alphas["M"]) / 2
    A2 = 0
    D = 0
    
  }else if(seg == "<nnxnp>"){
    A1 = 0
    A2 = (alphas["N"] - alphas["P"]) / 2
    D = 0
  }
  out = c(A1=A1,A2=A2,D=D)
  return(out)
}




################### SimulPhenosBiparentalCross ##########################
# This function is only useful to produce plot of response which are looking like plotResponseVsGclasses
# with real data

#'
#' Phenotypes are simulated for offsprings from a bi-parental cross between outbred parents (AB x CD):
#' y = g + e where y is the phenotyic value, g the genotypic value and e the error.
#' See Xie & Xu (1999) and Doligez et al (2013).
#' By definition:
#' For a marker segregating <abxcd> :
#' A1 = (alphaA - alphaB) / 2
#' A2 = (alphaC - alphaD) / 2
#' We fix alphaA = alphaC = 0 as constraints.
#' Thus: 
#' alphaB = - 2 * A1 
#' alphaD = - 2 * A2
#' Also by definition:
#' D = deltaAC - deltaAD - deltaBC + deltaBD
#' We also fix DeltaAD = DeltaBC = 0 (see also param constraints.deltas).
#' 
#' For a marker segregating as <efxeg>: (A=E B=F C=E D=G)
#' alphaE = 0
#' alphaF = - 2 * A1
#' alphaG = - 2 * A2
#' deltaEG = deltaEF = 0
#' deltaFG = D
#' 
#' For the cross <hkxhk> :
#' AlphaH = 0
#' AlphaK = - 2 * A1
#' DeltaHH = 0
#' DeltaHK = D / 2
#' DeltaKK = D / 2
#' 
#' For the cross <lmxll>
#' AlphaL = 0
#' AlphaM = - 2 * A1 
#' DeltaLM = 0
#' 
#' For the cross <nnxnp>
#' AlphaN = 0
#' AlphaP = -2 * A1
#' DeltaNP = 0
#' 
#' Additive genetic variance according to Xie & Xu :
#' sigma2A = A1^2 +A2^2
#' and also dominance genetic variance :
#' sigma2D = (1/16) * D^2
#' 
#' 
#' 
#' @param clgeno data frame of marker genotypes with individuals in columns (colnames = geno, first column contains the marker segregation), markers in rows (rownames = marker name)
#' @param qtl a specific rowname from clgeno which will be the QTL (only one for the moment)
#' @param Gij mean response value on genotypic class,  given in the funtion MeanResponseOnGenoclass 
#' @param A1 half of the additive effect of the parent 1, given in the funtion MeanResponseOnGenoclass 
#' @param A2 half of the additive effect of the parent 2, given in the funtion MeanResponseOnGenoclass
#' @param D (half of?) the dominance effect
#' @param constraints.deltas by default DeltaAC = DeltaBD = D / 2, 
#' vector of weights for each Delta, the sum is equal to 1 
#' @param h2 heritability 
#' @param plot.type whether "beanplot" or "boxplot"
#' @param ylim default = NULL, limits of the y axis
  
#' @return a list of several elements :
#' - Delta : the dominance effect
#' - Alpha : half the additive effect of each parent
#' - data : a data frame containing the phenotype in column pheno, 
#' the genotypic effect in column g, 
#' the error term in column e,
#' the genotypic class at the QTL in column geno.qtl

SimulPhenosBiparentalCross = function (clgeno, marker, Gij, A1, A2, 
                                       D,constraints.deltas = c(0.5,0.5),
                                       h2,
                                       plot.type = "beanplot", ylim=NULL, logy = "") {

  stopifnot(sum(constraints.deltas) == 1)
  
  n = ncol(clgeno) - 2 # nb individuals
  geno.marker = clgeno[rownames(clgeno) == marker,3:ncol(clgeno)]
  row.names(geno.marker ) <- NULL
  colnames(geno.marker) <- NULL
  
  dat = data.frame(geno = colnames(clgeno)[-c(1,2)],
                   geno.qtl = 0,
                   g = 0,
                   pheno = 0)
  geno.marker = t(as.matrix(geno.marker))
  dat$geno.qtl <- geno.marker
  
  sigma2A = A1^2 + A2^2
  sigma2D = (1/16) * D^2
  sigma2G = sigma2A + sigma2D
  sigma2E = (sigma2G *(1-h2))/h2
  
  
##### AB x CD  
  if(clgeno[rownames(clgeno) == marker,"seg"] == "<abxcd>"){ # 4 alleles 

    AlphaA = 0  
    AlphaC = 0
    AlphaB = -2 * A1 
    AlphaD = -2 * A2
    
    DeltaAC = constraints.deltas[1] * D
    DeltaAD = 0
    DeltaBC = 0
    DeltaBD = constraints.deltas[2] * D
    Deltas = list(AC = DeltaAC, AD = DeltaAD, BC = DeltaBC, BD = DeltaBD)
    Alphas = list(A = AlphaA, B = AlphaB, C = AlphaC, D = AlphaD)
  
  
  for (i in 1:nrow (dat)){
    if(dat$geno.qtl[i] == 'ac'){
      dat$g[i] = mean(Gij) + AlphaA + AlphaC + DeltaAC
    } else if (dat$geno.qtl[i] == 'ad'){
      dat$g[i] = mean(Gij) + AlphaA + AlphaD + DeltaAD
    } else if (dat$geno.qtl[i] == 'bc'){
      dat$g[i] = mean(Gij) + AlphaB + AlphaC + DeltaBC
    } else if (dat$geno.qtl[i] == 'bd'){
      dat$g[i] = mean(Gij) + AlphaB + AlphaD + DeltaBD
     }
   } 
  }
  
##### EF x EG
  if (clgeno[rownames(clgeno) == marker, "seg"] == "<efxeg>"){ # 3 alleles
    AlphaE = 0
    AlphaF = - 2 * A1
    AlphaG = -2 * A2
    DeltaEF = 0
    DeltaEG = 0
    DeltaFG = D
    Deltas = list(EF = DeltaEF, EG = DeltaEG, FG = DeltaFG)
    Alphas = list(E = AlphaE, F = AlphaF, G = AlphaG)
    
    for (i in 1:nrow(dat)){
      if(dat$geno.qtl[i] == 'ee'){
        dat$g[i] = mean(Gij) + 2*AlphaE
      } else if (dat$geno.qtl[i] == 'ef'){
        dat$g[i] = mean(Gij) + AlphaE + AlphaF + DeltaEF
      } else if (dat$geno.qtl[i] == 'eg'){
        dat$g[i] = mean(Gij) + AlphaE + AlphaG + DeltaEG
      } else if (dat$geno.qtl[i] == 'fg'){
        dat$g[i] = mean(Gij) + AlphaF + AlphaG + DeltaFG
      }
    }
  }
  
  
##### HK x HK
  if(clgeno[rownames(clgeno) == marker,'seg'] == "<hkxhk>"){
    AlphaH = 0 # default
    AlphaK = - 2*A1
    DeltaHK = 0
    Deltas = list(HK = DeltaHK)
    Alphas = list(H = AlphaH, K = AlphaK)
    
    for (i in 1:nrow(dat)){
      if(dat$geno.qtl[i] == 'hh'){
        dat$g[i] = mean(Gij) + 2*AlphaH
      } else if (dat$geno.qtl[i] == 'hk'){
        dat$g[i] = mean(Gij) + AlphaH + AlphaK + DeltaHK
      } else if (dat$geno.qtl[i] == 'kk'){
        dat$g[i] = mean(Gij) + 2*AlphaK
      } 
    }
  }
  
##### LM x LL
  if(clgeno[rownames(clgeno) == marker,'seg'] == "<lmxll>"){
    AlphaL = 0 # default
    AlphaM = -2*A1
    DeltaLM = 0
    Deltas = list(LM = DeltaLM)
    Alphas = list(L = AlphaL, M = AlphaM)
    
    for (i in 1:nrow(dat)){
      if(dat$geno.qtl[i] == 'll'){
        dat$g[i] = mean(Gij) + 2*AlphaL
      } else if (dat$geno.qtl[i] == 'lm'){
        dat$g[i] = mean(Gij) + AlphaL + AlphaM + DeltaLM
      } 
    }
  }
  
##### NN x NP
  if(clgeno[rownames(clgeno) == marker,'seg'] == "<nnxnp>"){
    AlphaN = 0 # default
    AlphaP = -2 * A2
    DeltaNP = 0
    Deltas = list(NP = DeltaNP)
    Alphas = list(N = AlphaN, P = AlphaP)
    
    for (i in 1:nrow(dat)){
      if(dat$geno.qtl[i] == 'nn'){
        dat$g[i] = mean(Gij) + 2*AlphaN
      } else if (dat$geno.qtl[i] == 'np'){
        dat$g[i] = mean(Gij) + AlphaN + AlphaP + DeltaNP
      } 
    }
  }
  dat$geno.qtl[dat$geno.qtl == '--'] <- NA
  e = rnorm(n = n, mean = 0, sd = sqrt(sigma2E))

  # y = g + e
  dat$pheno = dat$g + e
  
  if (plot.type == 'boxplot') {
    boxplot(pheno ~ geno.qtl,
            data = dat,
      ylim = ylim,
      xlab = "genotypic class",
      ylab = "simulated blup",
      main = paste0("H2 = ", h2, "-", marker),
      las = 1,
      varwidth = TRUE)
  } else if (plot.type == 'beanplot') {
    beanplot(pheno ~ geno.qtl,
             data = dat,
      ylim = ylim,
      xlab = "genotypic class",
      ylab = "simulated blup",
      main = paste0("H2 = ", h2, "-", marker),
      las = 1, log = logy, bw="nrd0", 
      what=c(1,1,1,0), col=c("grey",0,0,"red"), border=NA)
  }
  
  out = list(Delta = Deltas, Alpha = Alphas,
             # y = y, g = g, e = e,
             data = dat)
  return(out)
}



############################## Calculate Beta  ########################################


#' CalcBeta
#'
#' @param seg segregation type
#' @param marker marker name with an effect on phenotype
#' @param X design matrix of genotypic markers
#' @param config 2 dimension vector which indicates the presence (1) or the absence (0) of a QTL
#' @param Alpha additive value of the allele, if NULL, genetic effects are drawn in a normal distrbution
#' @param Delta Dominance value, if NULL, genetic effects are drawn in a normal distrbutio
#' @param sigma2A additive genetic variance
#' @param sigma2D dominance genetic variance
#' @param VbA matrix of variance covariance of additive effect
#' @param VbD matrix of variance covariance of dominance effect
#'
#' @return a df Betas which contains 3 columns, B ; B.1 ; B.2, some of them are null according to the configuration

CalcBeta = function(seg, marker, X, config, Alpha = NULL, Delta = NULL, 
                    sigma2A, sigma2D, VbA, VbD){
  requireNamespace("MASS")
  is.qtl = config == 1

  randomB = FALSE
  if(all(is.null(Alpha), is.null(Delta)))
    randomB = TRUE
  
  p = dim(X)[2]# nb marker per individual
  n = dim(X)[1]# nb observation per modality
  name.mk = colnames(X)
  Betas = data.frame(B=rep(0,p), B.1=rep(0,p), B.2=rep(0,p))
  row.names(Betas) = colnames(X)
  
  ## ABxCD
  if(seg =="<abxcd>"){
    mk1 = paste0(marker,".a") ;idx1 = which(mk1 == rownames(Betas))
    mk2 = paste0(marker, ".b") ;idx2 = which(mk2 == rownames(Betas))
    mk3 = paste0(marker, ".c") ;idx3 = which(mk3 == rownames(Betas))
    mk4 = paste0(marker, ".d") ;idx4 = which(mk4 == rownames(Betas))
    mk5 = paste0(marker, ".ac") ;idx5 = which(mk5 == rownames(Betas))
    mk6 = paste0(marker,".ad") ;idx6 = which(mk6 == rownames(Betas))
    mk7 = paste0(marker, ".bc") ;idx7 = which(mk7 == rownames(Betas))
    mk8 = paste0(marker, ".bd") ;idx8 = which(mk8 == rownames(Betas))
    if (randomB == FALSE){
      Betas$B[idx1] = Alpha$A
      Betas$B[idx2] = Alpha$B
      Betas$B[idx3] = Alpha$C
      Betas$B[idx4] = Alpha$D
      Betas$B[idx5] = Alpha$A + Alpha$C + Delta$AC
      Betas$B[idx6] = Alpha$A + Alpha$D + Delta$AD
      Betas$B[idx7] = Alpha$B + Alpha$C + Delta$BC
      Betas$B[idx8] = Alpha$B + Alpha$D + Delta$BD
    }else if(randomB == TRUE){
      if (!all(config == 1)){
        draw = rnorm(n = 4, mean = 0, sd = sqrt(sigma2A[is.qtl]))
        j = 1
        for(i in idx1:idx4){ # additive statistical effect
          Betas$B[i] = draw[j]
          j = j + 1 
        }
        k = 1
        draw = rnorm(n = 4, mean = 0, sd = sqrt(sigma2D[is.qtl]))
        for(i in idx5:idx8){ # dominance statistical effect
          Betas$B[i] = draw[k]
          k = k + 1 
        }

      }else if(all(config == 1)){ # config = c(1,1)
        draw = MASS::mvrnorm(n = 4, mu = rep(0,2), Sigma = VbA) 
        j = 1
        for(i in idx1 : idx4){ # additive statistical effect 
          Betas$B.1[i] = draw[j,1] # value in the 1st condition
          Betas$B.2[i] = draw[j,2] # value in the 2nd condition
          j = j + 1
        }
        draw = MASS::mvrnorm(n = 4, mu = rep(0,2), Sigma = VbD)
        k = 1
        for(i in idx5:idx8){ # dominance statistical effect
          Betas$B.1[i] = draw[k,1]
          Betas$B.2[i] = draw[k,2]
          k = k + 1 
        }
      }
    }
    
    ## EFxEG
  }else if(seg =="<efxeg>"){
    mk1 = paste0(marker,".e") ;idx1 = which(mk1 == rownames(Betas))
    mk2 = paste0(marker, ".f") ;idx2 = which(mk2 == rownames(Betas))
    mk3 = paste0(marker, ".g") ;idx3 = which(mk3 == rownames(Betas))
    mk4 = paste0(marker, ".ee") ;idx4 = which(mk4 == rownames(Betas))
    mk5 = paste0(marker, ".ef") ;idx5 = which(mk5 == rownames(Betas))
    mk6 = paste0(marker,".eg") ;idx6 = which(mk6 == rownames(Betas))
    mk7 = paste0(marker, ".fg") ;idx7 = which(mk7 == rownames(Betas))
    if(randomB == FALSE){
      Betas$B[idx1] = Alpha$E
      Betas$B[idx2] = Alpha$F
      Betas$B[idx3] = Alpha$G
      Betas$B[idx4] = 2*Alpha$E
      Betas$B[idx5] = Alpha$E + Alpha$F + Delta$EF
      Betas$B[idx6] = Alpha$E + Alpha$G + Delta$EG
      Betas$B[idx7] = Alpha$F + Alpha$G + Delta$FG
    }else if(randomB == TRUE){
      if (!all(config == 1)){
        draw = rnorm(n = 3, mean = 0, sd = sqrt(sigma2A[is.qtl]))
        j = 1
        for(i in idx1:idx3){ # additive statistical effect
          Betas$B[i] = draw[j]
          j = j + 1 
        }
        k = 1
        draw = rnorm(n = 4, mean = 0, sd = sqrt(sigma2D[is.qtl]))
        for(i in idx4:idx7){ # dominance statistical effect
          Betas$B[i] = draw[k]
          k = k + 1 
        }
        
      }else if(all(config == 1)){
        draw = MASS::mvrnorm(n = 3, mu = rep(0,2), Sigma = VbA) 
        j = 1
        for(i in idx1 : idx3){
          Betas$B.1[i] = draw[j,1] # value in the 1st condition
          Betas$B.2[i] = draw[j,2] # value in the 2nd condition
          j = j + 1
        }
        draw = MASS::mvrnorm(n = 4, mu = rep(0,2), Sigma = VbD)
        k = 1
        for(i in idx4:idx7){ # dominance statistical effect
          Betas$B.1[i] = draw[k,1]
          Betas$B.2[i] = draw[k,2]
          k = k + 1 
        }
      }
    }

    
    ## HKxHK
  }else if(seg =="<hkxhk>"){
    mk1 = paste0(marker,".h") ;idx1 = which(mk1 == rownames(Betas))
    mk2 = paste0(marker, ".k") ;idx2 = which(mk2 == rownames(Betas))
    mk3 = paste0(marker, ".hh") ;idx3 = which(mk3 == rownames(Betas))
    mk4 = paste0(marker, ".hk") ;idx4 = which(mk4 == rownames(Betas))
    mk5 = paste0(marker, ".kk") ;idx5 = which(mk5 == rownames(Betas))
    if(randomB == FALSE){
      Betas$B[idx1] = Alpha$H
      Betas$B[idx2] = Alpha$K
      Betas$B[idx3] = 2*Alpha$H
      Betas$B[idx4] = Alpha$H + Alpha$K + Delta$HK
      Betas$B[idx5] = 2*Alpha$K
      
    } else if(randomB == TRUE){
      if (!all(config == 1)){
        draw = rnorm(n = 2, mean = 0, sd = sqrt(sigma2A[is.qtl]))
        j = 1
        for(i in idx1:idx2){ # additive statistical effect
          Betas$B[i] = draw[j]
          j = j + 1 
        }
        k = 1
        draw = rnorm(n = 3, mean = 0, sd = sqrt(sigma2D[is.qtl]))
        for(i in idx3:idx5){ # dominance statistical effect
          Betas$B[i] = draw[k]
          k = k + 1 
        }
        
      }else if(all(config == 1)){
        draw = MASS::mvrnorm(n = 2, mu = rep(0,2), Sigma = VbA) 
        j = 1
        for(i in idx1 : idx2){
          Betas$B.1[i] = draw[j,1] # value in the 1st condition
          Betas$B.2[i] = draw[j,2] # value in the 2nd condition
          j = j + 1
        }
        draw = MASS::mvrnorm(n = 3, mu = rep(0,2), Sigma = VbD)
        k = 1
        for(i in idx3:idx5){ # dominance statistical effect
          Betas$B.1[i] = draw[k,1]
          Betas$B.2[i] = draw[k,2]
          k = k + 1 
        }
      }
    }
    
    ## LMxLL
  }else if (seg =="<lmxll>"){
    mk1 = paste0(marker,".l") ;idx1 = which(mk1 == rownames(Betas))
    mk2 = paste0(marker, ".m") ;idx2 = which(mk2 == rownames(Betas))
    mk3 = paste0(marker, ".ll") ;idx3 = which(mk3 == rownames(Betas))
    mk4 = paste0(marker, ".lm") ;idx4 = which(mk4 == rownames(Betas))
    if(randomB == FALSE){
      Betas$B[idx1] = Alpha$L
      Betas$B[idx2] = Alpha$M
      Betas$B[idx3] = 2*Alpha$L
      Betas$B[idx4] = Alpha$L + Alpha$M + Delta$LM
    }else if(randomB == TRUE){
      if (!all(config == 1)){
        draw = rnorm(n = 2, mean = 0, sd = sqrt(sigma2A[is.qtl]))
        j = 1
        for(i in idx1:idx2){ # additive statistical effect
          Betas$B[i] = draw[j]
          j = j + 1 
        }
        k = 1
        draw = rnorm(n = 2, mean = 0, sd = sqrt(sigma2D[is.qtl]))
        for(i in idx3:idx4){ # dominance statistical effect
          Betas$B[i] = draw[k]
          k = k + 1 
        }
        
      }else if(all(config == 1)){
        draw = MASS::mvrnorm(n = 2, mu = rep(0,2), Sigma = VbA) 
        j = 1
        for(i in idx1 : idx2){
          Betas$B.1[i] = draw[j,1] # value in the 1st condition
          Betas$B.2[i] = draw[j,2] # value in the 2nd condition
          j = j + 1
        }
        draw = MASS::mvrnorm(n = 3, mu = rep(0,2), Sigma = VbD)
        k = 1
        for(i in idx3:idx4){ # dominance statistical effect
          Betas$B.1[i] = draw[k,1]
          Betas$B.2[i] = draw[k,2]
          k = k + 1 
        }
      }
    }
    
    ## NNxNP
  }else if (seg =="<nnxnp>"){
    mk1 = paste0(marker,".n") ;idx1 = which(mk1 == rownames(Betas))
    mk2 = paste0(marker, ".p") ;idx2 = which(mk2 == rownames(Betas))
    mk3 = paste0(marker, ".nn") ;idx3 = which(mk3 == rownames(Betas))
    mk4 = paste0(marker, ".np") ;idx4 = which(mk4 == rownames(Betas))
    if(randomB == FALSE){
      Betas$B[idx1] = Alpha$N
      Betas$B[idx2] = Alpha$P
      Betas$B[idx3] = 2*Alpha$N
      Betas$B[idx4] = Alpha$N + Alpha$P + Delta$NP
    }else if(randomB == TRUE){
      if (!all(config == 1)){
        draw = rnorm(n = 2, mean = 0, sd = sqrt(sigma2A[is.qtl]))
        j = 1
        for(i in idx1:idx2){ # additive statistical effect
          Betas$B[i] = draw[j]
          j = j + 1 
        }
        k = 1
        draw = rnorm(n = 2, mean = 0, sd = sqrt(sigma2D[is.qtl]))
        for(i in idx3:idx4){ # dominance statistical effect
          Betas$B[i] = draw[k]
          k = k + 1 
        }
        
      }else if(all(config == 1)){
        draw = MASS::mvrnorm(n = 2, mu = rep(0,2), Sigma = VbA) 
        j = 1
        for(i in idx1 : idx2){
          Betas$B.1[i] = draw[j,1] # value in the 1st condition
          Betas$B.2[i] = draw[j,2] # value in the 2nd condition
          j = j + 1
        }
        draw = MASS::mvrnorm(n = 4, mu = rep(0,2), Sigma = VbD)
        k = 1
        for(i in idx3:idx4){ # dominance statistical effect
          Betas$B.1[i] = draw[k,1]
          Betas$B.2[i] = draw[k,2]
          k = k + 1 
        }
      }
    }
  }
 out = Betas
 return(out)
}



######################## SimulResponse  #############################################

#' This function is aiming to simulate genotypic values and response(s) 
#' in a segregating population where only 1 marker has an effect. 
#'In a multivariate case, the multiple linear regression besomes :
#'     Y = XB +  E 
#'Y is a matrix of NxT dimension with N observations and T responses
#'
#'X is the design matrix of markers (NxP) with P predictors, here markers
#'
#'B is the marker effects matrix. 2 choices: 
#'    - B is calculated from Alphas and Deltas from the previous
#'function 
#'    - B is drawn from a Matrix Normal distribution : B ~ MN(0,UbVb)
#'       We assume that Ub = Id (no covariance between markers)
#'       Vb is a TxT variance covariance matrix between responses
#'       VbA is the matrix with additive genetic effect based on the statistical model and
#'       VbD is the matrix with dominance statistical effect
#'       VbA =  | sigma²1A             rhoB.sigma1A.sigma2A |
#'             | rhoB.sigma1A.sigma2A   sigma²2A           |
#'             rhoB is the covariance between the responses due to genetics
#'             
#'E is the matrix of residual error (NxT)
#'    E ~ MN(0,UeVe)
#'    - we assume Ue = Id : no covariance between individuals due to error because we have 
#'    a good experimental design
#'    - Ve =  | sigma²1e             rhoE.sigma1e.sigma2e |
#'            | rhoE.sigma1e.sigma2e   sigma²2 e          |
#'            With rhoE the correlation between responses not due to genetic variation
#' Sigma²1e and Sigma²2e are calculated from Xie & Xu formula with A1, A2 and D
#'
#'
#' @param marker which marker has an impact on the phenotype
#' @param X design matrix additively coded in 0,1,2; the number of individuals
#' genotypes should be n. Rownames(X) should correspond to the marker names.
#' @param seg is the segregation observed at the marker, written in this format : "<mmxmm>" with m letters
#' according to the updated Joinmap format, for example for a 4 alleles segregation : "<abxcd>"
#' @param config a 2 dimension vector, if c(1), assume that there is one condition with a QTL,
#' if config = c(1,0), there is a QTL in the condition 1 but not 2, inverse for c(0,1),
#' c(0,0) = no QTL and c(1,1) = QTL in both conditions.
#' The dimension of config indicates the dimension of the response : 1 or 2 columns
#' @param h2 heritability (0<h2<1), amount of noise to add to simulate phenotypes based on real genetic effect: if h2 = 1, there is no noise.
#' Dimension of h2 = dimension of the response (and of config). If there is no QTL, h2 = 0
#' @param mean.response of the response; same dimension as response (1 or 2 columns). By default, mean.response = 0
#' @param rhoE correlation between the 2 conditions (not due to genetic variation) in the variance of the error term
#' @param Alpha half of the additive effect of each parent, calculated in MeanResponseOnGenoClass function
#' If Alpha = NULL, the beta's are drawn from an univariate or multivariate normal distribution with a mean of 0
#' and a standard deviation of sigma2B
#' @param Delta dominance effect, calculated in MeanResponseOnGenoClass function.
#' If Alpha = NULL, Delta must be NULL too
#' @param sigma2A additive genetic variance (in the statistical model) in the condition where there is a QTL, same dimension as response. 
#' If there is no QTL, sigma2A must be equal to 0
#' @param sigma2D dominance genetic variance (in the statistical model), must be equal to 0 if there is no QTL, 
#' can be 0 in the response where there  is a QTL
#' @param rhoB correlation between the 2 conditions in the variance of the B term, only used if config = c(1,1)
#' @return list of 3 elements : 
#' - y which is the response matrix (phenotype), with 1 or 2 columns according to config
#' - B which is the marker effect, with 1 or 2 columns according to config
#' - X the genotypic matrix 
#' @param sigma2.E is the error variance. In config [11] or [1], the error variance is defined by h2 and sigma2G only. If 
#' config is [01] or [10], one error variance is not defined, we have to add it separately, it is sigma2.E.
SimulResponse = function(marker, X, seg, config=c(0), h2=c(0), mean.response=c(0),
                         rhoE=0,
                         Alpha=NULL, Delta=NULL, sigma2A=c(0.1), sigma2D=c(0.1), rhoB=0, sigma2.E=0){
  requireNamespace("MASS")
  stopifnot(is.character(marker), # marker in right class
            ! is.null(colnames(X)), # marker names in colnames X present
            ! is.null(rownames(X)), # individual names in colnames X
            all(config %in% c(0,1)), # format of config
            all(h2 >= 0), 
            all(h2 <= 1),
            length(h2) == length(config),
            length(mean.response) == length(config),
            rhoE >= -1,
            rhoE <= 1,
            length(sigma2.E) == 1,
            all(sigma2A >= 0), # variances are positive
            all(sigma2D >= 0),
            xor(all(is.null(Alpha), is.null(Delta)), # if Alpha is NULL, Delta must be NULL too
                all(! is.null(Alpha), ! is.null(Delta))))
  if (all(all(config == 1) & length(config) == 2)){
    stopifnot(rhoB >= -1,
    rhoB <= 1)
  }
  
  possible.seg = c("<abxcd>","<efxeg>","<hkxhk>",'<lmxll>','<nnxnp>')
  stopifnot(is.seg = seg %in% possible.seg) # segregation corresponds to the updated Joinmap format
  markers = 0
  for(i in 1:ncol(X)){
    markers[i] = unlist(strsplit(colnames(X)[i], ".", fixed = TRUE))[1]
  }
  stopifnot(marker %in% markers)
  
  randomB = FALSE
  if(all(is.null(Alpha), is.null(Delta))){
    randomB = TRUE
    stopifnot(length(sigma2A) == length(config))
    stopifnot(length(sigma2D) == length(config))
  }
  
  
  is.qtl = logical(0)
  is.qtl[1] = config[1] == 1 ; is.qtl[2] = config[2] == 1
  if (length(config) == 1){
    stopifnot(h2[is.qtl][1] != 0)
  } else {
    if(any(is.qtl)){
      stopifnot(all(h2[is.qtl] != 0),
                all(h2[!is.qtl] == 0))
      
      if(randomB == TRUE){
        stopifnot(all(sigma2A[is.qtl]!=0), # must be !=0 if there is a QTL whereas sigma2D can be 0
                  all(sigma2A[!is.qtl] == 0),
                  all(sigma2D[!is.qtl] == 0))
      }
    } else {
      stopifnot(all(h2 ==0))
    }
  }

  p = dim(X)[2]# nb marker per individual
  n = dim(X)[1]# nb observation per modality
  
  name.mk = colnames(X)
  B = matrix(0,ncol = 1, nrow=p)
  row.names(B) = name.mk
  if(all(all(config == 1) & length(config) == 2)){
  # Vb matrix for additive effects
    VbA = as.matrix(cbind(c(sigma2A[1], rhoB*sqrt(sigma2A[1]*sigma2A[2])),
                       c(rhoB*sqrt(sigma2A[1]*sigma2A[2]), sigma2A[2])))
  # Vb matrix for dominance effects
    VbD = as.matrix(cbind(c(sigma2D[1], rhoB*sqrt(sigma2D[1]*sigma2D[2])),
                        c(rhoB*sqrt(sigma2D[1]*sigma2D[2]), sigma2D[2])))
  }else {
    VbA = 0
    VbD = 0
  }
  # Use the sub function CalcBeta to calculate the Betas according to the segregation and if randomB = TRUE or FALSE
  Betas = CalcBeta(seg = seg, marker=marker, X=X, config=config, Alpha=Alpha, Delta=Delta, 
                   VbA=VbA, VbD=VbD, sigma2A=sigma2A, sigma2D=sigma2D)
  Betas = as.data.frame(Betas)
  B0 = matrix(0, nrow=p)# cond with no impact on phenotype
  
 # sigma2G is used to calculate sigma2E with the heritatability value
    ## if there is no QTL, sigma2G = 0 <=> there is no genetic variance
  sigma2G = c(0) ; if(all(config == 0)){sigma2G = c(0,0)}
  sigma2E=c(0)
  if (randomB == FALSE){
    sigma2A = A1^2 + A2^2
    sigma2D = (1/16) * D^2
    sigma2G[is.qtl] = sigma2A + sigma2D
  }else{
    sigma2G = sigma2A + sigma2D
  }
  
  # Different configuration according to the configuration (QTL on which condition)
  ##### 1 response
  if (all(length(config) == 1 & config == 1)){ 
    # Variance of errors according to the heritability formula
    sigma2E = (sigma2G *(1-h2))/h2 # for 1 condition
    # Error simulation for 1 response
    E = as.matrix(rnorm(n = n, mean = mean.response, sd = (sqrt(sigma2E[1]))))
    B = Betas$B
    y = X %*% B + E
    
    ##### QTL in both conditions
  } else if (all(config == 1)){ 
    
    if(randomB == FALSE){
      B = cbind(Betas$B,Betas$B) # no covariance between conditions
      
    } else if(randomB == TRUE){
      
      if (all(config == 1) & length(config) == 2){
        B = cbind(Betas$B.1,Betas$B.2) # with covariances between conditions
      ## For 2 conditions, there is some covariance between condition and error
        sigma2E[1] = (sigma2G[1] *(1-h2[1]))/h2[1]
        sigma2E[2] = (sigma2G[2] *(1-h2[2]))/h2[2]

    
      ## Variance covariance matrix of error
      Ve = as.matrix(cbind(c(sigma2E[1], rhoE*sqrt(sigma2E[1]*sigma2E[2])),
                         c(rhoE*sqrt(sigma2E[1]*sigma2E[2]), sigma2E[2]))) 
    
      colnames(B) = c("Condition 1","Condition 2")
      E = MASS::mvrnorm(n = n, mu = mean.response, Sigma = Ve) # for config !=NA
      y = X %*% B + E # y has 2 columns, 1 for each condition
      }
    }
    
 ###### No QTL   
  } else if (all(config == 0)){
    ## For 2 conditions, there is some covariance between condition and error
    sigma2E[1] = 0
    sigma2E[2] = 0
    
    ## Variance covariance matrix of error
    Ve = as.matrix(cbind(c(sigma2E[1], rhoE*sqrt(sigma2E[1]*sigma2E[2])),
                         c(rhoE*sqrt(sigma2E[1]*sigma2E[2]), sigma2E[2])))  
    
    B = cbind(B0,B0)
    colnames(B) = c("Condition 1","Condition 2")
    E = MASS::mvrnorm(n = n, mu = mean.response, Sigma = Ve ) # for config !=NA
    y = X %*% B + E # y has 2 columns, 1 for each condition

    ###### QTL in 1 condition on the 2
  } else if (all(config[1] == 0 & config[2] == 1) |
             all(config[1] == 1 & config[2] == 0)) { # QTL in 1 condition 
    ## For 2 conditions, there is some covariance between condition and error
    sigma2E[is.qtl] = (sigma2G[is.qtl]*(1-h2[is.qtl]))/h2[is.qtl]
    # the response where there is no QTL, sigma2G and h2 are not defined, so we have to define sigma2E directly
    sigma2E[!is.qtl] = sigma2.E

    
    ## Variance covariance matrix of error
    Ve = as.matrix(cbind(c(sigma2E[1], rhoE*sqrt(sigma2E[1]*sigma2E[2])),
                         c(rhoE*sqrt(sigma2E[1]*sigma2E[2]), sigma2E[2]))) 
    
    if (is.qtl[1]){ # QTL in the 1st condition, config 10
      B.1 = Betas$B
      B.2 = B0
    }else { # QTL in the 2nd condition, config 01
      B.1 = B0
      B.2 = Betas$B
    }
    B = cbind(B.1,B.2)
    colnames(B) = c("Condition 1","Condition 2")
    E = MASS::mvrnorm(n = n, mu = mean.response, Sigma = Ve ) # for config !=NA
    y = X %*% B + E # y has 2 columns, 1 for each condition
    
  } 
  B <- as.matrix(B)
  rownames(B) = colnames(X)
  out = list(y=y, X=X, B=B)
  return(out)
}



############################## QTL detection ###################



#' univQTL
#'
#' @param dat.doligez2013 a 4way cross class object
#' @param response a data frame with 1 column which contains the response to analyze, rownames must
#' be individual names
#' @param nperm number of permutation to apply to find the LOD threshold, by default 100 
#' @param qtl.dist confidence interval to apply on both sides of the marker
#' @param method by default "em"
#' @param plot logical, to plot or not the result of the scanone function
#' @param ylim limit of the plot
#'
#' @return logical, TRUE if the QTL was above the threshold, FALSE else
#'
univQTL <- function (dat.doligez2013, response, nperm=100,qtl.dist=5,
                     method="em", plot=TRUE, ylim=c(0,15), main){
  
  out <- data.frame(0)
  out$chromosome <- NA
  out$position <- NA
  out$LOD <- NA
  requireNamespace("qtl")
  stopifnot(!is.null(rownames(response)))
  dat.doligez2013$pheno[["resp"]] <- NA
  dat.doligez2013 <- calc.genoprob(dat.doligez2013, step=1, map.function="kosambi")

  for(i in 1:nrow(dat.doligez2013$pheno)){
    ind <- as.character(dat.doligez2013$pheno$indiv[i])
    if(ind %in% rownames(response)){
      dat.doligez2013$pheno[["resp"]][i] <- response[ind,]
    }
  }
  qtl.em.perm <- qtl::scanone(dat.doligez2013, pheno.col="resp", method=method, n.perm = nperm)
  threshold.LOD <- summary(qtl.em.perm)[1] # find the correct threshold
  qtl.em <- qtl::scanone(dat.doligez2013, pheno.col="resp", method=method)
  if(plot){
    plot(qtl.em, ylim = ylim, main=main)
      abline(h=threshold.LOD) 
  }

  chr <- substr(qtl.em[which(qtl.em$lod > threshold.LOD),"chr"], start=2, stop=3)
  pos <- qtl.em[which(qtl.em$lod > threshold.LOD),"pos"]
  LOD <- qtl.em[which(qtl.em$lod > threshold.LOD), "lod"]
  nb_marker <- length(chr)
  out2 <- data.frame(chromosome=0, position=0, LOD=0)

  if(nb_marker != 0){
    out <- data.frame(chromosome=rep(0,nb_marker), position=rep(0,nb_marker), LOD=rep(0,nb_marker))
    for(i in 1:nb_marker){
      out$chromosome[i] <- chr[i]
      out$position[i] <- pos[i]
      out$LOD[i] <- LOD[i]
    }
    # For each LOD score superior than threshold, a new line is added in df out, we only keep
    # 1 marker per chromosome with the highest LOD score contained in out2 data frame
    out$chromosome <- as.factor(out$chromosome)
    nb_chr <- nlevels(out$chromosome)
    for(j in 1: nb_chr){
      chr <- levels(out$chromosome)[j]
      out_chr <- subset(out, out$chromosome == chr)
      if(nb_chr != 0){
        max_lod <- suppressWarnings(max(out_chr$LOD))
        out_chr <- out_chr[which(out_chr$LOD == max_lod),]
        out2 <- rbind(out2, out_chr)
      }
    }
  }
  out2 <- out2[-which(out2$chromosome == 0),]
  out2
  return(out2)
  }


glmnet.multivar <- function(x, y, q, type = c("conservative", "anticonservative"), ...) {
  if (type == "conservative") 
    fit <- suppressWarnings(glmnet::glmnet(x, y, pmax = q, ...))
  if (type == "anticonservative") 
    fit <- glmnet::glmnet(x, y, dfmax = q - 1, ...)
  selected <- predict(fit, type = "nonzero")
  selected <- selected[[length(selected)]]
  ret <- logical(ncol(x))
  ret[selected] <- TRUE
  cf <- fit$beta[[1]]
  sequence <- as.matrix(cf != 0)
  return(list(selected = ret, path = sequence))
}




#' glmnetmultivarQTL
#'
#' @param X design matrix of genotypic markers
#' @param Y data frame of response
#' @param marker marker name (character) which has an effect
#' @param map genetic map, 2 columns V1 and V2, one for the marker and the other for the position
#' @param choose.alpha value for parameter alpha in elastic net
#' @param PFER per family error rate for stability selection
#' @param B number of subsampling to apply for stability selection
#' @param cutoff after which probability do we select predictors after stability selection
#' @param plot logical, plot or not of stability paths
#' @param nb.cores result (-1) of the function DetectCores
#'
#' @return

glmnetmultivarQTL <- function(X, Y, marker, map, choose.alpha=0.1, 
                      PFER=1, B=500, cutoff=0.85, plot=TRUE, nb.cores){
  
  requireNamespace("glmnet")
  requireNamespace("stabs")
  
  # if(choose.alpha == FALSE){
  #   seq.alpha <- seq(0, 1, len=10)
  #   nfold <- 5
  #   nrep  <- 10 
  #   min.cv <- do.call(rbind, pbmclapply(1:nrep, function(i) {
  #     foldid <- sample(rep(seq(nfold), length=nrow(X)))
  #      sapply(seq.alpha, function(Falpha) {
  #        min(cv.glmnet(x=X, y=Y, family="mgaussian", alpha=Falpha,
  #                      standardize=FALSE, foldid=foldid)$cvm)
  #      })}, mc.cores=nb.cores)
  #    )
  #    dplot <- data.frame(mean = colMeans(min.cv),
  #                        sd   = sd(min.cv) / sqrt(nrep),
  #                        alpha = seq.alpha)
  #    ggplot(dplot, aes(x=alpha, y=mean, ymin=mean-sd, ymax=mean+sd)) +
  #      geom_smooth(stat="identity") + theme_bw() +
  #      geom_vline(xintercept=0.1)
  #    # alpha = min(dplot$mean),alpha
  # }
  alpha=choose.alpha 
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  q <- max(5, min(2 * round(nrow(X)/log(ncol(X))/10)*10, ncol(X)))
  system.time(
    stabsel.menet <- 
      stabsel(X, Y, q=q, PFER=1, B=500, 
                fitfun=glmnet.multivar,
                args.fitfun=list(alpha=alpha, family="mgaussian",
                                 standardize=FALSE),
                mc.cores=nb.cores, mc.preschedule=TRUE))
  if(plot){
    plot(stabsel.menet, type="path", main=paste("Stability paths","glmnet"))
      abline(h=cutoff, col="blue")
    plot(stabsel.menet, type="maxsel", main=paste("Maximum selection frequencies", "glmnet"))
      abline(v=cutoff, col="blue")
  }
  sel <- as.matrix(sort(stabsel(stabsel.menet, cutoff= cutoff)$selected, decreasing = TRUE))
  sel.mrks.glmnet <- data.frame(marker=0, allele=0,
                                linkage_group=0, position=0, probability=0, idx=0) 
  if(length(sel)!=0){ # at least 1 marker has been detected
    for(i in 1:nrow(sel)){
      if (is.na(unlist(strsplit(rownames(sel)[i], ".", fixed = TRUE))[3])){#no composed marker name
        sel.mrks.glmnet[i,2] <- unlist(strsplit(rownames(sel)[i], ".", fixed = TRUE))[2] # allele
        sel.mrks.glmnet[i,1] <- unlist(strsplit(rownames(sel)[i], ".", fixed = TRUE))[1] # marker 
      }else{
        sel.mrks.glmnet[i,2] <- unlist(strsplit(rownames(sel)[i], ".", fixed = TRUE))[3] # allele
        sel.mrks.glmnet[i,1] <- paste0(unlist(strsplit(rownames(sel)[i], ".", fixed = TRUE))[1],"-", # marker part1
                                       unlist(strsplit(rownames(sel)[i], ".", fixed = TRUE))[2])# marker part2
      }
      sel.mrks.glmnet[i,3] <- substr(mrk2lg[sel.mrks.glmnet[i,1]],start=2, stop=3) # linkage group
      if(!sel.mrks.glmnet[i,1] %in% map$V1){
        sel.mrks.glmnet[i,4] <- NA # too complicated marker name
        
      }else {
        sel.mrks.glmnet[i,4] <- as.numeric(map[map$V1 == sel.mrks.glmnet[i,1],2]) # position
      }
      sel.mrks.glmnet[i,5] <-  stabsel.menet$max[sel][i] # probability of selection
    }
  }
 sel.mrks.glmnet <- sel.mrks.glmnet[order(sel.mrks.glmnet$probability, decreasing=TRUE),]

 return(sel.mrks.glmnet)
}




#' fus2modQTL
#'
#' @param X design matrix for genotypic markers
#' @param Y data frame of response (multivariate)
#' @param marker marker which has an effect on phenotype
#' @param map genetic map to localize markers
#' @param group name of 2 group responses, by default it is "R1" and "R2"
#' @param nrep 2 dimension vector for number repetition for cross validation and stability selection
#' @param plot logical, plot or not of stability path
#' @param cutoff after which probability do we select predictors after stability selection
#' @param lambda value of lambda for lasso parameter, default is null
#' @param nb.cores result (-1) of the function DetectCores
#'
#' @return

fus2modQTL <- function(X, Y, marker, map, group=c("R1","R2"), nrep=c(5,100), 
                       plot=c(FALSE,TRUE),cutoff=0.85,lambda=NULL,nb.cores){
  requireNamespace("Fus2mod")
  requireNamespace("stabs")
  nb_ind <- dim(Y)[1]
  y = as.matrix(c(Y[,1],Y[,2]))
  group <- c(rep(group[1],nb_ind),rep(group[2],nb_ind))
  regressors <- rbind(X,X) # in this function, we need to duplicate the design matrix of markers
  X <- model.matrix(~group + group:regressors -1)
  system.time(
    CV <- Fus2mod::cv.fl2(response = y, regressors = regressors, group = group, mina = 0.1, nfold = 5, nrep=nrep[1],
                          nb.cores = nb.cores,plot = plot[1]))
  ddl <- CV[CV$Criterion =="ddl",]
  a <- ddl[which.min(ddl$mean), "a"]
  
  stab <- Fus2mod::stab.fl2_fixa(response = y, regressors, group, a,
                                 lambda = lambda,  nrep= nrep[2],
                                 nb.cores = nb.cores, plot = plot[2])
  plot(stab, type="maxsel", main=paste("Max sel frequencies","Fus2mod"," - cutoff =", cutoff))
  abline(v=cutoff, col="blue")
  plot(stab, type="path", main=paste("Stability paths", "Fus2mod", " - cutoff =", cutoff))
  abline(h=cutoff, col="blue")
  
  ## Extract markers with a selection probability superior to the cutoff ##
  
  sel <- names(stabsel(stab, cutoff=cutoff)$selected)[-c(1,2)]
  sel.mrks.fus <- data.frame(marker=0, allele=0, response=0, 
                             linkage_group=0, position=0, probability=0, idx=0)
  if (length(sel) !=0){
    for(i in 1:length(sel)){
      ### predictor has been detected with an effect on 1 condition
      if(substr(sel[i], start=1, stop=5) == "group"){ 
        mark.all <- unlist(strsplit(sel[i], "regressors"))[2] # marker with allele
        ##### no composed marker name
        if (is.na(unlist(strsplit((mark.all)[i], ".", fixed=TRUE))[3])){
          sel.mrks.fus[i,1] <- unlist(strsplit((mark.all), ".", fixed=TRUE))[1] # marker
          sel.mrks.fus[i,2] <- unlist(strsplit((mark.all), ".", fixed=TRUE))[2] # allele
          ##### composed marker name
        }else{
          sel.mrks.fus[i,2] <- unlist(strsplit((mark.all)[i], ".", fixed=TRUE))[3] # allele
          sel.mrks.fus[i,1] <- paste0(unlist(strsplit((mark.all)[i], ".", fixed=TRUE))[1],# marker part 1
                                      ".",
                                      unlist(strsplit((mark.all)[i], ".", fixed=TRUE))[2])# marker part 2
        }
        sel.mrks.fus[i,3] <- substr(sel[i], start=6, stop=7) # extract group response
        
        ### predictor has been detected with an effect on both conditions
      }else{
        #### no composed marker name
        if (is.na(unlist(strsplit((sel)[i], ".", fixed=TRUE))[3])){
          sel.mrks.fus[i,1] <- unlist(strsplit((sel)[i], ".", fixed=TRUE))[1] # marker
          sel.mrks.fus[i,2] <- unlist(strsplit((sel)[i], ".", fixed=TRUE))[2] # allele
          #### composed marker name
        }else{
          sel.mrks.fus[i,2] <- unlist(strsplit((sel)[i], ".", fixed=TRUE))[3] # allele
          sel.mrks.fus[i,1] <- paste0(unlist(strsplit((sel)[i], ".", fixed=TRUE))[1],# marker part 1
                                      "-",
                                      unlist(strsplit((sel)[i], ".", fixed=TRUE))[2])# marker part 2
        }
        sel.mrks.fus[i,3] <- "both"
      }
      sel.mrks.fus[i,4] <- substr(mrk2lg[sel.mrks.fus[i,1]],start=2, stop=3) # linkage group
      ### extract marker position
      if(sel.mrks.fus[i,1] %in% map$V1){
        sel.mrks.fus[i,5] <- as.numeric(map[map$V1 == sel.mrks.fus[i,1],2]) # position
      }else {
        sel.mrks.fus[i,5] <- NA # too complicated marker name
      }
      sel.mrks.fus[i,6] <- stab$max[sel][i] # probablity of selection
    }
  }

  
  
  
  ## Now select the 5 first markers with a selection probability of 0.55 ##
  
  sel.first <- names(stabsel(stab, cutoff=0.55)$selected)[-c(1,2)][1:5]
  sel.first.mrks.fus <- data.frame(marker=0, allele=0, response=0, 
                             linkage_group=0, position=0, probability=0, idx=0)
  if (all(!is.na(sel.first))& length(sel.first )!=0){
    for(i in 1:length(sel.first)){
      ### predictor has been detected with an effect on 1 condition
      if(substr(sel.first[i], start=1, stop=5) == "group"){
        mark.all <- unlist(strsplit(sel.first[i], "regressors"))[2] # marker with allele
        #### no composed marker name
        if (is.na(unlist(strsplit((mark.all)[i], ".", fixed=TRUE))[3])){
          sel.first.mrks.fus[i,1] <- unlist(strsplit((mark.all), ".", fixed=TRUE))[1] # marker
          sel.first.mrks.fus[i,2] <- unlist(strsplit((mark.all), ".", fixed=TRUE))[2] # allele
          #### composed marker name
        }else{
          sel.first.mrks.fus[i,2] <- unlist(strsplit((mark.all)[i], ".", fixed=TRUE))[3] # allele
          sel.first.mrks.fus[i,1] <- paste0(unlist(strsplit((mark.all)[i], ".", fixed=TRUE))[1],# marker part 1
                                      ".",
                                      unlist(strsplit((mark.all)[i], ".", fixed=TRUE))[2])# marker part 2
        }
        sel.first.mrks.fus[i,3] <- substr(sel.first[i], start=6, stop=7) # extract group response
        
        ### predictor has been detected with an effect on both conditions 
      }else { 
        #### no composed marker name
        if (is.na(unlist(strsplit((sel.first)[i], ".", fixed=TRUE))[3])){
          sel.first.mrks.fus[i,1] <- unlist(strsplit((sel.first)[i], ".", fixed=TRUE))[1] # marker
          sel.first.mrks.fus[i,2] <- unlist(strsplit((sel.first)[i], ".", fixed=TRUE))[2] # allele
          #### composed marker name
        }else{
          sel.first.mrks.fus[i,2] <- unlist(strsplit((sel.first)[i], ".", fixed=TRUE))[3] # allele
          sel.first.mrks.fus[i,1] <- paste0(unlist(strsplit((sel.first)[i], ".", fixed=TRUE))[1],# marker part 1
                                      "-",
                                      unlist(strsplit((sel.first)[i], ".", fixed=TRUE))[2])# marker part 2
        }
        sel.first.mrks.fus[i,3] <- "both"
      }
      sel.first.mrks.fus[i,4] <- substr(mrk2lg[sel.first.mrks.fus[i,1]],start=2, stop=3) # linkage group
      ### extract marker position
      if(sel.first.mrks.fus[i,1] %in% map$V1){
        sel.first.mrks.fus[i,5] <- as.numeric(map[map$V1 == sel.first.mrks.fus[i,1],2]) # position
      }else {
        sel.first.mrks.fus[i,5] <- NA # too complicated marker name
      }
      sel.first.mrks.fus[i,6] <- stab$max[sel.first][i] # probablity of selection
    } # end for
  } # end if
  
  ## Order the results by selection probability
  sel.first.mrks.fus <- sel.first.mrks.fus[order(sel.first.mrks.fus$probability, decreasing=TRUE),]
  sel.mrks.fus <- sel.mrks.fus[order(sel.mrks.fus$probability, decreasing=TRUE),]
  
  ## Out
  out = list(sel.mrks.fus, sel.first.mrks.fus)
  return(out)
}

############# Convert stat to functionnal effect ####################
# Create a function which convert statistical genetic effects (Ealpha) to functionnal statistical effects (Ef)

# For that, according to Yang et al 2008 Functionnal and statistical genetic effects with multiple alleles,
# we computeTalpha such that Ealpha = Talpha x Ef(0)

S2F = function(N=cbind(c(2,1,0,1,0,0,1,0,0,0),c(0,1,2,0,1,0,0,1,0,0),c(0,0,0,1,1,2,0,0,1,0),c(0,0,0,0,0,0,1,1,1,2)), 
               p=c(p11, p12, p22, p13, p23, p33, p14, p24, p34, p44), 
               alpha=c(G11, G12, G22, G13, G23, G33, G14, G24, G34, G44) ){
  G = solve(solve(t(N)*p*diag(10)*N)*t(N)*p*diag(10)*(diag(10)-t(p)))*alpha
  
  
  alphaG = solve(t(N)*p*diag(10)*N)*t(N)*p*diag(10)*(diag(10)-t(p))*G
  
}




