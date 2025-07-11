#' Adjusted p values for simultaneous testing in the two-sample problem.
#'
#' This function runs a number of two sample tests using Rcpp and parallel computing and then finds the correct p value for the combined tests.
#'
#' For details consult vignette("R2sample","R2sample")
#' 
#' @param  x  a vector of numbers if data is continuous or of counts  if data is discrete, or a list with the data.
#' @param  y a vector of numbers if data is continuous or of counts  if data is discrete.
#' @param  vals =NA, a vector of numbers, the values of a discrete random variable. NA if data is continuous data.
#' @param  TS routine to calculate test statistics for non-chi-square tests
#' @param  TSextra additional info passed to TS, if necessary
#' @param  wx A numeric vector of weights of x.
#' @param  wy A numeric vector of weights of y.
#' @param  B =c(5000, 1000), number of simulation runs for permutation test
#' @param  nbins =c(50,10), number of bins for chi square tests.
#' @param  minexpcount = 5, minimum required expected counts for chi-square tests
#' @param  samplingmethod ="independence" or "Binomial" for discrete data
#' @param  rnull routine for parametric bootstrap
#' @param  SuppressMessages = FALSE print informative messages?
#' @param  doMethods ="all" a vector of codes for the methods to include. If "all", all methods are used.
#' @return A list of two numeric vectors, the test statistics and the p values. 
#' @examples
#'  x=rnorm(100)
#'  y=rt(200, 4)
#'  R2sample::twosample_test_adjusted_pvalue(x, y, B=c(500, 500))
#'  vals=1:5
#'  x=table(c(1:5, sample(1:5, size=100, replace=TRUE)))-1
#'  y=table(c(1:5, sample(1:5, size=100, replace=TRUE, prob=c(1,1,3,1,1))))-1
#'  R2sample::twosample_test_adjusted_pvalue(x, y, vals, B=c(500, 500))
#' @export
twosample_test_adjusted_pvalue=function(x, y, vals=NA, TS, TSextra, wx=rep(1, length(x)),
                        wy=rep(1, length(y)), B=c(5000, 1000), nbins=c(50,10),
                        minexpcount=5, samplingmethod="independence", 
                        rnull, SuppressMessages=FALSE, doMethods) {
    default.methods = list(cont=c("ES small", "ZA", "ZK", "Wassp1","Kuiper"), 
                           disc=c("small", "ZA", "Kuiper", "Wassp1"))
    all.methods = list(cont=c("KS","Kuiper","CvM","AD","LR","ZA","ZK","ZC","Wassp1",
                             "ES large","ES small","EP large","EP small"),
                       disc=c("KS","Kuiper","CvM","AD","LR","ZA","Wassp1","large","small"))                                          
    if(is.list(x)) {
      y=x$y
      if("vals"%in%names(x)) vals=x$vals
      if("wx"%in%names(x)) wx=x$wx
      if("wy"%in%names(x)) wy=x$wy
      x=x$x
    }
    samplingmethod=ifelse(samplingmethod=="independence", 1, 2)
    if(missing(TSextra)) TSextra = list(samplingmethod=samplingmethod)
    else TSextra=c(TSextra, samplingmethod=samplingmethod)
    if(!missing(rnull)) TSextra=c(TSextra, rnull=rnull)
    Continuous=ifelse(any(is.na(vals)), TRUE, FALSE)
    CustomTS=ifelse(missing(TS), FALSE, TRUE)
    WithWeights=ifelse(all(c(wx,wy)==1), FALSE, TRUE)
    if(length(B)==1) B=c(B, B)
    if(missing(doMethods)) {
      if(Continuous) doMethods=default.methods$cont
      else doMethods=default.methods$disc
    }
    if(doMethods[1]=="all") {
      if(Continuous) doMethods=all.methods$cont
      else doMethods=all.methods$disc
    }
    if(test_methods(doMethods,Continuous, FALSE, WithWeights))
      return(NULL)
# what methods are to be run?  
    if(!CustomTS) { # do included methods
      if(Continuous) {
        if(!WithWeights) { # No weights
          typeTS = 1
          TS=TS_cont
          dta=list(x=x, y=y)
          tmp=TS(x,y)
        }
        else { # weighted data
          typeTS = 2
          TS=TSw_cont
          dta=list(x=x, y=y, wx=wx, wy=wy)
          tmp=TS(x, y, wx, wy)
          if(min(wx,wy)<0.01 | max(wx,wy)>5)
            if(!SuppressMessages) message("Some of the weights are either exceptionally small or large. In either case the tests can be unreliable!")
        }
      }
      else {
        typeTS=5
        TS=TS_disc
        dta=list(x=x, y=y, vals=vals)
        TSextra=c(TSextra, list(adw=weights(dta)))
        tmp=TS(x, y, vals, TSextra$adw)
      }  
    }  
    else { # do user-supplied tests
      if(substr(deparse(TS)[2], 1, 5)==".Call") {
        if(maxProcessor>1) {
          if(!SuppressMessages) message("Parallel Programming is not possible if custom TS is written in C++. Switching to single processor")  
          maxProcessor=1
        }  
      }
      if(Continuous) {
        if(length(formals(TS))==2) {
          typeTS=3
          dta=list(x=x, y=y)
          tmp=TS(x, y)
        }
        else {
          typeTS=4
          dta=list(x=x, y=y)
          tmp=TS(x, y, TSextra)
        }
      }
      else {
        if(length(formals(TS))==3) {
          typeTS=6
          dta=list(x=x, y=y, vals=vals)
          tmp=TS(x, y, vals)
        }  
        else {
          typeTS=7
          dta=list(x=x, y=y, vals=vals)
          tmp=TS(x, y, vals, TSextra)
        }
      }
      if(is.null(names(tmp))) {
        message("output of TS routine has to be a named vector!")
        return(NULL)
      }  
    }
    
#   Find matrix of values of test statistics for simulated data
    num_tests=length(tmp)
    A=matrix(0, B[1], num_tests)
    for(i in 1:B[1]) {
       simdta=gen_sim_data(dta, TSextra)
       A[i, ]=calcTS(simdta, TS, typeTS, TSextra)
    }  
    pvalsTS=matrix(0, B[2]+1, num_tests)
    colnames(pvalsTS)=names(tmp)
    pvalsCHI=matrix(0, B[2]+1, ifelse(Continuous, 4, 2))
    if(Continuous) colnames(pvalsCHI)=c("ES large","ES small","EP large","EP small")
    else colnames(pvalsCHI)=c("large","small")
    for(i in 1:(B[2]+1)) {
      if(i==1) {
        tmp=calcTS(dta, TS, typeTS, TSextra)
        if(!CustomTS) 
          pvalsCHI[i, ]=chi_test(dta, nbins=nbins, minexpcount=minexpcount, typeTS=typeTS)$p.values
      }  
      else {
        simdta=gen_sim_data(dta, TSextra)
        tmp=calcTS(simdta, TS, typeTS, TSextra)
        if(!CustomTS) 
          pvalsCHI[i, ]=chi_test(simdta, nbins=nbins, minexpcount=minexpcount, typeTS=typeTS)$p.values
      }  
      for(j in 1:num_tests) pvalsTS[i,j]=pvalsTS[i,j]+sum(tmp[j]>A[,j])/B[1]    
    }
    if(CustomTS)  pvalsCHI=NULL
    pvals=cbind(pvalsTS, pvalsCHI) 
    pvals=pvals[,doMethods,drop=FALSE]
    minp_x=min(pvals[1, ])
    minp_sim=apply(pvals[-1, ,drop=FALSE], 1, min)
    z=seq(0, 1, length=250)
    y=z
    for(i in 1:250) y[i]=sum(minp_sim<=z[i])/length(minp_sim)
    I=c(1:250)[z>minp_x][1]-1
    slope=(y[I+1]-y[I])/(z[I+1]-z[I])
    minp_adj=round(y[I]+slope*(minp_x-z[I]),4)
    message("p values of individual tests:")
    for(i in 1:ncol(pvals)) message(paste(colnames(pvals)[i],": ", pvals[1,i]))
    message(paste0("adjusted p value of combined tests: ", minp_adj))
    
}
