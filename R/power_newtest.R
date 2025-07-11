#' Power for tests with p values
#' 
#' This function estimates the power of test routines that calculate p value(s)
#'
#' @param TS routine to calculate test statistics.
#' @param rxy routine that generates data.
#' @param xparam values of first parameter under the alternative hypothesis. 
#' @param yparam values of second parameter under the alternative hypothesis.  
#' @param TSextra list passed to TS.
#' @param alpha =0.05  type I error.
#' @param B = 1000 number of simulation runs to estimate the power.
#' @return A matrix of power values

power_newtest <- function(TS, f, xparam, yparam, TSextra, alpha=0.05, B=1000) {
     Continuous=ifelse("vals" %in% names(formals(TS)), FALSE, TRUE) 
     if(missing(yparam)) dta=f(xparam[1])
     else dta=f(xparam[1], yparam[1])
     if(Continuous) {
        if(length(formals(TS))==2) tmp=TS(dta$x, dta$y)
        else tmp=TS(dta$x, dta$y, TSextra)
     }
     else {
       if(length(formals(TS))==3) tmp=TS(dta$x, dta$y, dta$vals)
        else tmp=TS(dta$x, dta$y, dta$vals, TSextra)
     } 
     pwr=matrix(0, length(xparam), length(TS))
     for(j in seq_along(xparam)) {
        for(i in 1:B) {
          if(missing(yparam)) dta=f(xparam[j])
          else dta=f(xparam[j], yparam[j])
          if(Continuous) {
            if(length(formals(TS))==2) tmp=TS(dta$x, dta$y)
             else tmp=TS(dta$x, dta$y, TSextra)
          }
          else {
            if(length(formals(TS))==3) tmp=TS(dta$x, dta$y, dta$vals)
            else tmp=TS(dta$x, dta$y, dta$vals, TSextra)
          }
          pwr[j, ]=pwr[j, ]+ifelse(tmp<alpha,1,0) 
        }
     }
     pwr/B
}
