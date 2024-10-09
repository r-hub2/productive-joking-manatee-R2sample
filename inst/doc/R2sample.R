## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(R2sample)

## -----------------------------------------------------------------------------
set.seed(123)

## -----------------------------------------------------------------------------
x1 = rnorm(10)
y1 = rnorm(12)
twosample_test(x1, y1, B=500, maxProcessor = 2)

## -----------------------------------------------------------------------------
twosample_test_adjusted_pvalue(x1, y1, B=c(500,500))

## -----------------------------------------------------------------------------
myTS1 = function(x, y) {
   out = c(0, 0)
   out[1] = abs(sd(x) - sd(y))
   out[2] = abs(mean(x)/sd(x) - mean(y)/sd(y))
   names(out) = c("std", "std t test")
   out
}

## -----------------------------------------------------------------------------
twosample_test(x1, y1, TS=myTS1, B=500, maxProcessor = 2)

## -----------------------------------------------------------------------------
x2 = table(c(0:5,rbinom(1000, 5, 0.5)))-1
y2 = table(c(0:5,rbinom(1200, 5, 0.55)))-1
rbind(x2, y2)
twosample_test(x2, y2, vals=0:5, B=500, maxProcessor = 2)$p.values

## -----------------------------------------------------------------------------
twosample_test(x2, y2, vals=0:5, TS=myTS2, B=500)

## -----------------------------------------------------------------------------
x=rnorm(10)
y=rnorm(12)
twosample_test(x, y, B=500, maxProcessor = 2, doMethods=c("KS","AD"))

## ----eval=FALSE---------------------------------------------------------------
#  run_shiny()

## ----eval=FALSE---------------------------------------------------------------
#  plot_power(twosample_power(
#    f=function(df) list(x=rnorm(100), y=rt(200, df)), df=1:10)
#  )

## ----eval=FALSE---------------------------------------------------------------
#  plot_power(twosample_power(
#    function(p) {
#      vals=0:10
#      x=table(c(vals, rbinom(100, 10, 0.5))) - 1 #Make sure each vector has length 11
#      y=table(c(vals, rbinom(120, 10, p))) - 1  #and names 1-10
#      vals=vals[x+y>0] # only vals with at least one count
#      list(x=x[as.character(vals)],y=y[as.character(vals)], vals=vals)
#    },
#    p=seq(0.5, 0.6, length=5)
#  ), "p")

## ----eval=FALSE---------------------------------------------------------------
#  plot_power(twosample_power(
#    f=function(n) list(x=rnorm(n), y=rnorm(n, 1)),
#    n=10*1:10), "n")

