library(smcure)
data(e1684)
library(Rcpp)
source("highCoxCure.r")
times <- e1684$FAILTIME
delta <- e1684$FAILCENS
z <- cbind(e1684$TRT, e1684$SEX, e1684$AGE)
z[is.na(z)] = mean(z, na.rm=TRUE)
n <- length(times)
fit <- highCoxCure(times, delta, z, fold=2,
            nlambda1 = 15, nlambda2 = 15, lambda.min = 0.01,
            lambda3 = c(0, 0.001, 0.01, 0.05, 0.1),	
            nu=4, xi=1/sqrt(n), eps=10^-4, maxitCD=1, maxitEM=10^3)
fit$bet.hat.gcv
fit$gam.hat.gcv