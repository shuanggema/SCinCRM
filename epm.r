# interface function for EPM
sourceCpp("trial.cpp")
source("standardize.r")
epm = function(times, delta, x, x1, nlambda1, nlambda2, lambda.min,
                lambda3, nu, xi, eps, maxitCD, maxitEM, maxDF = sum(delta)/4)
{
    std    = standardize(x)
    z      = std[[1]]
    center = std[[2]]
    scales = std[[3]]
	
	std1    = standardize(x1)
    z1      = std1[[1]]
    center1 = std1[[2]]
    scales1 = std1[[3]]
	
    z1 = cbind(1, z1)
    fit = EPM(times, delta, z, z1, nlambda1, nlambda2, lambda.min,
                lambda3, nu, xi, eps, maxitCD, maxitEM, maxDF)
	beta  = fit$beta/scales
	gamma = array(0, dim=c(ncol(z1), nlambda1, nlambda2))
    gamma[-1, , ] =	fit$gamma[-1,,]/scales1
	for (i in 1:nlambda1) 
	    for (j in 1:nlambda2)
	        gamma[1, i, j] =	fit$gamma[1, i, j]	- t(center1/scales1) %*% fit$gamma[-1, i, j]
			
			
	beta[abs(beta) <= 10^-2] = 0		
	gamma[abs(gamma) <= 10^-2] = 0
    list(gamma   = gamma,
	     beta    = beta, 
		 lambda1 = fit$lambda1,
		 lambda2 = fit$lambda2,
		 loglik1 = fit$loglik1,
         loglik2 = fit$loglik2, 	
         loglik  = fit$loglik,	
         W       = fit$W,
         S0      = fit$S0,
         eta     = fit$eta,
         pi      = fit$pi,
         converge= fit$converge,
         iterEM  = fit$it,
         S       = fit$S  		 
		 )

 }
 
