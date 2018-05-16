# cox cure  
source("epm.r")
source("LogLik.r")
highCoxCure = function(times, delta, z, fold=3,
            nlambda1 = 15, nlambda2 = 15, lambda.min = 0.01,
            lambda3 = c(0, 0.001, 0.01, 0.05, 0.1),	
            nu=4, xi=1/sqrt(n), eps=10^-4, maxitCD=1, maxitEM=10^3)
{
    n = nrow(z)
	p = ncol(z)
    L = length(lambda3)
    bet.bic = bet.gcv = matrix(0, p, L)
	gam.bic = gam.gcv = matrix(0, p+1, L)
	opt.bic = opt.gcv = matrix(0, 2, L)
	BIC = LogLik = numeric(L)
	index = sample(1:fold, n, replace=TRUE) 
	# make sure the samples with the largest death time is in the training set
	# or the observied full log-likelihood on testing set will be -inf.
	index[which(times == max(times[delta==1]))] = fold
	times_te = times[index==1]; delta_te=delta[index==1]; z_te=z[index==1,]
	times_tr = times[index!=1]; delta_tr=delta[index!=1]; z_tr=z[index!=1,]
	tic = proc.time()
    for (l in 1:L)
    {

        fit = epm(times_tr, delta_tr, z_tr, z_tr, nlambda1, nlambda2,
           		lambda.min, lambda3[l], 
				nu, xi, eps, maxitCD, maxitEM, maxDF=sum(delta_tr)/4)
        # bic
		dfs = matrix(0, nlambda1, nlambda2)
        for (j in 1:nlambda2) 
            for (i in 1:nlambda1) 
			    dfs[i, j] = sum(fit$beta[,i,j]!=0) +
				                sum(fit$gamma[,i,j]!=0)
		 
        loglik.train = matrix(0, nrow=nlambda1, nlambda2)
		for (j in 1:nlambda2) {
            for (i in 1:nlambda1) {
			   NAs = sum(is.na(c(fit$beta[,i,j],fit$gamma[,i,j])))
			   if ( all(fit$W[,i,j]==0) | (NAs != 0) ) loglik.train[i, j] = -Inf 
               else	 loglik.train[i, j] =  fullloglik(times_tr, delta_tr, z_tr, fit$beta[,i,j], 
										        fit$gamma[,i,j], fit$W[,i,j], 
			                                    times_tr, delta_tr, z_tr)  								
			}									
		}										
 
        #bic = -2*fit$loglik + log(n)*dfs
		bic = -2*loglik.train + log(n)*dfs
		
		# likelihood based on test data.
		loglik.test = matrix(0, nrow=nlambda1, nlambda2)
		
		for (j in 1:nlambda2) {
            for (i in 1:nlambda1) {
			   NAs = sum(is.na(c(fit$beta[,i,j],fit$gamma[,i,j])))
			   if ( all(fit$W[,i,j]==0)  | (NAs != 0) ) loglik.test[i, j] = -Inf  
			   else loglik.test[i, j] =  fullloglik(times_tr, delta_tr, z_tr, fit$beta[,i,j], 
										       fit$gamma[,i,j], fit$W[,i,j], 
			                                   times_te, delta_te, z_te)  
			}								   
		}
        BIC[l] = min(bic, na.rm = TRUE)		
        opt.bic[,l] = which(bic == min(bic, na.rm = TRUE), arr.ind=TRUE)[1,]
        bet.bic[,l] = fit$beta[,opt.bic[1,l], opt.bic[2,l]]
		gam.bic[,l] = fit$gamma[,opt.bic[1,l], opt.bic[2,l]]
		
		LogLik[l] = max(loglik.test, na.rm = TRUE)
		opt.gcv[,l] = which(loglik.test == max(loglik.test, na.rm = TRUE), arr.ind=TRUE)[1,]
        bet.gcv[,l] = fit$beta[,opt.gcv[1,l], opt.gcv[2,l]]
		gam.gcv[,l] = fit$gamma[,opt.gcv[1,l], opt.gcv[2,l]]
    }
	toc = proc.time()
	list(
	    bet.mcp.bic = bet.bic[,1],
	    gam.mcp.bic = gam.bic[,1],
	    bet.mcp.gcv = bet.gcv[,1],
	    gam.mcp.gcv = gam.gcv[,1],
	    bet.hat.bic  = bet.bic[,which.min(BIC)],
	    gam.hat.bic  = gam.bic[,which.min(BIC)],
	    bet.hat.gcv  = bet.gcv[,which.max(LogLik)],
	    gam.hat.gcv  = gam.gcv[,which.max(LogLik)],
		
		opt.bic      = opt.bic,     # local of optimal tuning for lambda1 and lambda2, under each Lambda3
	    opt.gcv      = opt.gcv,
	    BIC          = BIC,         # minimal BIC under each Lambda3
	    LogLik       = LogLik,      # maximal loglik under each Lambda3.
	    lambda1      = fit$lambda1,
	    lambda2      = fit$lambda2,
	    lambda3      = lambda3)
}
	