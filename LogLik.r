# observed full likelihood function


fullloglik = function(times, status, z, bet, gam, w, tnew, statusnew, znew)  {
    expeta = exp(expsafe(gam[1] + c(znew %*% gam[-1])))
    pi = expeta/(1 + expeta)
			
    S0 = survbaseline(times, status, z, bet, w)
    surv0 = stepfun(sort(times, decreasing = FALSE),
                    c(1, S0[order(times)]))
    surv = surv0(tnew)^exp(expsafe(c(znew %*% bet)))

    h = hazardbaseline(times, status, z, bet, w)
    hazard0 = stepfun(h[,1], c(exp(-200), h[,3]))
    hazard = hazard0(tnew) * exp(expsafe(c(znew %*% bet)))

    sum(log( (pi*hazard*surv)^statusnew * ((1-pi)+pi*surv)^(1-statusnew) ))
}
