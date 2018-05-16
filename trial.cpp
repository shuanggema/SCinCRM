#include <algorithm>     // for count
# include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
 
using namespace Rcpp;
using namespace arma;


// [[Rcpp::export]]
colvec expsafe(colvec x) {
	for (int i = 0; i < x.n_rows; i++) {
		 if (x(i) < -150.0) x(i) = -150.0;
	     if (x(i) >  20.0)  x(i) = 20.0;
    }
    return x;		 
}


// abc thresholding rule 
double threshold(double a, double b , double c)
{  
  double val;
  if (fabs(b) <= c) val = 0;
  else {
	    if (b > c) val = (b - c)/a;
	    else val = (b + c)/a;
  }
  //if (abs(val) < 10^-6) val = 0;
  return val;
}


// [[Rcpp::export]]
mat hazardbaseline (colvec  time,// output: ordered deathpoint, deathcount, baseline hazard.
                     colvec  status,
					 mat     z,
					 colvec  beta,
					 colvec  w) {  					 
    uvec ids = find(status == 1);	
    colvec deathtime = sort(unique(time.elem(ids)));
    colvec phi = exp(expsafe(z * beta));
   
    int k = deathtime.n_rows;
    colvec event(k);
    colvec hazard(k);
    for (int it = 0; it < k; it++) {
	    // # of death
	    event(it) = std::count(time.begin(), time.end(), deathtime(it)); 
		// sum_l [w_l exp(z\beta)]
		uvec riskset = find(time >= deathtime(it));
		hazard(it) = event(it)/sum( phi.elem(riskset) % w.elem(riskset) );
    }
     mat res = zeros<mat>(k, 3);
    res.col(0) = deathtime;
    res.col(1) = event;
    res.col(2) = hazard;
    return res;
}


// [[Rcpp::export]]
colvec survbaseline (colvec  time,        //Nelson-Aalen
                     colvec  status,
					 mat     z,
					 colvec  beta,
					 colvec  w) {  					 
   uvec ids = find(status == 1);	
   colvec deathtime = sort(unique(time.elem(ids)));
   colvec coxexp = exp(expsafe(z * beta));
   
   int n = z.n_rows;
   int k = deathtime.n_rows;
   colvec event(k);
   colvec hazard(k);
   for (int it = 0; it < k; it++) {
	    // # of death
	    event(it) = std::count(time.begin(), time.end(), deathtime(it)); 
		// sum_l [w_l exp(z\beta)]
		uvec riskset = find(time >= deathtime(it));
		hazard(it) = event(it)/sum( coxexp.elem(riskset) % w.elem(riskset) );
   }
   
   colvec Hazard(n);
   for (int it = 0; it <n; it++) {
	    // died just before t
	    uvec deathset = find(deathtime <= time(it));
	    Hazard(it) =  sum(hazard.elem(deathset));
		if(time(it) > max(deathtime)) Hazard(it) = R_PosInf;
		if(time(it) < min(deathtime)) Hazard(it) = 0;
   }
   colvec S0 = exp(expsafe(-Hazard));
   return S0;
}


/*
// [[Rcpp::export]]
colvec survbaseline (colvec  time,    //Kaplan-Meier
                     colvec  status,
					 mat     z,
					 colvec  beta,
					 colvec  w) {  					 
   uvec ids = find(status == 1);	
   colvec deathtime = sort(unique(time.elem(ids)));
   colvec coxexp = exp(expsafe(z * beta));
   
   int n = z.n_rows;
   int k = deathtime.n_rows;
   colvec alpha(k);

   for (int it = 0; it < k; it++) {
	   // sum_l [w_l exp(z\beta)]
		uvec riskset = find(time >= deathtime(it));
		// index of the patient died at time tao_j.
		uvec j = find(time == deathtime(it));
		vec risk = coxexp.elem(riskset) % w.elem(riskset);
	    double a = 1.0 - coxexp(j(0)) /sum(risk);
		double b = 1.0/coxexp(j(0));
		alpha(it) = pow(a, b);
   }
   
   colvec S0(n);
   for (int it = 0; it <n; it++) {
	    // died just before t
	    uvec deathset = find(deathtime <= time(it));
	    S0(it) =  prod(alpha.elem(deathset));
		if(time(it) >= max(deathtime)) S0(it) = 0;
		if(time(it) < min(deathtime)) S0(it) = 1;
   }
   return S0;
}
*/

// [[Rcpp::export]]
double coxloglik (colvec  time,  // full loglik for cox model.
                  colvec  status,
				  mat     z,
				  colvec  beta,
				  colvec  w,
				  colvec  S) { 
	vec eta = expsafe(z*beta);
	mat res = hazardbaseline(time, status, z, beta, w); 
	double ell2part = 0;
	for (int l=0; l < z.n_rows; l++)  // to aviod log0.
		if(w(l) != 0) ell2part = ell2part + w(l) * log(S(l));
			 
	//uvec timeindex = sort_index(time);
	uvec ind = find(status == 1);
	//uvec eventindex = timeindex.elem(ind);
	//double ell2 =  sum(res.col(1) % log(res.col(2) % eta.elem(eventindex))) + ell2part;
	double ell2 =  sum(res.col(1) % log(res.col(2))) +  sum(eta.elem(ind)) + ell2part;
	return ell2 ;
}

// [[Rcpp::export]]
double parloglik (colvec  time,  // parital loglik for cox model.
                  colvec  status,
				  mat     z,
				  colvec  beta,
				  colvec  w) { 
	uvec ids = find(status == 1);
	vec eta = expsafe(z*beta);
    colvec deathtime = sort(unique(time.elem(ids)));
    colvec phi = exp(expsafe(z * beta));

    int k = deathtime.n_rows;
    double ell2part = 0;
    for (int it = 0; it < k; it++) {
	    // # of death
	    int event = std::count(time.begin(), time.end(), deathtime(it)); 
		uvec riskset = find(time >= deathtime(it));
		ell2part = ell2part + event * log(sum( phi.elem(riskset) % w.elem(riskset) ));
    }
	double ell2 =  sum(eta.elem(ids)) - ell2part;
	return ell2 ;
}


// [[Rcpp::export]]
colvec loglikDerv (colvec time, 
                   colvec status,
 				   mat    z,
				   colvec beta, 
				   colvec w,
				   int    i       
				   )  {
    uvec ids = find(status == 1);	
    colvec deathtime  = sort(unique(time.elem(ids)));
    colvec wExpEta    = w % exp(expsafe(z * beta));
	colvec wExpEtaz   = wExpEta % z.col(i);
	colvec wExpEtaz2  = wExpEtaz % z.col(i);;
	
    int k = deathtime.n_rows;
	
	colvec x(k);
    colvec zi = z.col(i); 
	for (int j = 0; j < k; j++) {
	     // died at time t
         uvec diedset  = find(time == deathtime(j));
	     x(j) = sum(zi.elem(diedset));
	}
   
	
	colvec event(k);
	colvec sumcoxexp(k);
	colvec sumcoxexpz(k);
	colvec sumcoxexpz2(k);
	
    event(0) = std::count(time.begin(), time.end(), deathtime(0)); 
    uvec riskset  = find(time >= deathtime(0));	
    sumcoxexp(0)  = sum( wExpEta.elem(riskset) );
    sumcoxexpz(0) = sum( wExpEtaz.elem(riskset) );
	sumcoxexpz2(0) = sum( wExpEtaz2.elem(riskset) );
	
	
    for (int j = 1; j < k; j++) {
	     // # of death
	    event(j) = std::count(time.begin(), time.end(), deathtime(j)); 
        uvec riskinterval = find((time >= deathtime(j-1)) % (time < deathtime(j))); 		 
        sumcoxexp(j)   = sumcoxexp(j-1) - sum( wExpEta.elem(riskinterval) ); 	
        sumcoxexpz(j)  = sumcoxexpz(j-1) - sum( wExpEtaz.elem(riskinterval) ); 	
        sumcoxexpz2(j) = sumcoxexpz2(j-1) - sum( wExpEtaz2.elem(riskinterval) ); 		
    }	

	colvec der(2); 
    der(0) = sum(x - event % sumcoxexpz / sumcoxexp);
	der(1) = sum(event % (pow((sumcoxexpz)/(sumcoxexp), 2.0) - sumcoxexpz2/(sumcoxexp)));
	return der;
}



// [[Rcpp::export]]
colvec cdUpdate(colvec time,
				colvec status,
                mat    z,
				mat    z1,      // (1, z) 
				colvec beta0,
				colvec gamma0,
                colvec w,    
				double lambda1,
				double lambda2,
				double lambda3,
				double nu,
				double xi,
				double eps,
				int    maxit
                ){    
					
	const int n = z.n_rows;
	const int p = z.n_cols;	 
	const int q = z1.n_cols;	
	
	vec    theta0   = join_vert(gamma0, beta0);
	vec    theta    = theta0;
	vec    gamma    = gamma0;
	vec    beta     = beta0;
	vec    pi(n), eta(n), expeta(n);
	
 	double converge = 1000, der1, der2;
	double lambda = lambda1 + lambda2 + lambda3;
	double a=0, b=0, c=0;
	int    it = 0;

	while (it < maxit &  converge > eps) {
		
		// update gamma[0]
		expeta = exp(expsafe(z1 * gamma));
		pi  = expeta/(1+expeta);
		der1 = sum(w - pi);
		der2 = -n/4;// -sum( pi % (1-pi) % z1.col(0) % z1.col(0)); //
		gamma(0) = gamma(0) - der1/der2;
		
	    for (int i=1; i < q; i++){
			expeta = exp(expsafe(z1 * gamma));
			pi  = expeta/(1+expeta);
			der1 = sum((w - pi) % z1.col(i));
		    der2 =-n/4;//-sum( pi % (1-pi) % z1.col(i) % z1.col(i) );// 
			
			if (lambda == 0) {     // when z and z1 are different covariates, we need to get this term independently.
				gamma(i) = gamma(i) - der1/der2;		
			}
			else {
				if (xi != 0) {
				    a = -2.0/n*der2 - (fabs(gamma(i)) < nu*lambda1)/nu + lambda3/pow(fabs(gamma(i)) + xi, 2.0);
			        b = 2.0/n*(-der2*gamma(i)  + der1) + lambda3*beta(i-1)/(fabs(beta(i-1)) + xi)/(fabs(gamma(i)) + xi);
		        } else  {
				    a = -2.0/n*der2 - (fabs(gamma(i)) < nu*lambda1)/nu + lambda3;
			        b = 2.0/n*(-der2*gamma(i)  + der1) + lambda3*beta(i-1);  
			    }
			    c = lambda1 * (fabs(gamma(i)) < nu*lambda1);
				gamma(i) = threshold(a, b, c);
			}
	        	
		}
		
		// update beta		 
        for (int i=0; i < p; i++) {
			colvec der = loglikDerv(time, status, z, beta, w, i);  
			if (lambda == 0 ) {
				beta(i) = beta(i) - der(0)/der(1);
			} else { 
			    if (xi != 0) { 
			        a = -2.0/n * der(1) - (fabs(beta(i)) < nu*lambda2)/nu + lambda3/pow(fabs(beta(i)) + xi, 2.0);	
			        b = 2.0/n * (-der(1)*beta(i) + der(0)) + lambda3*gamma(i+1)/(fabs(beta(i)) + xi)/(fabs(gamma(i+1)) + xi); 
			    }
		        else {
				    a = -2.0/n * der(1) - (fabs(beta(i)) < nu*lambda2)/nu + lambda3;	
			        b = 2.0/n * (-der(1)*beta(i) + der(0)) + lambda3*gamma(i+1); 
			    }
			    c = lambda2 * (fabs(beta(i)) < nu*lambda2); 
			    beta(i) = threshold(a, b, c); 
			}
		}
		 
		theta = join_vert(gamma, beta);	 
	    converge = sum(square(theta - theta0));
	    it++;
    }	 
	  
	return theta;
}

// [[Rcpp::export]]
vec grid(colvec time,   
         colvec status,
		 mat    z,
		 mat    z1){ // calculate the lambda.max which produces zero estimates.
		
		int n = z.n_rows;
        int p = z.n_cols;
        int q = z1.n_cols;		
		vec w = status;
        vec beta(p, fill::zeros); 			
		vec derell1(q, fill::zeros);	  
		vec derell2(p, fill::zeros);	
        vec lambdamax(2);
		
		double pi = 0.5;
		double converge = 1000;
		colvec eta = exp(expsafe(z * beta));;
		colvec S(n, fill::zeros);
		colvec S0(n, fill::zeros);
		colvec wOld(n);
		int iter = 0;
		while ((converge > 0.001) & (iter < 100)) {
			wOld = w;
		    S0 = survbaseline(time, status, z, beta, w);  
		    for (int ii=0; ii < n; ii++) {
			    S(ii)  = pow(S0(ii), eta(ii));
		    }	    	
	     	w = status + (1-status) % (pi * S)/((1-pi) + pi * S);
			converge = sum(pow(w - wOld, 2));
			iter++;
		}
		
		for  (int i=0; i < q; i++) {
			derell1(i) = sum((w - pi) % z1.col(i));
		}
		
		for (int i=0; i < p; i++) {
	        vec der = loglikDerv(time, status, z, beta, w, i);
			derell2(i) = der(0); 
		}
		lambdamax(0) = max(abs(derell1))*2.0/n;
		lambdamax(1) = max(abs(derell2))*2.0/n;
		
		return lambdamax;	  
}

// [[Rcpp::export]]
List EPM (colvec time,
          colvec status,
		  mat    z,  // cox 
		  mat    z1, // logistic, c(1, z)
		  int    nlambda1,
		  int    nlambda2,
		  double lambdaminratio,
		  double lambda3,
		  double nu,
		  double xi,
		  double eps = 10^-4,
		  int    maxitCD = 1,
		  int    maxitEM = 10,
		  int    maxDF = 30
         ) {
// inputs
    const int n = z.n_rows;
    const int p = z.n_cols;
	const int q = z1.n_cols;

	// lambda sequences 
    vec lambdamax = grid(time, status, z, z1);
	vec lambda1 = log(linspace(exp(lambdamax(0)), exp(lambdamax(0)*lambdaminratio), nlambda1));
	vec lambda2 = log(linspace(exp(lambdamax(1)), exp(lambdamax(1)*lambdaminratio), nlambda2));
	
// containers
    cube gamma(q, nlambda1, nlambda2); // initialize gamma to 0
    gamma.fill(0.0); 
    cube beta(p, nlambda1, nlambda2);  // initialize betas to 0
    beta.fill(0.0);  
	colvec eta(n, fill::zeros);
	colvec pi(n, fill::zeros);    // uncure probability
	colvec S(n, fill::zeros);     // suvival function
	mat loglik1(nlambda1, nlambda2, fill::zeros);
	mat loglik2(nlambda1, nlambda2, fill::zeros);
    mat df(nlambda1, nlambda2, fill::zeros);
	
// algorithm
	colvec w  = status;
	cube W(n, lambda1.n_rows, lambda2.n_rows);
	W.fill(0.0);
	colvec S0 = survbaseline(time, status, z, beta.slice(0).col(0), w);  
	colvec S0new(n, fill::zeros); 
	vec theta(p+q, fill::zeros), thetanew(p+q, fill::zeros);

	double converge;
	uword i, j;
	int it, df_beta, df_gamma;
	colvec bet, gam;
    for (j = 0; j < lambda2.n_rows; j++) {
	    for (i = 0; i < lambda1.n_rows; i++) {
            converge = 1000.0; 
		    it = 0; 
			df_beta = 0;
			df_gamma = 0;
			
			if (i ==0 & j == 0) {
				bet = beta.slice(j).col(i);
				gam = gamma.slice(j).col(i);
			    } else if (i == 0 & j !=0 ) {
					bet = beta.slice(j-1).col(i);
					gam = gamma.slice(j-1).col(i);
			    }  else { 
					bet = beta.slice(j).col(i-1);
					gam = gamma.slice(j).col(i-1);
			}
            while ( (it < maxitEM) & (converge > eps) ) {
				vec ezgamma = exp(expsafe(z1*gam));
		        pi = ezgamma/(ezgamma+1.0);
		        eta = exp(expsafe(z * bet));
		        for (int ii=0; ii < n; ii++) {
			        S(ii)  = pow(S0(ii), eta(ii));
			    }
		        // E step	
		        w = status + (1-status) % (pi % S)/((1-pi) + pi % S);
		
	            // M step, update gamma and beta --- coordinate descent 		   
		        thetanew = cdUpdate(time, status, z, z1, bet, 
			                    gam, w,
		                        lambda1(i), lambda2(j), lambda3,
								nu, xi, eps, maxitCD);
		        S0new    = survbaseline(time, status, z, theta(span(q, p+q-1)), w);  
		        converge = sum(pow(theta - thetanew, 2)) + sum(pow(S0 - S0new, 2));
		   
		        S0 = S0new;
		        theta = thetanew;
			    gamma.slice(j).col(i) = theta(span(0, q-1));				
			    beta.slice(j).col(i) = theta(span(q, p+q-1)); 
				bet = beta.slice(j).col(i);
				gam = gamma.slice(j).col(i);				
	            it++;
            }
			W.slice(j).col(i) = w;
			loglik1(i, j) = sum(w % log(pi) + (1-w) % log(1-pi));
			loglik2(i, j) = parloglik (time, status, z, theta(span(q, p+q-1)), w);
			//loglik2(i, j) = coxloglik(time, status, z, theta(span(q, p+q-1)), w, S);
			
			df_gamma = sum(gamma.slice(j).col(i) != 0);
			df_beta  = sum(beta.slice(j).col(i) != 0);
			df(i, j) = df_gamma + df_beta;
			if (df(i, j) > maxDF) break;
	    }
		if (df(0, j) > maxDF) break;
    }	
// returns
    List ret ;
    ret["gamma"] = gamma;
	ret["beta"] = beta;
	ret["lambda1"] = lambda1;
	ret["lambda2"] = lambda2;
	ret["loglik1"]  = loglik1;
	ret["loglik2"]  = loglik2;
	ret["loglik"]  = loglik1+loglik2;
	ret["W"] = W;
	ret["S0"]  = S0;
	ret["eta"]  = eta;
	ret["pi"]     = pi;
	ret["converge"]   = converge;
	ret["iterEM"]   = it;
	ret["S"] = S;
    return(ret) ;
}