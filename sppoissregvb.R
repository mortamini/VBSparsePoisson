sppoissregvb <- function(X, y, init, prior = "CS", 
		eps = 1e-6, maxiter = 100){
	if(! (prior %in% c("CS","Laplace","Bernulli"))) 
		stop("prior must be one of 'CS','Laplace' or 'Bernulli'")
	mu_beta = init$mu_beta
  	sigma_beta = init$sigma_beta
  	D_beta = sigma_beta + mu_beta %*% t(mu_beta)
	pst = init$pst
	p = ncol(X)
	p0 = max(0.01,min(sum(pst)/p,0.99))
	n = length(y)
	if(nrow(X) != n) stop("number of samples in X and y are not equal!")
	if(length(mu_beta) != p | nrow(sigma_beta) != p)
		stop("mu_beta or sigma_beta in init are not consistent with samples!")
	Pstar=diag(pst)
	sigmoid <- function(v) 1/(1+exp(-v))
	tr <- function(Mat) sum(diag(Mat))
	alphast = rep(1,p)
	betast = rep(1,p)
	alpha0 = 1
  	beta0 = (1 - p0)/p0
	elambda = 0.01
	if(prior == "Bernulli"){
		a0 = 1e-2
		b0 = 1e-2
		a_beta = rep(a0 + 0.5,p)
	} else if(prior == "CS"){
  		c = 1e-2
  		A = 1e-2
  		Esigm2 = 1
  		Eam = 1
    	alpha_sigma2 = 0.5 * (p - 1)
	} else if(prior == "Laplace"){
  		nu = 0.0001
  		delta = 0.01
  		elambda = 0.01
  		ebeta2 = diag(D_beta[-1,-1])
  		astar = elambda
  		bstar = ebeta2
		tau = .tau_update(astar,bstar)
		Etaum1 = tau$Etaum1
		Etau = tau$Etau
  		A = 1e-2
  		Eam1 = 0.01 
  		Etau0m1 = (0.5 * D_beta[1,1] + Eam1)^(-1)
  		pst = rep(1,p)
	}
	Pstar=diag(pst)
    omega = pst %*% t(pst) + Pstar * (diag(p) - Pstar)
	kesi = X %*% diag(pst) %*% mu_beta
	#####kesi = sqrt(diag(X %*% (D_beta * omega) %*% t(X)))
	Mkesi = exp(kesi)*(1-kesi)
	diff = 10 * eps 
  	ELBO = -Inf
  	cntr = 0
  	while(diff > eps & cntr < maxiter){
		cntr = cntr + 1
    		SX = matrix(0,p,p)
    		for(i in 1:n){
      		SX = SX + exp(kesi[i]) * X[i,] %*% t(X[i,])
    		}
		if(prior == "Bernulli"){
			b_beta= b0 + 0.5*diag(D_beta)
    		omega = pst %*% t(pst) + Pstar * (diag(p) - Pstar)
			GM = diag(a_beta/b_beta)
		} else if(prior == "CS"){
			omega = matrix(1,p,p)
			GM = Esigm2 * (diag(pst)+c^(-1)*(diag(p) - diag(pst)))
			Pstar = diag(p)
		} else if(prior == "Laplace"){
			omega = matrix(1,p,p)
			GM = diag(c(Etau0m1,Etaum1))
			Pstar = diag(p)
		}
		sigma_beta = ginv(SX * omega + GM)
  		mu_beta = sigma_beta %*% Pstar %*% t(X) %*% (y-Mkesi)
  		D_beta = sigma_beta + mu_beta %*% t(mu_beta)
		if(prior == "Bernulli"){
    		pst = as.vector(sigmoid( 
				t(y-Mkesi) %*% X * t(mu_beta)+ 
				digamma(alphast) - digamma(betast) -
           		0.5 * diag(SX * D_beta) - 
				0.5 * rowSums(SX * D_beta) * pst +
				diag(SX * D_beta) * pst))
			pst[1] = 1
    		Pstar = diag(pst)
		} else if(prior == "CS"){
			beta_sigma2 = 0.5 * tr(diag(pst) %*% D_beta) + 0.5 * c^(-1) * tr(diag(1-pst) %*% D_beta )+ Eam
    		Esigm2 = alpha_sigma2/beta_sigma2
    		Eam = (1/A + Esigm2)^(-1)
   			pst = as.vector(sigmoid(
				digamma(alphast)- digamma(betast)-
				0.5*Esigm2*diag(D_beta)*(1-c^(-1))))
			pst[1] = 1
    		Pstar = diag(pst)
		} else if(prior == "Laplace"){
    		elambda = (p + nu)/(delta + 0.5 * sum(Etau))
  			ebeta2 = diag(D_beta[-1,-1])
  			astar = elambda
  			bstar = ebeta2
			tau = .tau_update(astar,bstar)
			Etaum1 = tau$Etaum1
			Etau = tau$Etau
			Elogtau = tau$Elogtau
  			Etau0m1 = (0.5 * D_beta[1,1] + Eam1)^(-1)
			Elogtau0 = log(0.5 * D_beta[1,1] + Eam1) - digamma(1)
  			Eam1 = (A^(-1) + Etau0m1)^(-1)
			Eloga = log(A^(-1) + Etau0m1) - digamma(1)
			pst = rep(1,p)
    		Pstar = diag(pst)
		}
			alphast = pst + alpha0
    		betast = beta0 - pst + 1
			Pstar=diag(pst)
    		omega = pst %*% t(pst) + Pstar * (diag(p) - Pstar)
			esi = X %*% diag(pst) %*% mu_beta
			#####kesi = sqrt(diag(X %*% (D_beta * omega) %*% t(X)))
    		Mkesi = exp(kesi) * (1-kesi)
		if(prior == "Bernulli"){
    		elbo = t(y - Mkesi) %*% X %*% Pstar %*% mu_beta - 
				0.5 * tr(D_beta %*% (SX * omega)) -
          		sum( exp(kesi) +  kesi - 0.5 * kesi^2) + 
				0.5 * sum(digamma(a_beta) - log(b_beta))  -
          		0.5 * tr(D_beta %*% diag(a_beta/b_beta)) + 
				sum(pst * log(pst+1e-12) + (1- pst) * log(1- pst+1e-12)) +
          		0.5 * logdet(sigma_beta) + sum(lgamma(alphast) + lgamma(betast)) + 
				0.5 * sum(log(b_beta)) - sum(a_beta*log(b_beta)) + b0 * sum(a_beta/b_beta)
		} else if(prior == "Laplace"){
    		elbo = - t(Mkesi) %*% (1 +X %*% mu_beta) -
           		sum( exp(kesi) * kesi^2/2) -
          		0.5 * tr(SX %*% D_beta) + t(y)%*% X %*% mu_beta -
				0.5 * Elogtau0  - 0.5 * Etau0m1 *  D_beta[1,1] -
				delta * elambda  - 2 * log(A^(-1) + Etau0m1) -
				0.5 * log(0.5 * D_beta[1,1] + Eam1) - 
				Eam1 * Etau0m1 - A^(-1)/(A^(-1) + Etau0m1) + 
				0.5 * logdet(sigma_beta) -
				0.25 * elambda * sum(1/diag(D_beta[-1,-1])) + 
				sum(log(sapply(sqrt(astar*bstar), function(z){BesselK(0.5, z)})))
		} else if(prior == "CS"){
		    elbo = - t(Mkesi) %*% (1+X %*% mu_beta) - 
				sum(exp(kesi) * kesi^2/2) -
           		0.5 * tr(SX %*% D_beta) + t(y)%*% X %*% mu_beta -
           		0.5 * (1-c^(-1))* Esigm2 * sum(Pstar * D_beta) -
	     		0.5 * c^(-1) * Esigm2 * tr(D_beta) + 
           		0.5 * (p+3) * (digamma(alpha_sigma2) - log(beta_sigma2)) -
           		1.5 * (digamma(alpha_sigma2) - log(beta_sigma2)) -
           		0.5 * Esigm2 / (A^(-1) + Esigm2) +
           		sum(pst * (digamma(alphast) - digamma(betast))) +
           		(alpha0 -1) * sum(digamma(alphast) - digamma(alphast+betast)) + 
				(beta0 -1) *sum(digamma(betast) - digamma(alphast+betast)) - 
				log(A^(-1) + Esigm2) - A^(-1) / (A^(-1) + Esigm2) +
	     		0.5 * logdet(Sigma_beta) - 
           		alpha_sigma2 * (log(beta_sigma2) - 1) +lgamma(alpha_sigma2) -
           		(alpha_sigma2 - 1) * (digamma(alpha_sigma2) - log(beta_sigma2)) -
           		sum(Pstar * log(Pstar+1e-5) + (1-Pstar) * log(1-Pstar+1e-5))
		}
	if(elbo == -Inf) elbo = -1e300
    	ELBO =c(ELBO,elbo)
    	diff = abs(ELBO[cntr +1] - ELBO[cntr]) 		
	}
	if(is.null(alphast)) alphast = 1
	if(is.null(betast)) betast = 1
	if(prior == "Bernulli"){
		sigmapars = list(a_beta = a_beta, b_beta = b_beta)
	} else if(prior == "Laplace"){
		sigmapars = list(astar = astar, bstar = bstar, Etau0m1 = Etau0m1)
	} else if(prior == "CS"){
		sigmapars = list(alpha_sigma2 = alpha_sigma2, beta_sigma2 = beta_sigma2)
	}
  	output = list(mu_beta = mu_beta ,
                sigma_beta = sigma_beta,
                pst = pst, 
				 alphast = alphast,
				 betast = betast, 
				 sigmapars = sigmapars,
				elambda = elambda, ELBO = ELBO) 
  	output 
}