library(mvtnorm)
generate_samples <- function(n = 30, p = 5, 
		beta , mu_x = 0, sigma_x = 1, 
		rhox = 0.3, link = "log",family="Poisson",r=3){
	Sigma = diag(p-1)
	Sigma = sigma_x * rhox^(abs(row(Sigma)-col(Sigma)))
	mu = rep(mu_x,p-1)
	X = rmvnorm(n, mu, Sigma)
	X = cbind(1,X)
	if(family == "Poisson"){
		if(link == "log"){
			lambda = exp(X %*% beta)
		}else if(link == "unity"){
			lambda = X %*% beta		
		}else{
			stop("Unknown link!")
		}
		y = rpois(n,lambda)
	} else if(family == "nbinom"){
		y = rnbinom(n,r,prob = 1/(1+exp(X %*% beta)))
	}else{
		stop("Unknown family!")
	}
	list(X = X, y = y)
}