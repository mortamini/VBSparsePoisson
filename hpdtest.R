hpdtest <- function(model,j){
	tt = seq(model$mu_beta[j] - 5*sqrt(abs(model$sigma_beta[j,j])),
			model$mu_beta[j] + 5*sqrt(abs(model$sigma_beta[j,j])),length.out=100)
			ss = dnorm(tt,model$mu_beta[j],sqrt(abs(model$sigma_beta[j,j])))
			object = list(x = tt, y = ss)
			hhpp <- hpd(object)
			lower = hhpp[1]
			upper = hhpp[2]
	pst = 1 - 1 * (lower < 0 & upper > 0) 
	pst
}