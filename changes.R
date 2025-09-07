changes = function(init){
	inits = list()
	p = length(init$mu_beta)
	errors = seq(0.1,0.5,0.1)
	direction = 2 * rbinom(p,1,0.5) - 1
	for(i in 1:5){
		inits[[i]] = init
		inits[[i]]$mu_beta = init$mu_beta * 
			(1 + direction * errors[i])
	}
	pst = init$pst
	k0 = sum(pst == 0)
	k1 = sum(pst == 1)
	J = min(c(k1,k0,3))
	if(J > 0){
		for(j in 1:J){
			inits[[j+5]] = init
			inds0 = which(pst == 0)
			inds1 = which(pst == 1)
			sam0 = sample(inds0,j)
			sam1 = sample(inds1,j)
			inits[[j+5]]$pst[sam0] = 1
			inits[[j+5]]$pst[sam1] = 0
		}
	}
	list(inits = inits, err = errors, J = J)
}