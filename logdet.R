logdet <- function(Sig){
	s <- svd(Sig)
	d <- s$d
	det <- 1
	i <- 0
	while(det > 1e-300 & i < length(d)){
		i = i + 1
		det = det * d[i]
	}
	if(det == Inf) det = 1e300
	#det = det / d[i]
	log(det) 
}