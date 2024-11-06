predict <- function(fit,Xtest,method){
	nt = nrow(Xtest)
	yhat = c()
	for(i in 1:nt){
		x0 = Xtest[i,]
		y0 = ytest[i]
		mu_beta0 = fit$mu_beta 
		sigma_beta0 = fit$sigma_beta 
		if(method == "CS" | method == "Laplace"){
			mean = as.vector(t(x0)%*%mu_beta0)
			sd = sqrt(as.vector(t(x0)%*%sigma_beta0%*%(x0)))
		}else if(method == "Bernoulli"){
			Gamma = diag(fit$pst)
			mean = as.vector(t(x0)%*%Gamma%*%mu_beta0)
			sd = sqrt(as.vector(t(x0)%*%Gamma%*%sigma_beta0%*%Gamma%*%(x0)))
		}
		yj = round(exp(mean-10*sd):exp(mean+10*sd)+1)
		if(length(yj)>50) yj = round((exp(mean)-25):(exp(mean)+25))
		pyj = sapply(yj,function(y00){
			integrate(function(x){exp(-x+y00*log(x)-lgamma(y00+1))*dlnorm(x,mean,sd)},
			min(yj),max(yj))$value})
		yhat = c(yhat,yj[which.max(pyj)])
	}
	yhat
}