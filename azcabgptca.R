if(!("generate_samples" %in% ls())){
	rm(list=ls())
	setwd(choose.dir())
	source("logdet.R")
	source("utility.R")
	source("sppoissregvb.R")
	source("generate_samples.R")
	source("predict.R")
	source("hpd.R")
	source("hpdplot.R")
	source("predict_mcmc.R")
}
if(!all(c("package:grpreg","package:glmnet","package:MASS") %in% search())){
	library(MASS)
	library(glmnet)
	library(grpreg)
	library(GPBayes)
	library(COUNT)
	library(R2OpenBUGS)
	library(ggplot2)
	library(gsl)
	library("rjags")
	library(dplyr)
}
#
data(azcabgptca)
head(azcabgptca)
y = azcabgptca[,5]
ggplot(as.data.frame(y), aes(x=y,y=..ncount..))+
  geom_histogram(color="darkblue", fill="lightblue",position = "stack")+
  ggtitle("Histogram of los for azcabgptca data")

X = as.matrix(azcabgptca[,-5])
hist(y)
n = nrow(X)
p = ncol(X)
set = matrix(0,5,10)
rownames(set) <- c("LASSO","SCAD","Bernulli-VB","CS-VB","LAPLACE-VB")
for(iter in 1:10){
	cat("Iteration",iter,"\n")
	y = azcabgptca[,5]
	X = as.matrix(azcabgptca[,-5])
	n = nrow(X)
	p = ncol(X)
	samp = sample(1:n,trunc(n*0.8))
	ytest = y[-samp]
	Xtest = X[-samp,]
	y  <- y[samp]
	X = X[samp,]
	X = cbind(1,X)
	Xtest = cbind(1,Xtest)
	#------------------------------------
	fit1 <- glmnet(X[,-1],y,family="poisson")
	betahat1 = coef(fit1)
	tLL <- fit1$nulldev - deviance(fit1)
	k <- fit1$df
	#n <- fit1$nobs
	AICc <- -tLL+2*k+2*k*(k+1)/(n-k-1)
	betahat1 = mu_beta = mu_0 = betahat1[,which.min(AICc)]
	pst1 = 1*(mu_0 !=0)
	yhat1 = exp(X %*% mu_beta)
	yhattest1 = exp(Xtest %*% mu_beta)
	#------------------------------------
	fit2 <- suppressWarnings(grpreg(X[,-1], y, penalty="grSCAD",
	family="poisson"))
	betahat3 = coef(fit2)
	tLL <- fit2$loss
	k <- fit2$df
	AICc <- tLL+2*k+2*k*(k+1)/(n-k-1)
	betahat2 = mu_1 = betahat3[,which.min(AICc)]
	pst2 = 1*(mu_1 !=0)
	yhat2 = exp(X %*% mu_1)
	yhattest2 = exp(Xtest %*% mu_1)
	#------------------------------------
	p = ncol(X)
	Sigma_beta = diag(rep(1e-10,p))
	pst0 = pst1
		####c(1,rep(1,p-1))
	astar <- bstar <- 1
	Esigm2 <- 0.1
	alphastar = rep(1,p)
	betastar = rep(1,p)
	alpha0 = beta0 = 1
	init = list(mu_beta = mu_beta,sigma_beta = Sigma_beta,pst = pst0,
		astar = astar, bstar = bstar, Esigm2 = Esigm2, alphastar = alphastar,
		betastar = betastar)
	fit3 <- sppoissregvb(X,y,init,prior="Bernulli")
	pst3 = 1*(fit3$pst > 0.5)
	betahat3 = fit3$mu_beta * fit3$pst
	yhat3 = exp(X %*% betahat3)
	yhattest3 = predict(fit3,Xtest,method = "Bernoulli")
	#------------------------------------
	c = 1e-2
	fit4 <- sppoissregvb(X,y,init,prior="CS")
	p0 = min(0.9,max(0.1,sum(pst1)/p))
	pst4 = 1 * ((fit4$mu_beta)^2>(-2*c/(1+c)*(0.5*log(c)+log(p0/(1-p0)))*fit4$sigmapars$beta_sigma2/(fit4$sigmapars$alpha_sigma2 - 1)))
	####pst4 = 1*((fit4$pst)/sqrt(fit4$pst*(1-fit4$pst))>qnorm(1-0.05/p))
	betahat4 = fit4$mu_beta * pst4
	yhat4 = exp(X %*% betahat4)
	yhattest4 = predict(fit4,Xtest,method = "CS")
	#------------------------------------
	fit5 <- sppoissregvb(X,y,init,prior="Laplace")
	risk <- function(gam){
		tauhat <- diag(fit5$sigma_beta)
		risk <- mean(tauhat) * (p - 2 * 
			sum(abs(fit5$mu_beta)/sqrt(tauhat)<gam)+
			sum((pmax(abs(fit5$mu_beta)/sqrt(tauhat),gam))^2))
		risk
	}
	a = c(0,fit5$mu_beta - 1e-5)
	b = c(fit5$mu_beta + 1e-5,0)
	gamseq = sort(c((a + b)/2,fit5$mu_beta))
	risks = sapply(gamseq,risk)
	gamopt = gamseq[which.min(risks[risks>0])]
	pst5 = 1 * (abs(fit5$mu_beta)>gamopt)
	pst5[1] = 1
	betahat5 = fit5$mu_beta * pst5
	yhat5 = exp(X %*% betahat5)
	yhattest5 = predict(fit5,Xtest,method = "Laplace")
	#------------------------------------
	set[1,iter] = mean((ytest-yhattest1)^2)/var(ytest)
	set[2,iter] = mean((ytest-yhattest2)^2)/var(ytest)
	set[3,iter] = mean((ytest-yhattest3)^2)/var(ytest)
	set[4,iter] = mean((ytest-yhattest4)^2)/var(ytest)
	set[5,iter] = mean((ytest-yhattest5)^2)/var(ytest)
}
round(rowMeans(set),3)
round(apply(set,1,sd),3)


