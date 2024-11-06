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
}
if(!all(c("package:grpreg","package:glmnet","package:MASS") %in% search())){
	library(MASS)
	library(glmnet)
	library(grpreg)
	library(GPBayes)
	library(COUNT)
	library(ggplot2)
	library(dplyr)
}
#
data(fishing)
fishing = fishing[,-c(1,5,6)]
head(fishing)
fishing = na.omit(fishing)
y = fishing[,1]
ggplot(as.data.frame(y), aes(x=y,y=..ncount..))+
  geom_histogram(color="darkblue", fill="lightblue",position = "stack")+
  ggtitle("Histogram of totabund for fishing data")
X = as.matrix(fishing[,-1])
y = y - min(y)
hist(y)
X = t((t(X) - colMeans(X))/apply(X,2,sd))
n = nrow(X)
p = ncol(X)
set = matrix(0,5,10)
rownames(set) <- c("LASSO","SCAD","Bernulli","CS","LAPLACE")
for(iter in 1:10){
	cat("Iteration",iter,"\n")
	y = fishing[,1]
	X = as.matrix(fishing[,-1])
	y = y - min(y)
	X = t((t(X) - colMeans(X))/apply(X,2,sd))
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
	pst0 = c(1,rep(1,p-1))
	astar <- bstar <- 1
	Esigm2 <- 0.1
	alphastar = rep(1,p)
	betastar = rep(1,p)
	alpha0 = beta0 = 1
	init = list(mu_beta = mu_beta,sigma_beta = Sigma_beta,pst = pst0,
		astar = astar, bstar = bstar, Esigm2 = Esigm2, alphastar = alphastar,
		betastar = betastar)
	fit3 <- sppoissregvb(X,y,init,prior="Bernulli")
	betahat3 = fit3$mu_beta * fit3$pst
	pst3 = fit3$pst
	yhat3 = exp(X %*% betahat3)
	#yhattest3 = exp(Xtest %*% betahat3)
	yhattest3 = predict(fit3,Xtest,method = "Bernoulli")
	#------------------------------------
	fit4 <- sppoissregvb(X,y,init,prior="CS")
	pst4 = 1*((fit4$pst)/sqrt(fit4$pst*(1-fit4$pst))>qnorm(1-0.05/p))
	betahat4 = fit4$mu_beta * pst4
	yhat4 = exp(X %*% betahat4)
	#yhattest4 = exp(Xtest %*% betahat4)
	yhattest4 = predict(fit4,Xtest,method = "CS")
	#------------------------------------
	fit5 <- sppoissregvb(X,y,init,prior="Laplace")
	pst5 = fit5$mu_beta/diag(fit5$sigma_beta) > qnorm(1-0.05/p)
	pst5[1] = 1
	betahat5 = fit5$mu_beta * pst5
	yhat5 = exp(X %*% betahat5)
	#yhattest5 = exp(Xtest %*% betahat5)
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
#==============================================
	y = fishing[,1]
	X = as.matrix(fishing[,-1])
	y = y - min(y)
	X = t((t(X) - colMeans(X))/apply(X,2,sd))
	n = nrow(X)
	p = ncol(X)
	X = cbind(1,X)
	#------------------------------------
	fit1 <- glmnet(X[,-1],y,family="poisson")
	betahat1 = coef(fit1)
	tLL <- fit1$nulldev - deviance(fit1)
	k <- fit1$df
	#n <- fit1$nobs
	AICc <- -tLL+2*k+2*k*(k+1)/(n-k-1)
	betahat1 = mu_beta = mu_0 = betahat1[,which.min(AICc)]
	#------------------------------------
	p = ncol(X)
	Sigma_beta = diag(rep(1e-10,p))
	pst0 = c(1,rep(1,p-1))
	astar <- bstar <- 1
	Esigm2 <- 0.1
	alphastar = rep(1,p)
	betastar = rep(1,p)
	alpha0 = beta0 = 1
	fit = list()
	init = list(mu_beta = mu_beta,sigma_beta = Sigma_beta,pst = pst0,
		astar = astar, bstar = bstar, Esigm2 = Esigm2, alphastar = alphastar,
		betastar = betastar)
	fit[[1]] <- sppoissregvb(X,y,init,prior="Bernulli")
	#------------------------------------
	fit[[2]] <- sppoissregvb(X,y,init,prior="CS")
	#------------------------------------
	fit[[3]] <- sppoissregvb(X,y,init,prior="Laplace")
	#------------------------------------
	lower = upper = matrix(0,p,3)
	for(j in 1:p){
		for(h in 1:3){
			tt = seq(fit[[h]]$mu_beta[j] - 3*sqrt(fit[[h]]$sigma_beta[j,j]),
			fit[[h]]$mu_beta[j] + 3*sqrt(fit[[h]]$sigma_beta[j,j]),length.out=100)
			ss = dnorm(tt,fit[[h]]$mu_beta[j],sqrt(fit[[h]]$sigma_beta[j,j]))
			object = list(x = tt, y = ss)
			hhpp <- hpd(object)
			lower[j,h] = hhpp[1]
			upper[j,h] = hhpp[2]
		}
	}

j = 4
h = 3
	func <- function(x) dnorm(x,fit[[h]]$mu_beta[j],
		sqrt(fit[[h]]$sigma_beta[j,j]))
	hpddata = data.frame(x = tt, y = ss)
	gg_density <- ggplot(data=hpddata, aes(x=x)) + 
  	geom_function(fun = func, xlim = c(lower[j,h]-0.5,
		upper[j,h]+0.5))
	gg_density %>% 
	plot_credible_interval(lower[j,h], upper[j,h])



