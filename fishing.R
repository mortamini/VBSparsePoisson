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
set = syst = matrix(0,8,10)
rownames(set) <- rownames(syst) <- c("LASSO","SCAD","Bernulli-VB","CS-VB","LAPLACE-VB","Bernulli-MCMC","CS-MCMC","LAPLACE-MCMC")
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
	syst[1,iter] = system.time(fit1 <- glmnet(X[,-1],y,
	family="poisson"))[[3]]
	betahat1 = coef(fit1)
	tLL <- fit1$nulldev - deviance(fit1)
	k <- fit1$df
	n <- fit1$nobs
	AICc <- -tLL+2*k+2*k*(k+1)/(n-k-1)
	betahat1 = mu_beta = mu_0 = betahat1[,which.min(AICc)]
	pst1 = 1*(mu_0 !=0)
	yhat1 = exp(X %*% mu_beta)
	yhattest1 = exp(Xtest %*% mu_beta)
	#------------------------------------
	syst[2,iter] = system.time(fit2 <- suppressWarnings(grpreg(X[,-1], y, penalty="grSCAD",
	family="poisson")))[[3]]
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
	syst[3,iter] = system.time(fit3 <- sppoissregvb(X,y,init,prior="Bernulli"))[[3]]
	pst3 = 1*(fit3$pst > 0.5)
	betahat3 = fit3$mu_beta * fit3$pst
	yhat3 = exp(X %*% betahat3)
	yhattest3 = predict(fit3,Xtest,method = "Bernoulli")
	#------------------------------------
	cv3fold <- function(gam, pri){
		nts=trunc(n/3)
		CVf=c()
		for(f in 1:3){
			tests=nts*(f-1)+1:nts
			ts.x=X[tests,]
			ts.y=y[tests]
			tr.y=y[-tests]
			tr.x=X[-tests,]
			fit <- tryCatch({
				sppoissregvb(tr.x,tr.y,init,prior=pri, 
				eps = 1e-5, maxiter = 100)},error=function(e){NULL})
			betahat = tryCatch({fit$mu_beta},error=function(e){NULL})
			betahat[-1] = tryCatch({betahat[-1] * (abs(betahat[-1]) > gam)},error=function(e){NULL})
			yhat = tryCatch({trunc(exp(ts.x %*% betahat))},error=function(e){NULL})
			CVf = c(CVf,tryCatch({mean(yhat - ts.y)^2},error=function(e){NULL}))
		}
		mean(CVf)
	}
	#------------------------------------
	AICf <- function(gam, prin){
		bethat = tryCatch({betahat[[prin]]},error=function(e){NULL})
		bethat[-1] = tryCatch({bethat[-1] * (abs(bethat[-1]) > gam)},error=function(e){NULL})
		yhat = tryCatch({trunc(exp(X %*% bethat))},error=function(e){NULL})
		sse = tryCatch({mean(yhat - y)^2},error=function(e){NULL})
		k = sum(bethat[-1] != 0)
		AIC <- sse + 2*k
		AIC
	}
	#------------------------------------
	c = 1e-2
	syst[4,iter] = system.time(fit4 <- sppoissregvb(X,y,init,prior="CS"))[[3]]
	p0 = min(0.9,max(0.1,sum(pst1)/p))
	a = c(0,abs(fit4$mu_beta[-1]) - 1e-5)
	b = c(abs(fit4$mu_beta[-1]) + 1e-5,0)
	gamseq = sort((a + b)/2)
	if(riskmethod == "cv3fold"){
		risks = sapply(gamseq,cv3fold,"CS")
	} else {
		risks = sapply(gamseq,AICf,4)
	}
	gamopt = gamseq[which.min(risks)]
	pst4 = c(1,1 * (abs(fit4$mu_beta[-1])> gamopt))
	betahat4 = fit4$mu_beta * pst4
	yhat4 = exp(X %*% betahat4)
	yhattest4 = predict(fit4,Xtest,method = "CS")
	#------------------------------------
	syst[5,iter] = system.time(fit5 <- sppoissregvb(X,y,init,prior="Laplace"))[[3]]
	a = c(0,abs(fit5$mu_beta[-1]) - 1e-5)
	b = c(abs(fit5$mu_beta[-1]) + 1e-5,0)
	gamseq = sort((a + b)/2)
	if(riskmethod == "cv3fold"){
		risks = sapply(gamseq,cv3fold,"Laplace")
	} else {
		risks = sapply(gamseq,AICf,5)
	}
	gamopt = gamseq[which.min(risks)]
	pst5 = c(1,1 * (abs(fit5$mu_beta[-1])> gamopt))
	betahat5 = fit5$mu_beta * pst5
	yhat5 = exp(X %*% betahat5)
	yhattest5 = predict(fit5,Xtest,method = "Laplace")
	#------------------------------------
	#------------------------------------
	#----------MCMC samplings------------
	#--------for last iteration----------
	#------------------------------------
	#------------------------------------
	X = X[,-1]
	p = p - 1
	n = 117
	data = list("n","p","y","X")
	inits <- function(){
	list(intercept = 0, 
	beta = rep(0.1,p), gamma = rep(1,p),
	pi = rep(0.5,p), sigma= 1)
	}
	mod_string="model {
	for(i in 1:n){
	for(j in 1:p){
	vec[i,j] <- X[i,j] * gamma[j] * beta[j]
	}
	power[i] <- intercept + sum(vec[i,])
	lambda[i] <- exp(power[i])
	y[i] ~ dpois(lambda[i])
	}
	intercept ~ dnorm(0,alpha0)
	for(j in 1:p){
	gamma[j] ~ dbin(pi[j],1)
	pi[j] ~ dbeta(theta1,theta2)
	beta[j] ~ dnorm(0,alpha[j])
	alpha[j] ~ dgamma(a0,b0)
	}
	alpha0 ~ dgamma(a0,b0)
	theta1 <- 1.5
	theta2 <- 1.5
	a0 <- 1.0E-6
	b0 <- 1.0E-6
	}
	"
	data = list(n=n,p=p,y=y,X=X)
	params=c("intercept", "beta" , "gamma","pi", 
	"alpha0", "alpha")
	inits=function(){
  		inits=list(intercept = 0, 
		beta = rep(0.1,p), gamma = rep(1,p),
		pi = rep(0.5,p), alpha0= 1, alpha = rep(1,p))
	}
	syst[6,iter] = system.time(mod<-jags.model(textConnection(mod_string),data=data,
	inits=inits,n.chain=1))[[3]]
	library(coda)
	syst[6,iter] = syst[6,iter] + system.time(update(mod,1000))[[3]]
	syst[6,iter] = syst[6,iter] + system.time(mod_sim1<-coda.samples(model = 
	mod,variable.names = params,n.iter = 5000))[[3]]
	yhattest6 = trunc(exp(Xtest %*%diag(c(1,colMeans(mod_sim1[[1]][,c(8:10)]))) %*% colMeans(mod_sim1[[1]][,c(11,5:7)])))
#------------------------------------
	inits <- function(){
	list(intercept = init$mu_beta[1], 
	beta = init$mu_beta[-1], z = init$pst[-1],
	pi = init$pst[-1], sigma = 1, ai = 1,
	eta1 = init$mu_beta[-1], eta2 = init$mu_beta[-1])
	}
	data = list(n=n,p=p,y=y,X=X)
	mod_string="model {
	for(i in 1:n){
	for(j in 1:p){
	vec[i,j] <- X[i,j] * beta[j]
	}
	power[i] <- intercept + sum(vec[i,])
	lambda[i] <- exp(power[i])
	y[i] ~ dpois(lambda[i])
	}
	intercept ~ dnorm(0,sigma)
	for(j in 1:p){
	z[j] ~ dbin(pi[j],1)
	pi[j] ~ dbeta(alpha0,beta0)
	eta1[j] ~ dnorm(0,sigma) 
	eta2[j] ~ dnorm(0,sigma2)
	beta[j] <- z[j] * eta1[j] + (1-z[j]) * eta2[j] 
	}
	sigma ~ dgamma(0.5,ai)
	sigma2 <- sigma/c
	ai ~ dgamma(0.5,1)
	alpha0 <- 1.5
	beta0 <- 1.5
	c <- 1.0E-2
	}
	"
	params=c("intercept", "z","pi", 
	"sigma", "ai", "eta1", "eta2")
	inits <- function(){
	inits = list(intercept = init$mu_beta[1], 
	z = init$pst[-1],
	pi = rep(0.5,3), sigma = 1, ai = 1,
	eta1 = init$mu_beta[-1], eta2 = init$mu_beta[-1])
	}
	syst[7,iter] = system.time(mod<-jags.model(textConnection(mod_string),data=data,
	inits=inits,n.chain=1))[[3]]
	library(coda)
	syst[7,iter] = syst[7,iter] + system.time(update(mod,1000))[[3]]
	syst[7,iter] = syst[7,iter] + system.time(mod_sim1<-coda.samples(model = 
	mod,variable.names = params,n.iter = 5000))[[3]]
	data = list(n=n,p=p,y=y,X=X,p0=p0)
	fit7 <- 
	bugs(data,inits,model.file = "cs.txt",
	parameters = c("intercept", "beta" , "z", "pi", 
	"sigma","ai","eta1","eta2"),
	n.chains = 1, n.iter = 10000,n.burnin=5000,
	,n.thin=10)
	yhattest7 = trunc(exp(Xtest %*% colMeans(fit7$sims.matrix[,1:4])))
#------------------------------------
	inits <- function(){
	list(intercept = init$mu_beta[1], 
	beta = init$mu_beta[-1], sigma = 1, 
	tau = rep(1,p), 
	lam2 = 0.5, ai = 1)
	}
	mod_string="model {
	for(i in 1:n){
	for(j in 1:p){
	vec[i,j] <- X[i,j] * beta[j]
	}
	power[i] <- intercept + sum(vec[i,])
	lambda[i] <- exp(power[i])
	y[i] ~ dpois(lambda[i])
	}
	intercept ~ dnorm(0,sigma)
	for(j in 1:p){
	beta[j] ~ dnorm(0,sigma2[j])
	tau[j] ~ dexp(lam)
	sigma2[j] <- 1/tau[j]
	}
	lam <- lam2/2
	lam2 ~ dgamma(nu,delta)
	sigma ~ dgamma(0.5,ai)
	ai ~ dgamma(0.5,1)
	alpha0 <- 1.5
	beta0 <- 1.5
	nu <- 1.0E-4
	delta <- 1.0E-2
	}
	"
	data = list(n=n,p=p,y=y,X=X)
	params=c("intercept", "beta" , "sigma","tau", 
	"lam2", "ai")
	inits=function(){
  		inits=list(intercept = init$mu_beta[1], 
	beta = init$mu_beta[-1], sigma = 1, 
	tau = rep(1,p), 
	lam2 = 0.5, ai = 1)
	}
	syst[8,iter] = tryCatch({system.time(mod<-jags.model(textConnection(mod_string),data=data,
	inits=inits,n.chain=1))[[3]]},error=function(e){syst[8,1]})
	syst[8,iter] = syst[8,iter] + system.time(update(mod,1000))[[3]]
	syst[8,iter] = syst[8,iter] + system.time(mod_sim3<-coda.samples(model = 
	mod,variable.names = params,n.iter = 5000))[[3]]
	yhattest8 = trunc(exp(Xtest %*% colMeans(mod_sim3[[1]][,c(5,2:4)])))
	#------------------------------------
	set[1,iter] = mean((ytest-yhattest1)^2)/var(ytest)
	set[2,iter] = mean((ytest-yhattest2)^2)/var(ytest)
	set[3,iter] = mean((ytest-yhattest3)^2)/var(ytest)
	set[4,iter] = mean((ytest-yhattest4)^2)/var(ytest)
	set[5,iter] = mean((ytest-yhattest5)^2)/var(ytest)
	set[6,iter] = mean((ytest-yhattest6)^2)/var(ytest)
	set[7,iter] = mean((ytest-yhattest7)^2)/var(ytest)
	set[8,iter] = mean((ytest-yhattest8)^2)/var(ytest)
}
round(rowMeans(set),3)
round(apply(set,1,sd),3)
syst[3,] <- syst[3,]/syst[6,]
syst[4,] <- syst[4,]/syst[7,]
syst[5,] <- syst[5,]/syst[8,]
rownames(syst)[3:5] <- 
c("Bernulli-VB/MCMC", 
"CS-VB/MCMC", "LAPLACE-VB/MCMC")
#------------------------------------
#------------------------------------
#------------------------------------
#---------------Box plots------------
#------------------------------------
#------------------------------------
pdf(file = "boxfishing.pdf",width=23,height=7)
dseb = stack(as.data.frame(t(set[,1:(iter-1)])))
names(dseb)[2] <- "Method"
myplot <- ggplot(dseb,aes(x = Method, y = values, color = Method)) + 
geom_boxplot(fatten = 2, notch=FALSE,outlier.shape = 16) +
scale_y_continuous(limits = quantile(dseb$values, c(0.1, 0.9), na.rm = T)) + 
ylim(0.1, 0.3) + xlab('') + 
ylab('Test relative error')
suppressWarnings(print(myplot))
dev.off()
#------------------------------------
pdf(file = "boxsyst.pdf",width=23,height=7)
dseb = stack(as.data.frame(t(syst[3:5,1:(iter-1)])))
names(dseb)[2] <- "Method"
myplot <- ggplot(dseb,aes(x = Method, y = values, color = Method)) + 
geom_boxplot(fatten = 2, notch=FALSE,outlier.shape = 16) +
scale_y_continuous(limits = quantile(dseb$values, c(0.1, 0.9), na.rm = T)) + 
ylim(0.0, 0.05) + xlab('') + 
ylab('Relative computation time')
suppressWarnings(print(myplot))
dev.off()
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
  	geom_function(fun = func, xlim = c(lower[j,h]-0.1,
		upper[j,h]+0.1))
	gg_density %>% 
	plot_credible_interval(lower[j,h], upper[j,h])



