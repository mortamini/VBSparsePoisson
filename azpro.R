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
	library(gsl)
	library(dplyr)
	library(mpath)
}
#
riskmethod = "aic"
#
data(azpro)
head(azpro)
y = azpro[,1]
ggplot(as.data.frame(y), aes(x=y,y=..ncount..))+
  geom_histogram(color="darkblue", fill="lightblue",position = "stack")+
  ggtitle("Histogram of los for azpro data")

X = as.matrix(azpro[,-1])
#
df = data.frame(cbind(X,y))
formul = as.formula(paste("y~",paste(colnames(X),collapse="+")))
model_poisson <- glm(formul, data=df,family = poisson())
AIC(model_poisson)
model_nb <- glm.nb(formul, data=df)
AIC(model_nb)
sum(residuals(model_nb)^2)
sum(residuals(model_poisson)^2)
#
library(AER)
dispersiontest(model_poisson)
#
hist(y)
n = nrow(X)
p = ncol(X)
set = matrix(0,6,10)
rownames(set) <- c("LASSO","SCAD",
"Bernulli-VB","CS-VB","LAPLACE-VB", "LASSO-NB")
pst =array(NA,dim = c(p+1,10,5))
#
for(iter in 1:10){
	cat("Iteration",iter,"\n")
	y = azpro[,1]
	X = as.matrix(azpro[,-1])
	n = nrow(X)
	p = ncol(X)
	samp = sample(1:n,trunc(n*0.8))
	ytest = y[-samp]
	Xtest = X[-samp,]
	y  <- y[samp]
	X = X[samp,]
	X = cbind(1,X)
	Xtest = cbind(1,Xtest)
	dframe = data.frame(cbind(X[,-1],y))
	#------------------------------------
	fit0 <- glmregNB(y ~ .,alpha=1, data = dframe)
	betahat0 = coef(fit0)
	tLL <- fit0$nulldev - deviance(fit0)
	k <- fit0$df
	n <- nrow(X)
	AICc <- -tLL+2*k+2*k*(k+1)/(n-k-1)
	betahat0 = betahat0[,which.min(AICc)]
	pst0 = 1*(betahat0 !=0)
	yhat0 = exp(X %*% betahat0)
	yhattest0 = exp(Xtest %*% betahat0)
	#------------------------------------
	fit1 <- glmnet(X[,-1],y,family="poisson")
	betahat1 = coef(fit1)
	tLL <- fit1$nulldev - deviance(fit1)
	k <- fit1$df
	n <- fit1$nobs
	AICc <- -tLL+2*k+2*k*(k+1)/(n-k-1)
	betahat1 = mu_beta = mu_0 = betahat1[,which.min(AICc)]
	pst1 = 1*(mu_0 !=0)
	yhat1 = exp(X %*% mu_beta)
	yhattest1 = exp(Xtest %*% mu_beta)
	pst[,iter,1] = pst1
	#------------------------------------
	fit2 <- suppressWarnings(grpreg(X[,-1], y, penalty="grSCAD",
	family="poisson"))
	betahat3 = coef(fit2)
	####tLL <- fit2$loss
	tLL <- fit2$deviance
	k <- fit2$df
	AICc <- tLL+2*k+2*k*(k+1)/(n-k-1)
	betahat2 = mu_1 = betahat3[,which.min(AICc)]
	pst2 = 1*(mu_1 !=0)
	yhat2 = exp(X %*% mu_1)
	yhattest2 = exp(Xtest %*% mu_1)
	pst[,iter,2] = pst2
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
	pst[,iter,3] = pst3
	betahat3 = fit3$mu_beta * fit3$pst
	yhat3 = exp(X %*% betahat3)
	#yhattest3 = exp(Xtest %*% betahat3)
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
	AICf <- function(gam, bet){
		bethat = tryCatch({bet},error=function(e){NULL})
		bethat[-1] = tryCatch({bethat[-1] * (abs(bethat[-1]) > gam)},error=function(e){NULL})
		yhat = tryCatch({trunc(exp(X %*% bethat))},error=function(e){NULL})
		sse = tryCatch({mean(yhat - y)^2},error=function(e){NULL})
		k = sum(bethat[-1] != 0)
		AIC <- sse + 2*k
		AIC
	}
	#------------------------------------
	c = 1e-2
	fit4 <- sppoissregvb(X,y,init,prior="CS")
	p0 = min(0.9,max(0.1,sum(pst1)/p))
	a = c(0,abs(fit4$mu_beta[-1]) - 1e-5)
	b = c(abs(fit4$mu_beta[-1]) + 1e-5,0)
	gamseq = sort((a + b)/2)
	if(riskmethod == "cv3fold"){
		risks = sapply(gamseq,cv3fold,"CS")
	} else {
		risks = sapply(gamseq,AICf,fit4$mu_beta)
	}
	gamopt = gamseq[which.min(risks)]
	pst4 = c(1,1 * (abs(fit4$mu_beta[-1])> gamopt))
	pst[,iter,4] = pst4
	betahat4 = fit4$mu_beta * pst4
	yhat4 = exp(X %*% betahat4)
	yhattest4 = predict(fit4,Xtest,method = "CS")
	#------------------------------------
	fit5 <- sppoissregvb(X,y,init,prior="Laplace")
	a = c(0,abs(fit5$mu_beta[-1]) - 1e-5)
	b = c(abs(fit5$mu_beta[-1]) + 1e-5,0)
	gamseq = sort((a + b)/2)
	if(riskmethod == "cv3fold"){
		risks = sapply(gamseq,cv3fold,"Laplace")
	} else {
		risks = sapply(gamseq,AICf,fit5$mu_beta)
	}
	gamopt = gamseq[which.min(risks)]
	pst5 = c(1,1 * (abs(fit5$mu_beta[-1])> gamopt))
	pst[,iter,5] = pst5
	betahat5 = fit5$mu_beta * pst5
	yhat5 = exp(X %*% betahat5)
	yhattest5 = predict(fit5,Xtest,method = "Laplace")
	#------------------------------------
	set[1,iter] = mean((ytest-yhattest1)^2)/var(ytest)
	set[2,iter] = mean((ytest-yhattest2)^2)/var(ytest)
	set[3,iter] = mean((ytest-yhattest3)^2)/var(ytest)
	set[4,iter] = mean((ytest-yhattest4)^2)/var(ytest)
	set[5,iter] = mean((ytest-yhattest5)^2)/var(ytest)
	set[6,iter] = mean((ytest-yhattest0)^2)/var(ytest)
}
round(apply(set,1,mean),3)
round(apply(set,1,sd),3)
pstperc = apply(pst,c(1,3),mean)
rownames(pstperc) = names(pst1)
colnames(pstperc) = rownames(set[1:5,])
t(pstperc)




