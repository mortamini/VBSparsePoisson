if(!("generate_samples" %in% ls())){
	rm(list=ls())
	setwd(choose.dir())
	source("logdet.R")
	source("utility.R")
	source("sppoissregvb.R")
	source("kdensity.R")
	source("hpd.R")
	source("hpdplot.R")
	source("predict.R")
	source("generate_samples.R")
}
if(!("package:dplyr" %in% search())){
	library(MASS)
	library(glmnet)
	library(grpreg)
	library(GPBayes)
	library(R2OpenBUGS)
	library(ggplot2)
	library(gsl)
	library("rjags")
	library(dplyr)
	library(tidyr)
}
#------------------------------------
Iterations = 1000
seb = se = set = FPR = FNR = syst = matrix(0,8,Iterations)
rownames(seb) <- rownames(se) <-
rownames(set) <- rownames(FPR) <-
rownames(FNR) <- rownames(syst) <- c("LASSO","SCAD","Bernoulli-VB","CS-VB",
"LAPLACE-VB","Bernoulli-MCMC","CS-MCMC",
"LAPLACE-MCMC")
riskmethod = "aic"
#------------------------------------
fit = betahat = pst = yhat = yhattest =  list()
accuracyl = accuracy2l = accuracy3l = coverage = 
accuracy4l = accuracy5l = accuracy6l = list()
p = 10
n = 100
mu_beta = 0.7
sigma_beta = 0.5
beta = rnorm(p, mu_beta, sigma_beta)
z = c(1, c(0,1,0,0,0,1,0,1,0))
realbeta = beta * z
cover = matrix(0,p,6)
rownames(cover) <- paste0("beta",0:(p-1))
colnames(cover) <- c("Bernoulli-VB","CS-VB",
"LAPLACE-VB","Bernoulli-MCMC","CS-MCMC",
"LAPLACE-MCMC")
for(iter in 1:Iterations){
	p = 10
	n = 100
	cat("Iteration",iter,"\n")
	sam = generate_samples(n=n, p=p, beta = realbeta, 
		mu_x = 0.1, sigma_x = 1)
	samp = sample(1:n,trunc(0.8*n))
	y  <- sam$y[samp]
	X = sam$X[samp,]
	ytest = sam$y[-samp]
	Xtest = sam$X[-samp,]
	realp = 1*(realbeta !=0)
	n <- length(y)
	#------------------------------------
	cat("LASSO \n")
	p = ncol(X)
	syst[1,iter] = system.time(fit[[1]] <- glmnet(X[,-1],y,family="poisson"))[[3]]
	betahat1 = coef(fit[[1]])
	tLL <- fit[[1]]$nulldev - deviance(fit[[1]])
	k <- fit[[1]]$df
	n <- fit[[1]]$nobs
	AICc <- -tLL+2*k+2*k*(k+1)/(n-k-1)
	betahat[[1]] = mu_beta = mu_0 = betahat1[,which.min(AICc)]
	pst[[1]] = 1*(mu_0 !=0)
	yhat[[1]] = exp(X %*% mu_beta)
	yhattest[[1]] = exp(Xtest %*% mu_beta)
	#------------------------------------
	cat("SCAD \n")
	syst[2,iter] = system.time(fit[[2]] <- suppressWarnings(grpreg(X[,-1], y, penalty="grSCAD",
	family="poisson")))[[3]]
	betahat3 = coef(fit[[2]])
	tLL <- fit[[2]]$loss
	k <- fit[[2]]$df
	AICc <- tLL+2*k+2*k*(k+1)/(n-k-1)
	betahat[[2]] = mu_1 = betahat3[,which.min(AICc)]
	pst[[2]] = 1*(mu_1 !=0)
	yhat[[2]] = exp(X %*% mu_1)
	yhattest[[2]] = exp(Xtest %*% mu_1)
	#------------------------------------
	cat("Bernoulli-VB \n")
	Sigma_beta = diag(rep(1e-5,p))
	pst0 = pst[[1]]
	astar <- bstar <- 1
	Esigm2 <- 1
	alphastar = rep(1,p)
	betastar = rep(1,p)
	alpha0 = beta0 = 1
	init = list(mu_beta = mu_beta,sigma_beta = Sigma_beta,pst = pst0,
		astar = astar, bstar = bstar, Esigm2 = Esigm2, alphastar = alphastar,
		betastar = betastar)
	syst[3,iter] = system.time(fit[[3]] <- tryCatch({
	sppoissregvb(X,y,init,prior="Bernulli", 
		eps = 1e-5, maxiter = 100)},error=function(e){NULL}))[[3]]
	betahat[[3]] = tryCatch({fit[[3]]$mu_beta * fit[[3]]$pst},error=function(e){NULL})
	pst[[3]] = tryCatch({1*(fit[[3]]$pst > 0.5)},error=function(e){NULL})
	yhat[[3]] = tryCatch({trunc(exp(X %*% betahat[[3]]))},error=function(e){NULL})
	yhattest[[3]] = tryCatch({predict(fit[[3]],Xtest,method ="Bernoulli")},error=function(e){NULL})
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
	cat("CS-VB \n")
	syst[4,iter] = system.time(fit[[4]] <- tryCatch({
	sppoissregvb(X,y,init,prior="CS", 
		eps = 1e-12, maxiter = 100)},error=function(e){NULL}))[[3]]
	p0 = max(0.01,min(sum(pst[[1]])/p,0.99))
	a = c(0,abs(fit[[4]]$mu_beta[-1]) - 1e-5)
	b = c(abs(fit[[4]]$mu_beta[-1]) + 1e-5,0)
	gamseq = sort((a + b)/2)
	if(riskmethod == "cv3fold"){
		risks = sapply(gamseq,cv3fold,"CS")
	} else {
		risks = sapply(gamseq,AICf,4)
	}
	gamopt = gamseq[which.min(risks)]
	pst[[4]] = c(1,1 * (abs(fit[[4]]$mu_beta[-1])> gamopt))
	betahat[[4]] = tryCatch({fit[[4]]$mu_beta * pst[[4]]},error=function(e){NULL})
	yhat[[4]] = tryCatch({trunc(exp(X %*% betahat[[4]]))},error=function(e){NULL})
	yhattest[[4]] = tryCatch({predict(fit[[4]],Xtest,method ="CS")},error=function(e){NULL})
	#------------------------------------
	cat("Laplace-VB \n")
	syst[5,iter] = system.time(fit[[5]] <- tryCatch({
	sppoissregvb(X,y,init,prior="Laplace", 
		eps = 1e-5, maxiter = 100)},error=function(e){NULL}))[[3]]
	a = c(0,abs(fit[[5]]$mu_beta[-1]) - 1e-5)
	b = c(abs(fit[[5]]$mu_beta[-1]) + 1e-5,0)
	gamseq = sort((a + b)/2)
	if(riskmethod == "cv3fold"){
		risks = sapply(gamseq,cv3fold,"Laplace")
	} else {
		risks = sapply(gamseq,AICf,5)
	}
	gamopt = gamseq[which.min(risks)]
	pst[[5]] = c(1,1 * (abs(fit[[5]]$mu_beta[-1])> gamopt))
	betahat[[5]] = tryCatch({fit[[5]]$mu_beta * pst[[5]]},error=function(e){NULL})
	yhat[[5]] = tryCatch({trunc(exp(X %*% betahat[[5]]))},error=function(e){NULL})
	yhattest[[5]] = tryCatch({predict(fit[[5]],Xtest,method ="Laplace")},error=function(e){NULL})
	#------------------------------------
	#------------------------------------
	#------------------------------------
	#----------MCMC samplings------------
	#------------------------------------
	#------------------------------------
	cat("Bernoulli-MCMC \n")
	XM = X[,-1]
	pM = p - 1
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
	data = list(n=n,p=pM,y=y,X=XM)
	params=c("intercept", "beta" , "gamma","pi", 
	"alpha0", "alpha")
	inits=function(){
	  inits=list(intercept = 0, 
		beta = rep(0.1,pM), gamma = rep(1,pM),
		pi = rep(0.5,pM), alpha0= 1, alpha = rep(1,pM))
	}
	syst[6,iter] = system.time(mod<-jags.model(textConnection(mod_string),data=data,
	inits=inits,n.chain=1))[[3]]
	library(coda)
	syst[6,iter] = syst[6,iter] + system.time(update(mod,5000))[[3]]
	syst[6,iter] = syst[6,iter] + system.time(mod_sim1<-coda.samples(model = 
	mod,variable.names = params,n.iter = 10000, thin = 10))[[3]]
	parsims = list()
	parsims[[1]] = as.matrix(mod_sim1)
	parsims[[1]] = cbind(parsims[[1]][,29],
	parsims[[1]][,c(11:19)],
	parsims[[1]][,c(20:28)],
	parsims[[1]][,30:38],
	parsims[[1]][,10],
	parsims[[1]][,1:9])
	colnames(parsims[[1]])[c(1,29)] <- c("intercept","alpha[0]")
	pst[[6]] = 1 * (c(1,colMeans(parsims[[1]][,11:19])) > 0.5)
	betahat[[6]] = colMeans(parsims[[1]][,1:10]) * pst[[6]]
	yhat[[6]] = trunc(exp(X %*% diag(pst[[6]]) %*% betahat[[6]]))	
	yhattest[[6]] = trunc(exp(Xtest %*% diag(pst[[6]]) %*% betahat[[6]]))	
	#------------------------------------
	cat("CS-MCMC \n")
	data = list(n=n,p=pM,y=y,X=XM,p0=p0)
	inits <- function(){
	list(intercept = init$mu_beta[1], 
	beta = init$mu_beta[-1], z = init$pst[-1],
	pi = init$pst[-1], sigma = 1, ai = 1,
	eta1 = init$mu_beta[-1], eta2 = init$mu_beta[-1])
	}
	if(iter == 1){timelim = Inf}else{timelim = syst[7,1]*1.2}
	syst[7,iter] = system.time(fit[[7]] <- tryCatch({
		setTimeLimit(elapsed = timelim)
		bugs(data,inits,model.file = "cs.txt",
		parameters = c("intercept", "beta" , "z", "pi", 
		"sigma","ai","eta1","eta2"),
		n.chains = 1, n.iter = 10000,n.burnin=5000,
		n.thin=10,debug = F)
	}, error = function(e) {
  		if (grepl("time limit", e$message)) {
    		message(paste("Execution timed out after", timelim, "seconds"))
    		return(NULL)
  		} else {
    		return(NULL)
  		}
	}, finally = {
  		setTimeLimit(elapsed = Inf)  # Reset time limit
	}))[[3]]
	pars = tryCatch({fit[[7]]$summary[,1]}, error = function(e) {return(NULL)})
	pst[[7]] = tryCatch({c(1,fit[[7]]$summary[11:19,3])}, error = function(e) {return(NULL)})
	parsims[[2]] = tryCatch({fit[[7]]$sims.matrix[(1:5000) %% 10 == 1, ]}, error = function(e) {return(NULL)})
	betahat[[7]] = tryCatch({colMeans(parsims[[2]][,1:10])}, error = function(e) {return(NULL)})
	yhat[[7]] = tryCatch({trunc(exp(X %*% betahat[[7]]))}, error = function(e) {return(NULL)})	
	yhattest[[7]] = tryCatch({trunc(exp(Xtest %*% betahat[[7]]))}, error = function(e) {return(NULL)})	
	#------------------------------------
	cat("Laplace-MCMC \n")
	inits <- function(){
	list(intercept = init$mu_beta[1], 
	beta = init$mu_beta[-1], sigma = 1, 
	tau = rep(1,pM), 
	lam2 = 1, ai = 1)
	}
	if(iter == 1){timelim = Inf}else{timelim = syst[8,1]*1.2}
	syst[8,iter] = system.time(fit[[8]] <- tryCatch({
  		setTimeLimit(elapsed = timelim)  # 60 seconds
		bugs(data,inits,model.file = "laplace.txt",
		parameters = c("intercept", "beta" , "sigma",
		"tau", "lam2", "ai"),
		n.chains = 1, n.iter = 10000,n.burnin=5000,
		n.thin=10,debug = F)
	}, error = function(e) {
  		if (grepl("time limit", e$message)) {
    		message(paste("Execution timed out after", timelim, "seconds"))
    		return(NULL)
  		} else {
    		return(NULL)
  		}
	}, finally = {
  		setTimeLimit(elapsed = Inf)  # Reset time limit
	}))[[3]]
	pars = tryCatch({fit[[8]]$summary[,1]}, error = function(e) {return(NULL)})
	parsims[[3]] = tryCatch({fit[[8]]$sims.matrix[(1:5000) %% 10 == 1, ]}, error = function(e) {return(NULL)})
	betahat[[8]] = tryCatch({colMeans(parsims[[3]][,1:10])}, error = function(e) {return(NULL)})
	cv3foldmcmc <- function(gam){
		nts=trunc(n/3)
		CVf=c()
		for(f in 1:3){
			cat("fold",f,"\n")
			tests=nts*(f-1)+1:nts
			ts.x=XM[tests,]
			ts.y=y[tests]
			tr.y=y[-tests]
			tr.x=XM[-tests,]
			data = list(n=n-nts,p=pM,y=tr.y,X=tr.x,p0=p0)
			fit <- tryCatch({
					bugs(data,inits,model.file = "laplace.txt",
					parameters = c("intercept", "beta" , "sigma",
					"tau", "lam2", "ai"),
					n.chains = 1, n.iter = 5000,n.burnin=2500,
					,n.thin=10,debug = F)},error=function(e){NULL})
			parsims = fit$sims.matrix[(1:2500) %% 10 == 1, ]
			betahat = colMeans(parsims[,1:5])
			betahat[-1] = tryCatch({betahat[-1] * (abs(betahat[-1]) > gam)},error=function(e){NULL})
			yhat = tryCatch({trunc(exp(cbind(1,ts.x) %*% betahat))},error=function(e){NULL})
			CVf = c(CVf,tryCatch({mean(yhat - ts.y)^2},error=function(e){NULL}))
		}
		mean(CVf)
	}
	a = tryCatch({c(0,abs(betahat[[8]][-1]) - 1e-5)}, error = function(e) {return(NULL)})
	b = tryCatch({c(abs(betahat[[8]][-1]) + 1e-5,0)}, error = function(e) {return(NULL)})
	gamseq = tryCatch({sort((a + b)/2)}, error = function(e) {return(NULL)})
	if(riskmethod == "cv3fold"){
		risks = tryCatch({sapply(gamseq,cv3foldmcmc)}, error = function(e) {return(NULL)})
	} else {
		risks = tryCatch({sapply(gamseq,AICf,8)}, error = function(e) {return(NULL)})
	}
	gamopt = tryCatch({gamseq[which.min(risks[risks>0])]}, error = function(e) {return(NULL)})
	pst[[8]] = tryCatch({c(1,1 * (abs(betahat[[8]][-1]) > gamopt))}, error = function(e) {return(NULL)})
	betahat[[8]] = tryCatch({betahat[[8]] * pst[[8]]},error=function(e){NULL})
	yhat[[8]] = tryCatch({trunc(exp(X %*% betahat[[8]]))}, error = function(e) {return(NULL)})
	yhattest[[8]] = tryCatch({trunc(exp(Xtest  %*% betahat[[8]]))}, error = function(e) {return(NULL)})	
	#------------------------------------
	#------------------------------------
	#------------------------------------
	#------------------------------------
	#------------------------------------
	#------------------------------------
	mcmcnames = c("Bernoulli", "CS", "Laplace")
	#------------------------------------
	#------------------------------------
	#----compute accuracies of betas-----
	#------------------------------------
	#------------------------------------
	cat("accuracies of betas \n")
	accuracy = matrix(NA,3,pM+1)
	colnames(accuracy) <- colnames(parsims[[1]])[1:(pM+1)]
	rownames(accuracy) <- mcmcnames
	accuracyl[[iter]] = accuracy
	if(length(parsims) == 3){
	for(j in 1:(pM+1)){
	if(iter == 1) pdf(file = paste0("dens",j,".pdf"),width=14,height=7)
	if(iter == 1) par(mfrow=c(1,3))
	limsx = limsy = matrix(0,2,3)
	for(h in 3:5){
		tt = seq(fit[[h]]$mu_beta[j] - 3*sqrt(fit[[h]]$sigma_beta[j,j]),
		fit[[h]]$mu_beta[j] + 3*sqrt(fit[[h]]$sigma_beta[j,j]),length.out=100)
		ss = dnorm(tt,fit[[h]]$mu_beta[j],sqrt(fit[[h]]$sigma_beta[j,j]))
		limsx[1,h-2] = fit[[h]]$mu_beta[j]-1.5*sqrt(fit[[h]]$sigma_beta[j,j])
		limsx[2,h-2] = fit[[h]]$mu_beta[j]+1.5*sqrt(fit[[h]]$sigma_beta[j,j])
		limsy[1,h-2] = 0
		limsy[2,h-2] = max(ss)
	}
	xlims = c(min(tt,min(limsx[1,])),max(tt,max(limsx[2,])))
	tt = seq(xlims[1],xlims[2],length.out=100)
	col = c("red","green","blue")
	for(h in 3:5){
		ylims = c(0,max(density(parsims[[h-2]][,j])$y,min(limsy[2,])))
		if(iter == 1) plot(dd <- density(parsims[[h-2]][,j]),main = paste0(mcmcnames[h-2],"-prior posterior-density estimate for ",colnames(parsims[[h-2]])[j],"=",round(realbeta[j],3)),
		xlim = xlims,
		ylim = ylims,xlab = "parameter",ylab = "posterior estimate")
		qfunc = function(x){ 
			a = dnorm(x,fit[[h]]$mu_beta[j],sqrt(fit[[h]]$sigma_beta[j,j]))
			if(h == 3){
				if(pst[[3]][j] == 1){
					return(a)
				} else {
					if(x == 0){
						return(1)
					} else {
						return(0)
					}
				}
			} else {
				return(a)
			}
		}
		post = kdensity(parsims[[h-2]][,j])$fhat
		ss = sapply(tt,qfunc)
		if(iter == 1) lines(ss~tt,col="red")
		if(iter == 1){
			if(sum(ss) == 1 | sum(ss) == 0) lines(c(0,0),c(0,max(dd$y)),lwd=3,col="red")
			legend(xlims[1],ylims[2],lty = rep(1,4),legend = c("MCMC","VB"),col=c("black","red"))
		}
		distt = function(theta) sapply(theta,function(th){
		abs(post(th)-qfunc(th))})
		accuracy[h-2,j] = tryCatch({100*(1-0.5*integrate(distt,-Inf,Inf,
			subdivisions=2000)$value)},error=function(e){NA})
	}
	if(iter == 1) dev.off()
	}
	accuracyl[[iter]] = round(accuracy,2)
	}
	#------------------------------------
	#------------------------------------
	#-compute accuracies of Z and gamma--
	#------------------------------------
	#------------------------------------
	cat("accuracies of Z and gamma \n")
	accuracy2 = matrix(NA,2,pM)
	colnames(accuracy2) <- colnames(parsims[[1]])[(pM+2):(2*pM+1)]
	rownames(accuracy2) <- mcmcnames[1:2]
	accuracy2l[[iter]] = accuracy2
	if(length(parsims) == 3){
	cntr = 0
	for(j in (pM+2):(2*pM+1)){
	for(h in 3:(1+length(parsims))){
		cntr = cntr + 1
		if(iter == 1) pdf(file = paste0("bar",cntr,".pdf"),width=7,height=7)
		tb = table(parsims[[h-2]][,j])/sum(table(parsims[[h-2]][,j]))
		if(length(tb) == 1){
			if(names(tb) == "1") tb = c(0,1) else tb = c(1,0)
		}
		number <- c(dbinom(0:1,1,fit[[h]]$pst[j-p]),tb)
		type <- c(rep("VB",2),rep("MCMC",2))
		circle <- data.frame(number,type)
		if(iter == 1) print(ggplot(circle, aes(factor(rep(0:1,2)),number, fill = type)) +
	  		geom_bar(stat="identity", position = "dodge") +
	  		labs(title=paste("Barplots for Bernulli posterior for",colnames(parsims[[h-2]])[j],"using prior",mcmcnames[h-2])) +
	  		xlab("values")+
	  		ylab("probability"))
		if(iter == 1) dev.off()
		accuracy2[h-2,j-10] = 100*(1-0.5*(abs(number[1] - number[3]) +abs(number[2] - number[4])))
	}
	}
	accuracy2l[[iter]] = round(accuracy2,2)
	}
	#------------------------------------
	#------------------------------------
	#------compute accuracies of pi------
	#------------------------------------
	#------------------------------------
	cat("accuracies of pi \n")
	accuracy3 = matrix(NA,2,pM)
	colnames(accuracy3) <- colnames(parsims[[1]])[(2*pM+2):(3*pM+1)]
	rownames(accuracy3) <- mcmcnames[1:2]
	accuracy3l[[iter]] = accuracy3
	if(length(parsims) == 3){
	for(j in (2*pM+2):(3*pM+1)){
	if(iter == 1) pdf(file = paste0("dens",j,".pdf"),width=14,height=7)
	if(iter == 1) par(mfrow=c(1,2))
	limsx = limsy = matrix(0,2,2)
	for(h in 3:4){
		tt = seq(0,1,length.out=100)
		ss = dbeta(tt,fit[[h]]$alphast[j-2*pM],fit[[h]]$betast[j-2*pM])
		limsy[1,h-2] = 0
		limsy[2,h-2] = max(ss)
	}
	xlims = c(0,1)
	for(h in 3:4){
		ylims = c(0,max(density(parsims[[h-2]][,j])$y,max(limsy[2,][is.finite(limsy[2,])])))
		if(iter == 1) plot(density(parsims[[h-2]][,j]),main = paste0(mcmcnames[h-2],"-prior posterior-density estimate for ",colnames(parsims[[h-2]])[j]),
		xlim = c(0,1),
		ylim = ylims,xlab = "parameter",ylab = "posterior estimate")
		ss = dbeta(tt,fit[[h]]$alphast[j-2*pM-1],fit[[h]]$betast[j-2*pM-1])
		if(iter == 1) lines(ss~tt,col="red")
		if(iter == 1) legend(xlims[1],ylims[2],lty = rep(1,2),legend = c("MCMC","VB"),col=c("black","red"))
		post = kdensity(parsims[[h-2]][,j])$fhat
		qfunc = function(x) dbeta(x,fit[[h]]$alphast[j-2*pM],fit[[h]]$betast[j-2*pM])
		distt = function(theta) sapply(theta,function(th){
		abs(post(th)-qfunc(th))})
		accuracy3[h-2,j-2*pM-1] = tryCatch({100*(1-0.5*
		integrate(distt,0,1)$value)},error=function(e){NA})
	}
	if(iter == 1) dev.off()
	}
	accuracy3l[[iter]] = round(accuracy3,2)
	}
	#------------------------------------
	#------------------------------------
	#--compute accuracies of variances---
	#------------------------------------
	#------------------------------------
	cat("accuracies of vars \n")
	dinvgamma <- function(x,alpha,beta) beta^alpha / gamma(alpha) * 
	(1/x)^(alpha+1) * exp(-beta/x) 
	dgeninvgauss <- function(x,a,b) (a/b)^0.25 / BesselK(0.5, sqrt(a*b))/2 * 
	x ^(-0.5) * exp(-(a*x + b/x)/2)
	#------------------------------------
	accuracy4 = matrix(NA,1,pM+1)
	colnames(accuracy4) <- colnames(parsims[[1]])[(3*pM+2):(4*pM+2)]
	rownames(accuracy4) <- mcmcnames[1]
	accuracy4l[[iter]] = accuracy4
	if(length(parsims) == 3){
	for(j in (3*pM+2):(4*pM+2)){
	if(iter == 1) pdf(file = paste0("dens",j,".pdf"),width=14,height=7)
	if(iter == 1) par(mfrow=c(1,1))
	limsx = limsy = matrix(0,2,1)
	h=3
	dd = density(1/parsims[[h-2]][,j])
	dd$y = dd$y/max(dd$y)
	tt = seq(0.00001,mean(dd$x),length.out=100)
	ss = dinvgamma(tt,fit[[h]]$sigmapars$a_beta[j-3*pM-1],fit[[h]]$sigmapars$b_beta[j-3*pM-1])
	limsy[1,h-2] = 0
	limsy[2,h-2] = max(ss)
	#
	xlims = c(0,3)
	h=3
	ylims = c(0,1)
	if(iter == 1) plot(dd,main = paste0(mcmcnames[h-2],"-prior (normalized) posterior-density estimate for ",colnames(parsims[[h-2]])[j]),
	xlim = c(0,mean(dd$x)),
	ylim = ylims,xlab = "parameter",ylab = "posterior estimate")
	ss = dinvgamma(tt,fit[[h]]$sigmapars$a_beta[j-3*pM-1],fit[[h]]$sigmapars$b_beta[j-3*pM-1])
	ss = ss/max(ss)
	if(iter == 1) lines(ss~tt,col="red")
	if(iter == 1) legend(xlims[1],ylims[2],lty = rep(1,2),legend = c("MCMC","VB"),col=c("black","red"))
	post = kdensity(1/parsims[[h-2]][,j])$fhat
	qfunc = function(x) dinvgamma(x,fit[[h]]$sigmapars$a_beta[j-3*pM-1],fit[[h]]$sigmapars$b_beta[j-3*pM-1])
	distt = function(theta) sapply(theta,function(th){
	abs(post(th)-qfunc(th))})
	accuracy4[h-2,j-3*pM-1] = tryCatch({100*(1-0.5*
		integrate(distt,0,Inf,subdivisions=2000,
		rel.tol =1e-3)$value)},error=function(e){NA})
	if(iter == 1) dev.off()
	}
	accuracy4l[[iter]] = round(accuracy4,2)
	}
	#------------------------------------
	#------------------------------------
	#--compute accuracies of variances---
	#------------------------------------
	#------------------------------------
	accuracy5 = matrix(NA,1,1)
	colnames(accuracy5) <- "sigma"
	rownames(accuracy5) <- mcmcnames[2]
	accuracy5l[[iter]] = accuracy5
	j=3*pM+2
	if(iter == 1) pdf(file = paste0("dens",19,".pdf"),width=14,height=7)
	if(iter == 1) par(mfrow=c(1,1))
	limsx = limsy = matrix(0,2,1)
	h=4
	if(length(parsims) == 3){
	dd = density(1/parsims[[h-2]][,j])
	dd$y = dd$y/max(dd$y)
	tt = seq(0.00001,mean(dd$x),length.out=100)
	ss = dinvgamma(tt,fit[[h]]$sigmapars$alpha_sigma2,fit[[h]]$sigmapars$beta_sigma2)
	limsy[1,h-3] = 0
	limsy[2,h-3] = max(ss)
	#
	xlims = c(0,3)
	h=4
	ylims = c(0,1)
	if(iter == 1) plot(dd,main = paste0(mcmcnames[h-2],"-prior (normalized) posterior-density estimate for ",colnames(parsims[[h-2]])[j]),
	xlim = c(0,mean(dd$x)),
	ylim = ylims,xlab = "parameter",ylab = "posterior estimate")
	ss = dinvgamma(tt,fit[[h]]$sigmapars$alpha_sigma2,fit[[h]]$sigmapars$beta_sigma2)
	ss = ss/max(ss)
	if(iter == 1) lines(ss~tt,col="red")
	if(iter == 1) legend(xlims[1],ylims[2],lty = rep(1,2),legend = c("MCMC","VB"),col=c("black","red"))
	post = kdensity(1/parsims[[h-2]][,j])$fhat
	qfunc = function(x) dinvgamma(x,fit[[h]]$sigmapars$alpha_sigma2,fit[[h]]$sigmapars$beta_sigma2)
	distt = function(theta) sapply(theta,function(th){
	abs(post(th)-qfunc(th))})
	accuracy5[h-3,j-3*p-1] = tryCatch({100*(1-0.5*
		integrate(distt,0,Inf,
		subdivisions=2000)$value)},error=function(e){NA})
	if(iter == 1) dev.off()
	#
	accuracy5l[[iter]] = round(accuracy5,2)
	}
	#------------------------------------
	#------------------------------------
	#--compute accuracies of variances---
	#------------------------------------
	#------------------------------------
	accuracy6 = matrix(NA,1,pM+1)
	colnames(accuracy6) <- paste0("tau[",0:pM,"]")
	rownames(accuracy6) <- mcmcnames[3]
	accuracy6l[[iter]] = accuracy6
	if(length(parsims) == 3){
	for(j in (pM+2):(2*pM+2)){
	if(iter == 1) pdf(file = paste0("dens",j+14,".pdf"),width=14,height=7)
	if(iter == 1) par(mfrow=c(1,1))
	limsx = limsy = matrix(0,2,1)
	h=5
	if(j == (pM+2)){
		dd = density(1/parsims[[h-2]][,j])
	}else{
		dd = density(parsims[[h-2]][,j])
	}
	dd$y = dd$y/max(dd$y)
	tt = seq(0.00001,mean(dd$x),length.out=100)
	if(j == (pM+2)){
		ss = dinvgamma(tt,1,(fit[[h]]$sigmapars$Etau0m1)^(-1))
	} else {
		aa = fit[[h]]$sigmapars$astar
		bb = fit[[h]]$sigmapars$bstar[j-pM-2]
		ss = dgeninvgauss(tt,aa,bb)
	}
	limsy[1,h-4] = 0
	limsy[2,h-4] = max(ss)
	#
	xlims = c(0,3)
	h=5
	ylims = c(0,1)
	if(iter == 1) plot(dd,main = paste0(mcmcnames[h-2],"-prior (normalized) posterior-density estimate for ",colnames(parsims[[h-2]])[j]),
	xlim = c(0,mean(dd$x)),
	ylim = ylims,xlab = "parameter",ylab = "posterior estimate")
	if(j == (pM+2)){
		ss = dinvgamma(tt,1,(fit[[h]]$sigmapars$Etau0m1)^(-1))
	} else {
		aa = fit[[h]]$sigmapars$astar
		bb = fit[[h]]$sigmapars$bstar[j-pM-2]
		ss = dgeninvgauss(tt,aa,bb)
	}
	ss = ss/max(ss)
	if(iter == 1) lines(ss~tt,col="red")
	if(iter == 1) legend(xlims[1],ylims[2],lty = rep(1,2),legend = c("MCMC","VB"),col=c("black","red"))
	if(j == (pM+2)){
		post = kdensity(1/parsims[[h-2]][,j])$fhat
	}else{
		post = kdensity(parsims[[h-2]][,j])$fhat
	}
	qfunc = function(x){
	if(j == (pM+2)){
		ss = dinvgamma(x,1,(fit[[h]]$sigmapars$Etau0m1)^(-1))
	} else {
		aa = fit[[h]]$sigmapars$astar
		bb = fit[[h]]$sigmapars$bstar[j-pM-2]
		ss = dgeninvgauss(x,aa,bb)
	}
	ss
	}
	distt = function(theta) sapply(theta,function(th){
	abs(post(th)-qfunc(th))})
	accuracy6[h-4,j-pM-1] = tryCatch({100*(1-0.5*
	integrate(distt,0.000001,Inf,
	subdivisions=2000)$value)},error=function(e){NA})
	if(iter == 1) dev.off()
	}
	accuracy6l[[iter]] = round(accuracy6,2)
	}
	#------------------------------------
	#------------------------------------
	#
	cat("criteria \n")
	for(j in 1:8){
		seb[j,iter] = tryCatch({mean((betahat[[j]]-realbeta)^2)/mean(realbeta^2)},error=function(e){NA})
		#
		se[j,iter] = tryCatch({mean((y-yhat[[j]])^2)/var(y)},error=function(e){NA})
		#
		set[j,iter] = tryCatch({mean((ytest-yhattest[[j]])^2)/var(ytest)},error=function(e){NA})
		#
		err1 = tryCatch({sum(realp == 0 & pst[[j]] !=0)/sum(realp == 0)},error=function(e){NA})
		if(is.nan(err1)) err1 = 0
		err2 = tryCatch({sum(realp != 0 & pst[[j]] ==0)/sum(realp != 0)},error=function(e){NA})
		#
		FPR[j,iter] = err1
		#
		FNR[j,iter] = err2
	}
	cat("coverage probs \n")
	lower = upper = matrix(0,p,6)
	for(j in 1:p){
		for(h in 3:8){
			if(h<=5){
				tt = seq(fit[[h]]$mu_beta[j] - 3*sqrt(abs(fit[[h]]$sigma_beta[j,j])),
				fit[[h]]$mu_beta[j] + 3*sqrt(abs(fit[[h]]$sigma_beta[j,j])),length.out=100)
				ss = dnorm(tt,fit[[h]]$mu_beta[j],sqrt(abs(fit[[h]]$sigma_beta[j,j])))
				object = list(x = tt, y = ss)
			} else {
				object = tryCatch({density(parsims[[h-5]][,j])},error=function(e){NA})
			}
			hhpp <- tryCatch({hpd(object)},error=function(e){NA})
			lower[j,h-2] = tryCatch({hhpp[1]},error=function(e){NA})
			upper[j,h-2] = tryCatch({hhpp[2]},error=function(e){NA})
			cover[j,h-2] = tryCatch({(realbeta[j] >= lower[j,h-2])*(realbeta[j] <= upper[j,h-2])},error=function(e){NA})
		}
	}
	coverage[[iter]] = cover
}
#------------------------------------
load(file.choose())
#------------------------------------
#------------------------------------
systr = syst
for(i in 1:8){
	systr[i,] = syst[i,]/syst[7,]
}
#------------------------------------
#------------------------------------
#-----------trimmed means------------
#------------------------------------
#------------------------------------
round(apply(seb,1,function(x){mean(x,trim = 0.25, na.rm = T)}),5)
round(apply(se,1,function(x){mean(x,trim = 0.25, na.rm = T)}),5)
round(apply(set,1,function(x){mean(x,trim = 0.25, na.rm = T)}),5)
round(apply(FPR,1,function(x){mean(x,trim = 0.25, na.rm = T)}),5)
round(apply(FNR,1,function(x){mean(x,trim = 0.25, na.rm = T)}),5)
round(apply(systr,1,function(x){mean(x,trim = 0.25, na.rm = T)}),5)
#------------------------------------
#------------------------------------
#----------------Accuracies----------
#------------------------------------
#------------------------------------
boxplot_mat <- function(matrix_list, ymin, ymax){
	n_rows <- nrow(matrix_list[[1]])
	n_cols <- ncol(matrix_list[[1]])
	df <- do.call(rbind, lapply(1:length(matrix_list), function(iter) {
  		mat <- matrix_list[[iter]]
  		data.frame(
    		Iteration = iter,
    		Row = rep(1:n_rows, times = n_cols),
    		Col = rep(1:n_cols, each = n_rows),
			rownames = rownames(matrix_list[[1]]),
			colnames = colnames(matrix_list[[1]]),
    		Value = as.vector(mat)
  		)
	}))
	df$Position <- interaction(df$rownames[df$Row], df$colnames[df$Col], sep = "-")
	df$Position <- gsub("CS-gamma", "CS-Z", df$Position)
	# Create the boxplot
	ggplot(df, aes(x = Position, y = Value)) +
  		geom_boxplot() + ylim(ymin, ymax) +
  		theme_bw() + 
  		theme(axis.text.x = element_text(size = 12, angle = 90, hjust = 1)) + 
  		theme(axis.text.y = element_text(size = 12)) +
  		labs(x = "Method-Parameter", y = "Accuracy")
}
#------------------------------------
pdf(file = "box_acc1.pdf",width=14,height=7)
boxplot_mat(accuracyl,25,100)
dev.off()
#------------------------------------
pdf(file = "box_acc2.pdf",width=14,height=7)
boxplot_mat(accuracy2l,25,100)
dev.off()
#------------------------------------
pdf(file = "box_acc3.pdf",width=14,height=7)
boxplot_mat(accuracy3l,25,100)
dev.off()
#------------------------------------
pdf(file = "box_acc4.pdf",width=14,height=7)
boxplot_mat(accuracy4l,25,100)
dev.off()
#------------------------------------
pdf(file = "box_acc5.pdf",width=14,height=7)
boxplot_mat(accuracy5l,25,100)
dev.off()
#------------------------------------
accuracy6l <- Map(function(x){colnames(x)[1] <- "tau[0]";x},accuracy6l)
pdf(file = "box_acc6.pdf",width=14,height=7)
boxplot_mat(accuracy6l,25,100)
dev.off()
#------------------------------------
#------------------------------------
#------------------------------------
#------------------------------------
nna1 <- Reduce("+",Map(function(x){sum(is.na(x))},
accuracyl))/length(accuracyl[[1]])
accnna <- Map(function(x){x[is.na(x)] <- 0;x},
accuracyl)
round((A1 <- Reduce("+",accnna)/(Iterations-nna1)),2)
round(sqrt((Reduce("+",Map(function(x){(x-A1)^2},
accnna))-nna1*A1^2)/(Iterations-nna1)),2)
#------------------------------------
nna2 <- Reduce("+",Map(function(x){sum(is.na(x))},
accuracy2l))/length(accuracy2l[[1]])
accnna <- Map(function(x){x[is.na(x)] <- 0;x},
accuracy2l)
round((A2 <- Reduce("+",accnna)/(Iterations-nna2)),2)
round(sqrt((Reduce("+",Map(function(x){(x-A2)^2},
accnna))-nna2*A2^2)/(Iterations-nna2)),2)
#------------------------------------
nna3 <- Reduce("+",Map(function(x){sum(is.na(x))},
accuracy3l))/length(accuracy3l[[1]])
accnna <- Map(function(x){x[is.na(x)] <- 0;x},
accuracy3l)
round((A3 <- Reduce("+",accnna)/(Iterations-nna3)),2)
round(sqrt((Reduce("+",Map(function(x){(x-A3)^2},
accnna))-nna3*A3^2)/(Iterations-nna3)),2)
#------------------------------------
nna4 <- Reduce("+",Map(function(x){sum(is.na(x))},
accuracy4l))/length(accuracy4l[[1]])
accnna <- Map(function(x){x[is.na(x)] <- 0;x},
accuracy4l)
round((A4 <- Reduce("+",accnna)/(Iterations-nna4)),2)
round(sqrt((Reduce("+",Map(function(x){(x-A4)^2},
accnna))-nna4*A4^2)/(Iterations-nna4)),2)
#------------------------------------
nna5 <- Reduce("+",Map(function(x){sum(is.na(x))},
accuracy5l))/length(accuracy5l[[1]])
accnna <- Map(function(x){x[is.na(x)] <- 0;x},
accuracy5l)
round((A5 <- Reduce("+",accnna)/(Iterations-nna5)),2)
round(sqrt((Reduce("+",Map(function(x){(x-A5)^2},
accnna))-nna5*A5^2)/(Iterations-nna5)),2)
#------------------------------------
nna6 <- Reduce("+",Map(function(x){sum(is.na(x))},
accuracy6l))/length(accuracy6l[[1]])
accnna <- Map(function(x){x[is.na(x)] <- 0;x},
accuracy6l)
round((A6 <- Reduce("+",accnna)/(Iterations-nna6)),2)
round(sqrt((Reduce("+",Map(function(x){(x-A6)^2},
accnna))-nna6*A6^2)/(Iterations-nna6)),2)
#------------------------------------
#------------------------------------
#------------------------------------
#--coverage probabilities for HPD----
#------------------------------------
#------------------------------------
nna <- Reduce("+",Map(function(x){is.na(x)},
coverage))
accnna <- Map(function(x){x[is.na(x)] <- 0;x},
coverage)
round((CM <- Reduce("+",accnna)/(Iterations-nna)),2)
round(sqrt((Reduce("+",Map(function(x){(x-CM)^2},
accnna))-nna*CM^2)/(Iterations-nna)),2)
#------------------------------------
CL = Map(colMeans,coverage)
covers = do.call(rbind, CL)
pdf(file = "box17r.pdf",width=14,height=7)
dseb = stack(as.data.frame(covers))
names(dseb)[2] <- "Method"
myplot <- ggplot(dseb,aes(x = Method, y = values, color = Method)) + 
geom_boxplot(fatten = 2, notch=FALSE,outlier.shape = 16) +
scale_y_continuous(limits = quantile(dseb$values, c(0.1, 0.9), na.rm = T)) + 
ylim(0.5, 1) + xlab('') + 
ylab('Coverage Probabilities')
suppressWarnings(print(myplot))
dev.off()
#------------------------------------
#------------------------------------
#------------------------------------
#---------------Box plots------------
#------------------------------------
#------------------------------------
pdf(file = "box11r.pdf",width=14,height=7)
dseb = stack(as.data.frame(t(seb[,1:(iter-1)])))
names(dseb)[2] <- "Method"
myplot <- ggplot(dseb,aes(x = Method, y = values, color = Method)) + 
geom_boxplot(fatten = 2, notch=FALSE,outlier.shape = 16) +
scale_y_continuous(limits = quantile(dseb$values, c(0.1, 0.9), na.rm = T)) + 
ylim(0.001, 0.06) + xlab('') +
theme_bw() + 
theme(legend.position="none") + 
theme(axis.text.x = element_text(size = 12)) + 
theme(axis.text.y = element_text(size = 12)) +
theme(axis.title.x = element_text(size = 12)) + 
theme(axis.title.y = element_text(size = 12)) + 
ylab('Coefficient estimate relative error (CRE)')
suppressWarnings(print(myplot))
dev.off()
#------------------------------------
pdf(file = "box12r.pdf",width=14,height=7)
dseb = stack(as.data.frame(t(se[,1:(iter-1)])))
names(dseb)[2] <- "Method"
myplot <- ggplot(dseb,aes(x = Method, y = values, color = Method)) + 
geom_boxplot(fatten = 2, notch=FALSE,outlier.shape = 16) +
scale_y_continuous(limits = quantile(dseb$values, c(0.1, 0.9), na.rm = T)) + 
ylim(0.0, 0.3) + xlab('') +
theme_bw() + 
theme(legend.position="none") + 
theme(axis.text.x = element_text(size = 12)) + 
theme(axis.text.y = element_text(size = 12)) +
theme(axis.title.x = element_text(size = 12)) + 
theme(axis.title.y = element_text(size = 12)) + 
ylab('Train relative error (TRRE)')
suppressWarnings(print(myplot))
dev.off()
#------------------------------------
pdf(file = "box13r.pdf",width=14,height=7)
dseb = stack(as.data.frame(t(set[,1:(iter-1)])))
names(dseb)[2] <- "Method"
myplot <- ggplot(dseb,aes(x = Method, y = values, color = Method)) + 
geom_boxplot(fatten = 2, notch=FALSE,outlier.shape = 16) +
scale_y_continuous(limits = quantile(dseb$values, c(0.1, 0.9), na.rm = T)) + 
ylim(0.0, 0.5) + xlab('') +
theme_bw() + 
theme(legend.position="none") + 
theme(axis.text.x = element_text(size = 12)) + 
theme(axis.text.y = element_text(size = 12)) +
theme(axis.title.x = element_text(size = 12)) + 
theme(axis.title.y = element_text(size = 12)) + 
ylab('Test relative error (TSRE)')
suppressWarnings(print(myplot))
dev.off()
#------------------------------------
pdf(file = "box14r.pdf",width=14,height=7)
dseb = stack(as.data.frame(t(FPR[,1:(iter-1)])))
names(dseb)[2] <- "Method"
myplot <- ggplot(dseb,aes(x = Method, y = values, color = Method)) + 
geom_boxplot(fatten = 2, notch=FALSE,outlier.shape = 16) +
scale_y_continuous(limits = quantile(dseb$values, c(0.1, 0.9), na.rm = T)) +
ylim(0, 0.9) + xlab('') + 
theme_bw() + 
theme(legend.position="none") + 
theme(axis.text.x = element_text(size = 12)) + 
theme(axis.text.y = element_text(size = 12)) +
theme(axis.title.x = element_text(size = 12)) + 
theme(axis.title.y = element_text(size = 12)) + 
ylab('False positive rate (FPR)')
suppressWarnings(print(myplot))
dev.off()
#------------------------------------
pdf(file = "box15r.pdf",width=14,height=7)
dseb = stack(as.data.frame(t(FNR[,1:(iter-1)])))
names(dseb)[2] <- "Method"
myplot <- ggplot(dseb,aes(x = Method, y = values, color = Method)) + 
geom_boxplot(fatten = 2, notch=FALSE,outlier.shape = 16) +
#scale_y_continuous(limits = quantile(dseb$values, c(0.1, 0.995), na.rm = T)) +
ylim(0, 0.75) + xlab('') + 
theme_bw() + 
theme(legend.position="none") + 
theme(axis.text.x = element_text(size = 12)) + 
theme(axis.text.y = element_text(size = 12)) +
theme(axis.title.x = element_text(size = 12)) + 
theme(axis.title.y = element_text(size = 12)) + 
ylab('False negative rate (FNR)')
suppressWarnings(print(myplot))
dev.off()
#------------------------------------
pdf(file = "box16r.pdf",width=14,height=7)
dseb = stack(as.data.frame(t(systr[1:5,1:(iter-1)])))
names(dseb)[2] <- "Method"
myplot <- ggplot(dseb,aes(x = Method, y = values, color = Method)) + 
geom_boxplot(fatten = 2, notch=FALSE,outlier.shape = NA) +
scale_y_continuous(limits = quantile(dseb$values, c(0.1, 0.9), na.rm = T)) +
ylim(0, 0.002) + xlab('') + 
theme_bw() + 
theme(legend.position="none") + 
theme(axis.text.x = element_text(size = 12)) + 
theme(axis.text.y = element_text(size = 12)) +
theme(axis.title.x = element_text(size = 12)) + 
theme(axis.title.y = element_text(size = 12)) + 
ylab('Computation time relative to CS-MCMC')
suppressWarnings(print(myplot))
dev.off()
#------------------------------------



