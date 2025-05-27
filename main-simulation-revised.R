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
}
#------------------------------------
Iterations = 100
seb = se = set = FPR = FNR = syst = matrix(0,5,Iterations)
rownames(seb) <- rownames(se) <-
rownames(set) <- rownames(FPR) <-
rownames(FNR) <- rownames(syst) <- c("LASSO","SCAD","Bernoulli-VB","CS-VB",
"LAPLACE-VB")
#------------------------------------
fit = betahat = pst = yhat = yhattest =  list()
p = 5
cover = matrix(0,p,3)
for(iter in 1:Iterations){
	n = 100
	sam = generate_samples(n=n, p=p, p0=0.3)
	samp = sample(1:n,trunc(0.8*n))
	y  <- sam$y[samp]
	X = sam$X[samp,]
	ytest = sam$y[-samp]
	Xtest = sam$X[-samp,]
	realbeta = sam$beta
	realp = 1*(realbeta !=0)
	n <- length(y)
	#------------------------------------
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
	Sigma_beta = diag(rep(1e-5,p))
	pst0 = pst[[1]]
		###c(1,rep(1,p-1))
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
	###pst[[3]] = tryCatch({fit[[3]]$pst},error=function(e){NULL})
	pst[[3]] = tryCatch({1*(fit[[3]]$pst > 0.5)},error=function(e){NULL})
	###pst[[3]] = sapply(1:p, function(j){hpdtest(fit[[3]],j)})
	yhat[[3]] = tryCatch({exp(X %*% betahat[[3]])},error=function(e){NULL})
	yhattest[[3]] = tryCatch({predict(fit[[3]],Xtest,method ="Bernoulli")},error=function(e){NULL})
	#------------------------------------
	c = 1e-2
	syst[4,iter] = system.time(fit[[4]] <- tryCatch({
	sppoissregvb(X,y,init,prior="CS", 
		eps = 1e-12, maxiter = 100)},error=function(e){NULL}))[[3]]
	####pst[[4]] = tryCatch({1*((fit[[4]]$pst)/sqrt(fit[[4]]$pst*(1-fit[[4]]$pst))>qnorm(1-0.05/p))},error=function(e){NULL})
	####pst[[4]] = sapply(1:p, function(j){hpdtest(fit[[4]],j)})
	###pst[[4]] = tryCatch({abs(fit[[4]]$mu_beta)/sqrt(diag(fit[[4]]$sigma_beta)) > qnorm(1-0.025/p)},error=function(e){NULL})
	p0 = max(0.01,min(sum(pst[[1]])/p,0.99))
	pst[[4]] = 1 * ((fit[[4]]$mu_beta)^2>(-2*c/(1+c)*(0.5*log(c)+log(p0/(1-p0)))*fit[[4]]$sigmapars$beta_sigma2/(fit[[4]]$sigmapars$alpha_sigma2 - 1)))
	betahat[[4]] = tryCatch({fit[[4]]$mu_beta * pst[[4]]},error=function(e){NULL})
	yhat[[4]] = tryCatch({exp(X %*% betahat[[4]])},error=function(e){NULL})
	yhattest[[4]] = tryCatch({predict(fit[[4]],Xtest,method ="CS")},error=function(e){NULL})
	#------------------------------------
	syst[5,iter] = system.time(fit[[5]] <- tryCatch({
	sppoissregvb(X,y,init,prior="Laplace", 
		eps = 1e-5, maxiter = 100)},error=function(e){NULL}))[[3]]
	#####pst[[5]] = rep(1,p)
	####pst[[5]] = tryCatch({abs(fit[[5]]$mu_beta)/sqrt(diag(fit[[5]]$sigma_beta)) > qnorm(1-0.025/p)},error=function(e){NULL})
	####pst[[5]] = sapply(1:p, function(j){hpdtest(fit[[5]],j)})
	risk <- function(gam){
		###yhat = exp(X %*% fit[[5]]$mu_beta)
		###tauhat <- mean((yhat - y)^2)/n
		tauhat <- diag(fit[[5]]$sigma_beta)
		####tauhat <- 1
		risk <- mean(tauhat) * (p - 2 * 
			sum(abs(fit[[5]]$mu_beta)/sqrt(tauhat)<gam)+
			sum((pmax(abs(fit[[5]]$mu_beta)/sqrt(tauhat),gam))^2))
		risk
	}
	a = c(0,fit[[5]]$mu_beta - 1e-5)
	b = c(fit[[5]]$mu_beta + 1e-5,0)
	gamseq = sort(c((a + b)/2,fit[[5]]$mu_beta))
	risks = sapply(gamseq,risk)
	gamopt = gamseq[which.min(risks[risks>0])]
	pst[[5]] = 1 * (abs(fit[[5]]$mu_beta)> gamopt)
	pst[[5]][1] = 1
	betahat[[5]] = tryCatch({fit[[5]]$mu_beta * pst[[5]]},error=function(e){NULL})
	yhat[[5]] = tryCatch({exp(X %*% betahat[[5]])},error=function(e){NULL})
	yhattest[[5]] = tryCatch({predict(fit[[5]],Xtest,method ="Laplace")},error=function(e){NULL})
	#------------------------------------
	#
	for(j in 1:5){
		seb[j,iter] = tryCatch({mean((betahat[[j]]-realbeta)^2)/mean(realbeta^2)},error=function(e){NA})
		#
		se[j,iter] = tryCatch({mean((y-yhat[[j]])^2)/var(y)},error=function(e){NA})
		#
		set[j,iter] = tryCatch({mean((ytest-yhattest[[j]])^2)/var(ytest)},error=function(e){NA})
		#
		err1 = sum(realp == 0 & pst[[j]] !=0)/sum(realp == 0)
		if(is.nan(err1)) err1 = 0
		err2 = sum(realp != 0 & pst[[j]] ==0)/sum(realp != 0)
		#
		FPR[j,iter] = err1
		#
		FNR[j,iter] = err2
	}
	lower = upper = matrix(0,p,3)
	for(j in 1:p){
		for(h in 3:5){
			tt = seq(fit[[h]]$mu_beta[j] - 3*sqrt(abs(fit[[h]]$sigma_beta[j,j])),
			fit[[h]]$mu_beta[j] + 3*sqrt(abs(fit[[h]]$sigma_beta[j,j])),length.out=100)
			ss = dnorm(tt,fit[[h]]$mu_beta[j],sqrt(abs(fit[[h]]$sigma_beta[j,j])))
			object = list(x = tt, y = ss)
			hhpp <- hpd(object)
			lower[j,h-2] = hhpp[1]
			upper[j,h-2] = hhpp[2]
			cover[j,h-2] = cover[j,h-2] + 
				(realbeta[j] >= lower[j,h-2])*(realbeta[j] <= upper[j,h-2])/Iterations 
		}
	}
	#-----------------------------------------
	cat("Iteration",iter,"\n")
}
#------------------------------------
#------------------------------------
#-----------trimmed means------------
#------------------------------------
#------------------------------------
round(apply(seb,1,function(x){mean(x,trim = 0.25)}),3)
round(apply(se,1,function(x){mean(x,trim = 0.25)}),3)
round(apply(set,1,function(x){mean(x,trim = 0.25)}),3)
round(apply(FPR,1,function(x){mean(x,trim = 0.25)}),3)
round(apply(FNR,1,function(x){mean(x,trim = 0.25)}),3)
round(apply(syst,1,function(x){mean(x,trim = 0.25)}),3)
#------------------------------------
#------------------------------------
#--coverage probabilities for HPD----
#------------------------------------
#------------------------------------
cover
#------------------------------------
#------------------------------------
#----------MCMC samplings------------
#--------for last iteration----------
#------------------------------------
#------------------------------------
X = X[,-1]
p = p - 1
data = list("n","p","y","X")
inits <- function(){
	list(intercept = 0, 
	beta = rep(0.1,p), gamma = rep(1,p),
	pi = rep(0.5,p), sigma= 1)
}
#syst6 = system.time(fit[[6]] <- 
#	bugs(data,inits,model.file = "bernulli.txt",
#	parameters = c("intercept", "beta" , "gamma","pi", "sigma"),
#	n.chains = 1, n.iter = 10000,n.burnin=5000,
#	n.thin=10,debug = T))[[3]]
#	pars = fit[[6]]$summary[,1]
#	pst[[6]] = c(1,fit[[6]]$summary[(p+2):(2*p+1),3])
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
mod<-jags.model(textConnection(mod_string),data=data,
inits=inits,n.chain=1)
library(coda)
update(mod,1000)
mod_sim1<-coda.samples(model = 
mod,variable.names = params,n.iter = 5000)
#------------------------------------
data = list(n=n,p=p,y=y,X=X,p0=p0)
inits <- function(){
	list(intercept = init$mu_beta[1], 
	beta = init$mu_beta[-1], z = init$pst[-1],
	pi = init$pst[-1], sigma = 1, ai = 1,
	eta1 = init$mu_beta[-1], eta2 = init$mu_beta[-1])
}
syst7 = system.time(fit[[7]] <- 
	bugs(data,inits,model.file = "cs.txt",
	parameters = c("intercept", "beta" , "z", "pi", 
	"sigma","ai","eta1","eta2"),
	n.chains = 1, n.iter = 10000,n.burnin=5000,
	,n.thin=10,debug = T))[[3]]
	pars = fit[[7]]$summary[,1]
	pst[[7]] = c(1,fit[[7]]$summary[(p+2):(2*p+1),3])
#------------------------------------
inits <- function(){
	list(intercept = init$mu_beta[1], 
	beta = init$mu_beta[-1], sigma = 1, 
	tau = rep(1,p), 
	lam2 = 1, ai = 1)
}
syst8 = system.time(fit[[8]] <- 
	bugs(data,inits,model.file = "laplace.txt",
	parameters = c("intercept", "beta" , "sigma",
	"tau", "lam2", "ai"),
	n.chains = 1, n.iter = 10000,n.burnin=5000,
	,n.thin=10,debug = T))[[3]]
	pars = fit[[8]]$summary[,1]
	pst[[8]] = c(1,fit[[8]]$summary[(p+2):(2*p+1),3])
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
ylim(0.1, 2) + xlab('') + 
ylab('Coefficient estimate relative error')
suppressWarnings(print(myplot))
dev.off()
#------------------------------------
pdf(file = "box12r.pdf",width=14,height=7)
dseb = stack(as.data.frame(t(se[,1:(iter-1)])))
names(dseb)[2] <- "Method"
myplot <- ggplot(dseb,aes(x = Method, y = values, color = Method)) + 
geom_boxplot(fatten = 2, notch=FALSE,outlier.shape = 16) +
scale_y_continuous(limits = quantile(dseb$values, c(0.1, 0.9), na.rm = T)) + 
ylim(0.0, 0.5) + xlab('') +
ylab('Train relative error')
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
ylab('Test relative error')
suppressWarnings(print(myplot))
dev.off()
#------------------------------------
pdf(file = "box14r.pdf",width=14,height=7)
dseb = stack(as.data.frame(t(FPR[,1:(iter-1)])))
names(dseb)[2] <- "Method"
myplot <- ggplot(dseb,aes(x = Method, y = values, color = Method)) + 
geom_boxplot(fatten = 2, notch=FALSE,outlier.shape = 16) +
scale_y_continuous(limits = quantile(dseb$values, c(0.1, 0.9), na.rm = T)) +
ylim(0, 1) + xlab('') + ylab('False positive rate')
suppressWarnings(print(myplot))
dev.off()
#------------------------------------
pdf(file = "box15r.pdf",width=14,height=7)
dseb = stack(as.data.frame(t(FNR[,1:(iter-1)])))
names(dseb)[2] <- "Method"
myplot <- ggplot(dseb,aes(x = Method, y = values, color = Method)) + 
geom_boxplot(fatten = 2, notch=FALSE,outlier.shape = 16) +
#scale_y_continuous(limits = quantile(dseb$values, c(0.1, 0.995), na.rm = T)) +
ylim(0, 1) + xlab('') + ylab('False negative rate')
suppressWarnings(print(myplot))
dev.off()
#------------------------------------
pdf(file = "box16r.pdf",width=14,height=7)
dseb = stack(as.data.frame(t(syst[,1:(iter-1)])))
names(dseb)[2] <- "Method"
myplot <- ggplot(dseb,aes(x = Method, y = values, color = Method)) + 
geom_boxplot(fatten = 2, notch=FALSE,outlier.shape = NA) +
scale_y_continuous(limits = quantile(dseb$values, c(0.1, 0.9), na.rm = T)) +
ylim(0, 0.3) + xlab('') + ylab('Computation time')
suppressWarnings(print(myplot))
dev.off()
#------------------------------------
colnames(cover) <- c("Bernoulli-VB","CS-VB",
"LAPLACE-VB")
pdf(file = "box17r.pdf",width=14,height=7)
dseb = stack(as.data.frame(cover))
names(dseb)[2] <- "Method"
myplot <- ggplot(dseb,aes(x = Method, y = values, color = Method)) + 
geom_boxplot(fatten = 2, notch=FALSE,outlier.shape = NA) +
scale_y_continuous(limits = quantile(dseb$values, c(0.1, 0.9))) + 
ylim(0.8, 1) + xlab('') + 
ylab('Coverage Probabilities')
suppressWarnings(print(myplot))
dev.off()
#------------------------------------
parsims = list()
parsims[[1]] = as.matrix(mod_sim1)
parsims[[1]] = cbind(parsims[[1]][,14],
parsims[[1]][,c(6:9)],
parsims[[1]][,c(10:13)],
parsims[[1]][,15:18],
parsims[[1]][,5],
parsims[[1]][,1:4])
colnames(parsims[[1]])[c(1,14)] <- c("intercept","alpha[0]")
parsims[[2]] = fit[[7]]$sims.matrix[(1:5000) %% 10 == 1, ]
parsims[[3]] = fit[[8]]$sims.matrix[(1:5000) %% 10 == 1, ]
mcmcnames = c("Bernoulli", "CS", "Laplace")
#------------------------------------
#------------------------------------
#----plot HPDs one by one -----------
#------------------------------------
#------------------------------------
j = 2
h = 3
	func <- function(x) dnorm(x,fit[[h]]$mu_beta[j],
		sqrt(fit[[h]]$sigma_beta[j,j]))
	hpddata = data.frame(x = tt, y = ss)
	gg_density <- ggplot(data=hpddata, aes(x=x)) + 
  	geom_function(fun = func, xlim = c(lower[j,h-2]-0.5,
		upper[j,h-2]+0.5))
	gg_density %>% 
	plot_credible_interval(lower[j,h-2], upper[j,h-2])
#------------------------------------
#------------------------------------
#----compute accuracies of betas-----
#------------------------------------
#------------------------------------
accuracy = matrix(0,3,p+1)
for(j in 1:(p+1)){
pdf(file = paste0("dens",j,".pdf"),width=14,height=7)
par(mfrow=c(1,3))
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
	plot(dd <- density(parsims[[h-2]][,j]),main = paste0(mcmcnames[h-2],"-prior posterior-density estimate for ",colnames(parsims[[h-2]])[j],"=",round(realbeta[j],3)),
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
	##ss = dnorm(tt,fit[[h]]$mu_beta[j],sqrt(fit[[h]]$sigma_beta[j,j]))
	post = kdensity(parsims[[h-2]][,j])$fhat
	ss = sapply(tt,qfunc)
	lines(ss~tt,col="red")
	if(sum(ss) == 1 | sum(ss) == 0) lines(c(0,0),c(0,max(dd$y)),lwd=3,col="red")
	legend(xlims[1],ylims[2],lty = rep(1,4),legend = c("MCMC","VB"),col=c("black","red"))
	distt = function(theta) sapply(theta,function(th){
	abs(post(th)-qfunc(th))})
	accuracy[h-2,j] = 100*(1-0.5*integrate(distt,-Inf,Inf,subdivisions=2000)$value)
}
dev.off()
}
colnames(accuracy) <- colnames(parsims[[1]])[1:(p+1)]
rownames(accuracy) <- mcmcnames
round(accuracy,2)
#------------------------------------
#------------------------------------
#-compute accuracies of Z and gamma--
#------------------------------------
#------------------------------------
accuracy2 = matrix(0,2,p)
cntr = 0
for(j in (p+2):(2*p+1)){
for(h in 3:4){
	cntr = cntr + 1
	pdf(file = paste0("bar",cntr,".pdf"),width=7,height=7)
	tb = table(parsims[[h-2]][,j])/sum(table(parsims[[h-2]][,j]))
	if(length(tb) == 1){
		if(names(tb) == "1") tb = c(0,1) else tb = c(1,0)
	}
	number <- c(dbinom(0:1,1,fit[[h]]$pst[j-p]),tb)
	type <- c(rep("VB",2),rep("MCMC",2))
	circle <- data.frame(number,type)
	print(ggplot(circle, aes(factor(rep(0:1,2)),number, fill = type)) +
  		geom_bar(stat="identity", position = "dodge") +
  		labs(title=paste("Barplots for Bernulli posterior for",colnames(parsims[[h-2]])[j],"using prior",mcmcnames[h-2])) +
  		xlab("values")+
  		ylab("probability"))
	dev.off()
	accuracy2[h-2,j-5] = 100*(1-0.5*(abs(number[1] - number[3]) +abs(number[2] - number[4])))
}
}
colnames(accuracy2) <- colnames(parsims[[1]])[(p+2):(2*p+1)]
rownames(accuracy2) <- mcmcnames[1:2]
round(accuracy2,2)
#------------------------------------
#------------------------------------
#------compute accuracies of pi------
#------------------------------------
#------------------------------------
accuracy3 = matrix(0,2,p)
for(j in (2*p+2):(3*p+1)){
pdf(file = paste0("dens",j,".pdf"),width=14,height=7)
par(mfrow=c(1,2))
limsx = limsy = matrix(0,2,2)
for(h in 3:4){
	tt = seq(0,1,length.out=100)
	ss = dbeta(tt,fit[[h]]$alphast[j-2*p],fit[[h]]$betast[j-2*p])
	limsy[1,h-2] = 0
	limsy[2,h-2] = max(ss)
}
xlims = c(0,1)
for(h in 3:4){
	ylims = c(0,max(density(parsims[[h-2]][,j])$y,max(limsy[2,][is.finite(limsy[2,])])))
	plot(density(parsims[[h-2]][,j]),main = paste0(mcmcnames[h-2],"-prior posterior-density estimate for ",colnames(parsims[[h-2]])[j]),
	xlim = c(0,1),
	ylim = ylims,xlab = "parameter",ylab = "posterior estimate")
	ss = dbeta(tt,fit[[h]]$alphast[j-2*p-1],fit[[h]]$betast[j-2*p-1])
	lines(ss~tt,col="red")
	legend(xlims[1],ylims[2],lty = rep(1,2),legend = c("MCMC","VB"),col=c("black","red"))
	post = kdensity(parsims[[h-2]][,j])$fhat
	qfunc = function(x) dbeta(x,fit[[h]]$alphast[j-2*p],fit[[h]]$betast[j-2*p])
	distt = function(theta) sapply(theta,function(th){
	abs(post(th)-qfunc(th))})
	accuracy3[h-2,j-2*p-1] = 100*(1-0.5*integrate(distt,0,1)$value)
}
dev.off()
}
colnames(accuracy3) <- colnames(parsims[[1]])[(2*p+2):(3*p+1)]
rownames(accuracy3) <- mcmcnames[1:2]
round(accuracy3,2)
#------------------------------------
#------------------------------------
#--compute accuracies of variances---
#------------------------------------
#------------------------------------
dinvgamma <- function(x,alpha,beta) beta^alpha / gamma(alpha) * 
(1/x)^(alpha+1) * exp(-beta/x) 
dgeninvgauss <- function(x,a,b) (a/b)^0.25 / BesselK(0.5, sqrt(a*b))/2 * 
x ^(-0.5) * exp(-(a*x + b/x)/2)
#------------------------------------
accuracy4 = matrix(0,1,p+1)
for(j in (3*p+2):(3*p+6)){
pdf(file = paste0("dens",j,".pdf"),width=14,height=7)
par(mfrow=c(1,1))
limsx = limsy = matrix(0,2,1)
h=3
dd = density(1/parsims[[h-2]][,j])
dd$y = dd$y/max(dd$y)
tt = seq(0.00001,mean(dd$x),length.out=100)
ss = dinvgamma(tt,fit[[h]]$sigmapars$a_beta[j-3*p-1],fit[[h]]$sigmapars$b_beta[j-3*p-1])
limsy[1,h-2] = 0
limsy[2,h-2] = max(ss)
#
xlims = c(0,3)
h=3
ylims = c(0,1)
plot(dd,main = paste0(mcmcnames[h-2],"-prior (normalized) posterior-density estimate for ",colnames(parsims[[h-2]])[j]),
xlim = c(0,mean(dd$x)),
ylim = ylims,xlab = "parameter",ylab = "posterior estimate")
ss = dinvgamma(tt,fit[[h]]$sigmapars$a_beta[j-3*p-1],fit[[h]]$sigmapars$b_beta[j-3*p-1])
ss = ss/max(ss)
lines(ss~tt,col="red")
legend(xlims[1],ylims[2],lty = rep(1,2),legend = c("MCMC","VB"),col=c("black","red"))
post = kdensity(1/parsims[[h-2]][,j])$fhat
qfunc = function(x) dinvgamma(x,fit[[h]]$sigmapars$a_beta[j-3*p-1],fit[[h]]$sigmapars$b_beta[j-3*p-1])
distt = function(theta) sapply(theta,function(th){
abs(post(th)-qfunc(th))})
accuracy4[h-2,j-3*p-1] = 100*(1-0.5*integrate(distt,0,Inf,subdivisions=2000,rel.tol =1e-3)$value)
dev.off()
}
colnames(accuracy4) <- colnames(parsims[[1]])[(3*p+2):(3*p+6)]
rownames(accuracy4) <- mcmcnames[1]
round(accuracy4,2)
#------------------------------------
#------------------------------------
#--compute accuracies of variances---
#------------------------------------
#------------------------------------
accuracy5 = matrix(0,1,1)
j=3*p+2
pdf(file = paste0("dens",19,".pdf"),width=14,height=7)
par(mfrow=c(1,1))
limsx = limsy = matrix(0,2,1)
h=4
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
plot(dd,main = paste0(mcmcnames[h-2],"-prior (normalized) posterior-density estimate for ",colnames(parsims[[h-2]])[j]),
xlim = c(0,mean(dd$x)),
ylim = ylims,xlab = "parameter",ylab = "posterior estimate")
ss = dinvgamma(tt,fit[[h]]$sigmapars$alpha_sigma2,fit[[h]]$sigmapars$beta_sigma2)
ss = ss/max(ss)
lines(ss~tt,col="red")
legend(xlims[1],ylims[2],lty = rep(1,2),legend = c("MCMC","VB"),col=c("black","red"))
post = kdensity(1/parsims[[h-2]][,j])$fhat
qfunc = function(x) dinvgamma(x,fit[[h]]$sigmapars$alpha_sigma2,fit[[h]]$sigmapars$beta_sigma2)
distt = function(theta) sapply(theta,function(th){
abs(post(th)-qfunc(th))})
accuracy5[h-3,j-3*p-1] = 100*(1-0.5*integrate(distt,0,Inf,subdivisions=2000)$value)
dev.off()
#
colnames(accuracy5) <- colnames(parsims[[2]])[3*p+2]
rownames(accuracy5) <- mcmcnames[2]
round(accuracy5,2)
#------------------------------------
#------------------------------------
#--compute accuracies of variances---
#------------------------------------
#------------------------------------
accuracy6 = matrix(0,1,p+1)
for(j in (p+2):(p+6)){
pdf(file = paste0("dens",j+14,".pdf"),width=14,height=7)
par(mfrow=c(1,1))
limsx = limsy = matrix(0,2,1)
h=5
if(j == (p+2)){
	dd = density(1/parsims[[h-2]][,j])
}else{
	dd = density(parsims[[h-2]][,j])
}
dd$y = dd$y/max(dd$y)
tt = seq(0.00001,mean(dd$x),length.out=100)
if(j == (p+2)){
	ss = dinvgamma(tt,2,fit[[h]]$sigmapars$Etau0m1)
} else {
	aa = fit[[h]]$sigmapars$lambdastar[j-p-2]/fit[[h]]$sigmapars$mustar[j-p-2]^2
	bb = fit[[h]]$sigmapars$lambdastar[j-p-2]
	ss = dgeninvgauss(tt,aa,bb)
}
limsy[1,h-4] = 0
limsy[2,h-4] = max(ss)
#
xlims = c(0,3)
h=5
ylims = c(0,1)
plot(dd,main = paste0(mcmcnames[h-2],"-prior (normalized) posterior-density estimate for ",colnames(parsims[[h-2]])[j]),
xlim = c(0,mean(dd$x)),
ylim = ylims,xlab = "parameter",ylab = "posterior estimate")
if(j == (p+2)){
	ss = dinvgamma(tt,2,fit[[h]]$sigmapars$Etau0m1)
} else {
	aa = fit[[h]]$sigmapars$lambdastar[j-p-2]/fit[[h]]$sigmapars$mustar[j-p-2]^2
	bb = fit[[h]]$sigmapars$lambdastar[j-p-2]
	ss = dgeninvgauss(tt,aa,bb)
}
ss = ss/max(ss)
lines(ss~tt,col="red")
legend(xlims[1],ylims[2],lty = rep(1,2),legend = c("MCMC","VB"),col=c("black","red"))
if(j == (p+2)){
	post = kdensity(1/parsims[[h-2]][,j])$fhat
}else{
	post = kdensity(parsims[[h-2]][,j])$fhat
}
qfunc = function(x){
if(j == (p+2)){
	ss = dinvgamma(x,2,fit[[h]]$sigmapars$Etau0m1)
} else {
	aa = fit[[h]]$sigmapars$lambdastar[j-p-2]/fit[[h]]$sigmapars$mustar[j-p-2]^2
	bb = fit[[h]]$sigmapars$lambdastar[j-p-2]
	ss = dgeninvgauss(x,aa,bb)
}
ss
}
distt = function(theta) sapply(theta,function(th){
abs(post(th)-qfunc(th))})
accuracy6[h-4,j-p-1] = 100*(1-0.5*integrate(distt,0,Inf,subdivisions=2000)$value)
dev.off()
}
colnames(accuracy6) <- colnames(parsims[[3]])[(p+2):(p+6)]
rownames(accuracy6) <- mcmcnames[3]
round(accuracy6,2)
#------------------------------------
#------------------------------------




