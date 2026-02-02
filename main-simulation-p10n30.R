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
Iterations = 100
seb = se = set = FPR = FNR = syst = matrix(0,5,Iterations)
rownames(seb) <- rownames(se) <-
rownames(set) <- rownames(FPR) <-
rownames(FNR) <- rownames(syst) <- c("LASSO","SCAD","Bernoulli-VB","CS-VB",
"LAPLACE-VB")
riskmethod = "aic"
#------------------------------------
fit = betahat = pst = yhat = yhattest = 
coverage = list()
p = 10
n = 30
mu_beta = 0.7
sigma_beta = 0.5
beta = rnorm(p, mu_beta, sigma_beta)
z = c(1, c(0,1,0,0,0,1,0,1,0))
realbeta = beta * z
cover = matrix(0,p,3)
rownames(cover) <- paste0("beta",0:(p-1))
colnames(cover) <- c("Bernoulli-VB","CS-VB",
"LAPLACE-VB")
for(iter in 1:Iterations){
	p = 10
	n = 30
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
	####tLL <- fit[[2]]$loss
	tLL <- fit[[2]]$deviance
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
	p0 = 0.5
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
	c = 1e-2
	syst[4,iter] = system.time(fit[[4]] <- tryCatch({
	sppoissregvb(X,y,init,prior="CS", 
		eps = 1e-12, maxiter = 100)},error=function(e){NULL}))[[3]]
	#####p0 = max(0.01,min(sum(pst[[1]])/p,0.99))
	p0 = 0.5
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
	#------------------------------------
	#
	cat("criteria \n")
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
	cat("coverage probs \n")
	lower = upper = matrix(0,p,6)
	for(j in 1:p){
		for(h in 3:5){
			tt = seq(fit[[h]]$mu_beta[j] - 3*sqrt(abs(fit[[h]]$sigma_beta[j,j])),
			fit[[h]]$mu_beta[j] + 3*sqrt(abs(fit[[h]]$sigma_beta[j,j])),length.out=100)
			ss = dnorm(tt,fit[[h]]$mu_beta[j],sqrt(abs(fit[[h]]$sigma_beta[j,j])))
			object = list(x = tt, y = ss)
			hhpp <- hpd(object)
			lower[j,h-2] = hhpp[1]
			upper[j,h-2] = hhpp[2]
			cover[j,h-2] = (realbeta[j] >= lower[j,h-2])*(realbeta[j] <= upper[j,h-2])
		}
	}
	coverage[[iter]] = cover
}
#------------------------------------
save(seb,se,set,FPR,FNR,syst,coverage,
file = "results2-new.Rdata")
#------------------------------------
load(file.choose())
#------------------------------------
#------------------------------------
systr = syst
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
#------------------------------------
#------------------------------------
#--coverage probabilities for HPD----
#------------------------------------
#------------------------------------
(CM <- Reduce("+",coverage)/Iterations)
round(sqrt(Reduce("+",Map(function(x){(x-CM)^2},
coverage))/Iterations),2)
#------------------------------------
#------------------------------------
#------------------------------------
#---------------Box plots------------
#------------------------------------
#------------------------------------
pdf(file = "box21r.pdf",width=14,height=7)
dseb = stack(as.data.frame(t(seb[,1:(iter-1)])))
names(dseb)[2] <- "Method"
myplot <- ggplot(dseb,aes(x = Method, y = values, color = Method)) + 
geom_boxplot(fatten = 2, notch=FALSE,outlier.shape = 16) +
scale_y_continuous(limits = quantile(dseb$values, c(0.1, 0.9), na.rm = T)) + 
ylim(0.001, 0.075) + xlab('') + 
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
pdf(file = "box22r.pdf",width=14,height=7)
dseb = stack(as.data.frame(t(se[,1:(iter-1)])))
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
ylab('Train relative error (TRRE)')
suppressWarnings(print(myplot))
dev.off()
#------------------------------------
pdf(file = "box23r.pdf",width=14,height=7)
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
pdf(file = "box24r.pdf",width=14,height=7)
dseb = stack(as.data.frame(t(FPR[,1:(iter-1)])))
names(dseb)[2] <- "Method"
myplot <- ggplot(dseb,aes(x = Method, y = values, color = Method)) + 
geom_boxplot(fatten = 2, notch=FALSE,outlier.shape = 16) +
scale_y_continuous(limits = quantile(dseb$values, c(0.1, 0.9), na.rm = T)) +
ylim(0, 1) + xlab('') + 
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
pdf(file = "box25r.pdf",width=14,height=7)
dseb = stack(as.data.frame(t(FNR[,1:(iter-1)])))
names(dseb)[2] <- "Method"
myplot <- ggplot(dseb,aes(x = Method, y = values, color = Method)) + 
geom_boxplot(fatten = 2, notch=FALSE,outlier.shape = 16) +
#scale_y_continuous(limits = quantile(dseb$values, c(0.1, 0.995), na.rm = T)) +
ylim(0, 0.9) + xlab('') + 
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
pdf(file = "box26r.pdf",width=14,height=7)
dseb = stack(as.data.frame(t(systr[,1:(iter-1)])))
names(dseb)[2] <- "Method"
myplot <- ggplot(dseb,aes(x = Method, y = values, color = Method)) + 
geom_boxplot(fatten = 2, notch=FALSE,outlier.shape = NA) +
scale_y_continuous(limits = quantile(dseb$values, c(0.1, 0.9), na.rm = T)) +
ylim(0, 0.1) + xlab('') + 
theme_bw() + 
theme(legend.position="none") + 
theme(axis.text.x = element_text(size = 12)) + 
theme(axis.text.y = element_text(size = 12)) +
theme(axis.title.x = element_text(size = 12)) + 
theme(axis.title.y = element_text(size = 12)) + 
ylab('Computation time')
suppressWarnings(print(myplot))
dev.off()
#------------------------------------
