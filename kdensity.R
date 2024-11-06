#-----------my density function---------
boxcar<-function(t) 0.5*(abs(t)<=1) 
gaussian<-function(t) dnorm(t)
epanechnikhov<-function(t) 0.75*(1-t^2)*(abs(t)<=1)
tricube<-function(t) 70/81*(1-abs(t)^3)^3*(abs(t)<=1)
#
kdensity<-function(x,kernel="gaussian",bandwidth=NULL){
	ker<-tryCatch({match.fun(kernel)},error=function(e){
		stop("kernel is not known")
	})
	if(missing(bandwidth)| is.null(bandwidth)){
		s=sd(x)
		Q=quantile(x,0.75)-quantile(x,0.25)
		sigma.hat=min(s,Q/1.34)
		h=1.06*sigma.hat/length(x)^(1/5)
	}else{
		h<-bandwidth
	}
	fhat<-function(t) mean(ker((t-x)/h))/h
	out = list(fhat=fhat,x=x,bandwidth=h)
	class(out)<-"kdensity"
	return(out)
}
