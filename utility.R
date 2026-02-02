fprime <- function(x){
  sqrt(pi/2/x) * expint_E1(2*x) * exp(x)
  #####sqrt(pi/2/x) * myexpint_E1(2*x) * exp(x)
}

.tau_update = function(astar,bstar){
		B1 = sapply(sqrt(astar*bstar), 
			function(z){BesselK(1.5, z)})
			#####function(z){mybesselK(1.5, z)})
		B2 = sapply(sqrt(astar*bstar), 
			function(z){BesselK(0.5, z)})
			#####function(z){mybesselK(0.5, z)})
		grad = sapply(sqrt(astar*bstar), 
			function(z){fprime(z)})/B2
  		Etaum1 = B1/B2 * sqrt(astar/bstar) - 
				1/bstar
  		Etau = B1/B2 * sqrt(bstar/astar)
		Elogtau = log(sqrt(bstar/astar)) + grad 
	list(Etaum1 = Etaum1, Etau = Etau, Elogtau = Elogtau)
}
readkeygraph <- function(prompt)
{
    getGraphicsEvent(prompt = prompt, 
                 onMouseDown = NULL, onMouseMove = NULL,
                 onMouseUp = NULL, onKeybd = onKeybd,
                 consolePrompt = "[click on graph then follow top prompt to continue]")
    Sys.sleep(0.01)
    return(keyPressed)
}

onKeybd <- function(key)
{
    keyPressed <<- key
}

mybesselK <- function(a, t){
	out <- besselK(a,t)
	if(!is.finite(out) | is.nan(out)) out <- 1e100
	out
}


myexpint_E1 <- function(a){
	out <- tryCatch({expint_E1(a)},error = function(e){1e-10})
	if(!is.finite(out) | is.nan(out)) out <- 1e-10
	out
}


gelman.diag2 <- function (x, confidence = 0.95,
autoburnin = FALSE){
    x <- as.mcmc.list(x)
    if (nchain(x) < 2) 
        stop("You need at least two chains")
    if (autoburnin && start(x) < end(x)/2) 
        x <- window(x, start = end(x)/2 + 1)
    Niter <- niter(x[[1]])
    Nchain <- nchain(x)
    Nvar <- nvar(x[[1]])
    xnames <- varnames(x[[1]])
    x <- lapply(x, as.matrix)
    S2 <- array(sapply(x, var, simplify = TRUE), dim = c(Nvar, 
        Nvar, Nchain))
    W <- apply(S2, c(1, 2), mean)
    xbar <- matrix(sapply(x, apply, 2, mean, simplify = TRUE), 
        nrow = Nvar, ncol = Nchain)
    B <- Niter * var(t(xbar))
    w <- diag(W)
    b <- diag(B)
    s2 <- matrix(apply(S2, 3, diag), nrow = Nvar, ncol = Nchain)
    muhat <- apply(xbar, 1, mean)
    var.w <- apply(s2, 1, var)/Nchain
    var.b <- (2 * b^2)/(Nchain - 1)
    cov.wb <- (Niter/Nchain) * diag(var(t(s2), t(xbar^2)) - 2 * 
        muhat * var(t(s2), t(xbar)))
    V <- (Niter - 1) * w/Niter + (1 + 1/Nchain) * b/Niter
    var.V <- ((Niter - 1)^2 * var.w + (1 + 1/Nchain)^2 * var.b + 
        2 * (Niter - 1) * (1 + 1/Nchain) * cov.wb)/Niter^2
    df.V <- (2 * V^2)/var.V
    df.V[is.nan(df.V)] <- 1
    df.adj <- (df.V + 3)/(df.V + 1)
    B.df <- Nchain - 1
    W.df <- (2 * w^2)/var.w
    W.df[is.nan(W.df)] <- 1
    R2.fixed <- (Niter - 1)/Niter
    R2.random <- (1 + 1/Nchain) * (1/Niter) * (b/w)
    R2.random[is.nan(R2.random)] <- 0
    R2.estimate <- R2.fixed + R2.random
    R2.upper <- R2.fixed + qf((1 + confidence)/2, B.df, W.df) * 
        R2.random
    psrf <- cbind(sqrt(df.adj * R2.estimate), sqrt(df.adj * R2.upper))
    dimnames(psrf) <- list(xnames, c("Point est.", "Upper C.I."))
    out <- list(psrf = psrf)
    class(out) <- "gelman.diag"
    out
}
