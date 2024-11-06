fprime <- function(x){
  sqrt(pi/2/x) * expint_E1(2*x) * exp(x)
}
.tau_update = function(astar,bstar){
		B1 = sapply(sqrt(astar*bstar), 
			function(z){BesselK(1.5, z)})
		B2 = sapply(sqrt(astar*bstar), 
			function(z){BesselK(0.5, z)})
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