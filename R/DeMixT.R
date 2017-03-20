###DeMixT running function for a pipeline
DeMixT <- function(inputdata, groupid, niter = 10, ninteg1 = 50,ninteg2 = 50, filter.option = 1, filter.criteria1 = c(0.5,0.5), filter.criteria2 = c(250,250), 
filter.criteria3 = 0.25, if.filter = FALSE){
##Begin Step 1
	print("Step 1 is started for estimating proportions.")
	res.S1 <- DeMixT.S1(inputdata, groupid, niter = niter, ninteg = ninteg1, filter.option = filter.option, filter.criteria1 = filter.criteria1, filter.criteria2 = filter.criteria2, 
			  filter.criteria3 = filter.criteria3, if.filter = if.filter)
              givenpi1 <- res.S1$pi 
			  givenpi <- c(t(givenpi1))# expand the piout matrix to a vector
    pi_iteration <- res.S1$pi_iteration
	
##Begin Step2
	print("Step 2 is started for deconvolution of expressions.")
	res.S2 <- DeMixT.S2(inputdata, groupid, givenpi, ninteg = ninteg2,  filter.option = filter.option)
	return(list(pi = givenpi1, decovExpr = res.S2$decovExpr, decovMu = res.S2$decovMu, decovSigma = res.S2$decovSigma, pi_iteration = pi_iteration))

}











