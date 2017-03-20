###DeMixT step II function
DeMixT.S2 <- function(inputdata, groupid, givenpi,ninteg = 50,  filter.option = 1){
	sample.id <- colnames(inputdata[, groupid == 3])
    if(is.null(row.names(inputdata))) row.names(inputdata) <- 1:nrow(inputdata)
	if(filter.option == 1){
		inputdata <- inputdata[apply(inputdata, 1, FUN = function(x) sum(x <= 0) == 0), ]
        gene.id <- row.names(inputdata)

	}else if(filter.option == 2){
	    inputdata <- inputdata + 1
		gene.id <- row.names(inputdata)
	}else{
		stop("the argument filter.option can only be 1 or 2")
	}
	
	res <- DeMixT.Kernel(inputdata, groupid, nhavepi = 1, givenpi=givenpi, givenpiT = rep(0,sum(groupid==3)), niter=1, ninteg=ninteg)

    decovExpr <- res$decovExpr; colnames(decovExpr) = sample.id; row.names(decovExpr) = gene.id
	decovMu <- res$decovMu; colnames(decovMu) = c('Mu'); row.names(decovMu) = gene.id
	decovSigma <- res$decovSigma; colnames(decovSigma) = c('Sigma'); row.names(decovSigma) = gene.id
	
	return(list(decovExpr = decovExpr, decovMu = decovMu, decovSigma = decovSigma))
	
}



