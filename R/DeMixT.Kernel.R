###main function
DeMixT.Kernel <- function(inputdata, groupid, nhavepi, givenpi, givenpiT, niter, ninteg){
	               if(!is.matrix(inputdata)) stop("argument inputdata must be a matrix");
	               if(!is.vector(groupid)) stop("argument groupid must be a vector"); 
				   if((FALSE %in% (unique(groupid) %in% c(1,2,3))) == TRUE) stop("argument groupid must be a vector of 1, 2, 3")
                   if((ninteg %% 1 != 0)||(ninteg <= 0)) stop(print("ninteg can only be an positive integer"))
				   if((niter %% 1 != 0)||(niter <= 0)) stop(print("niter can only be an positive integer"))
	               nsub <- as.integer(dim(inputdata)[2])   ## number of all the samples(normal + mixed)
				   wgenes <- as.integer(dim(inputdata)[1]) ## numebr of genes
	               intx=sum(groupid == 3);                ## number of mixed samples
				   intn=nsub-intx;                        ## number of normal samples
	               nCid = length(unique(groupid))-1       ## number of known components
				   groupid<-(as.array(groupid))
                   ######################process for reorganizing groupid###############################
	               groupid1 = which(groupid == 1)
	               if(nCid > 1) groupid2 = which(groupid == 2)
				   groupidT = which(groupid == 3)
	
	               sample.id = colnames(inputdata[,groupidT])
	               gene.id = row.names(inputdata)
	
				   if(nCid == 1){
		               inputdata = cbind(inputdata[, groupid1], inputdata[, groupidT])
					   groupid = c(rep(1, length(groupid1)), rep(3, length(groupidT)))
				   }else if(nCid == 2){
					   inputdata = cbind(inputdata[, groupid1], inputdata[, groupid2], inputdata[, groupidT])
					   groupid = c(rep(1, length(groupid1)), rep(2, length(groupid2)), rep(3, length(groupidT)))
				   }
					   
	               dataarray1 <- (as.array(matrix(inputdata, nrow = 1, byrow = F))) ## expression profile of the normal sample as an array
                   ## Set given proportions
	               if(nhavepi==1){
					   if(!is.vector(givenpi)) stop("argument option must be a vector if pi is known");
					   givenpi=(as.array(givenpi))
					   givenpi3=(as.array(rep(0, intx)))
				   }else if(nhavepi==2){
					   givenpi=(as.array(rep(0, nCid*intx)))
					   givenpi3 = as.array(givenpiT)    
				   }else if(nhavepi==0){
					   givenpi=(as.array(rep(0, nCid*intx)))
					   givenpi3=(as.array(rep(0, intx)))
				   }else{
				       stop("nhavepi argument must be set 0, 1 or 2")
				   }
	
	              if(nCid==1){
					  givenpi1 = givenpi[1:intx]
					  givenpi2 = (as.array(rep(0, intx)))
				  }else{
					  givenpi1 = givenpi[1:intx]
					  givenpi2 = givenpi[(intx+1):(2*intx)]
				  }
	
	            rres <- .C("Tdemix", dataarray1, as.integer(groupid), as.integer(nsub), as.integer(wgenes), as.integer(nhavepi), givenpi1, givenpi2, givenpi3, as.integer(nCid), as.integer(niter), as.integer(ninteg), rep(0, 2*intx),  rep(0, intx*wgenes), rep(0, niter*wgenes), rep(0, niter*wgenes), rep(0, niter*intx), rep(0, niter*intx))
	            outcome1<-matrix(rres[[12]], ncol=intx, nrow=2, byrow=T)
	            outcome2<-matrix(rres[[13]], ncol=(intx), nrow=wgenes, byrow = T)
	            outcome3<-matrix(rres[[14]], ncol=niter, nrow=wgenes,byrow=T)
	            outcome4<-matrix(rres[[15]], ncol=niter, nrow=wgenes,byrow=T)
	            outcome5<-matrix(rres[[16]], ncol=niter, nrow=intx,byrow=T)
	            outcome6<-matrix(rres[[17]], ncol=niter, nrow=intx,byrow=T)
	
	            print(paste0('Totally ', round(100*sum(outcome4[,niter]>0.99)/wgenes,3),'% genes estimated touch the optimization bound'))
			    return(list(pi = outcome1, decovExpr = outcome2, decovMu = outcome3, decovSigma = outcome4, pi1 = outcome5, pi2= outcome6))
}

