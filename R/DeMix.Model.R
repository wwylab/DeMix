#DeMix.Model, which is the main function 
############################################################
DeMix.Model <- function(inputdata, groupid, givenpi,
                        NormalizationStep, NormalizationMethod,
                        nStage,niter=1500, seeds = c(629555906.00, 921927245.00, 1265635378.00), pgenes=c(1,2,3,4,5,6,7,8,9,10)) 
{
  #Define all the parameters
  # niter=1500 iteration is defaulted as 1500
  
  nsub=length(groupid)
  intx=sum(groupid)
  intn=nsub-intx
  conc=0.8
  fc=1.2
  ninteg = 30 
  core.num <- 1

  #nthread <- core.num 
  cbit <- 64;
							
  col.id <- colnames(inputdata)[groupid==1]

  ##if((ninteg %% 1 != 0)||(ninteg <= 0)) stop(print("ninteg can only be an positive integer"))
  ##if((nthread %% 1 != 0)||(nthread <= 0)) stop(print("nthread can only be an positive integer"))
  if (.Machine$sizeof.pointer < 8) cbit = 32
  
  #Determine whether we need to do the normalization or not, which is optional for users.
  if (NormalizationStep == 0){
      print("No need to do the normalization for your input data")
      newt <- inputdata
  }
    
  if (NormalizationStep == 1){
    newt <- DeMix.Normalization(inputdata,NormalizationMethod,groupid)
  }
    
  
  #Next is to filter the data
  inputdata1 <- DeMix.Filter(newt, groupid=c(rep(0,intn),rep(1,intx)), zero_filter=TRUE, conc=0.8, fc=1.2)
							
  if (sum(inputdata1 <= 0) > 0) stop(print("the input data matrix can only be positive values"))

  inputdata1 <- as.matrix(inputdata1)
  
  #Process data based on clean data
  
  nsub <- as.integer(dim(inputdata1)[2])
  wgenes <- as.integer(dim(inputdata1)[1])
  
  rnan1<-inputdata1[, groupid==0] #denote normal as 0
  rnat1<-inputdata1[, groupid==1] #denote tumor as 1
  ovsn1=((apply(rnan1,1,  sd)^2)-apply(rnan1,1,  mean)+1) / (apply(rnan1,1,  mean)+1)^2
  ovst1=((apply(rnat1,1,  sd)^2)-apply(rnat1,1,  mean)+1) / (apply(rnat1,1,  mean)+1)^2
  
  
  newovs<-c(abs(ovsn1), abs(ovst1))
  
  if(!is.vector(pgenes) || !is.numeric(pgenes)){
    pgenes=seq(1,10,1)
  } else {
    if(length(pgenes)<10)
      pgenes=c(pgenes, seq(1,(10-length(pgenes)),1))
    else
      pgenes=pgenes[1:10]
    
    if(max(pgenes)>wgenes)
      pgenes=seq(1,10,1)
    
    
  }
  
  
  
 if (givenpi == 1) {
   print("nStage can only be either 1 or 3 with unknown pi")
   #Stage 1 with the givenpi=1
   groupid<-(as.array(groupid))
   # if(ninteg<10) ninteg=10;
   tmp <- 1;
   #Run the main C code
 
   rres1 <- .C("Bdemix", inputdata1, as.integer(groupid), as.integer(nsub),
             as.integer(wgenes),  as.integer(cbit), as.integer(tmp), givenpi,
             as.integer(ninteg), as.integer(niter), newovs, as.integer(pgenes),  rep(0, intx*3), rep(0, intx*3),
             rep(0, nsub*wgenes), rep(0, nsub*wgenes), rep(0, 2*wgenes),seeds, rep(0, niter*intx), rep(0, niter*10))


   pi_out<-rbind(matrix(rres1[[12]], ncol=intx, nrow=3, byrow=T), matrix(rres1[[13]], ncol=intx, nrow=3, byrow=T))
   rownames(pi_out) <- c('Posterior Mean (Poisson)', '95% CI (Poisson)', '95% CI (Poisson)', 
                      'Posterior Mean (NB)', '95% CI (NB)', '95% CI (NB)')
   colnames(pi_out) <- col.id
   mu<-matrix(rres1[[16]],ncol=2, nrow=wgenes, byrow=T)
   colnames(mu) <- c('normal', 'tumor')

   decon_normal<-matrix(rres1[[15]], ncol=(nsub) , nrow=wgenes, byrow = T)

   decon_tumor<-matrix(rres1[[14]], ncol=(nsub) , nrow=wgenes, byrow = T)

   post_pi<-matrix(rres1[18], ncol=(nsub), nrow=niter, byrow = T)
   post_mu<-matrix(rres1[19], ncol=(nsub), nrow=niter, byrow = T)

   #Determine the stage number to decide to continue or stop
   out2 = pi_out
 
   #if (nStage!=1||nStage!=3) stop(print("Stage number must be 1 or 2"))

   if(nStage==1) 
     return(list(pi=pi_out,decovMu=mu, decovExpr_normal=decon_normal, decovExpr_tumor=decon_tumor,
                 post_pi=post_pi,post_mu=post_mu))
 	 

   #Go to run C code again when stage is 3
   if(nStage==3) {
      givenpi <- (as.array(out2[4,]))
      inputdata3 <- newt[apply(newt, 1, FUN = function(x) sum(x <= 0) == 0), ]
      rnan3<-inputdata3[, groupid==0] #denote normal as 0
      rnat3<-inputdata3[, groupid==1] #denote tumor as 1
      ovsn3=((apply(rnan3,1,  sd)^2)-apply(rnan3,1,  mean)+1) / (apply(rnan3,1,  mean)+1)^2
      ovst3=((apply(rnat3,1,  sd)^2)-apply(rnat3,1,  mean)+1) / (apply(rnat3,1,  mean)+1)^2
      newovs<-c(abs(ovsn3), abs(ovst3))
      
      row.id <- row.names(newt)[apply(newt, 1, FUN = function(x) sum(x <= 0) == 0)]

      wgenes <- as.integer(dim(inputdata3)[1])
      tmp <- 2
      rres3 <- .C("Bdemix", inputdata3,  as.integer(groupid), as.integer(nsub),
               as.integer(wgenes),  as.integer(cbit), as.integer(tmp), givenpi,
               as.integer(ninteg), as.integer(niter), newovs, as.integer(pgenes),  
               rep(0, intx*3), rep(0, intx*3), rep(0, nsub*wgenes), rep(0, nsub*wgenes), 
               rep(0, 2*wgenes),seeds, rep(0, niter*intx), rep(0, niter*10))

      mu<-matrix(rres3[[16]], ncol=2, nrow=wgenes, byrow=T)
      colnames(mu) <- c('normal', 'tumor')

      decon_normal            <- matrix(rres3[[15]], ncol=(nsub) , nrow=wgenes, byrow = T)
      decon_tumor             <- matrix(rres3[[14]], ncol=(nsub) , nrow=wgenes, byrow = T)
      row.names(decon_tumor)  <- row.id 
      row.names(decon_normal) <- row.id
      row.names(mu)           <- row.id
							

      decon_tumor <- decon_tumor[, groupid == 1]							
      decon_normal <- decon_normal[, groupid == 1]	
      colnames(decon_tumor) = col.id; colnames(decon_normal) = col.id;
      post_pi<-matrix(rres3[18], ncol=(nsub), nrow=niter, byrow = T)
      post_mu<-matrix(rres3[19], ncol=(nsub), nrow=niter, byrow = T)
			
      return(list(decovExpr_normal = decon_normal, decovExpr_tumor = decon_tumor, pi = out2, 
                  decovMu = mu,post_pi=post_pi,post_mu=post_mu)) 
      }
 } else{   
     simu_givenpi <- (as.array(givenpi))
     nStage <- 2
     inputdata2 <- newt[apply(newt, 1, FUN = function(x) sum(x <= 0) == 0), ]
     rnan2<-inputdata2[, groupid==0] #denote normal as 0
     rnat2<-inputdata2[, groupid==1] #denote tumor as 1
     ovsn2=((apply(rnan2,1,  sd)^2)-apply(rnan2,1,  mean)+1) / (apply(rnan2,1,  mean)+1)^2
     ovst2=((apply(rnat2,1,  sd)^2)-apply(rnat2,1,  mean)+1) / (apply(rnat2,1,  mean)+1)^2
     newovs<-c(abs(ovsn2), abs(ovst2))
     
     row.id <- row.names(newt)[apply(newt, 1, FUN = function(x) sum(x <= 0) == 0)]
 
     wgenes <- as.integer(dim(inputdata2)[1])
     if(wgenes > 1000) {print("The number of genes has exceeded 1000, please expect long runtime for your data. We suggest parallel computing in order to reduce the runtime, e.g., by dividing the genes into smaller groups, such as 1000 genes per group. You may find an example in the help file of DeMix.Model() by ‘?DeMix.Model’.")}
     tmp <- 2
     rres2 <- .C("Bdemix", inputdata2,  as.integer(groupid), as.integer(nsub),
             as.integer(wgenes),  as.integer(cbit), as.integer(tmp), simu_givenpi,
             as.integer(ninteg), as.integer(niter), newovs, as.integer(pgenes),  rep(0, intx*3), rep(0, intx*3),
             rep(0, nsub*wgenes), rep(0, nsub*wgenes), rep(0, 2*wgenes),seeds, rep(0, niter*intx), rep(0, niter*10))
 
     mu<-matrix(rres2[[16]], ncol=2, nrow=wgenes, byrow=T)
     colnames(mu) <- c('normal', 'tumor')
 
     decon_normal <- matrix(rres2[[15]], ncol=(nsub) , nrow=wgenes, byrow = T)
     decon_tumor <- matrix(rres2[[14]], ncol=(nsub) , nrow=wgenes, byrow = T)
     row.names(decon_tumor) = row.id; row.names(decon_normal) = row.id;row.names(mu) = row.id
 
 
     decon_tumor <- decon_tumor[, groupid == 1]							
     decon_normal <- decon_normal[, groupid == 1]	
     colnames(decon_tumor) = col.id; colnames(decon_normal) = col.id;
     post_pi<-matrix(rres2[18], ncol=(nsub), nrow=niter, byrow = T)
     post_mu<-matrix(rres2[19], ncol=(nsub), nrow=niter, byrow = T)
 
     return(list(decovExpr_normal = decon_normal, decovExpr_tumor = decon_tumor, pi = simu_givenpi, 
                 decovMu = mu,post_pi=post_pi,post_decovMu=post_mu)) 
 }
}
