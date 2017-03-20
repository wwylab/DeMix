/***************************************************************************************************************
* DeMixT Latest Version
* By Zeya Wang
* Initial Date 12.13.2014
****************************************************************************************************************/

  #include <stdio.h>
  #include <stdlib.h>
  #include <math.h>
  #include <time.h>
  #include <string.h>
  #include <float.h> /* DBL_EPSILON */
  #include "DeMixTH.h"




void Tdemix(double *data, int *nGroup, int *nsamp, int *ngenes, int *npi, double *fixpi1, double *fixpi2, double *fixpi3,
            int *nCid, int *niter, int *ninteg, double *output1, double *output3, double *output5, double *output7,
            double *output9, double *output11)
{
  //nCid = 1, we have just one stroma component; =2, we have 2
  int i, j, k, l;
  int iteration;

  double *st1_mu, *st2_mu, *st1_sig, *st2_sig;
  double *st1_mu_2, *st2_mu_2, *st1_sig_2, *st2_sig_2;
  double **stroma1, **stroma2;
  double **stroma1_2, **stroma2_2;

  clock_t start_t, end_t;
  opt_seq = 1;


  //get value from input parameter
  nS=*nsamp;            // Number of Samples
  nG=*ngenes;           // Number of Genes
  nHavepi =  *npi;      // Have pi or not: 0: Not given pi1 and pi2; 1: given pi1 and pi2; 2: given piT
  Cid = *nCid; //indicator for number of components
  integ = *ninteg;
  //printf("integ is %d\n", integ);

  iteration = *niter; // number of iterations
  printf("iteration is %d\n", iteration);



  // 1 for component 1, 2 for component 2, 3 for mixed sample
  fNorm1=0;
  fNorm2=0;
  for(j=0;j<nS;j++)
  {
    if(nGroup[j]==1) fNorm1++;
    if(nGroup[j]==2) fNorm2++;
    printf("%d : %d \n", j, nGroup[j]);
  }

  // Determine the number of tumor sample
  fNorm = fNorm1 + fNorm2;
  intx=nS-fNorm;




  nmle = (Cid+1)*intx;
  // Random number generation. *ct can have 32 or 64 depending on the machine
  //M = exp(-*ct*log(2.0));

  p = (PARAMM *)calloc(1,sizeof(PARAMM));

  // FD, Parameter initialized
  FD = calloc(nS ,sizeof(double *));


  for(j=0;j<nS;j++) FD[j]= calloc(nG, sizeof(double));


  printf("Setting is over\n");

  // Data transfer
  initialSet(p);
  load_data(data);

  printf("Loading is over\n");

  tmppi1 =calloc(iteration ,sizeof(double *));  // save pi1 estimates in each iteration
  tmppi2 =calloc(iteration ,sizeof(double *));  // save pi2 estimates in each iteration

  for(j=0;j<iteration;j++)
  {
    tmppi1[j]= calloc(intx, sizeof(double ));
    tmppi2[j]= calloc(intx, sizeof(double ));
  }
  printf("Loading1 is over\n");

  //calculate mean and sd of normal samples
  stroma1 = calloc(nG, sizeof(double *));
  stroma2 = calloc(nG, sizeof(double *));
  st1_mu = calloc(nG ,sizeof(double));
  st2_mu = calloc(nG ,sizeof(double));
  st1_sig = calloc(nG ,sizeof(double));
  st2_sig = calloc(nG ,sizeof(double));
    
  stroma1_2 = calloc(nG, sizeof(double *));
  stroma2_2 = calloc(nG, sizeof(double *));
  st1_mu_2 = calloc(nG ,sizeof(double));
  st2_mu_2 = calloc(nG ,sizeof(double));
  st1_sig_2 = calloc(nG ,sizeof(double));
  st2_sig_2 = calloc(nG ,sizeof(double));

  printf("Loading2 is over\n");

  for(j=0;j<nG;j++) stroma1[j]= calloc(fNorm1, sizeof(double));
  for(j=0;j<nG;j++) stroma2[j]= calloc(fNorm2, sizeof(double));
  for(j=0;j<nG;j++) stroma1_2[j]= calloc(fNorm1, sizeof(double));
  for(j=0;j<nG;j++) stroma2_2[j]= calloc(fNorm2, sizeof(double));

  printf("Loading3 is over\n");

if(opt_seq == 1){
  for(i=0;i<nG;i++)
  {
    for(j=0;j<fNorm1;j++) stroma1[i][j] = log2(FD[j][i]);
    if(Cid == 2){
      for(j=0;j<fNorm2;j++) stroma2[i][j] = log2(FD[j+fNorm1][i]);
    }
  }
    for(i=0;i<nG;i++)
    {
        st1_mu[i] = mean(stroma1[i], fNorm1);
        st1_sig[i] = sd(stroma1[i], fNorm1);
        if(Cid == 2){
            st2_mu[i] = mean(stroma2[i], fNorm2);
            st2_sig[i] = sd(stroma2[i], fNorm2);
        }
    }
}else{
    for(i=0;i<nG;i++)
    {
        for(j=0;j<fNorm1;j++) stroma1_2[i][j] = FD[j][i];
        if(Cid == 2){
            for(j=0;j<fNorm2;j++) stroma2_2[i][j] = FD[j+fNorm1][i];
        }
    }
    for(i=0;i<nG;i++)
    {
        st1_mu_2[i] = mean(stroma1_2[i], fNorm1);
        st1_sig_2[i] = sd(stroma1_2[i], fNorm1);
        if(Cid == 2){
            st2_mu_2[i] = mean(stroma2_2[i], fNorm2);
            st2_sig_2[i] = sd(stroma2_2[i], fNorm2);
        }
    }
}
    
  for(i=0;i<nG;i++) free(stroma1[i]);
  for(i=0;i<nG;i++) free(stroma2[i]);
  for(i=0;i<nG;i++) free(stroma1_2[i]);
  for(i=0;i<nG;i++) free(stroma2_2[i]);

  free(stroma1);
  free(stroma2);
  free(stroma1_2);
  free(stroma2_2);
    
    
  for(j=0;j<nG;j++)
  {
 
    if(opt_seq == 1){
       p->Navg1[j]=st1_mu[j];
       if(Cid == 2) p->Navg2[j]=st2_mu[j];
       p->Tavg[j]=7.0;
       p->Nsigma1[j]= pow(st1_sig[j], 2.0);
       if(Cid == 2) p->Nsigma2[j]= pow(st2_sig[j], 2.0);
       p->Tsigma[j]=0.01;
    }else{
       p->Navg1[j]=st1_mu_2[j];
       if(Cid == 2) p->Navg2[j]=st2_mu_2[j];
       p->Tavg[j]=1000.0;
       p->Nsigma1[j]= (pow(st1_sig_2[j], 2.0) - p->Navg1[j])/pow(p->Navg1[j], 2.0);
       if(Cid == 2) p->Nsigma2[j]= (pow(st2_sig_2[j], 2.0) - p->Navg2[j])/pow(p->Navg2[j], 2.0);
       p->Tsigma[j] = 0.5;
    }
  }


  if(nHavepi==1)
  {
    for(k=0;k<intx;k++)
    {
      p->pi1[k]= fixpi1[k];
      if(Cid == 2) p->pi2[k]= fixpi2[k];
    }
  }else if(nHavepi == 2){
    for(k=0;k<intx;k++) p->piT[k] = fixpi3[k];
  }

  printf("Loading is over\n");

  free(st1_mu);
  free(st2_mu);
  free(st1_sig);
  free(st2_sig);

  //begin our iteration

//  for(j=0;j<nG;j++)
  //  {
  //    printf("nsigma1 for %d is %lf\n",j, p->Nsigma1[j]);
  //    printf("nsigma2 for %d is %lf\n",j, p->Nsigma2[j]);
  //    printf("navg1 for %d is %lf\n",j, p->Navg1[j]);
  //    printf("navg2 for %d is %lf\n",j, p->Navg2[j]);
//    }


  CD = calloc(intx ,sizeof(double *));
  for(j=0;j<intx;j++) CD[j]= calloc(nG, sizeof(double));
  printf("Iteration is starting\n");

  avgparT = calloc(iteration, sizeof(double *));
  sigparT = calloc(iteration, sizeof(double *));
  for(j = 0; j < iteration; j++) avgparT[j] = calloc(nG, sizeof(double));
  for(j = 0; j < iteration; j++) sigparT[j] = calloc(nG, sizeof(double));


  for(i=0;i<iteration;i++)
  {

    //updating pi value
    if(nHavepi==0)
    {
      for(l=0;l<intx;l++) getpi(l, Cid);
    }else if(nHavepi == 2){
      for(l=0;l<intx;l++) getpiT(l); // this is to estimate pi1 and pi2 given piT in the three component case

    }

    for(j=0;j<intx;j++)
    {
      tmppi1[i][j] =p->pi1[j];
      printf("%15.3f \t",  p->pi1[j]);
      if(Cid == 2)
      {
        tmppi2[i][j] =p->pi2[j];
        printf(" %15.3f \t",  p->pi2[j]);
      }
      printf("\n");
    }
    printf("iteration %d Updating Purity\n", i+1);

   



    if(i==0) start_t=clock();

    for(j=0;j<nG;j++)
    {
      gettumor(j, Cid);
      avgparT[i][j] = p->Tavg[j];
      sigparT[i][j] = p->Tsigma[j];
      // printf("Getting Iteration Tumor Parameters  %d:%d\n", i, j);
    }

      /*for(j=0;j<nG;j++)
      {
          printf("sigmat is %lf\n", p->Tsigma[j]);
          printf("avgt is %lf\n", p->Tavg[j]);
      }*/

    printf("iteration %d Updating Parameters\n", i+1);


    // Total time expectation depending on the first five iterations
    if(i==1)
    {
      end_t=clock();
      printf("================================================================\n");
      printf("Deconvolution is estimated to be finished in %lf hours\n", (double)(end_t-start_t)/CLOCKS_PER_SEC*((double)iteration/3600.0));
      printf("================================================================\n");
    }
  }

  //final pi value

  if(nHavepi==0)
  {
    for(l=0;l<intx;l++) getpi(l, Cid);
  }else if(nHavepi == 2){
    for(l=0;l<intx;l++) getpiT(l); // this is to estimate pi1 and pi2 given piT in the three component case
  }



  // deconvolve mle

  for(j=0;j<nG;j++)
  {
    for(k=0;k<intx;k++)
    {
      getmle(k, j, Cid);
      CD[k][j] = p->mle[j][k]; //deconvolved value
    }

  }

  fflush(NULL);
  saveFiles(output1, output3, output5, output7, output9, output11, iteration);

  printf("Step 4: Save files done : \n");
  //free variable
  free(p->Navg1);
  free(p->Navg2);
  free(p->Tavg);
  free(p->Nsigma1);
  free(p->Nsigma2);
  free(p->Tsigma);
  free(p->pi1);
  free(p->pi2);

  for(i=0;i<nG;i++)
  {
    free(p->mle[i]);
  }

  free(p->mle);
  free(p);
}




/************************************************************************************
  * Function to load data and save data
************************************************************************************/
  // Data transfer
static void load_data(double *mat1)
{
  int j,k;

  for(k=0;k<nS;k++)
  {
    for(j=0;j<nG;j++)
    {
      FD[k][j]=  mat1[nG*k+j];
    }
  }
  printf("There are  %d normals and %d tumors\n", fNorm,intx);

}




static void initialSet(PARAMM *qq)
{
  int i;
  qq->Navg1 =calloc(nG ,sizeof(double ));
  qq->Navg2 =calloc(nG ,sizeof(double ));
  qq->Tavg =calloc(nG ,sizeof(double ));
  qq->Nsigma1 =calloc(nG ,sizeof(double ));
  qq->Nsigma2 =calloc(nG ,sizeof(double ));
  qq->Tsigma = calloc(nG ,sizeof(double ));
  qq->pi1 =calloc(intx ,sizeof(double ));
  qq->pi2 =calloc(intx ,sizeof(double ));
  qq->piT = calloc(intx ,sizeof(double ));

  qq->mle =calloc(nG ,sizeof(double *));

  for(i=0;i<nG;i++) qq->mle[i]=calloc(nmle, sizeof(double ));

}







// Deconvolved expressions
static void saveFiles(double *put1, double *put3, double *put5, double *put7, double *put9, double *put11, int iteration)
{
  int i, j;
  //pi1 and pi2 save
  if(Cid == 2) //tow components
{
  for(i=0;i<intx;i++)
  {
    put1[i+intx*0] = p->pi1[i];
    put1[i+intx*1] = p->pi2[i];
  }
}else{
  for(i=0;i<intx;i++) put1[i+intx*0] = p->pi1[i];
}
  // pure tumor expression save.
  for(j=0;j<nG;j++)
  {
    for(i=0;i<intx;i++)
    {
      put3[j*intx+i] =  CD[i][j];
    }
  }
  //mu and sigma save for T
  for(j=0;j<nG;j++)
  {
    for(i=0;i<iteration;i++)
    {
      put5[j*iteration+i] = avgparT[i][j];
    }
  }

  for(j=0;j<nG;j++)
  {
    for(i=0;i<iteration;i++)
    {
      put7[j*iteration+i] = sigparT[i][j];
    }
  }

  for(j=0;j<intx;j++)
  {
    for(i=0;i<iteration;i++)
    {
      put9[j*iteration+i] = tmppi1[i][j];
    }
  }


  for(j=0;j<intx;j++)
  {
    for(i=0;i<iteration;i++)
    {
      put11[j*iteration+i] = tmppi2[i][j];
    }
  }
}



/************************************************************************************
  * get tumor function
*************************************************************************************/


  /*Tumor Component Update*/
static void gettumor(int genes, int h) // option h = 1 for one component; h = 2 for two
{
    double mu, sigma,upp;
    if(opt_seq == 1){
        upp = 16;
    }else{
        upp = 1500;
    }
    
  if(h == 1)
  {
    mu = tmin_y(2, upp, genes,h, mint, 0.001);
    
    sigma =  tmin_y2(0.0001, 1, genes, mu, tf_y, 0.001);
  }else{
    mu = tmin_y(2, upp, genes,h, mint, 0.001);
    sigma =  tmin_y2(0.0001, 1, genes, mu, tf_y2, 0.001);

  }
  p->Tavg[genes] = mu;
  p->Tsigma[genes] = sigma;

}

double tf_y(int genes, double mu, double sigma)
{
  double new_like;
  double tmp;
  int i;
  new_like=0.0;
  if(opt_seq == 1){ //microarray
     for (i=0; i<intx; i++){
         tmp = ft_y(FD[i+fNorm][genes],p->Navg1[genes], mu, p->Nsigma1[genes], sigma, p->pi1[i], 0);
         new_like+= tmp;
     }
  }else{
     for (i=0; i<intx; i++){
         tmp = ft_nb(FD[i+fNorm][genes],p->Navg1[genes], mu, p->Nsigma1[genes], sigma, p->pi1[i], 0);
         new_like+= tmp;
     }
  }
  return (new_like*(-1.0));
}


double tf_y2(int genes, double mu, double sigma)
{
  double new_like;
  double tmp;
  int i;
  new_like=0.0;
  if(opt_seq == 1){//microarray
      for (i=0; i<intx; i++){
          tmp = ft_y2(FD[i+fNorm][genes], p->Navg1[genes], p->Navg2[genes], mu, p->Nsigma1[genes], p->Nsigma2[genes], sigma, p->pi1[i], p->pi2[i]);
          new_like+= tmp;
      }
  }else{
      for (i=0; i<intx; i++){ 
          tmp = ft_nb2(FD[i+fNorm][genes], p->Navg1[genes], p->Navg2[genes], mu, p->Nsigma1[genes], p->Nsigma2[genes], sigma, p->pi1[i], p->pi2[i]);
          new_like+= tmp;
      }
  }
  return (new_like*(-1.0));
}


double mint(int genes, int h, double mu) // opt h =1:one component; h =2:two component
{
  double sigma;
  double tmp;
  if(h == 1)
  {
    sigma =  tmin_y2(0.0001, 1, genes, mu, tf_y, 0.001);
    tmp = tf_y(genes, mu, sigma);
  }else{
    sigma =  tmin_y2(0.0001, 1, genes, mu, tf_y2, 0.001);
    tmp = tf_y2(genes, mu, sigma);
  }

  return (tmp);
}



double tmin_y(double ax, double bx, int genes, int h, double (*f)(int, int, double), double tol)
{
    double a, b, d, e, p, q, r, u, v, w, x;
    double t2, fu, fv, fw, fx, xm, eps, tol1, tol3;
    const double c = (3. - sqrt(5.)) * .5;

    eps = DBL_EPSILON;
    tol1 = eps + 1.;
    eps = sqrt(eps);
    a = ax;
    b = bx;
    v = a + c * (b - a);
    w = v;
    x = v;
    d = 0.;
    e = 0.;
    fx = (*f)(genes, h, x);
    fv = fx;
    fw = fx;
    tol3 = tol / 3.;
    for(;;) {
      xm = (a + b) * .5;
      tol1 = eps * fabs(x) + tol3;
      t2 = tol1 * 2.;
      if (fabs(x - xm) <= t2 - (b - a) * .5) break;
      p = 0.;
      q = 0.;
      r = 0.;
      if (fabs(e) > tol1) {
                            r = (x - w) * (fx - fv);
                            q = (x - v) * (fx - fw);
                            p = (x - v) * q - (x - w) * r;
                            q = (q - r) * 2.;
                            if (q > 0.) p = -p; else q = -q;
                            r = e;
                            e = d;
      }

      if (fabs(p) >= fabs(q * .5 * r) ||
            p <= q * (a - x) || p >= q * (b - x)) {
                                                    if (x < xm) e = b - x; else e = a - x;
                                                    d = c * e;
      }
      else {
               d = p / q;
               u = x + d;
               if (u - a < t2 || b - u < t2) {
                 d = tol1;
                 if (x >= xm) d = -d;
               }
      }
        if (fabs(d) >= tol1)
          u = x + d;
        else if (d > 0.)
          u = x + tol1;
        else
          u = x - tol1;
      fu = (*f)(genes, h, u);
        if (fu <= fx) {
          if (u < x) b = x; else a = x;
          v = w;    w = x;   x = u;
          fv = fw; fw = fx; fx = fu;
        } else {
          if (u < x) a = u; else b = u;
          if (fu <= fw || w == x) {
            v = w; fv = fw;
            w = u; fw = fu;
          } else if (fu <= fv || v == x || v == w) {
            v = u; fv = fu;
          }
        }
    }
    return x;
}




double tmin_y2(double ax, double bx, int genes, double mu, double (*f)(int, double, double), double tol)
{

    double a, b, d, e, p, q, r, u, v, w, x;
    double t2, fu, fv, fw, fx, xm, eps, tol1, tol3;
    const double c = (3. - sqrt(5.)) * .5;

    eps = DBL_EPSILON;
    tol1 = eps + 1.;
    eps = sqrt(eps);
    a = ax;
    b = bx;
    v = a + c * (b - a);
    w = v;
    x = v;
    d = 0.;
    e = 0.;
    fx = (*f)(genes, mu, x);
    fv = fx;
    fw = fx;
    tol3 = tol / 3.;
    for(;;) {
      xm = (a + b) * .5;
      tol1 = eps * fabs(x) + tol3;
      t2 = tol1 * 2.;
      if (fabs(x - xm) <= t2 - (b - a) * .5) break;
      p = 0.;
      q = 0.;
      r = 0.;
      if (fabs(e) > tol1) {
                            r = (x - w) * (fx - fv);
                            q = (x - v) * (fx - fw);
                            p = (x - v) * q - (x - w) * r;
                            q = (q - r) * 2.;
                            if (q > 0.) p = -p; else q = -q;
                            r = e;
                            e = d;
      }
      if (fabs(p) >= fabs(q * .5 * r) ||
            p <= q * (a - x) || p >= q * (b - x)) {
                                                    if (x < xm) e = b - x; else e = a - x;
                                                    d = c * e;
      }
      else {
          d = p / q;
             u = x + d;
               if (u - a < t2 || b - u < t2) {
                 d = tol1;
                 if (x >= xm) d = -d;
               }
      }
        if (fabs(d) >= tol1)
          u = x + d;
        else if (d > 0.)
          u = x + tol1;
        else
          u = x - tol1;
        fu = (*f)(genes, mu, u);
        if (fu <= fx) {
          if (u < x) b = x; else a = x;
          v = w;    w = x;   x = u;
          fv = fw; fw = fx; fx = fu;
        } else {
          if (u < x) a = u; else b = u;
          if (fu <= fw || w == x) {
            v = w; fv = fw;
            w = u; fw = fu;
          } else if (fu <= fv || v == x || v == w) {
            v = u; fv = fu;
          }
        }
    }
    return x;

}
/************************************************************************************
  * negative binomial function
*************************************************************************************/
static double NB_nt(double nt, double y, double mung, double mutg, double neta, double teta, double pi1, double pi2)
{
	 double tmp = 0.0;
	 if(nt == 0.0||nt == y) return 0.0;
	 tmp +=  gammaln(nt/pi1 + 1.0/neta) - gammaln(nt/pi1 + 1.0) - gammaln(1.0/neta);
	 tmp += -1.0/neta*log(1.0 + mung*neta) + nt/pi1*log(mung*neta/(1 + mung*neta))-log(pi1);
	 tmp += gammaln((y - nt)/(1 - pi1 - pi2) + 1.0/teta) - gammaln((y - nt)/(1 - pi1 - pi2) + 1.0) - gammaln(1.0/teta);
     tmp += -1.0/teta*log(1.0 + mutg*teta) + (y - nt)/(1 - pi1 - pi2)*log(mutg*teta/(1 + mutg*teta)) - log(1 - pi1 - pi2);
     //printf("tmp is %lf\n", tmp);
     return exp(tmp);
}

 double ft_nb(double y, double mung, double mutg, double neta, double teta, double pi1, double pi2)
 {
	 int i;
	    double rtval;
	    double tmp;
	    double tmp_pos[integ];
	    rtval=0;
	    for(i=0;i<integ;i++) tmp_pos[i] = (y)/(double)integ*(i+0.5);
        for(i=0;i<integ;i++)
        {
        	tmp = NB_nt(tmp_pos[i], y, mung, mutg, neta, teta, pi1, pi2);
            rtval += tmp;
        }
        if(rtval <= 0) rtval=1e-313;
        //printf("rtval for ft_nb is %lf\n", rtval);

        return (log(rtval/(double)integ*y));
 }


 double ft_nb2(double y, double mung1, double mung2, double mutg, double neta1, double neta2, double teta, double pi1, double pi2)
 {
   int i;
   double rtval;
   double tmp;
   double res;
   double tmp_pos[integ];
   //printf("integ is %d\n", integ);

   rtval=0;
   for(i=0;i<integ;i++) tmp_pos[i] = (y)/(double)integ*(i+0.5);

   for(i=0;i<integ;i++)
   {
     res =y-tmp_pos[i];
     tmp =  gammaln(tmp_pos[i]/pi2 + 1.0/neta2) - gammaln(tmp_pos[i]/pi2 + 1.0) - gammaln(1.0/neta2);
     //printf("tmp for ft_nb2 is %lf\n", tmp);
     tmp += -1.0/neta2*log(1.0 + mung2*neta2) + tmp_pos[i]/pi2*log(mung2*neta2/(1 + mung2*neta2))-log(pi2);
     tmp += ft_nb(res, mung1, mutg, neta1, teta, pi1, pi2);
     tmp = exp(tmp);
     rtval +=tmp;
     //printf("tmp for ft_nb2 is %lf\n", tmp);

   }

   if(rtval<=0)
     rtval=1e-313;
   //printf("rtval for ft_nb2 is %lf\n", rtval);

   return (log(rtval/(double)integ*y));
 }

/************************************************************************************
  * getpi function
*************************************************************************************/



  /*optimize pi for updating*/
  //for calculating pi1


/* For update ratio calculation, we need to fix position at the same number*/
static  double ft_y(double y, double mung, double mutg, double sng, double stg, double pi1, double pi2)
{
    int i;
    double rtval;
    double tmp;
    double tmp_pos[integ];

    rtval=0;
    for(i=0;i<integ;i++) tmp_pos[i] = (y)/(double)integ*(i+0.5);

    for(i=0;i<integ;i++)
    {
      tmp= -0.5*log(sng)-0.5*log(stg) -log(y-tmp_pos[i])-log(tmp_pos[i]);
      tmp= tmp-0.5*pow(log2(tmp_pos[i])-(log2(pi1) + mung), 2.0) /sng;
      tmp= tmp-0.5*pow(log2(y - tmp_pos[i])-(log2(1-pi1-pi2) + mutg), 2.0)/stg;

      tmp= exp(tmp);

      rtval +=tmp;
    }

    if(rtval<=0)
      rtval=1e-313;

    return (log(rtval/(double)integ *y));
  }

  
    double ft_y_SC(double y, double mung, double mutg, double sng, double stg, double pi1, double pi2)
{
    int i;
    double rtval;
    double tmp;
    double tmp_pos[integ];

    rtval=0;
    for(i=0;i<integ;i++) tmp_pos[i] = (y)/(double)integ*(i+0.5);

    for(i=0;i<integ;i++)
    {
      tmp= -log(y-tmp_pos[i])-log(tmp_pos[i]);
      tmp= tmp-0.5*pow(log2(tmp_pos[i])-(log2(pi1) + mung), 2.0) /sng;
      tmp= tmp-0.5*pow(log2(y - tmp_pos[i])-(log2(1-pi1-pi2) + mutg), 2.0)/stg;

      tmp= exp(tmp);

      rtval +=tmp;
    }

    if(rtval<=0)
      rtval=1e-313;

    return (log(rtval/(double)integ *y));
  }

double pf_y(int samp, double pi1)
{
  double new_like;
  double tmp;
  int i;
  new_like=0.0;
  if(opt_seq == 1){//microarray
    for (i=0; i<nG; i++){
       tmp = ft_y(FD[samp+fNorm][i],p->Navg1[i], p->Tavg[i], p->Nsigma1[i], p->Tsigma[i], pi1, 0);
       new_like+= tmp;
    }
  }else{
    for (i=0; i<nG; i++){
        tmp = ft_nb(FD[samp+fNorm][i],p->Navg1[i], p->Tavg[i], p->Nsigma1[i], p->Tsigma[i], pi1, 0);
        new_like+= tmp;
      }
  }
  return (new_like*(-1.0));
}


/*Codes used for our functions pmin_y, pmin_y2, tmin_y, tmin_y2, fmin1, fmin2, which perform golden section
 search with parabolic interpolation, are modified from code of R base function,
 which references a Fortran code http://www.netlib.org/fmm/fmin.f*/

double pmin_y(double ax, double bx, int samp, double (*f)(int, double), double tol)
{

    double a, b, d, e, p, q, r, u, v, w, x;
    double t2, fu, fv, fw, fx, xm, eps, tol1, tol3;
    const double c = (3. - sqrt(5.)) * .5;

    eps = DBL_EPSILON;
    tol1 = eps + 1.;
    eps = sqrt(eps);
    a = ax;
    b = bx;
    v = a + c * (b - a);
    w = v;
    x = v;
    d = 0.;
    e = 0.;
    fx = (*f)(samp, x);
    fv = fx;
    fw = fx;
    tol3 = tol / 3.;
    for(;;) {
      xm = (a + b) * .5;
      tol1 = eps * fabs(x) + tol3;
      t2 = tol1 * 2.;
      if (fabs(x - xm) <= t2 - (b - a) * .5) break;
      p = 0.;
      q = 0.;
      r = 0.;
      if (fabs(e) > tol1) {
                            r = (x - w) * (fx - fv);
                            q = (x - v) * (fx - fw);
                            p = (x - v) * q - (x - w) * r;
                            q = (q - r) * 2.;
                            if (q > 0.) p = -p; else q = -q;
                            r = e;
                            e = d;
      }
      if (fabs(p) >= fabs(q * .5 * r) ||
            p <= q * (a - x) || p >= q * (b - x)) {
                                                    if (x < xm) e = b - x; else e = a - x;
                                                    d = c * e;
      }
      else {
             d = p / q;
             u = x + d;
               if (u - a < t2 || b - u < t2) {
                 d = tol1;
                 if (x >= xm) d = -d;
               }
      }
        if (fabs(d) >= tol1)
          u = x + d;
        else if (d > 0.)
          u = x + tol1;
        else
          u = x - tol1;
      fu = (*f)(samp, u);
        if (fu <= fx) {
          if (u < x) b = x; else a = x;
          v = w;    w = x;   x = u;
          fv = fw; fw = fx; fx = fu;
        } else {
          if (u < x) a = u; else b = u;
          if (fu <= fw || w == x) {
            v = w; fv = fw;
            w = u; fw = fu;
          } else if (fu <= fv || v == x || v == w) {
            v = u; fv = fu;
          }
        }
    }
    return x;

}


////another function to get pi
static void getpi(int samp, int h)  	// option h = 1 for 1 component, 2 for two component
{
  double pii1, pii2;
  double upp;
  pii1 = 0.0;
  pii2= 0.0;

  if(h == 1)
  {
    pii1 =  pmin_y(0.01, 0.99, samp, pf_y, 0.0001);
    p->pi1[samp] = pii1;

  }else{ //two component
         pii2 = pmin_y(0.01, 0.99, samp, minpi, 0.0001);
         upp = 1 - pii2;
         pii1 = pmin_y2(0.01, upp, samp, pii2, pf_y2, 0.0001);
         p->pi1[samp] = pii1;
         p->pi2[samp] = pii2;
  }
}


////another function to get pi given piT
static void getpiT(int samp)  	// for two component given piT
{
  double pii1, pii2, piiT;
  double upp;
  pii1 = 0.0;
  pii2= 0.0;
  piiT = p->piT[samp];
  //two component
  upp = 1 - piiT;
  pii1 = pmin_y2(0.01, upp, samp, piiT, pf_yT, 0.0001);
  pii2 = 1 - piiT - pii1;
  p->pi1[samp] = pii1;
  p->pi2[samp] = pii2;

}


// outer integration
double ft_y2(double y, double mung1, double mung2, double mutg, double sng1, double sng2, double stg, double pi1, double pi2)
{
  int i;
  double rtval;
  double tmp;
  double res;
  double tmp_pos[integ];

  rtval=0;
  for(i=0;i<integ;i++) tmp_pos[i] = (y)/(double)integ*(i+0.5);

  for(i=0;i<integ;i++)
  {
    res =y-tmp_pos[i];
    tmp = -log(tmp_pos[i])-0.5*pow(log2(tmp_pos[i])-(log2(pi2) + mung2), 2.0) /sng2;
    tmp += ft_y_SC(res, mung1, mutg, sng1, stg, pi1, pi2);
    tmp = exp(tmp);
    rtval +=tmp;
  }
  rtval=rtval/sqrt(sng1*sng2*stg);

  if(rtval<=0)
    rtval=1e-313;

  return (log(rtval/(double)integ*y));
}


double pf_y2(int samp, double pi1, double pi2)
{
  double new_like;
  double tmp;
  int i;
  new_like=0.0;
      if(opt_seq == 1){//microarray
        for (i=0; i<nG; i++){
            tmp = ft_y2(FD[samp+fNorm][i], p->Navg1[i], p->Navg2[i], p->Tavg[i], p->Nsigma1[i], p->Nsigma2[i], p->Tsigma[i], pi1, pi2);
            new_like+= tmp;
        }
      }else{
        for (i=0; i<nG; i++){
        tmp = ft_nb2(FD[samp+fNorm][i], p->Navg1[i], p->Navg2[i], p->Tavg[i], p->Nsigma1[i], p->Nsigma2[i], p->Tsigma[i], pi1, pi2);
        new_like+= tmp;
      }
      }
    
  return (new_like*(-1.0));
}


//the joint likelihood function given piiT
double pf_yT(int samp, double pi1, double piT)
{
  double new_like;
  double tmp;
  double pi2;
  int i;
  new_like=0.0;
  pi2 = 1 - piT -pi1;
  if(opt_seq == 1){//microarray
     for (i=0; i<nG; i++){
         tmp = ft_y2(FD[samp+fNorm][i], p->Navg1[i], p->Navg2[i], p->Tavg[i], p->Nsigma1[i], p->Nsigma2[i], p->Tsigma[i], pi1, pi2);
         new_like+= tmp;
     }
  }else{
     for (i=0; i<nG; i++){
         tmp = ft_nb2(FD[samp+fNorm][i], p->Navg1[i], p->Navg2[i], p->Tavg[i], p->Nsigma1[i], p->Nsigma2[i], p->Tsigma[i], pi1, pi2);
         new_like+= tmp;
    }
  }
  return (new_like*(-1.0));
}




double pmin_y2(double ax, double bx, int samp, double pi2, double (*f)(int, double, double), double tol)
{

    double a, b, d, e, p, q, r, u, v, w, x;
    double t2, fu, fv, fw, fx, xm, eps, tol1, tol3;
    const double c = (3. - sqrt(5.)) * .5;

    eps = DBL_EPSILON;
    tol1 = eps + 1.;
    eps = sqrt(eps);

    a = ax;
    b = bx;
    v = a + c * (b - a);
    w = v;
    x = v;
    d = 0.;
    e = 0.;
    fx = (*f)(samp, x, pi2);
    fv = fx;
    fw = fx;
    tol3 = tol / 3.;
    for(;;) {
      xm = (a + b) * .5;
      tol1 = eps * fabs(x) + tol3;
      t2 = tol1 * 2.;
        if (fabs(x - xm) <= t2 - (b - a) * .5) break;
      p = 0.;
      q = 0.;
      r = 0.;
      if (fabs(e) > tol1) {
                            r = (x - w) * (fx - fv);
                            q = (x - v) * (fx - fw);
                            p = (x - v) * q - (x - w) * r;
                            q = (q - r) * 2.;
                            if (q > 0.) p = -p; else q = -q;
                            r = e;
                            e = d;
      }

      if (fabs(p) >= fabs(q * .5 * r) ||
            p <= q * (a - x) || p >= q * (b - x)) {
                                                    if (x < xm) e = b - x; else e = a - x;
                                                    d = c * e;
      }
      else {
               d = p / q;
               u = x + d;
               if (u - a < t2 || b - u < t2) {
                 d = tol1;
                 if (x >= xm) d = -d;
               }
      }
        if (fabs(d) >= tol1)
          u = x + d;
        else if (d > 0.)
          u = x + tol1;
        else
          u = x - tol1;

        fu = (*f)(samp, u, pi2);

        if (fu <= fx) {
          if (u < x) b = x; else a = x;
          v = w;    w = x;   x = u;
          fv = fw; fw = fx; fx = fu;
        } else {
          if (u < x) a = u; else b = u;
          if (fu <= fw || w == x) {
            v = w; fv = fw;
            w = u; fw = fu;
          } else if (fu <= fv || v == x || v == w) {
            v = u; fv = fu;
          }
        }
    }
    return x;

}


//get pi1 function value given pi2

double minpi(int samp, double pi2)
{
  double pi1;
  double tmp;
  double upp = 1 - pi2;
  pi1 =  pmin_y2(0.01, upp, samp, pi2, pf_y2, 0.0001);
  tmp = pf_y2(samp, pi1, pi2);
  return (tmp);
}








/************************************************************************************
  *mle optimization function
************************************************************************************/
  //optimize function
void getmle(int opt, int opt2, int h)
{
  double tmp, tmp1, tmp2;
  double upp;
  //double tmpy;
  double y_sum, pi_sum;
  if(h == 1) //one component
{
  upp = (FD[opt+fNorm][opt2] - 1 + p->pi1[opt])/p->pi1[opt];
  y_sum = FD[opt+fNorm][opt2];
  tmp = fmin1(1, upp, opt2, opt, y_sum, 1, lf1, 0.0001);
  p->mle[opt2][opt+intx] = tmp;
  p->mle[opt2][opt] = (FD[opt+fNorm][opt2] - p->pi1[opt]*tmp )/(1 - p->pi1[opt]);
}else{
  upp = (FD[opt+fNorm][opt2] - 1 + p->pi1[opt] + p->pi2[opt])/p->pi2[opt];
 tmp2 = fmin2(1, upp, opt2, opt, lf2, 0.0001);
  p->mle[opt2][opt+2*intx] = tmp2;
  upp = (FD[opt+fNorm][opt2] - 1 + p->pi1[opt] + p->pi2[opt] - p->pi2[opt]*tmp2)/p->pi1[opt];
  y_sum = FD[opt+fNorm][opt2] - tmp2*p->pi2[opt];
  pi_sum = 1 - p->pi2[opt];
  tmp1 = fmin1(1, upp, opt2, opt, y_sum, pi_sum, lf1, 0.0001);
  p->mle[opt2][opt+intx] = tmp1;
  p->mle[opt2][opt] = (FD[opt+fNorm][opt2] - p->pi1[opt]*tmp1 - p->pi2[opt]*tmp2)/(1 - p->pi1[opt] - p->pi2[opt]);
}
}
// kernel optimize function for just component 1 normal and tumor

double lf1(int opt, int opt2,  double y_sum, double pi_sum, double nval1)
{
  double new_like;
  double tmp;
  new_like=0.0;

  new_like = -log(nval1) - 0.5*log(p->Nsigma1[opt2]) - 0.5*(log2(nval1) - p->Navg1[opt2])*(log2(nval1) - p->Navg1[opt2])/p->Nsigma1[opt2];
  tmp= (y_sum - p->pi1[opt]*nval1)/(pi_sum - p->pi1[opt]);
  new_like += -log(tmp) -0.5*log(p->Tsigma[opt2]) - 0.5*(log2(tmp) - p->Tavg[opt2])*(log2(tmp) - p->Tavg[opt2])/p->Tsigma[opt2];
  return (new_like*(-1.0));
}

double fmin1(double ax, double bx,int iG, int iS, double y_sum, double pi_sum, double (*f)(int, int, double, double, double), double tol)
{


    double a, b, d, e, p, q, r, u, v, w, x;
    double t2, fu, fv, fw, fx, xm, eps, tol1, tol3;
    const double c = (3. - sqrt(5.)) * .5;

    eps = DBL_EPSILON;
    tol1 = eps + 1.;
    eps = sqrt(eps);

    a = ax;
    b = bx;
    v = a + c * (b - a);
    w = v;
    x = v;

    d = 0.;
    e = 0.;
    fx = (*f)(iS, iG, y_sum, pi_sum, x);

    fv = fx;
    fw = fx;
    tol3 = tol / 3.;

    for(;;) {
      xm = (a + b) * .5;
      tol1 = eps * fabs(x) + tol3;
      t2 = tol1 * 2.;
        if (fabs(x - xm) <= t2 - (b - a) * .5) break;
      p = 0.;
      q = 0.;
      r = 0.;
      if (fabs(e) > tol1) {
                            r = (x - w) * (fx - fv);
                            q = (x - v) * (fx - fw);
                            p = (x - v) * q - (x - w) * r;
                            q = (q - r) * 2.;
                            if (q > 0.) p = -p; else q = -q;
                            r = e;
                            e = d;
      }

      if (fabs(p) >= fabs(q * .5 * r) ||
            p <= q * (a - x) || p >= q * (b - x)) {
                                                    if (x < xm) e = b - x; else e = a - x;
                                                    d = c * e;
      }
      else {
             d = p / q;
             u = x + d;
               if (u - a < t2 || b - u < t2) {
                 d = tol1;
                 if (x >= xm) d = -d;
               }
      }

        if (fabs(d) >= tol1)
          u = x + d;
      else if (d > 0.)
          u = x + tol1;
      else
          u = x - tol1;

      fu = (*f)(iS, iG, y_sum, pi_sum, u);

        if (fu <= fx) {
          if (u < x) b = x; else a = x;
          v = w;    w = x;   x = u;
          fv = fw; fw = fx; fx = fu;
        } else {
          if (u < x) a = u; else b = u;
          if (fu <= fw || w == x) {
            v = w; fv = fw;
            w = u; fw = fu;
          } else if (fu <= fv || v == x || v == w) {
            v = u; fv = fu;
          }
        }
    }
    return x;

}
// wrapped optimize function for component 2

double lf2(int opt, int opt2, double nval2)
{
  double tmp;
  double upp;
  double y_sum, pi_sum;
  double nval1;
  double new_like;
  new_like=0.0;
  upp = (FD[opt+fNorm][opt2] - 1 + p->pi1[opt] + p->pi2[opt] - p->pi2[opt]*nval2)/p->pi1[opt];
  y_sum = FD[opt+fNorm][opt2] - nval2*p->pi2[opt];
  pi_sum = 1 - p->pi2[opt];
  nval1 = fmin1(1, upp, opt2, opt, y_sum, pi_sum, lf1, 0.0001);

  new_like = -log(nval1) - 0.5*log(p->Nsigma1[opt2]) - 0.5*(log2(nval1) - p->Navg1[opt2])*(log2(nval1) - p->Navg1[opt2])/p->Nsigma1[opt2];
  new_like += -log(nval2) - 0.5*log(p->Nsigma2[opt2]) - 0.5*(log2(nval2) - p->Navg2[opt2])*(log2(nval2) - p->Navg2[opt2])/p->Nsigma2[opt2];
  tmp= (FD[opt+fNorm][opt2] - p->pi1[opt]*nval1 - p->pi2[opt]*nval2)/(1 - p->pi1[opt] - p->pi2[opt]);
  new_like += -log(tmp) -0.5*log(p->Tsigma[opt2]) - 0.5*(log2(tmp) - p->Tavg[opt2])*(log2(tmp) - p->Tavg[opt2])/p->Tsigma[opt2];
  return (new_like*(-1.0));
}


double fmin2(double ax, double bx,int iG, int iS, double (*f)(int, int, double), double tol)
{

    double a, b, d, e, p, q, r, u, v, w, x;
    double t2, fu, fv, fw, fx, xm, eps, tol1, tol3;
    const double c = (3. - sqrt(5.)) * .5;

    eps = DBL_EPSILON;

    tol1 = eps + 1.;
    eps = sqrt(eps);

    a = ax;
    b = bx;
    v = a + c * (b - a);
    w = v;
    x = v;

    d = 0.;
    e = 0.;
    fx = (*f)(iS, iG, x);
    fv = fx;
    fw = fx;
    tol3 = tol / 3.;

    for(;;) {
      xm = (a + b) * .5;
      tol1 = eps * fabs(x) + tol3;
      t2 = tol1 * 2.;

        if (fabs(x - xm) <= t2 - (b - a) * .5) break;
      p = 0.;
      q = 0.;
      r = 0.;
      if (fabs(e) > tol1) {
                            r = (x - w) * (fx - fv);
                            q = (x - v) * (fx - fw);
                            p = (x - v) * q - (x - w) * r;
                            q = (q - r) * 2.;
                            if (q > 0.) p = -p; else q = -q;
                            r = e;
                            e = d;
      }

      if (fabs(p) >= fabs(q * .5 * r) ||
            p <= q * (a - x) || p >= q * (b - x)) {
                                                    if (x < xm) e = b - x; else e = a - x;
                                                    d = c * e;
      }
      else {
               d = p / q;
             u = x + d;

               if (u - a < t2 || b - u < t2) {
                 d = tol1;
                 if (x >= xm) d = -d;
               }
      }
        if (fabs(d) >= tol1)
          u = x + d;
      else if (d > 0.)
        u = x + tol1;
      else
        u = x - tol1;

      fu = (*f)(iS, iG, u);

        if (fu <= fx) {
          if (u < x) b = x; else a = x;
          v = w;    w = x;   x = u;
          fv = fw; fw = fx; fx = fu;
        } else {
          if (u < x) a = u; else b = u;
          if (fu <= fw || w == x) {
            v = w; fv = fw;
            w = u; fw = fu;
          } else if (fu <= fv || v == x || v == w) {
            v = u; fv = fu;
          }
        }
    }
    return x;
}




