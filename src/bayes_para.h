#ifndef _Jaeil_H
#define _Jaeil_H

#define iteration 1500
#define burn 1000
#define subSample 10
#define subSample2 100
#define _MAX_BUF_    2000


// parameter definition
typedef struct param{
	double *Navg;
	double *Tavg;
  double *Favg;
	double *Nsigma;
	double *Tsigma;
	
	double *Neta;
	double *Teta;
	double *pi;
	
} PARAM;



// Function initialize
void initialSet(PARAM *qq);
void ReSet(PARAM *qq, double *obs, unsigned long *seed);
void load_data(double *mat1, int *group,  double *fixpi, unsigned long *seed);
void saveFiles( double *put3, double *put4);
void minimax(double **mat, int nrow, int ncol, double *rtmat);
void setmin();
void getlog2(unsigned long *seed, int s);
double NB_nt(double xt, double yt, int samp, int genes);
double NB_nt2(double xt, double yt, int samp, int genes);
long double NB_nt3(double xt, double yt, int samp, int genes);
long double Poisson_nt3(double yt, double xt, int samp, int genes);
double ft_y(int ntmp, double *tmppos, int samp, int genes);
double ft_lny(int ntmp, double *tmppos, double y, double mung, double mutg, double sng, double stg, double pi, double minx);


void Set_use();
double Fake_avg(int gene);



// Functions for Bayesian
double PoiGammaPi(int Sid, int Gid, double tmppi);
double tdensity(double posi, int samp, int genes);
double btranspi(int opt, double trpi);
double transpi(int opt, double trpi);
double transpi2(int opt, int opt2,  double trpi);
double transpi3(int opt, int opt2, double trpi);



// Expression deconvolution
void GetPures();
void getpi(int opt, unsigned long *seed);
void getxi(int opt, int genes, unsigned long *seed);

double getpi3(int samp, double tmppi);
double getpi2(int samp, unsigned long *seed);


void setnormalmean();


void gettumor(int genes, int opt, unsigned long *seed);		// opt 0 for Mu_n, 1 for sigma_n
void getnormal(int genes, int opt, unsigned long *seed);		// opt 0 for Mu_n, 1 for sigma_n

// Function Min
double golden(int v3, int ftopt);
double golden3(int v3, double (*fun)(int, double));

double golden2(int v3, int v4, int which, double (*fun)(int, int, double));
void shft2(double *a, double *b, double c);
void shft3(double *a, double *b, double *c, const double d);
double MutliLike(int opt, int opt2,  double trpi);

// Random Number generation function
double kiss(unsigned long *);
int runiform_n(int,unsigned long *);
void permute_sample(double *,int,unsigned long *); 
void permute_sample_int(int *,int,unsigned long *); 
long rmultinomial2(double *,long,unsigned long *);
double runif_atob(unsigned long *,double,double);
double rexp(double,unsigned long *);
double sgamma(double,unsigned long *);
double fsign(double,double);
double snorm(unsigned long *);
double snorm2(unsigned long *);
double rnorm(double,double,unsigned long *);
double rgamma(double,double,unsigned long *);
double rinverse_gamma(double,double,unsigned long *);
double fact_ln(int);
int rpois(double,unsigned long *);
double *rdirichlet(double *,int,unsigned long *);
double rbeta(double,double,unsigned long *);
double set_initphi(int sub);

// Math Function
int fact(int xx);
int choose(int x, int y);
double qnorm(double xx);
double phinorm(double mean, double var, double opt, double min, double max);
double gammaln(double xx);
double gamma(double xx);
double factln(int n);
double factrl(int n);


// Variable allocation, arguements for paramters
PARAM *p;
double **FD, ***CD, **NIG, *minG, *fake;
int *Dgroup;
int nUseO,nG, nS, nPer, WG,nVar, nLog, nSave, nPoi, nHavepi;
int nSpace, nOpt;
int Yespi;
int fNorm;			// Number of normal samples;

int nRef, nBoot, nCalc, intx;				// Number of reference genes
double RefRatio;			// Reference genes subsampling ratio
int *resam, *use, *caluse;
double *tmp_pos, *postpi;
double *like;


int nBoot;
int Fixpi;

int integ;        // Used for the integration



unsigned long *seed;
double M;


double gammaln(double xx)
{
	double x,y,tmp,ser;
	static double cof[14]={57.1562356658629235, -59.5979603554754912, 14.1360979747417471, -0.491913816097620199,
		0.339946499848118887e-4, 0.465236289270485756e-4, -0.983744753048795646e-4, 0.158088703224912494e-3,
		-0.210264441724104883e-3, 0.217439618115212643e-3, -0.164318106536763890e-3, 0.844182239838527433e-4,
		-0.261908384015814087e-4, 0.368991826595316234e-5};
	int j;
	
	y=x=xx;
	tmp=x+5.24218750000000000;
	tmp = (x+0.5)*log(tmp) -tmp;
	ser=0.999999999999997092;
	for (j=0;j<14;j++) ser += cof[j]/++y;
	return tmp+log(2.5066282746310005*ser/x);
}

double gamma(double xx)
{
	if(xx==0.0) return 99999.0;
	
	
	return exp(gammaln(xx));
	
}

double factln(int n)
{
	static int NTOP=2000;
	static double a[2000];
	static int init=1;
	int j;
	
	if(init)
	{
		init=0;
		for(j=0;j<NTOP;j++) 
			a[j]=gammaln(j+1.);
	}
	
	if (n <= 1) return 0.0;
	if (n < NTOP) return a[n];
	else return gammaln(n+0.01);
}


double factrl(int n)
{
	static double a[171];
	int j;
	static int init=1;
	
	
	if(init>0)
	{
		init=0;
		a[0]=1.;
		for(j=1;j<171;j++) a[j]=j*a[j-1];
	}
	
	
	if (n < 0 || n>170 ) printf("Negative factorial %d in routine factrl \n", n);
	return a[n];
}





double kiss(unsigned long *seed)
/* Generator proposed by Marsaglia and Zaman, 1993. See
   Robert and Casella (1999, pages 41-43) for details.  
   Watch out: the last line
        x = ((double) (*i+*j+*k)*exp(-32*log(2.0)));
   must be calibrated, depending on the precision 
   of the computer. */
{
  seed[1] = seed[1] ^ (seed[1]<<17);
  seed[2] = (seed[2] ^ (seed[2]<<18)) & 0x7FFFFFFF;
  seed[0] = 69069*seed[0]+23606797;
  seed[1] ^= (seed[1]>>15);
  seed[2] ^= (seed[2]>>13);
  return((seed[0]+seed[1]+seed[2])*M);
}

int runiform_n(int n,unsigned long *seed)
{
  int i;

  i = (int)floor(n*kiss(seed));
  if (i == n)
    return(i-1);
  else 
    return(i);
}

void permute_sample(double *v,int len,unsigned long *seed) 
{
  int i,j;
  double x;
  
  for (i=len;i>0;i--) {
    j = runiform_n(i,seed);
    x = v[j];
    v[j] = v[i-1];
    v[i-1] = x;
  }
}

void permute_sample_int(int *v,int len,unsigned long *seed) 
{
  int i,j;
  int x;
  
  for (i=len;i>0;i--) {
    j = runiform_n(i,seed);
    x = v[j];
    v[j] = v[i-1];
    v[i-1] = x;
  }
}

long rmultinomial(double *prob,long len,unsigned long *seed)
{
 long i;
 double y;

 y = kiss(seed);
 for (i=0;i<len;i++)
   if (y < prob[i])
     return i;
 return len-1;
}

long rmultinomial2(double *prob,long len,unsigned long *seed)
{
	long i;
	double y;
	double tmp=0;

	y = kiss(seed);

	tmp=prob[0];i=0;
	for(i=0;i<len;i++)
	{
		if(y<tmp)
			return i;
		else 
			tmp=tmp+prob[i+1];
	}
	return (len-1);
}




inline double runif_atob(unsigned long *seed,double a,double b)
{
 return ((b-a)*kiss(seed) + a);
}


double rexp(double beta,unsigned long *seed)
{
 return -log(kiss(seed))/beta;
}

double sgamma(double a,unsigned long* seed)
/*
**********************************************************************


     (STANDARD-)  G A M M A  DISTRIBUTION


**********************************************************************
**********************************************************************

               PARAMETER  A >= 1.0  !

**********************************************************************

     FOR DETAILS SEE:

               AHRENS, J.H. AND DIETER, U.
               GENERATING GAMMA VARIATES BY A
               MODIFIED REJECTION TECHNIQUE.
               COMM. ACM, 25,1 (JAN. 1982), 47 - 54.

     STEP NUMBERS CORRESPOND TO ALGORITHM 'GD' IN THE ABOVE PAPER
                                 (STRAIGHTFORWARD IMPLEMENTATION)

     Modified by Barry W. Brown, Feb 3, 1988 to use RANF instead of
     SUNIF.  The argument IR thus goes away.

**********************************************************************

               PARAMETER  0.0 < A < 1.0  !

**********************************************************************

     FOR DETAILS SEE:

               AHRENS, J.H. AND DIETER, U.
               COMPUTER METHODS FOR SAMPLING FROM GAMMA,
               BETA, POISSON AND BINOMIAL DISTRIBUTIONS.
               COMPUTING, 12 (1974), 223 - 246.

     (ADAPTED IMPLEMENTATION OF ALGORITHM 'GS' IN THE ABOVE PAPER)

**********************************************************************
     INPUT: A =PARAMETER (MEAN) OF THE STANDARD GAMMA DISTRIBUTION
     OUTPUT: SGAMMA = SAMPLE FROM THE GAMMA-(A)-DISTRIBUTION
     COEFFICIENTS Q(K) - FOR Q0 = SUM(Q(K)*A**(-K))
     COEFFICIENTS A(K) - FOR Q = Q0+(T*T/2)*SUM(A(K)*V**K)
     COEFFICIENTS E(K) - FOR EXP(Q)-1 = SUM(E(K)*Q**K)
     PREVIOUS A PRE-SET TO ZERO - AA IS A', AAA IS A"
     SQRT32 IS THE SQUAREROOT OF 32 = 5.656854249492380
*/
{
extern double fsign( double num, double sign );
static double q1 = 4.166669E-2;
static double q2 = 2.083148E-2;
static double q3 = 8.01191E-3;
static double q4 = 1.44121E-3;
static double q5 = -7.388E-5;
static double q6 = 2.4511E-4;
static double q7 = 2.424E-4;
static double a1 = 0.3333333;
static double a2 = -0.250003;
static double a3 = 0.2000062;
static double a4 = -0.1662921;
static double a5 = 0.1423657;
static double a6 = -0.1367177;
static double a7 = 0.1233795;
static double e1 = 1.0;
static double e2 = 0.4999897;
static double e3 = 0.166829;
static double e4 = 4.07753E-2;
static double e5 = 1.0293E-2;
static double aa = 0.0;
static double aaa = 0.0;
static double sqrt32 = 5.656854;
static double sgamma,s2,s,d,t,x,u,r,q0,b,si,c,v,q,e,w,p;

    if(a == aa) goto S10;
    if(a < 1.0) goto S120;
/*
     STEP  1:  RECALCULATIONS OF S2,S,D IF A HAS CHANGED
*/
    aa = a;
    s2 = a-0.5;
    s = sqrt(s2);
    d = sqrt32-12.0*s;
S10:
/*
     STEP  2:  T=STANDARD NORMAL DEVIATE,
               X=(S,1/2)-NORMAL DEVIATE.
               IMMEDIATE ACCEPTANCE (I)
*/
    t = snorm(seed);
    x = s+0.5*t;
    sgamma = x*x;
    if(t >= 0.0) return sgamma;
/*
     STEP  3:  U= 0,1 -UNIFORM SAMPLE. SQUEEZE ACCEPTANCE (S)
*/
    u = kiss(seed);
    if(d*u <= t*t*t) return sgamma;
/*
     STEP  4:  RECALCULATIONS OF Q0,B,SI,C IF NECESSARY
*/
    if(a == aaa) goto S40;
    aaa = a;
    r = 1.0/ a;
    q0 = ((((((q7*r+q6)*r+q5)*r+q4)*r+q3)*r+q2)*r+q1)*r;
/*
               APPROXIMATION DEPENDING ON SIZE OF PARAMETER A
               THE CONSTANTS IN THE EXPRESSIONS FOR B, SI AND
               C WERE ESTABLISHED BY NUMERICAL EXPERIMENTS
*/
    if(a <= 3.686) goto S30;
    if(a <= 13.022) goto S20;
/*
               CASE 3:  A .GT. 13.022
*/
    b = 1.77;
    si = 0.75;
    c = 0.1515/s;
    goto S40;
S20:
/*
               CASE 2:  3.686 .LT. A .LE. 13.022
*/
    b = 1.654+7.6E-3*s2;
    si = 1.68/s+0.275;
    c = 6.2E-2/s+2.4E-2;
    goto S40;
S30:
/*
               CASE 1:  A .LE. 3.686
*/
    b = 0.463+s+0.178*s2;
    si = 1.235;
    c = 0.195/s-7.9E-2+1.6E-1*s;
S40:
/*
     STEP  5:  NO QUOTIENT TEST IF X NOT POSITIVE
*/
    if(x <= 0.0) goto S70;
/*
     STEP  6:  CALCULATION OF V AND QUOTIENT Q
*/
    v = t/(s+s);
    if(fabs(v) <= 0.25) goto S50;
    q = q0-s*t+0.25*t*t+(s2+s2)*log(1.0+v);
    goto S60;
S50:
    q = q0+0.5*t*t*((((((a7*v+a6)*v+a5)*v+a4)*v+a3)*v+a2)*v+a1)*v;
S60:
/*
     STEP  7:  QUOTIENT ACCEPTANCE (Q)
*/
    if(log(1.0-u) <= q) return sgamma;
S70:
/*
     STEP  8:  E=STANDARD EXPONENTIAL DEVIATE
               U= 0,1 -UNIFORM DEVIATE
               T=(B,SI)-DOUBLE EXPONENTIAL (LAPLACE) SAMPLE
*/
    e = rexp(1.0,seed);
    u = kiss(seed);
    u += (u-1.0);
    t = b+fsign(si*e,u);
/*
     STEP  9:  REJECTION IF T .LT. TAU(1) = -.71874483771719
*/
    if(t < -0.7187449) goto S70;
/*
     STEP 10:  CALCULATION OF V AND QUOTIENT Q
*/
    v = t/(s+s);
    if(fabs(v) <= 0.25) goto S80;
    q = q0-s*t+0.25*t*t+(s2+s2)*log(1.0+v);
    goto S90;
S80:
    q = q0+0.5*t*t*((((((a7*v+a6)*v+a5)*v+a4)*v+a3)*v+a2)*v+a1)*v;
S90:
/*
     STEP 11:  HAT ACCEPTANCE (H) (IF Q NOT POSITIVE GO TO STEP 8)
*/
    if(q <= 0.0) goto S70;
    if(q <= 0.5) goto S100;
    w = exp(q)-1.0;
    goto S110;
S100:
    w = ((((e5*q+e4)*q+e3)*q+e2)*q+e1)*q;
S110:
/*
               IF T IS REJECTED, SAMPLE AGAIN AT STEP 8
*/
    if(c*fabs(u) > w*exp(e-0.5*t*t)) goto S70;
    x = s+0.5*t;
    sgamma = x*x;
    return sgamma;
S120:
/*
     ALTERNATE METHOD FOR PARAMETERS A BELOW 1  (.3678794=EXP(-1.))
*/
    aa = 0.0;
    b = 1.0+0.3678794*a;
S130:
    p = b*kiss(seed);
    if(p >= 1.0) goto S140;
    sgamma = exp(log(p)/ a);
    if(rexp(1.0,seed) < sgamma) goto S130;
    return sgamma;
S140:
    sgamma = -log((b-p)/ a);
    if(rexp(1.0,seed) < (1.0-a)*log(sgamma)) goto S130;
    return sgamma;
}

double fsign( double num, double sign )
/* Transfers sign of argument sign to argument num */
{
if ( ( sign>0.0f && num<0.0f ) || ( sign<0.0f && num>0.0f ) )
    return -num;
else return num;
}

double snorm(unsigned long *seed)
/*
 **********************************************************************
 
 
 (STANDARD-)  N O R M A L  DISTRIBUTION
 
 
 **********************************************************************
 **********************************************************************
 
 FOR DETAILS SEE:
 
 AHRENS, J.H. AND DIETER, U.
 EXTENSIONS OF FORSYTHE'S METHOD FOR RANDOM
 SAMPLING FROM THE NORMAL DISTRIBUTION.
 MATH. COMPUT., 27,124 (OCT. 1973), 927 - 937.
 
 ALL STATEMENT NUMBERS CORRESPOND TO THE STEPS OF ALGORITHM 'FL'
 (M=5) IN THE ABOVE PAPER     (SLIGHTLY MODIFIED IMPLEMENTATION)
 
 Modified by Barry W. Brown, Feb 3, 1988 to use RANF instead of
 SUNIF.  The argument IR thus goes away.
 
 **********************************************************************
 THE DEFINITIONS OF THE CONSTANTS A(K), D(K), T(K) AND
 H(K) ARE ACCORDING TO THE ABOVEMENTIONED ARTICLE
 */
{
static double a[33] = {0,
		0.00000000000000,0.03917608550309,0.07841241273311,0.11776987457909,
		0.15731068461017,0.19709908429430,0.23720210932878,0.27769043982157,
		0.31863936396437,0.36012989178957,0.40225006532172,0.44509652498551,
		0.48877641111466,0.53340970624127,0.57913216225555,0.62609901234641,
		0.67448975019607,0.72451438349236,0.77642176114792,0.83051087820539,
		0.88714655901887,0.94678175630104,1.00999016924958,1.07751556704027,
		1.15034938037600,1.22985875921658,1.31801089730353,1.41779713799625,
		1.53412054435253,1.67593972277344,1.86273186742164,2.15387469406144
	};
	static double d[32] = {0,
		0.67448975019607,0.47585963017993,0.38377116397654,0.32861132306910,
		0.29114282663980,0.26368432217502,0.24250845238097,0.22556744380930,
		0.21163416577204,0.19992426749317,0.18991075842246,0.18122518100691,
		0.17360140038056,0.16684190866667,0.16079672918053,0.15534971747692,
		0.15040938382813,0.14590257684509,0.14177003276856,0.13796317369537,
		0.13444176150074,0.13117215026483,0.12812596512583,0.12527909006226,
		0.12261088288608,0.12010355965651,0.11774170701949,0.11551189226063,
		0.11340234879117,0.11140272044119,0.10950385201710
	};
	static double t[32] = {0,
		0.00076738283767,0.00230687039764,0.00386061844387,0.00543845406707,
		0.00705069876857,0.00870839582019,0.01042356984914,0.01220953194966,
		0.01408124734637,0.01605578804548,0.01815290075142,0.02039573175398,
		0.02281176732513,0.02543407332319,0.02830295595118,0.03146822492920,
		0.03499233438388,0.03895482964836,0.04345878381672,0.04864034918076,
		0.05468333844273,0.06184222395816,0.07047982761667,0.08113194985866,
		0.09462443534514,0.11230007889456,0.13649799954975,0.17168856004707,
		0.22762405488269,0.33049802776911,0.58470309390507
	};
	static double h[33] = {0,
		0.03920617164634,0.03932704963665,0.03950999486086,0.03975702679515,
		0.04007092772490,0.04045532602655,0.04091480886081,0.04145507115859,
		0.04208311051344,0.04280748137995,0.04363862733472,0.04458931789605,
		0.04567522779560,0.04691571371696,0.04833486978119,0.04996298427702,
		0.05183858644724,0.05401138183398,0.05654656186515,0.05953130423884,
		0.06308488965373,0.06737503494905,0.07264543556657,0.07926471414968,
		0.08781922325338,0.09930398323927,0.09930398323927,0.14043438342816,
		0.18361418337460,0.27900163464163,0.70104742502766
	};

	static long i;
	static double snorm,u,s,ustar,aa,w,y,tt;
    u = 0; 
	
    while (u <= 1e-37) 
		u = kiss(seed);
    s = 0.0;
    if(u > 0.5) s = 1.0;
    u += (u-s);
    u = 32.0*u;
    i = (long)u;
    /*    if(i == 32) i = 31;*/

    if (i>31) i=31;
    if(i == 0) goto S100;
	/*
	 START CENTER
	 */
    ustar = u-(double)i;
    aa = *(a+i-1);
S40:
    if(ustar <= *(t+i-1)) goto S60;
    w = (ustar-*(t+i-1))**(h+i-1);
S50:
	/*
	 EXIT   (BOTH CASES)
	 */
    y = aa+w;
    snorm = y;
    if(s == 1.0) snorm = -y;
    return snorm;
S60:
	/*
	 CENTER CONTINUED
	 */
    u = kiss(seed);
    w = u*(*(a+i)-aa);
    tt = (0.5*w+aa)*w;
    goto S80;
S70:
    tt = u;
    ustar = kiss(seed);
S80:
    if(ustar > tt) goto S50;
    u = kiss(seed);
    if(ustar >= u) goto S70;
    ustar = kiss(seed);
    goto S40;
S100:
	/*
	 START TAIL
	 */
    i = 6;
    aa = *(a+31);
    goto S120;
S110:
    aa += *(d+i-1);
    i += 1;
	
S120:
    u += u;
    if(u < 1.0) goto S110;
    u -= 1.0;
S140:
    w = u**(d+i-1);
    tt = (0.5*w+aa)*w;
    goto S160;
S150:
    tt = u;
S160:
    ustar = kiss(seed);
    if(ustar > tt) goto S50;
    u = kiss(seed);
    if(ustar >= u) goto S150;
    u = kiss(seed);
    goto S140;
}

double snorm2(unsigned long *seed)
/*
 **********************************************************************
 
 
 (STANDARD-)  N O R M A L  DISTRIBUTION
 
 
 **********************************************************************
 **********************************************************************
 
 FOR DETAILS SEE:
 
 AHRENS, J.H. AND DIETER, U.
 EXTENSIONS OF FORSYTHE'S METHOD FOR RANDOM
 SAMPLING FROM THE NORMAL DISTRIBUTION.
 MATH. COMPUT., 27,124 (OCT. 1973), 927 - 937.
 
 ALL STATEMENT NUMBERS CORRESPOND TO THE STEPS OF ALGORITHM 'FL'
 (M=5) IN THE ABOVE PAPER     (SLIGHTLY MODIFIED IMPLEMENTATION)
 
 Modified by Barry W. Brown, Feb 3, 1988 to use RANF instead of
 SUNIF.  The argument IR thus goes away.
 
 **********************************************************************
 THE DEFINITIONS OF THE CONSTANTS A(K), D(K), T(K) AND
 H(K) ARE ACCORDING TO THE ABOVEMENTIONED ARTICLE
 */
{
	static double a[32] = {
		0.0,3.917609E-2,7.841241E-2,0.11777,0.1573107,0.1970991,0.2372021,0.2776904,
		0.3186394,0.36013,0.4022501,0.4450965,0.4887764,0.5334097,0.5791322,
		0.626099,0.6744898,0.7245144,0.7764218,0.8305109,0.8871466,0.9467818,
		1.00999,1.077516,1.150349,1.229859,1.318011,1.417797,1.534121,1.67594,
		1.862732,2.153875
	};
	static double d[31] = {
		0.0,0.0,0.0,0.0,0.0,0.2636843,0.2425085,0.2255674,0.2116342,0.1999243,
		0.1899108,0.1812252,0.1736014,0.1668419,0.1607967,0.1553497,0.1504094,
		0.1459026,0.14177,0.1379632,0.1344418,0.1311722,0.128126,0.1252791,
		0.1226109,0.1201036,0.1177417,0.1155119,0.1134023,0.1114027,0.1095039
	};
	static double t[31] = {
		7.673828E-4,2.30687E-3,3.860618E-3,5.438454E-3,7.0507E-3,8.708396E-3,
		1.042357E-2,1.220953E-2,1.408125E-2,1.605579E-2,1.81529E-2,2.039573E-2,
		2.281177E-2,2.543407E-2,2.830296E-2,3.146822E-2,3.499233E-2,3.895483E-2,
		4.345878E-2,4.864035E-2,5.468334E-2,6.184222E-2,7.047983E-2,8.113195E-2,
		9.462444E-2,0.1123001,0.136498,0.1716886,0.2276241,0.330498,0.5847031
	};
	static double h[31] = {
		3.920617E-2,3.932705E-2,3.951E-2,3.975703E-2,4.007093E-2,4.045533E-2,
		4.091481E-2,4.145507E-2,4.208311E-2,4.280748E-2,4.363863E-2,4.458932E-2,
		4.567523E-2,4.691571E-2,4.833487E-2,4.996298E-2,5.183859E-2,5.401138E-2,
		5.654656E-2,5.95313E-2,6.308489E-2,6.737503E-2,7.264544E-2,7.926471E-2,
		8.781922E-2,9.930398E-2,0.11556,0.1404344,0.1836142,0.2790016,0.7010474
	};
	int flag = 1;
	static long i;
	static double snorm,u,s,ustar,aa,w,y,tt;
    u = 0; 
	
    while (u <= 1e-37) 
		u = kiss(seed);
    s = 0.0;
    if(u > 0.5)
		s = 1.0;
    u += (u-s);
    u = 32.0*u;
    i = (long)u;
    /*    if(i == 32) i = 31;*/
    if (i>31) i=31;
    if(i == 0) { //start tail
		i = 6;
		aa = *(a+31);
		u += u;
		while(u < 1.0) {
			aa += *(d+i-1);
			i += 1;
			u += u;
		}
		u -= 1.0;
		w = u*(*(d+i-1));
		tt = (0.5*w+aa)*w;
		ustar = kiss(seed);
		while (!(ustar > tt)) {
			u = kiss(seed);
			if (!(ustar < u)) {
				tt = u;
				ustar = kiss(seed);
			}
			else {
				u = kiss(seed);
				w = u*(*(d+i-1));
				tt = (0.5*w + aa)*w;
				ustar = kiss(seed);
			}
		}
		y = aa + w;
		snorm = y;
		if (s == 1.0)
			snorm = -y;
		return (snorm);
	}
	else { //start center
//		printf("center ");
		ustar = u-(double)i;
		aa = *(a+i-1);
		while (ustar <= *(t+i-1)) {
			u = kiss(seed);
			w = u*(*(a+i)-aa);
			tt = (0.5*w + aa)*w;
			while (!(ustar > tt)) {
				u = kiss(seed);
				if (!(ustar < u)) {
					tt = u;
					ustar = kiss(seed);
				}
				else { // (ustar < 0)
					ustar = kiss(seed);
					flag = 0;
					break;
				}
			}
			if (flag) {
				y = aa+w;
				snorm = y;
				if(s == 1.0) 
					snorm = -y;
				return snorm;
			}
			else
				flag = 1;
		}
		w = (ustar-*(t+i-1))**(h+i-1);
		y = aa+w;
		snorm = y;
		if(s == 1.0) 
			snorm = -y;
		return snorm;
		
	}
}

double rnorm(double mean,double stdev,unsigned long *seed)
{
 return mean + stdev*snorm(seed);
}

double rgamma(double alpha,double beta,unsigned long *seed)
{
 return sgamma(alpha,seed)/beta;
}

double rinverse_gamma(double alpha,double beta,unsigned long *seed)
{
 return beta/sgamma(alpha,seed);
}

double fact_ln(int k)
{
  int i;
  double fact;

  if (k == 0) return 0;
  fact = 0.0;
  for (i=1;i<=k;i++)
    fact += log((double)i);
  return fact;
}

int rpois(double lambda,unsigned long *seed)
{
 int i;
 double U,P;

 U = kiss(seed);
 i = 0;

 P = exp(-lambda);
 while (P <= U) {
   i++;
   P += exp(-lambda + i * log(lambda) - fact_ln(i));
 }
 return i;
}

double *rdirichlet(double *alpha,int len,unsigned long *seed)
{
 int i;
 double *theta,denom=0.0;

 theta = (double *)calloc(len,sizeof(double));
 for (i=0;i<len;i++) {
   theta[i] = sgamma(alpha[i],seed);
   denom += theta[i];
 }
 for (i=0;i<len;i++)
   theta[i] /= denom;

 return theta;
}

double rbeta(double alpha,double beta,unsigned long *seed)
{
 double a,b;

 a = sgamma(alpha,seed);
 b = sgamma(beta,seed);

 return a/(a+b);
}





void adjust_acceptance(double x,double *X,double rate)
{
	double y;
	
	y = 1. + 1000.*(x-rate)*(x-rate)*(x-rate);
	if (y < .9)
		y = .9;
	if (y > 1.1)
		y = 1.1;
	*X *= y;


}


/*****************************************************************************************88
* Math Function
******************************************************************************************/

int fact(int xx)
{
	int i, rtv=1;
	
	if(xx <2 ) return 1;

	for(i=xx;i>0;i--)
	{
		rtv*=i;
	}
	return rtv;
}

int choose(int x, int y)
{

	return fact(x)/(fact(y)*fact(x-y));

}
// Density value
double phinorm(double mean, double var, double opt, double min, double max)
{
	double tmp;
	
	if(opt==0)
	{
		tmp= 1.0/sqrt(2*3.141591*var) * exp(-0.5/var *pow(mean, 2));
		return tmp;
	} else{
		tmp= 1.0/sqrt(2*3.141591*var) * exp(-0.5/var *pow(mean, 2));
		tmp= tmp/(qnorm(max)-qnorm(min));
		return tmp;
	}

}


double qnorm(double xx)
{
	if(xx>6.0) return 1;
	if(xx<-6.0) return 0;

	double b1=0.31938153;
	double b2=-0.356563782;
	double b3=1.781477937;
	double b4=-1.821255978;
	double b5=1.330274429;

	double p=0.2316419;
	double c2=0.3989423;

	double a=fabs(xx);
	double t= 1.0/(1.0+a*p);
	double b=c2*exp((-xx)*(xx/2.0));
	double n = ((((b5*t+b4)*t+b3)*t+b2)*t+b1)*t;
	
	n=1.0-b*n;
	if(xx<-0.0) n=1.0-n;

	return n;

}


#endif


