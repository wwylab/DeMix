
#ifndef Zeya_Wang
#define Zeya_Wang
#define max(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a > _b ? _a : _b; })
#define min(a,b) \
	({__typeof__ (a) _a = (a); \
      __typeof__ (b) _b = (b); \
      _a < _b ? _a : _b;})

// parameter definition
typedef struct paramm{
	//average parameter
	double *Navg1;
	double *Navg2;
	double *Tavg;
   //variance parameter
	double *Nsigma1;
	double *Nsigma2;
	double *Tsigma;
   //pi
	double *pi1;
	double *pi2;
	double *piT;
	//mle
	double **mle;
} PARAMM;



// Function initialize
static void initialSet(PARAMM *qq);
static void load_data(double *mat1);
static void saveFiles(double *put1, double *put3, double *put5, double *put7, double *put9, double *put11, int iteration);


//my new function


static void gettumor(int genes, int h); // option h = 1 for one component; h = 2 for two

double lf1(int opt, int opt2,  double y_sum, double pi_sum, double nval1);
double fmin1(double ax, double bx,int iG, int iS, double y_sum, double pi_sum, double (*f)(int, int, double, double, double), double tol);
double lf2(int opt, int opt2, double nval2);
double fmin2(double ax, double bx,int iG, int iS, double (*f)(int, int, double), double tol);
double pf1(double *y_sum, double pi_sum, double *nval1, double pii1, int *Gid1, int size);
double pmin1(double ax, double bx, double *y_sum, double pi_sum, double *nval1, int *Gid1, int size, double (*f)(double*, double, double*, double, int*, int), double tol);
double pf2(int opt, double *nval1, double *nval2, double pii2, int *Gid2, int size);
double pmin2(double ax, double bx,int iS, double *nval1, double *nval2, int *Gid2, int size, double (*f)(int, double*, double*, double, int*, int), double tol);
void getmle(int opt, int opt2, int h);
///new function to add
static double ft_y(double y, double mung, double mutg, double sng, double stg, double pi1, double pi2);
double ft_y_SC(double y, double mung, double mutg, double sng, double stg, double pi1, double pi2);
double ft_y2(double y, double mung1, double mung2, double mutg, double sng1, double sng2, double stg, double pi1, double pi2);
double pf_y(int samp, double pi1);
double pf_y2(int samp, double pi1, double pi2);
double pmin_y(double ax, double bx, int samp, double (*f)(int, double), double tol);
double pmin_y2(double ax, double bx, int samp, double pi2, double (*f)(int, double, double), double tol);
double minpi(int samp, double pi2);

static void getpi(int samp, int h);
static void getpiT(int samp);
static double gammaln(double xx);


double tf_y(int genes, double mu, double sigma);
double tf_y2(int genes, double mu, double sigma);
double mint(int genes, int h, double mu); // opt h =1:one component; h =2:two component
double tmin_y(double ax, double bx, int genes, int h, double (*f)(int, int, double), double tol);
double tmin_y2(double ax, double bx, int genes, double mu, double (*f)(int, double, double), double tol);
double pf_yT(int samp, double pi1, double piT);
static double NB_nt(double nt, double y, double mung, double mutg, double neta, double teta, double pi1, double pi2);
double ft_nb(double y, double mung, double mutg, double neta, double teta, double pi1, double pi2);
double ft_nb2(double y, double mung1, double mung2, double mutg, double neta1, double neta2, double teta, double pi1, double pi2);











///other functions
static double sum(double *array, int size);

// Random Number generation function

static double sd(double *, int);
static double mean(double *, int);



// Variable allocation, arguements for paramters
static PARAMM *p;
static double **FD, **CD;
static double **avgparN, **sigparN, **avgparT, **sigparT;
static double  **tmppi1, **tmppi2;

static int nP,nG, nS, nHavepi, Cid, nmle;
static double vap1, del1, vap1m, del1m;
static double pvap1, pdel1, pvap1m, pdel1m;

static int fNorm, fNorm1, fNorm2, intx;			// Number of normal samples;
static int integ; //number of integration bins
static int opt_seq; // indicator of sequencing technique
static unsigned long *seed;
static double M;




static double sum(double *array, int size)
{
	int i;
	double sum = 0.0;
	for (i=0;i<size;i++){
		sum += array[i];
	}
	return(sum);
}



static double mean(double *data, int n)
{
    return sum(data, n)/((double)n);
}



static double sd(double *data, int n)
{
    double mu=0.0, sum_deviation=0.0, tmp;
    int i;
    mu = mean(data, n);
    for(i=0; i<n;i++)
    {
    tmp = pow((data[i] - mu),2);
    sum_deviation += tmp;
    }
    return sqrt(sum_deviation/((double)n-1));
}

static double gammaln(double xx)
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


#endif





