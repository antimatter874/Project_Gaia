
#include <string.h>
#include <math.h>
// #include <stdio.h>
#include <sstream>
#include <iomanip>
#include <stdio.h>
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <random>
#include <thread>
#include <future>
//#include "c:/sm/mongo.h"


#define MAXTEFF 6000.0

//#define LAMOST
#define RAVE
#define OUTDIR "output/"
#define SUFFIX "distpem0t40bmr15pe10P6d10s10q40pm48pp00pf1uvs12000"

#define PAROFFSET 0.048		       //55
#define PARERRP  0.000     //40     //54
#define PARERRF  1.0        //1.08
#define RVOFFSET 0.0
#define OFFCONST 0.0
#define OFFCONE  0.0


#define PARERRMIN  0.000
#define PARERRMAX  0.040


#define PARQUALLIMIT (4.0)   //limit for accepting star par/parerr
#define PARQUALLIMITC 3.999  //must be smaller than PARQUALLIMIT

//#define SETERR
#define ERRORPLACE  0.0      //0.04        // 0.043
#define ERRFACT  1.0

#define PRIORD  0.15
#define PDREL  0.08

#define PFMAX 35
//#define REDOPRIOR   //specifically refit the prior for each subset
//#define REDOFULLPRIOR
#define SMAX 4.5



#define SAMB 12000     //30000
#define DSAM 4000   //10000


//#define USEU         //use exactly one of these three. Chooses which statistics to use, USEUV uses the full f, USEV uses f_U
//#define USEV         //uses f_V
#define USEUV


//#define USECARTVE   // cartesion velocity ellipsoid
#define USERADVE    // proper spherical velocity ellipsoid

#define ORDLIM 10.0    //original 1/par lim
#define ORDMLIM -0.3      //minimum value for 1/par


#define VARXYCOND  (dab < 10.0)                  //((fabs(varxyvepart)/(fabs(Tuw) + fabs(Tvw))) < 0.1)            // (varxyvepart
#define ADDCOND (1 == 1)                 // (fabs(zz) < 0.25)
#define ADDCOND2  (1 == 1)       //	 (starsort[ii+tt][26] > 5.0)                        //(fabs(starsort[ii+tt][24]) < 0.2)

//        double varxyvepart = Tuw*cos(galpha))*sin(2.0*gbeta) -0.5*(Tvw*sin(galpha);


#define PMERRLIM 1.0

#define BLIM   -500.0         // -500.0   //lower limit on latitude b
#define LLIM   -500.0       // -500.0   //lower limit on longitude l
#define LULIM  500.0      // 500.0    //upper limit on longitude l
#define ANGACCEPTl 100.0   //angle of acceptance in both directions, suggest 15.0
#define ANGACCEPTb 100.0   //angle of acceptance in both directions, suggest 15.0
#define VVRLIMIT  -10000.0     //minimum Vphi used
#define VVCOND (/*(fabs(eclat) < 25.0) && */(gmag > 0.1) && (bmag > 0.1) && (rmag > 0.1) && (gmag < 14.5) && ((bmag - rmag ) <= 1.5)  && /* ((bmag - rmag) <= 0.8) && */ ((dab < ORDLIM) && (dab > ORDMLIM) && (parerrread < PARERRMAX) && (parerrread > PARERRMIN)) && (RVerr < 10.0) && (fabs(gbr) > 0.1745329) && (fabs(RV) < 550.0) && (nvis > 5.0) && (excessnoise < 1.0)) //  && (fabs(zz) < 0.25))      // && (cos(glr) < 0.0))         //condition for doing statistics
#define VVCONDg ((dab < ORDLIM) && (dab > ORDMLIM) && (parerr < PARERRMAX) && (parerr > PARERRMIN))
#define WMAX   200.0
#define UUHLIM  1000000.0     // 10000000.0      //limit for deprojected U velocity, should not play a role...
#define INVERTT

//#define USEPAR 1  //useparallax




#define UCO  (ULSR - uustar)
#define VCO  (VSUN - vvstar)
#define WCO  (WLSR)

#define RVATERM  0.0                          //((WCO*sb + UCO*cb*cl) + VCO*cb*sl)
#define PMBTERM  ((1.0/(dab*facr))*(WCO*cb - UCO*sb*cl - VCO*sl*sb))
#define PMLTERM  ((1.0/(dab*facr))*(-UCO*sl + VCO*cl))

#define UCOR  (uustar)       //(ULSR)      //(uustar)         //(ULSR)           // (uustar)           // correct the hcorr velocity components
#define VCOR  (vvstar)           //(VSUN)           //(vvstar)                // (VSUN)       //(vvstar)
#define WCOR  (0.0)             //(WLSR)       // 0.0


//#define MOCK

//#define TESTVELBSTAT

#define VROT  231.0       //231.0

#define CREATENEWTABLE

#define LSTARN 2300000
#define RSTARN 600000
#define STARN 8000000


//#define PRIORALLRAVE	  // which selection function prior to use - all stars RAVE
//#define PRIORNFRAVE           // flagged RAVE stars
//#define PRIORLAMOST        //LAMOST
//#define PRIORGAIA
#define PRIORDR2     // PRIORDR2c
//#define PRIORDR2b
//#define NOSELECTION

#define ZOFFSET 0.02


//#define PRINTSTARS   //prints the stars

#define NTHREADS 40


#define TESTDISTPRIOR
#define PRIORLENGTH  1000  //integer

#define CORRERR              //correct f for the bias terms
#define SYSERR  0.09    // 0.09      // 0.09    //square of the relative uncertainty to each of the bias terms
#define SYSERRRV  0.01   // 0.01   // 0.01   //square for velocity dispersion term

//#define NOVELUSE


#define DIAL 1    // 0 1/par, 1 with prior
#define WWHLIM 200.0
#define Teffmin 3700.0
#define Teffmax 6000.0
#define loggmin -1.0
#define loggmax 4.0
#define fehmin -4.0
#define fehmax -0.0
#define distmin 0.5
#define distmax 25.0
#define ORDERING  1            // 21 par/parerr, 9 = ulennart, 15 = eclat    //9 = r, 15 = b, 16 = g examine distances ordered in the following parameter 9=feh  10=Teff    11=logg




#define PI 3.14159265358979
#define parsec 3.08567759756e+16
#define VC (4.7406*6.379*8.2 - 12.24)
#define RSUN 8.27
#define VSUN 250.0

#define VLSR 12.24        //12.24
#define ULSR 11.1       //11.1
#define WLSR 7.24


using namespace std;
using std::setprecision;


#define erfinv_a3 -0.140543331
#define erfinv_a2 0.914624893
#define erfinv_a1 -1.645349621
#define erfinv_a0 0.886226899

#define erfinv_b4 0.012229801
#define erfinv_b3 -0.329097515
#define erfinv_b2 1.442710462
#define erfinv_b1 -2.118377725
#define erfinv_b0 1

#define erfinv_c3 1.641345311
#define erfinv_c2 3.429567803
#define erfinv_c1 -1.62490649
#define erfinv_c0 -1.970840454

#define erfinv_d2 1.637067800
#define erfinv_d1 3.543889200
#define erfinv_d0 1

#ifdef WIN32
double erfinv (double x)
{
 it_fprintf (stderr, "undefined function erf()\n");
 return (NAN);
}
#else
double erfinv (double x)
{
 double x2, r, y;
 int  sign_x;

 if (x < -1 || x > 1)
   return NAN;

 if (x == 0)
   return 0;

 if (x > 0)
   sign_x = 1;
  else {
   sign_x = -1;
   x = -x;
 }

 if (x <= 0.7) {
   x2 = x * x;
   r = x * (((erfinv_a3 * x2 + erfinv_a2) * x2 + erfinv_a1) * x2 + erfinv_a0);
   r /= (((erfinv_b4 * x2 + erfinv_b3) * x2 + erfinv_b2) * x2 +
  erfinv_b1) * x2 + erfinv_b0;
 }
 else {
   y = sqrt (-log ((1 - x) / 2));
   r = (((erfinv_c3 * y + erfinv_c2) * y + erfinv_c1) * y + erfinv_c0);
   r /= ((erfinv_d2 * y + erfinv_d1) * y + erfinv_d0);
 }

  r = r * sign_x;
  x = x * sign_x;
  r -= (erf (r) - x) / (2 / sqrt (M_PI) * exp (-r * r));
  r -= (erf (r) - x) / (2 / sqrt (M_PI) * exp (-r * r));
  return r;
}
#endif

#undef erfinv_a3
#undef erfinv_a2
#undef erfinv_a1
#undef erfinv_a0

#undef erfinv_b4
#undef erfinv_b3
#undef erfinv_b2
#undef erfinv_b1
#undef erfinv_b0

#undef erfinv_c3
#undef erfinv_c2
#undef erfinv_c1
#undef erfinv_c0

#undef erfinv_d2
#undef erfinv_d1
#undef erfinv_d0



void quicksort(double arr[STARN], int arri[STARN], int left, int right) {
      int i = left, j = right;
      double tmp;
      int tmpi;
      double pivot = arr[(left + right) / 2];
      while (i <= j) {
            while (arr[i] < pivot) {
                  i++;
	    }
            while (arr[j] > pivot) {
                  j--;
	    }
            if (i <= j) {
                  tmp = arr[i];
                  arr[i] = arr[j];
                  arr[j] = tmp;
		  tmpi = arri[i];
		  arri[i] = arri[j];
		  arri[j] = tmpi;
                  i++;
                  j--;
            }
      };
//      cout << pivot << "  " << left << "  " << right << "\n";

      if (left < j) {
	quicksort(arr, arri, left, j);
      }
      if (i < right) {
	quicksort(arr, arri, i, right);
      }
}





double selfunction(double x);


double selfunction(double s) {
    double a = log(0.19);
    double b = 1.9;
   double d = 1.3/9.4;
   double e = 5.0;
//    double c = 6.2
#ifdef PRIORALLRAVE
    a = log(0.1135);
    b = 1.40218;
    d = 1.30174/9.43961;
    e = 5.06179;
#endif
#ifdef PRIORNFRAVE
    a = log(0.208731);
    b = 2.05523;
    d = 0.3286/5.75583;
    e = 1.5;
#endif
    double arg1 = (log(s) - a)*b;
    double x = exp(-(arg1*arg1*0.5)) + d*exp(-s*e);
//  double x = exp(-s/0.1);
//   if (s > 1.0) {
//  double arg11 = (log(1.0) - a)*b;
//   double x11 = exp(-(arg11*arg11*0.5)) + d*exp(-e);
//   double x12 = exp(-4.0*(s-1.0))*x11;
//   if (x > x12) {
//      x = x12;
//   }
//   }
#ifdef PRIORLAMOST
    b = 13.2056;   //a = 2593.38
    x = s*s*exp(-b*s);
#endif

#ifdef PRIORGAIA
    a = 24.0521;
    b = 5.50701;
    double c = 1.0;
    d = 47.5349;
    if (d*s < 300.0) {
	    x = a*exp(-b*s)*(1.0-c*exp(-d*s));                //x = a*exp(-b*s);
	    }
    else {
	x = a*exp(-b*s);
            }
#endif

#ifdef PRIORDR2
//fit to b > 10 deg, p/sp > 5

    a = 93.5385;
    b = 4.09583;
    double c = 0.0205;
    d = 0.28;
    double h = 5.19962;
    double j = 0.9965;
    double k2 = 3.7;
    double l = 2.5;
    double l2 = 2.85;

/*

    a = 102.237;    //          +/- 12.05        (11.79%)
b = 4.06928;  //          +/- 0.04178      (1.027%)
c               = 0.0170198;  //        +/- 0.002602     (15.29%)
j               = 0.966779;  //         +/- 0.03458      (3.577%)
h               = 4.22118;  //          +/- 0.5859       (13.88%)
l2              = 2.59594;  //          +/- 0.3373       (12.99%)
k2              = 1.98113;  //          +/- 1.933        (97.55%)
l               = 1.19384;  //          +/- 1.012        (84.77%)
d = 0.1;
*/

   a = 95.6598;
   b = 4.06271;
   c = 0.0151212;
   j = 0.971966;
   h = 5.33691;
   l2 = 2.64516;
   k2 = 3.07334;
   l = 2.5;
   d = 0.1;



a               = 96.9171;  //          +/- 1.868        (1.928%)
b               = 4.06452;  //          +/- 0.0104       (0.2559%)
c               = 0.0208475;  //        +/- 0.00235      (11.27%)
d               = 0.273609; //         +/- 0.0756       (27.63%)
l               = 1.79708;  //         +/- 0.4605       (25.62%)
l2              = 2.95633;  //          +/- 0.07223      (2.443%)
k2              = 2.63033; //          +/- 0.6542       (24.87%)
j               = 1.02999;  //          +/- 0.02321      (2.253%)
h               = 4.79127;  //          +/- 0.2468       (5.151%)


a               = 98.1719;   //        +/- 2.803        (2.855%)
b               = 4.105;       //     +/- 0.01726      (0.4205%)
c               = 0.0157616;    //    +/- 0.003893     (24.7%)
d               = 0.143119;     //    +/- 0.1486       (103.9%)
l               = 2.29157;      //    +/- 1.107        (48.33%)
l2              = 2.95772;      //    +/- 0.1004       (3.394%)
k2              = 2.22226;      //    +/- 0.8793       (39.57%)
j               = 0.966758;     //    +/- 0.05921      (6.124%)
h               = 4.45773;


a               = 97.917;    //          +/- 1.601        (1.635%)
b               = 4.10469;     //     +/- 0.01302      (0.3173%)
c               = 0.0136705;   //     +/- 0.002588     (18.93%)
d               = 0.0754458;  //      +/- 0.1095       (145.1%)
l               = 2.60013;   //       +/- 0.8208       (31.57%)
l2              = 2.9438;    //       +/- 0.07186      (2.441%)
k2              = 2.0213;    //       +/- 0.4322       (21.38%)
j               = 0.929942;  //       +/- 0.04592      (4.938%)
h               = 4.69577;   //       +/- 0.6698       (14.26%)



a               = 97.9594;  //          +/- 0.9051       (0.924%)
c               = 0.0159642;  //        +/- 0.000411     (2.575%)
h               = 4.91193;  //          +/- 0.3901       (7.943%)
j               = 0.956364; //         +/- 0.01573      (1.645%)
b               = 4.10228;  //          +/- 0.007057     (0.172%)
l2              = 2.97547; //           +/- 0.03825      (1.286%)
l               = 2.22279; //          +/- 0.3845       (17.3%)
k2              = 2.0704; //           +/- 0.165        (7.97%)
	 d = PRIORD;                      // 0.15;






    if (b*s < 300.0) {
     x =  a*(exp(-b*s) + c*exp(-d*s)*(1.0/(1.0 + exp(-h*(s-j)))))*((1.0/(0.5*PI + k2))*((atan(l*(l2 - s)) + k2)));
    }
    else {
	if (d*s < 300.0) {
	     x =  a*(c*exp(-d*s)*(1.0/(1.0 + exp(-h*(s-j)))))*((1.0/(0.5*PI + k2))*((atan(l*(l2 - s)) + k2)));
	}
	else {
		x = 0.0;
	}
    }
    if (s < 5.0) {
	x *= (1.0 - exp(-s/0.024));
    }


    if (s > 4.0) {
	    double sbs = s;
	    if (sbs > 10.0) {sbs = 10.0;}
  	    x *= exp(+(sbs-4.0)*PDREL);
    }


#endif

#ifdef PRIORDR2c
//fit to b > 10 deg, p/sp > 5

    a = 93.5385;
    b = 4.09583;
    double c = 0.0205;
    d = 0.28;
    double h = 5.19962;
    double j = 0.9965;
    double k2 = 3.7;
    double l = 2.5;
    double l2 = 2.85;
    if (s < 3.6) {
	    if (b*s < 300.0) {
	     x =  a*(exp(-b*s) + c*exp(-d*s)*(1.0/(1.0 + exp(-h*(s-j)))))*((1.0/(0.5*PI + k2))*((atan(l*(l2 - s)) + k2)));
	    }
	    else {
		if (d*s < 300.0) {
		     x =  a*(c*exp(-d*s)*(1.0/(1.0 + exp(-h*(s-j)))))*((1.0/(0.5*PI + k2))*((atan(l*(l2 - s)) + k2)));
		}
		else {
			x = 0.0;
		}
	    }
     }
     else {
	double s2 = 3.6;
	    if (b*s2 < 300.0) {
	     x =  a*(exp(-b*s2) + c*exp(-d*s2)*(1.0/(1.0 + exp(-h*(s2-j)))))*((1.0/(0.5*PI + k2))*((atan(l*(l2 - s2)) + k2)));
	    }
	    else {
		if (d*s2 < 300.0) {
		     x =  a*(c*exp(-d*s2)*(1.0/(1.0 + exp(-h*(s2-j)))))*((1.0/(0.5*PI + k2))*((atan(l*(l2 - s2)) + k2)));
		}
		else {
			x = 0.0;
		}
	    }
	double bb = 0.1;
	if (bb*(s-s2) < 300.0) {
		x *= exp(-bb*(s-s2));
	}
	else {
		x = 0.0;
	}
     }


/*






    a = 93.5385;    //          +/- 0.0739       (0.1049%)
    double c = 0.0205537;   //     +/- 0.002687     (5.905%)
    double h = 5.19962;   //          +/- 0.1306       (5.896%)
    double j = 0.9963458;  //           +/- 0.03769      (3.048%)
    b = 4.09583;    //          +/- 0.008424     (0.2024%)
    d = 0.333333;    //         +/- 0.01663      (2.624%)
    if (b*s < 300.0) {
	    x = a*(exp(-b*s) + c*exp(-d*s)*(1.0/(1.0 + exp(-h*(s-j)))));
      }
     else {
	if (d*s < 300.0) {
		x = a*(c*exp(-d*s)*(1.0/(1.0 + exp(-h*(s-j)))));
	}
	else {
		x = 0.0;
      }
}

*/


/*


    a = 1.0;
    b = 3.1499;
    double c = 0.012;
    d = 1.0/2.0;
    double l = 1.1;
    if (s < l) {
	x = a*exp(-b*s);
      }
    else {
	if (d*s < 300.0) {
		x = (a*(exp(-b*l))*exp(-d*(x-l)));
      }
	else {x = 0.0;}
}
*/

#endif


#ifdef PRIORDR2b
//fit to b > 10 deg, p/sp > 5
    a = 71.2707;    //          +/- 0.0739       (0.1049%)
    double c = 0.0676285;   //     +/- 0.002687     (5.905%)
    double h = 1.46219;   //          +/- 0.1306       (5.896%)
    double j = 1.55987;  //           +/- 0.03769      (3.048%)
    b = 4.240313;    //          +/- 0.008424     (0.2024%)
    d = 0.702286;    //         +/- 0.01663      (2.624%)
    double db = 0.4;
    if (s < 3.0) {
    if (b*s < 300.0) {
	    x = a*(exp(-b*s) + c*exp(-d*s)*(1.0/(1.0 + exp(-h*(s-j)))));
      }
     else {
	if (d*s < 300.0) {
		x = a*(c*exp(-d*s)*(1.0/(1.0 + exp(-h*(s-j)))));
	}
	else {
		x = 0.0;
      }
     }
    }
    else {
	double sb = s;
	s = 3.0;
    if (b*s < 300.0) {
	    x = a*(exp(-b*s) + c*exp(-d*s)*(1.0/(1.0 + exp(-h*(s-j)))));
      }
     else {
	if (d*s < 300.0) {
		x = a*(c*exp(-d*s)*(1.0/(1.0 + exp(-h*(s-j)))));
	}
	else {
		x = 0.0;
      }
     }
     x *= exp(-db*(sb-3.0));
     }


/*


    a = 1.0;
    b = 3.1499;
    double c = 0.012;
    d = 1.0/2.0;
    double l = 1.1;
    if (s < l) {
	x = a*exp(-b*s);
      }
    else {
	if (d*s < 300.0) {
		x = (a*(exp(-b*l))*exp(-d*(x-l)));
      }
	else {x = 0.0;}
}
*/

#endif



#ifdef NOSELECTION
    x = 1.0;
#endif
  return x;
}


double selfunctionre(double s, double *priorfit, double *priorfoot, int &itest) {
    double sele = selfunction(s);
    int i = itest;
    double s0 = priorfoot[i];
    while ((s < s0) && (i > 0)) {
        i--;
        s0 = priorfoot[i];
    }
    double s1 = priorfoot[i+1];
    while ((s > s1)) {
        i++;
        s1 = priorfoot[i+1];
    }
    itest = i;
    double fac = (s1 - s)/(s1-s0);
    return (sele*(priorfit[i] + fac*(priorfit[i+1] - priorfit[i])));
}




class rave {
private:
public:
	int *set;
	double **coords;
	double *RV; double *RVerr; double *RVsd;
	double *SN;
	double *Teff; double *Tefferr;
	double *logg; double *loggerr;
	double *met; double *meterr;
	double *fe; double *si; double *mg; double *ti; double *alpha; double *ni; double *al;
	double *chisq;
	int *algoconv;
	int *flagc0; int *flagc1; char *clustflag;
	double *TeffIR; double *TeffIRerr;
	double **TGAScoords;
	double **pm; double **pmerr;
	double *par; double *parerr;
	double *distbin; double *distbinerr;
        double *gmag;
        char **identifier;
        int *identifierlength;
        char **tgasidentifier;
        int *tgasidentifierlength;
	int *tgassum;
        rave() {
                set = new int[RSTARN];
                RV = new double[RSTARN];
		RVerr = new double[RSTARN];
		RVsd = new double[RSTARN];
		SN = new double[RSTARN];
		coords = new double*[RSTARN];
		Teff = new double[RSTARN]; Tefferr = new double[RSTARN];
		logg = new double[RSTARN]; loggerr = new double[RSTARN];
		met = new double[RSTARN]; meterr = new double[RSTARN];
		fe = new double[RSTARN]; si = new double[RSTARN]; mg = new double[RSTARN]; ti = new double[RSTARN];
		alpha = new double[RSTARN]; ni = new double[RSTARN]; al = new double[RSTARN];
		chisq = new double[RSTARN];
		algoconv = new int[RSTARN];
		flagc0 = new int[RSTARN]; flagc1 = new int[RSTARN]; clustflag = new char[RSTARN];
		TeffIR = new double[RSTARN]; TeffIRerr = new double[RSTARN];
		TGAScoords = new double*[RSTARN];
		pm = new double*[RSTARN]; pmerr = new double*[RSTARN];
		par = new double[RSTARN]; parerr = new double[RSTARN];
		distbin = new double[RSTARN]; distbinerr = new double[RSTARN];
                gmag = new double[RSTARN];
                identifierlength = new int[RSTARN];
                identifier = new char*[RSTARN];
                tgasidentifier = new char*[RSTARN];
		tgassum = new int[STARN];
                tgasidentifierlength = new int[RSTARN];
		for (int n = 0; n < RSTARN; n++) {
//                      cout << "n " << n << "\n";
                        set[n] = -1;
			tgassum[n] = 0;
			coords[n] = new double[4];
			TGAScoords[n] = new double[4];
			pm[n] = new double[4]; pmerr[n] = new double[4];
			RV[n] = -999999.99; RVerr[n] = -999999.99; RVsd[n] = -999999.99;
			SN[n] = -999999.99;
			Teff[n] = -999999.99; Tefferr[n] = -999999.99;
			logg[n] = -999999.99; loggerr[n] = -999999.99;
			met[n] = -999999.99; meterr[n] = -999999.99;
			fe[n] = -999999.99; si[n] = -999999.99; mg[n] = -999999.99; ti[n] = -999999.99;
			alpha[n] = -999999.99; ni[n] = -999999.99; al[n] = -999999.99;
			chisq[n] = -999999.99;
			algoconv[n] = 0;
			flagc0[n] = 0; flagc1[n] = 0; clustflag[n] = 0;
			TeffIR[n] = -999999.99; TeffIRerr[n] = -999999.99;
			par[n] = -999999.99; parerr[n] = -999999.99;
			distbin[n] = -999999.99; distbinerr[n] = -999999.99;
			gmag[n] = -999999.99;
			for (int ii = 0; ii < 4; ii++) {
				coords[n][ii] = -999999.99;
				TGAScoords[n][ii] = -999999.99;
				pm[n][ii] = -999999.99; pmerr[n][ii] = -999999.99;
			}
                        identifier[n] = new char[50];
			tgasidentifier[n] = new char[50];
			for (int ii = 0; ii < 50; ii++) {
                            identifier[n][ii] = ' ';
                            tgasidentifier[n][ii] = ' ';
                        }
                        identifierlength[n] = 0;
                        tgasidentifierlength[n] = 0;
			tgassum[n] = 0;
                }
	}
        int readdata();
	void resetrave(int ci);
	void copyrave(int ci, int ti);
	void copyotherrave(rave orave, int ci, int n);
	void sortrave(rave orave, int length);
};


void rave::resetrave(int n) {
//			cout << "n " << n << "\n";
	set[n] = -1;
	RV[n] = -999999.99; RVerr[n] = -999999.99; RVsd[n] = -999999.99;
	SN[n] = -999999.99;
	Teff[n] = -999999.99; Tefferr[n] = -999999.99;
	logg[n] = -999999.99; loggerr[n] = -999999.99;
	met[n] = -999999.99; meterr[n] = -999999.99;
	fe[n] = -999999.99; si[n] = -999999.99; mg[n] = -999999.99; ti[n] = -999999.99;
	alpha[n] = -999999.99; ni[n] = -999999.99; al[n] = -999999.99;
	chisq[n] = -999999.99;
	algoconv[n] = 0;
	flagc0[n] = 0; flagc1[n] = 0; clustflag[n] = 0;
	TeffIR[n] = -999999.99; TeffIRerr[n] = -999999.99;
	par[n] = -999999.99; parerr[n] = -999999.99;
	distbin[n] = -999999.99; distbinerr[n] = -999999.99;
	gmag[n] = -999999.99;
	for (int ii = 0; ii < 4; ii++) {
		coords[n][ii] = -999999.99;
		TGAScoords[n][ii] = -999999.99;
		pm[n][ii] = -999999.99; pmerr[n][ii] = -999999.99;
	}
	for (int ii = 0; ii < 50; ii++) {
	  identifier[n][ii] = ' ';
	}
	identifierlength[n] = 0;
	tgassum[n] = 0;
}


void rave::copyrave(int ci, int n) {
	set[n] = set[ci];
	RV[n] = RV[ci]; RVerr[n] = RVerr[ci]; RVsd[n] = RVsd[ci];
	SN[n] = SN[ci];
	Teff[n] = Teff[ci]; Tefferr[n] = Tefferr[ci];
	logg[n] = logg[ci]; loggerr[n] = loggerr[ci];
	met[n] = met[ci]; meterr[n] = meterr[ci];
	fe[n] = fe[ci]; si[n] = si[ci]; mg[n] = mg[ci]; ti[n] = ti[ci];
	alpha[n] = alpha[ci]; ni[n] = ni[ci]; al[n] = al[ci];
	chisq[n] = chisq[ci];
	algoconv[n] = algoconv[ci];
	flagc0[n] = flagc0[ci]; flagc1[n] = flagc1[ci]; clustflag[n] = clustflag[ci];
	TeffIR[n] = TeffIR[ci]; TeffIRerr[n] = TeffIRerr[ci];
	par[n] = par[ci]; parerr[n] = parerr[ci];
	distbin[n] = distbin[ci]; distbinerr[n] = distbinerr[ci];
	gmag[n] = gmag[ci];
	for (int ii = 0; ii < 4; ii++) {
		coords[n][ii] = coords[ci][ii];
		TGAScoords[n][ii] = TGAScoords[ci][ii];
                pm[n][ii] = pm[ci][ii]; pmerr[n][ii] = pmerr[ci][ii];
        }
        identifierlength[n] = identifierlength[ci];
        tgasidentifierlength[n] = tgasidentifierlength[ci];
        for (int ii = 0; ii < 50; ii++) {
          identifier[n][ii] = identifier[ci][ii];
          tgasidentifier[n][ii] = tgasidentifier[ci][ii];
        }
        tgassum[n] = tgassum[ci];
}

void rave::copyotherrave(rave orave, int ci, int n) {
	set[n] = orave.set[ci];
	RV[n] = orave.RV[ci]; RVerr[n] = orave.RVerr[ci]; RVsd[n] = orave.RVsd[ci];
	SN[n] = orave.SN[ci];
	Teff[n] = orave.Teff[ci]; Tefferr[n] = orave.Tefferr[ci];
	logg[n] = orave.logg[ci]; loggerr[n] = orave.loggerr[ci];
	met[n] = orave.met[ci]; meterr[n] = orave.meterr[ci];
	fe[n] = orave.fe[ci]; orave.si[n] = si[ci]; mg[n] = orave.mg[ci]; ti[n] = orave.ti[ci];
	alpha[n] = orave.alpha[ci]; ni[n] = orave.ni[ci]; al[n] = orave.al[ci];
	chisq[n] = orave.chisq[ci];
	algoconv[n] = orave.algoconv[ci];
	flagc0[n] = orave.flagc0[ci]; flagc1[n] = orave.flagc1[ci]; clustflag[n] = orave.clustflag[ci];
	TeffIR[n] = orave.TeffIR[ci]; TeffIRerr[n] = orave.TeffIRerr[ci];
	par[n] = orave.par[ci]; parerr[n] = orave.parerr[ci];
	distbin[n] = orave.distbin[ci]; distbinerr[n] = orave.distbinerr[ci];
	gmag[n] = orave.gmag[ci];
	for (int ii = 0; ii < 4; ii++) {
		coords[n][ii] = orave.coords[ci][ii];
		TGAScoords[n][ii] = orave.TGAScoords[ci][ii];
                pm[n][ii] = orave.pm[ci][ii]; pmerr[n][ii] = orave.pmerr[ci][ii];
        }
        identifierlength[n] = orave.identifierlength[ci];
        tgasidentifierlength[n] = orave.tgasidentifierlength[ci];
        for (int ii = 0; ii < 50; ii++) {
          identifier[n][ii] = orave.identifier[ci][ii];
          tgasidentifier[n][ii] = orave.tgasidentifier[ci][ii];
        }
        tgassum[n] = orave.tgassum[ci];
}



void rave::sortrave(rave inrave, int length) {
  double min = 1.0e+40;
  double max = -1.0e-40;
  int mini = -1;
  int maxi = -1;
  int minl = 0;
  double val[length];
  int vali[length];
  for (int i = 0; i < length; i++) {
    val[i] = inrave.TGAScoords[i][0];                   //.tgassum[i];
    vali[i] = i;
  }
  //  cout << "enter quicksort" << time(0) << "\n";
  quicksort(val, vali, 0, length-1);
  //  cout << "left quicksort " << time(0) << "\n";
  for (int i = 0; i < length; i++) {
    int ii = vali[i];
    copyotherrave(inrave, ii, i);
    //    cout << i << "  " << starz[i][n] << "  " << starz[i][3] << "  " << starz[i][4] << "\n";
  }
}

int rave::readdata() {
	FILE *in = fopen("RAVE_DR5.csv", "r");
	char buffer[16192];
	char **pointit = NULL;
	int counter = -1;
	int counter2 = -1;
	double radfac = PI/180.0;
	int duplicates = 0;
	int duplicates2 = 0;
	int ntgasn = 0;
	int *index = new int[RSTARN];
        int gaian = 0;
        int gaianu = 0;
	for (int ii = 0; ii < RSTARN; ii++) {
		index[ii] = ii;
	}
	while (fgets(buffer, 16192, in)) {
		int flag = 1;
		if (counter >= 0) {
			int nn = 0;
			while(buffer[nn] != ',') {nn++;}
			nn++; //RAVE_OBS_ID
			int ilen = 0;
			while(buffer[nn] != ',') {
			  identifier[counter][ilen] = buffer[nn];
			  ilen++;
			  nn++;
			}
			identifierlength[counter] = ilen;
			nn++; //HEALPix
			while(buffer[nn] != ',') {nn++;}
			nn++; //RAVEID
			while(buffer[nn] != ',') {nn++;}
			nn++; //RAdeg
			if (buffer[nn] != ',') {coords[counter][0] = strtod(&buffer[nn], pointit);}
			while(buffer[nn] != ',') {nn++;}
			nn++; //DECdeg
			if (buffer[nn] != ',') {coords[counter][1] = strtod(&buffer[nn], pointit);}
			while(buffer[nn] != ',') {nn++;}
			nn++; //glon
			if (buffer[nn] != ',') {coords[counter][2] = strtod(&buffer[nn], pointit);}
			while(buffer[nn] != ',') {nn++;}
			nn++; //glat
			if (buffer[nn] != ',') {coords[counter][3] = strtod(&buffer[nn], pointit);}
			coords[counter][2] *= radfac; coords[counter][3] *= radfac;
			while(buffer[nn] != ',') {nn++;}
			nn++; //RV
			if (buffer[nn] != ',') {RV[counter] = strtod(&buffer[nn], pointit);}
			while(buffer[nn] != ',') {nn++;}
			nn++; //RVerr
			if (buffer[nn] != ',') {RVerr[counter] = strtod(&buffer[nn], pointit);}
			if (RVerr[counter] < 1.0e-5) {RVerr[counter] = 1.0e+6;}
			while(buffer[nn] != ',') {nn++;}
			nn++; //sdRV
			if (buffer[nn] != ',') {RVsd[counter] = strtod(&buffer[nn], pointit);}
			while(buffer[nn] != ',') {nn++;}
			nn++; //MADRV
			while(buffer[nn] != ',') {nn++;}
			nn++; //SPARV
			while(buffer[nn] != ',') {nn++;}
			nn++; //SNR_K
			if (buffer[nn] != ',') {SN[counter] = strtod(&buffer[nn], pointit);}
			while(buffer[nn] != ',') {nn++;}
			nn++; //Teff_K
			if (buffer[nn] != ',') {Teff[counter] = strtod(&buffer[nn], pointit);}
			while(buffer[nn] != ',') {nn++;}
			nn++; //Teff_N_K
			while(buffer[nn] != ',') {nn++;}
			nn++; //eTeff_K
			if (buffer[nn] != ',') {Tefferr[counter] = strtod(&buffer[nn], pointit);}
			while(buffer[nn] != ',') {nn++;}
			nn++; //MAD_Teff_K
			while(buffer[nn] != ',') {nn++;}
			nn++; //std_Teff_K
			while(buffer[nn] != ',') {nn++;}
			nn++; //logg_K
			if (buffer[nn] != ',') {logg[counter] = strtod(&buffer[nn], pointit);}
			while(buffer[nn] != ',') {nn++;}
			nn++; //logg_N_K
			while(buffer[nn] != ',') {nn++;}
			nn++; //elogg_K
			if (buffer[nn] != ',') {loggerr[counter] = strtod(&buffer[nn], pointit);}
			while(buffer[nn] != ',') {nn++;}
			nn++; //MAD_logg_K
			while(buffer[nn] != ',') {nn++;}
			nn++; //StdDev_logg_K
			while(buffer[nn] != ',') {nn++;}
			nn++; //met_K
			if (buffer[nn] != ',') {met[counter] = strtod(&buffer[nn], pointit);}
			while(buffer[nn] != ',') {nn++;}
			nn++; //met_N_K
			while(buffer[nn] != ',') {nn++;}
			nn++; //emet_K
			if (buffer[nn] != ',') {meterr[counter] = strtod(&buffer[nn], pointit);}
			while(buffer[nn] != ',') {nn++;}
			nn++; //MAD_met_K
			while(buffer[nn] != ',') {nn++;}
			nn++; //StdDev_met_K
			while(buffer[nn] != ',') {nn++;}
			nn++; //chisq_K
			if (buffer[nn] != ',') {chisq[counter] = strtod(&buffer[nn], pointit);}
			while(buffer[nn] != ',') {nn++;}
			nn++; //ALGO_CONV
			algoconv[counter] = 0;
			if (buffer[nn] != '1') {algoconv[counter] = 1;}
			while(buffer[nn] != ',') {nn++;}
			nn++; //Teff_IR
			if (buffer[nn] != ',') {TeffIR[counter] = strtod(&buffer[nn], pointit);}
			while(buffer[nn] != ',') {nn++;}
			nn++; //eTeff_IR
			if (buffer[nn] != ',') {TeffIRerr[counter] = strtod(&buffer[nn], pointit);}
			while(buffer[nn] != ',') {nn++;}
			nn++; //IRdirect
			while(buffer[nn] != ',') {nn++;}
			nn++; //mgh
			if (buffer[nn] != ',') {mg[counter] = strtod(&buffer[nn], pointit);}
			while(buffer[nn] != ',') {nn++;}
			nn++; //mgn
			while(buffer[nn] != ',') {nn++;}
			nn++; //alh
			if (buffer[nn] != ',') {al[counter] = strtod(&buffer[nn], pointit);}
			while(buffer[nn] != ',') {nn++;}
			nn++; //aln
			while(buffer[nn] != ',') {nn++;}
			nn++; //sih
			if (buffer[nn] != ',') {si[counter] = strtod(&buffer[nn], pointit);}
			while(buffer[nn] != ',') {nn++;}
			nn++; //sin
			while(buffer[nn] != ',') {nn++;}
			nn++; //tih
			if (buffer[nn] != ',') {ti[counter] = strtod(&buffer[nn], pointit);}
			while(buffer[nn] != ',') {nn++;}
			nn++; //tin
			while(buffer[nn] != ',') {nn++;}
			nn++; //feh
			if (buffer[nn] != ',') {fe[counter] = strtod(&buffer[nn], pointit);}
			while(buffer[nn] != ',') {nn++;}
			nn++; //fen
			while(buffer[nn] != ',') {nn++;}
			nn++; //nih
			if (buffer[nn] != ',') {ni[counter] = strtod(&buffer[nn], pointit);}
			while(buffer[nn] != ',') {nn++;}
			nn++; //nin
			while(buffer[nn] != ',') {nn++;}
			nn++; //alpha
			if (buffer[nn] != ',') {alpha[counter] = strtod(&buffer[nn], pointit);}
			for (int ii = 0; ii < 4; ii++) { //chisqc, fracc, avS, binneyd
//				cout << "|" << buffer[nn];
				while(buffer[nn] != ',') {nn++;}
				nn++;
			} //binney dist
			if (buffer[nn] != ',') {distbin[counter] = strtod(&buffer[nn], pointit);}
			while(buffer[nn] != ',') {nn++;}
			nn++;  //binney dist err
			if (buffer[nn] != ',') {distbinerr[counter] = strtod(&buffer[nn], pointit);}
			for (int ii = 0; ii < 18; ii++) {
				while(buffer[nn] != ',') {nn++;}
				nn++;
			} //flags
			int hflagc0 = 0;
			int hflagc1 = 0;
//			cout << buffer[nn-2] << buffer[nn-1] << buffer[nn] << buffer[nn+1] << buffer[nn+2] << buffer[nn+3]<< buffer[nn+4] << "\n";
//			cout << "d" << distbin[counter] << " " << distbinerr[counter] << " " << fe[counter] << " " << mg[counter] << "|";
			for (int ii = 0; ii < 20; ii++) {
				while(buffer[nn] != ',') {nn++;}
				nn++;
//			cout << ii << "  " << buffer[nn] << buffer[nn+1] << buffer[nn+2] << buffer[nn+3]<< buffer[nn+4] << buffer[nn+5] << "\n";

				if (buffer[nn] == '\"') {nn++;}
//				cout << buffer[nn];
				if (buffer[nn] != 'n') {
					hflagc0++;
					if (buffer[nn] == 'b') {
						hflagc1++;
					}
				}
			}  //flags still
			for (int ii = 0; ii < 2; ii++) {
				while (buffer[nn] != ',') {nn++;}
				nn++;
			} //clustflag
			char hclustflag = buffer[nn];
			if (hclustflag == ',' || hclustflag == '\"') {hclustflag = '0';}
                        clustflag[counter] = hclustflag;
                        flagc0[counter] = hflagc0;
                        flagc1[counter] = hflagc1;
                        for (int ii = 0; ii < 2; ii++) {
                                while (buffer[nn] != ',') {nn++;}
                                nn++;
                        } //ID_TGAS_source                 RA_TGAS
                        int tgasilen = 0;
			tgassum[counter] = 0;
                        if (buffer[nn] != ',') {
                          while(buffer[nn] != ',') {
                          tgasidentifier[counter][tgasilen] = buffer[nn];
                          tgassum[counter] += (int)(buffer[nn]);
			  tgasilen++;
                          nn++;
			  }
                        }
                        tgasidentifierlength[counter] = tgasilen;
			if (tgasilen > 2) {
				gaian++;
				gaianu++;
			}
                        nn++; // tgasidflag
//                      if (buffer[nn] != ',') {TGAScoords[counter][0] = strtod(&buffer[nn], pointit);}
                        while (buffer[nn] != ',') {nn++;}
                        nn++; //RA_TGAS
                        if (buffer[nn] != ',') {TGAScoords[counter][0] = strtod(&buffer[nn], pointit);}
                        while (buffer[nn] != ',') {nn++;}
                        nn++; //DEC_Tgas
			if (buffer[nn] != ',') {TGAScoords[counter][1] = strtod(&buffer[nn], pointit);}
			while (buffer[nn] != ',') {nn++;}
			nn++; //pmra
			if (buffer[nn] != ',') {pm[counter][0] = strtod(&buffer[nn], pointit);}
			while (buffer[nn] != ',') {nn++;}
			nn++; //pmraerr
			if (buffer[nn] != ',') {pmerr[counter][0] = strtod(&buffer[nn], pointit);}
			while (buffer[nn] != ',') {nn++;}
			nn++; //pmdec
			if (buffer[nn] != ',') {pm[counter][1] = strtod(&buffer[nn], pointit);}
			while (buffer[nn] != ',') {nn++;}
			nn++; //pmdecerr
			if (buffer[nn] != ',') {pmerr[counter][1] = strtod(&buffer[nn], pointit);}
			while (buffer[nn] != ',') {nn++;}
			nn++; //pi_Tgas
			if (buffer[nn] != ',') {par[counter] = strtod(&buffer[nn], pointit);}
			while (buffer[nn] != ',') {nn++;}
			nn++; //pi_Tgas
			if (buffer[nn] != ',') {parerr[counter] = strtod(&buffer[nn], pointit);}
			while (buffer[nn] != ',') {nn++;}
			nn++; //gmag
			if (buffer[nn] != ',') {gmag[counter] = strtod(&buffer[nn], pointit);}
			if (counter >= 1) {
			int clength = counter/2;
			int cl0 = 0;
			int cl0l = 0;
			int cl0u = counter -1;
			int cl0m = (cl0u-cl0l)/2;
			while (cl0u > cl0l+1) {
				cl0m = (cl0u + cl0l)/2;
				if (coords[index[cl0m]][0] <= coords[index[counter]][0]) {
					cl0l = cl0m;
				}
				else {
					cl0u = cl0m;
				}
//				cout << cl0l << " " << cl0u << " " << coords[index[cl0l]][0] << " " << coords[index[counter]][0] << "\n";
			}
			cl0 = cl0l;
			int cflag = 0;
			int ccflag = 0;
			while (cflag == 0) {
				if (coords[index[cl0]][0] <= coords[index[counter]][0]) {
					cl0++;
					ccflag = 1;
				}
				else {
					cflag = 1;
				}
				if (cl0 == counter) {
					cflag = 2;
				}
			}
			cl0--;
			if (counter < 100) {
				for (int jii = 0; jii < counter; jii++) {
					cout << jii << " " << cl0 << " " << coords[index[jii]][0] << " " << coords[index[counter]][0] << "\n";
				}
			}
			int dupflag = 0;
			if (cl0 >= 0) {
		if ((fabs(coords[index[counter]][0] - coords[index[cl0]][0]) < 1.0e-3)) {
				int jj = cl0;
				if ((fabs(coords[index[counter]][1] - coords[index[cl0]][1]) < 1.0e-3) && (fabs(pm[index[jj]][0] - pm[index[counter]][0]) < 0.1) && (fabs(pm[index[jj]][1] - pm[index[counter]][1]) < 0.1)) {
					duplicates++;
					dupflag = 1;
					if (TGAScoords[index[counter]][0] > -1000.0) {
						duplicates2++;
					}
					if (RVerr[index[counter]] >= RVerr[index[cl0]]) {
						cout << "1 found duplicate " << cl0 << " " << counter << " " << RV[index[jj]] << " " << RV[index[counter]] << " " << RVerr[index[jj]] << " " << RVerr[index[counter]] << " " << coords[index[jj]][0] << " " << coords[index[jj]][1] << " " << coords[index[counter]][0] << " " << coords[index[counter]][1] << " " << clustflag[index[jj]] << " " << flagc0[index[jj]] << "\n";
						resetrave(index[counter]);
						counter--;
						if (tgasilen > 2) {
							gaianu--;
						}
					}
					else {
						cout << "2 found duplicate " << cl0 << " " << counter << " " << RV[index[jj]] << " " << RV[index[counter]] << " " << RVerr[index[jj]] << " " << RVerr[index[counter]] << " " << coords[index[jj]][0] << " " << coords[index[jj]][1] << " " << coords[index[counter]][0] << " " << coords[index[counter]][1] << " " << clustflag[index[jj]] << " " << flagc0[index[jj]] << "\n";
						copyrave(index[counter], index[cl0]);
						resetrave(index[counter]);
						counter--;
						if (tgasilen > 2) {
							gaianu--;
						}
					}
				}
				else {
					if (TGAScoords[index[counter]][0] > -1000.0) {
						ntgasn++;

					}
				}
			}
			} //end if cl0 > 0
			else {
				cl0 = 0;
			}
		if (cl0 < counter -1) {
		if ((fabs(coords[index[counter]][0] - coords[index[cl0]][0]) < 1.0e-3)) {
				int jj = cl0;
				if ((fabs(coords[index[counter]][1] - coords[index[cl0]][1]) < 1.0e-3) && (fabs(pm[index[jj]][0] - pm[index[counter]][0]) < 0.1) && (fabs(pm[index[jj]][1] - pm[index[counter]][1]) < 0.1)) {
					duplicates++;
					dupflag = 1;
					if (TGAScoords[index[counter]][0] > -1000.0) {
						duplicates2++;
					}
					if (RVerr[index[counter]] >= RVerr[index[cl0]]) {
						cout << "3 found duplicate " << cl0 << " " << counter << " " << RV[index[jj]] << " " << RV[index[counter]] << " " << RVerr[index[jj]] << " " << RVerr[index[counter]] << " " << coords[index[jj]][0] << " " << coords[index[jj]][1] << " " << coords[index[counter]][0] << " " << coords[index[counter]][1] << " " << clustflag[index[jj]] << " " << flagc0[index[jj]] << "\n";
						resetrave(index[counter]);
						counter--;
						if (tgasilen > 2) {
							gaianu--;
						}

					}
					else {
						cout << "4 found duplicate " << cl0 << " " << counter << " " << RV[jj] << " " << RV[counter] << " " << RVerr[jj] << " " << RVerr[counter] << " " << coords[jj][0] << " " << coords[jj][1] << " " << coords[counter][0] << " " << coords[counter][1] << " " <<clustflag[jj] << " " << flagc0[jj] << "\n";
						copyrave(index[counter], index[cl0]);
						resetrave(index[counter]);
						counter--;
						if (tgasilen > 2) {
							gaianu--;
						}

					}
				}
				else {
					if (TGAScoords[index[counter]][0] > -1000.0) {
						ntgasn++;
					}
				}
			}
		} //end if cl0
			if ((dupflag == 0) && (cflag < 2)) {
				if (ccflag == 1 && cl0 < counter) {
					cl0++;
				}
				index[counter+1] = index[counter];
				for (int jji = counter-1; jji >= cl0; jji--) {
					index[jji+1] = index[jji];
				}
				index[cl0] = index[counter+1];
				index[counter+1] = counter+1;
			}

			}  //end if counter > 1

		}
		counter2++;
		if (flag == 1) {
			counter++;
		}
	}
	cout << "rave readdata counter " << counter << " " << counter2 << " " << ntgasn << " " << duplicates <<  " " << duplicates2 << " Gaia source " << gaian << " " << gaianu << "\n";
	delete [] index;
	return counter;
}











class lamost {
private:
public:
	int *set;
	int *planid; int *spid; int *fiberid;
	double *radecobs[2];
	double *radec[2];
	double *sn[5];
	double *mag[7];
	bool *magflag[7];
	double *rv;
	double *rverr;
	double *feh; double *feherr; double *logg; double *loggerr; double *Teff; double *Tefferr;
	int *identifierlength;
	char **identifier;
	lamost() {
		set = new int[LSTARN];
		for (int i = 0; i < 7; i++) {
			mag[i] = new double[LSTARN];
			magflag[i] = new bool[LSTARN];
		}
		for (int i = 0; i < 5; i++) {
			sn[i] = new double[LSTARN];
		}
		for (int i = 0; i < 2; i++) {
			radecobs[i] = new double[LSTARN];
			radec[i] = new double[LSTARN];
		}
		planid = new int[LSTARN]; spid = new int[LSTARN]; fiberid = new int[LSTARN];
		rv = new double[LSTARN]; rverr = new double[LSTARN];
		feh = new double[LSTARN]; feherr = new double[LSTARN]; logg = new double[LSTARN]; loggerr = new double[LSTARN]; Teff = new double[LSTARN]; Tefferr = new double[LSTARN];
		identifierlength = new int[LSTARN];
		identifier = new char*[LSTARN];
	    for (int ll = 0; ll < LSTARN; ll++) {
		set[ll] = 0;
		for (int i = 0; i < 7; i++) {
			mag[i][ll] = 0.0;
			magflag[i][ll] = false;
		}
		for (int i = 0; i < 5; i++) {
			sn[i][ll] = 0.0;
		}
		for (int i = 0; i < 2; i++) {
			radecobs[i][ll] = 0.0;
			radec[i][ll] = 0.0;
		}
		planid[ll] = -1; spid[ll] = -1; fiberid[ll] = -1;
		rv[ll] = 0.0; rverr[ll] = 0.0;
		feh[ll] = 0.0; feherr[ll] = 0.0; logg[ll] = 0.0; loggerr[ll] = 0.0; Teff[ll] = 0.0; Tefferr[ll] = 0.0;
		identifierlength[ll] = 0;
		identifier[ll] = new char[50];
		for (int ii = 0; ii < 50; ii++) {
		  identifier[ll][ii] = ' ';
		}
	    }
	}
	int readdata();
};


int lamost::readdata() {
	FILE *in = fopen("dr2_stellar.csv", "r");
	char buffer[4096];
	char **pointit = NULL;
	int counter = -1;
	int counter2 = -1;
	while (fgets(buffer, 4095, in)) {
		int flag = 1;
		if (counter >= 0) {
			flag = 0;
			int nn = 0;
			while (buffer[nn] == '|') {
				nn++;
			} //obsid
			while (buffer[nn] != '|') {
				nn++;
			}
			nn++;  //designation
			int ilen = 0;
			while (buffer[nn] != '|') {
				identifier[counter2][ilen] = buffer[nn];
				ilen++;
				nn++;
			}
			identifierlength[counter2] = ilen;
			nn++;  //obsdate
			while (buffer[nn] != '|') {
				nn++;
			}
			nn++;  //lmjd
			while (buffer[nn] != '|') {
				nn++;
			}
			nn++;  //mjd
			while (buffer[nn] != '|') {
				nn++;
			}
			nn++;  //planid
			planid[counter2] = (int)(strtod(&buffer[nn], pointit));
			while (buffer[nn] != '|') {
				nn++;
			}
			nn++;  //spid
			spid[counter2] = (int)(strtod(&buffer[nn], pointit));
			while (buffer[nn] != '|') {
				nn++;
			}
			nn++;  //fiberid
			fiberid[counter2] = (int)(strtod(&buffer[nn], pointit));
			while (buffer[nn] != '|') {
				nn++;
			}
			nn++;  //RA_obs
/*
			if (counter < 10) {
				cout << "c" << counter << " ";
				for (int jj = 0; jj < 20; jj++) {
					cout << buffer[nn+jj];
				}
				cout << "\n";
			}
*/
			radecobs[0][counter2] = strtod(&buffer[nn], pointit);
			while (buffer[nn] != '|') {
				nn++;
			}
			nn++;  //DEC_obs
			radecobs[1][counter2] = strtod(&buffer[nn], pointit);
			while (buffer[nn] != '|') {
				nn++;
			}
			nn++;  //SN
			sn[0][counter2] = strtod(&buffer[nn], pointit);
			while (buffer[nn] != '|') {
				nn++;
			}
			nn++;  //SN
			sn[1][counter2] = strtod(&buffer[nn], pointit);
			while (buffer[nn] != '|') {
				nn++;
			}
			nn++;  //SN
			sn[2][counter2] = strtod(&buffer[nn], pointit);
			while (buffer[nn] != '|') {
				nn++;
			}
			nn++;  //SN
			sn[3][counter2] = strtod(&buffer[nn], pointit);
			while (buffer[nn] != '|') {
				nn++;
			}
			nn++;  //SN
			sn[4][counter2] = strtod(&buffer[nn], pointit);
			while (buffer[nn] != '|') {
				nn++;
			}
			nn++;  //objtype
			flag = 0;
/*
			if (counter < 10) {
				cout << "c" << counter << " ";
				for (int jj = 0; jj < 20; jj++) {
					cout << buffer[nn+jj];
				}
				cout << "\n";
			}
*/
			if (buffer[nn] == 'S' || buffer[nn] == 's') {
				flag = 1;
			}
			while (buffer[nn] != '|') {
				nn++;
			}
			nn++;  //class
			while (buffer[nn] != '|') {
				nn++;
			}
			nn++;  //subclass
			while (buffer[nn] != '|') {
				nn++;
			}
			nn++;  //z
			while (buffer[nn] != '|') {
				nn++;
			}
			nn++;  //zerr
			while (buffer[nn] != '|') {
				nn++;
			}
			nn++;  //magtype
			for (int j = 0; j < 7; j++) {
				while (buffer[nn] != '|') {
					nn++;
				}
				nn++;  //magtype
				mag[j][counter2] = strtod(&buffer[nn], pointit);
			}
			while (buffer[nn] != '|') {
				nn++;
			}
			nn++; //tsource	//JF_LEGAS_S
			while (buffer[nn] != '|') {
				nn++;
			}
			nn++; //fibertype ("Obj")?
			while (buffer[nn] != '|') {
				nn++;
			}
			nn++; //trom "-"
			while (buffer[nn] != '|') {
				nn++;
			}
			nn++; //tcomment "."
			while (buffer[nn] != '|') {
				nn++;
			}
			nn++; //offsets
			while (buffer[nn] != '|') {
				nn++;
			}
			nn++; //offset_v
			while (buffer[nn] != '|') {
				nn++;
			}
			nn++; //RA
			radec[0][counter2] = strtod(&buffer[nn], pointit);
			while (buffer[nn] != '|') {
				nn++;
			}
			nn++; //RA
			radec[1][counter2] = strtod(&buffer[nn], pointit);
			while (buffer[nn] != '|') {
				nn++;
			}
			nn++; //Teff
			Teff[counter2] = strtod(&buffer[nn], pointit);

			while (buffer[nn] != '|') {
				nn++;
			}
			nn++; //Tefferr
			Tefferr[counter2] = strtod(&buffer[nn], pointit);
			while (buffer[nn] != '|') {
				nn++;
			}
			nn++; //logg
			logg[counter2] = strtod(&buffer[nn], pointit);
			while (buffer[nn] != '|') {
				nn++;
			}
			nn++; //loggerr
			loggerr[counter2] = strtod(&buffer[nn], pointit);
			while (buffer[nn] != '|') {
				nn++;
			}
			nn++; //feh
			feh[counter2] = strtod(&buffer[nn], pointit);
			while (buffer[nn] != '|') {
				nn++;
			}
			nn++; //feherr
			feherr[counter2] = strtod(&buffer[nn], pointit);
			while (buffer[nn] != '|') {
				nn++;
			}
			nn++; //RV
			rv[counter2] = strtod(&buffer[nn], pointit);
			while (buffer[nn] != '|') {
				nn++;
			}
			nn++; //RVerr
			rverr[counter2] = strtod(&buffer[nn], pointit);
		}
		counter++;
		if (flag == 1) {
			counter2++;
		}
	}
	cout << "lamost.readdata counter " << counter << " " << counter2 << "\n";
	return counter;
}




class gaia {
private:
public:
	int set;
	double radec[4];
	double radecerr[4];
	double lb[4];
	double lberr[4];
	double par;
	double parerr;
	double errcorr[10];    //radec, rapar, rapmra, rapmdec, decpar, decpmra, decpmdex, parpmra, parpmdec, pmrapmdec
	char solid[100];
	char sourcid[100];
	int solidl;
	int sourcidl;
	int nobs[8];    // nobs, nobsac, nobsal, ngood, nbad
	double deltaq;
	double excessnoise;
	double excessnoisesig;
	int primaryflag;
	int relegationfactor;
        double weight[2];  //al, ac
        int photnobs;
        double phot[3]; // mean flux, flux error, mag
        char identifier[50];
        int identifierlength;
        int photvarflag;
        gaia() {
                for (int i = 0; i < 50; i++ ) {
                  identifier[i] = ' ';
                }
                identifierlength = 0;
                for (int i = 0; i < 10; i++) {
                        errcorr[i] = 0.0;
                }
		for (int i = 0; i < 8; i++) {
			nobs[i] = 0;
		}
		for (int i = 0; i < 50; i++) {
			solid[i] = ' ';
			sourcid[i] = ' ';
		}
		for (int i = 0; i < 4; i++) {
			set = 0;
			radec[i] = 0.0;
			radecerr[i] = 0.0;
			lb[i] = 0.0;
			lberr[i] = 0.0;
		}
		for (int i = 0; i < 3; i++) {
			phot[i] = 0.0;
		}
		for (int i = 0; i < 2; i++) {
			weight[i] = 0.0;
		}
		set = 0;
		par = 0.0;
		solidl = 0;
		sourcidl = 0;
		parerr = 0.0;
		deltaq = 0.0;
		excessnoise = 0.0;
		excessnoisesig = 0.0;
		primaryflag = -1;
		relegationfactor = -1.0;
		photnobs = -1;
		photvarflag = -1;
	}
	int readbuffer(char buffer[4096], char **pointit);
};


int gaia::readbuffer(char buffer[4096], char **pointit) {
	int nn = 0; //hip
	while (buffer[nn] != ',') {
		nn++;
	}
	nn++; //tycho2
	while (buffer[nn] != ',') {
		nn++;
	}
	nn++; //solution_id
	int mmi = 0;
	int mmj = 0;
	solidl = 0;
	sourcidl = 0;
	while (buffer[nn] != ',') {
		solid[mmi] = buffer[nn];
		if ((buffer[nn] != ' ') && (buffer[nn] != ',')) {mmi++; solidl++;}
                nn++;
        }
        nn++; //source_id
        while (buffer[nn] != ',') {
		sourcid[mmj] = buffer[nn];
		if ((buffer[nn] != ' ') && (buffer[nn] != ',')) {mmj++; sourcidl++;}
          nn++;
        }
        nn++; //random_index
        while (buffer[nn] != ',') {
                nn++;
	}
	nn++; //ref_epoch
	while (buffer[nn] != ',') {
		nn++;
	}
	nn++; //RA
	radec[0] = strtod(&buffer[nn], pointit);
	while (buffer[nn] != ',') {
		nn++;
	}
	nn++; //RAerr
	radecerr[0] = strtod(&buffer[nn], pointit);
	while (buffer[nn] != ',') {
		nn++;
	}
	nn++; //dec
	radec[1] = strtod(&buffer[nn], pointit);
	while (buffer[nn] != ',') {
		nn++;
	}
	nn++; //decerr
	radecerr[1] = strtod(&buffer[nn], pointit);
	while (buffer[nn] != ',') {
		nn++;
	}
	nn++; //par
	par = strtod(&buffer[nn], pointit);
	while (buffer[nn] != ',') {
		nn++;
	}
	nn++; //parerr
	parerr = strtod(&buffer[nn], pointit);
	while (buffer[nn] != ',') {
		nn++;
	}
	nn++; //pmra
	radec[2] = strtod(&buffer[nn], pointit);
	while (buffer[nn] != ',') {
		nn++;
	}
	nn++; //pmraerr
	radecerr[2] = strtod(&buffer[nn], pointit);
	while (buffer[nn] != ',') {
		nn++;
	}
	nn++; //pmdec
	radec[3] = strtod(&buffer[nn], pointit);
	while (buffer[nn] != ',') {
		nn++;
	}
	nn++; //pmdecerr
	radecerr[3] = strtod(&buffer[nn], pointit);
	for (int ji = 0; ji < 10; ji++) {
		while(buffer[nn] != ',') {
			nn++;
		}
		nn++;
		errcorr[ji] = strtod(&buffer[nn], pointit);
	}
	for (int ji = 0; ji < 6; ji++) {
		while(buffer[nn] != ',') {
			nn++;
		}
		nn++;
		nobs[ji] = (int)(strtod(&buffer[nn], pointit));
	}
	while (buffer[nn] != ',') {
		nn++;
	}
	nn++;
	deltaq = strtod(&buffer[nn], pointit);
	while (buffer[nn] != ',') {
		nn++;
	}
	nn++;
	excessnoise = strtod(&buffer[nn], pointit);
	while (buffer[nn] != ',') {
		nn++;
	}
	nn++;
	excessnoisesig = strtod(&buffer[nn], pointit);
	while (buffer[nn] != ',') {
		nn++;
	}
	nn++;
	primaryflag = (int)(strtod(&buffer[nn], pointit));
	while (buffer[nn] != ',') {
		nn++;
	}
	nn++;
	relegationfactor = strtod(&buffer[nn], pointit);
	while (buffer[nn] != ',') {
		nn++;
	}
	nn++;
	weight[0] = strtod(&buffer[nn], pointit);
	while (buffer[nn] != ',') {
		nn++;
	}
	nn++;
	weight[1] = strtod(&buffer[nn], pointit);
	while (buffer[nn] != ',') {
		nn++;
	}
	nn++; //priors used
	while (buffer[nn] != ',') {
		nn++;
	}
	nn++; //matched observations
	while (buffer[nn] != ',') {
		nn++;
	}
	nn++; //duplicated source
	while (buffer[nn] != ',') {
		nn++;
	}
	nn++; //scan-dir-strength k1
	while (buffer[nn] != ',') {
		nn++;
	}
	nn++; //scan-dir-strength k2
	while (buffer[nn] != ',') {
		nn++;
	}
	nn++; //scan-dir-strength k3
	while (buffer[nn] != ',') {
		nn++;
	}
	nn++; //scan-dir-strength k4
	while (buffer[nn] != ',') {
		nn++;
	}
	nn++; //scan-dir-mean k1
	while (buffer[nn] != ',') {
		nn++;
	}
	nn++; //scan-dir-mean k2
	while (buffer[nn] != ',') {
		nn++;
	}
	nn++; //scan-dir-mean k3
	while (buffer[nn] != ',') {
		nn++;
	}
	nn++; //scan-dir-mean k4
	while (buffer[nn] != ',') {
		nn++;
	}
	nn++; //photnobs
	photnobs = (int)(strtod(&buffer[nn], pointit));
	while (buffer[nn] != ',') {
		nn++;
	}
	nn++; //gmag mean flux
	phot[1] = strtod(&buffer[nn], pointit);
	while (buffer[nn] != ',') {
		nn++;
	}
	nn++; //gmag err
	phot[2] = strtod(&buffer[nn], pointit);
	while (buffer[nn] != ',') {
		nn++;
	}
	nn++; //gmag mean mag
	phot[0] = strtod(&buffer[nn], pointit);
	while (buffer[nn] != ',') {
		nn++;
	}
	nn++; //varflag
	while (buffer[nn] != ',') {
		nn++;
	}
	nn++; //l
	lb[0] = strtod(&buffer[nn], pointit);
	while (buffer[nn] != ',') {
		nn++;
	}
	nn++; //l
	lb[1] = strtod(&buffer[nn], pointit);
	return nn;
}



double propmotranslation(double ra, double dec, double mura, double mudec, double backback[4], int epoch) {
  double gl = 0.0;
  double gb = 0.0;
  double mugl = 0.0;
  double mugb = 0.0;
  double radians = PI/180.0;
  double rar = ra*radians;
  double decr = dec*radians;
  double alphag = 192.85948;            //193.04207;            //192.85948;
  double deltag = 27.12825;                 //27.04688;                   //27.12825;
  double lomega = 122.93192;            //32.93192;
  if (epoch == 2015) {
	alphag = 193.0420686;
	deltag = 27.0468781;
  }
  double alphagr = alphag*radians;
  double deltagr = deltag*radians;
  double lomegar = lomega*radians;
  double C1 = sin(deltagr)*cos(decr) - cos(deltagr)*sin(decr)*cos(rar-alphagr);
  double C2 = cos(deltagr)*sin(rar-alphagr);
  double cosb = sqrt(C1*C1 + C2*C2);
  mugl = (C1*mura + C2*mudec)/cosb;
  mugb = (-C2*mura + C1*mudec)/cosb;
  backback[2] = mugl;
  backback[3] = mugb;
  double sinb = cos(decr)*cos(deltagr)*cos(rar - alphagr) + sin(decr)*sin(deltagr);
  double sinl = cos(decr)*sin(rar - alphagr)/cosb;
  gb = asin(sinb);
  gl = -asin(sinl)+lomegar;
   double cosl = (sin(decr)*cos(deltagr)-cos(rar -alphagr)*sin(deltagr)*cos(decr))/cosb;
   double gl1;

   if (sinl >=0.)
     {
       gl1=acos(cosl);
     }
   if (sinl <0.)
     {
       gl1=2.*PI-acos(cosl);
     }
   gb = asin(sinb);
   gl = lomegar-gl1;
   if (gl <0.)
     {
       gl+=2.*PI;
     }
  backback[0] = gl;
  backback[1] = gb;
  return 1.0;
}





void readgaia(lamost la, int lan) {
        double facr = (parsec)*(PI)*0.001/(365.2425*24.0*3600.0*3600.0*180.0);
  double expos[200];
  FILE *fout = fopen("lamostout", "w");
  ostringstream ostr;
  double sume = 0.0;
  for (int ii = 0; ii < 200; ii++) {
	double dii = ((double)(ii-100))/25.0;
        expos[ii] = exp(-dii*dii/2.0);
	sume += expos[ii];
  }
  sume = 1.0/sume;
  for (int ii = 0; ii < 200; ii++) {
	expos[ii] *= sume;
  }
	int countthem = 0;
	for (int j = 0; j < 15; j++) {
	FILE *in;
	if (j == 0) {in = fopen("tgas/TgasSource_000-000-000.csv", "r");}
	if (j == 1) {in = fopen("tgas/TgasSource_000-000-001.csv", "r");}
	if (j == 2) {in = fopen("tgas/TgasSource_000-000-002.csv", "r");}
	if (j == 3) {in = fopen("tgas/TgasSource_000-000-003.csv", "r");}
	if (j == 4) {in = fopen("tgas/TgasSource_000-000-004.csv", "r");}
	if (j == 5) {in = fopen("tgas/TgasSource_000-000-005.csv", "r");}
	if (j == 6) {in = fopen("tgas/TgasSource_000-000-006.csv", "r");}
	if (j == 7) {in = fopen("tgas/TgasSource_000-000-007.csv", "r");}
	if (j == 8) {in = fopen("tgas/TgasSource_000-000-008.csv", "r");}
	if (j == 9) {in = fopen("tgas/TgasSource_000-000-009.csv", "r");}
	if (j == 10) {in = fopen("tgas/TgasSource_000-000-010.csv", "r");}
	if (j == 11) {in = fopen("tgas/TgasSource_000-000-011.csv", "r");}
	if (j == 12) {in = fopen("tgas/TgasSource_000-000-012.csv", "r");}
	if (j == 13) {in = fopen("tgas/TgasSource_000-000-013.csv", "r");}
	if (j == 14) {in = fopen("tgas/TgasSource_000-000-014.csv", "r");}
	if (j == 15) {in = fopen("tgas/TgasSource_000-000-015.csv", "r");}

	gaia maia;
	char buffer[4096];
	char **pointit = NULL;
	int linei = -1;
	while (fgets(buffer, 4095, in)) {

	    if (linei > -1) {
		int nn = maia.readbuffer(buffer, pointit);
		int mflag = 0;
		int jjh = 0;
		double dold = 1.0;
		double dra = 1.0;
		double ddec = 1.0;
		double dtot = 1.0;
		for (int jj = 0; jj < lan; jj++) {
			dra = (maia.radec[0] - la.radec[0][jj]) - 0.00000416666667*maia.radec[2];
			ddec = (maia.radec[1] - la.radec[1][jj]) - 0.00000416666667*maia.radec[3];
			dtot = sqrt(dra*dra + ddec*ddec);
			if (dtot < 0.0012) {
				if (dtot <= dold) {
//				cout << countthem << " " << linei << " " << jj << " " << maia.phot[0] << " " << la.mag[1][jj] << " " << maia.radec[0]-la.radec[0][jj] << " " << maia.radec[1] - la.radec[1][jj] <<  " " << (fabs((maia.radec[0] - la.radec[0][jj]) - 0.000004166666667*maia.radec[2])) << " " << (fabs((maia.radec[1] - la.radec[1][jj]) - 0.00000416666666667*maia.radec[3])) << "\n";
					mflag = 1;
					jjh = jj;
					dold = dtot;
				}
			}  // end if
		}  // end for

		if (mflag == 1) {
			dra = (maia.radec[0] - la.radec[0][jjh]) - 0.00000416666667*maia.radec[2];
			ddec = (maia.radec[1] - la.radec[1][jjh]) - 0.00000416666667*maia.radec[3];
			dtot = sqrt(dra*dra + ddec*ddec);

		double backback[4];
		double pie = maia.parerr;
		double pip = maia.par; // + erfinv(fsel)*pie;
		double RA_Tgas = maia.radec[0];
		double DEC_Tgas = maia.radec[1];
		double pmra_Tgas = maia.radec[2];
		double pmdec_Tgas = maia.radec[3];
		double pmraerr_Tgas = maia.radecerr[2];
		double pmdecerr_Tgas = maia.radecerr[3];
		double glr = maia.lb[0]*PI/180.0;
		double gbr = maia.lb[1]*PI/180.0;
		double RV = la.rv[jjh];
		double RVerr = la.rverr[jjh];
		double eRV = RVerr;
		for (int ii = 0; ii < 4; ii++) {
			backback[ii] = 0.0;
		}
		propmotranslation(RA_Tgas, DEC_Tgas, pmra_Tgas, pmdec_Tgas, backback, 2000);
//		cout << backback[0] << " " << maia.lb[0]*PI/180.0 << " "  << backback[1] << " " << maia.lb[1]*PI/180.0 << "\n";
		double dab = 0.0;
		double sumdab = 0.0;
		for (int i = 0; i < 200; i++) {
			double di = ((double)(i-100))/25.0;
			double pip0 = pip + pie*di;
			if (pip0 > 1.0e-5) {
				double dab0 = 1.0/pip0;
				double dfac = dab0*dab0*exp(-dab0/1.0)*expos[i];
				sumdab += dfac;
				dab += dfac*dab0;
//				cout << i << " "  << dab0 << " " << pip << " " << pie << " " << pip0 << " " << dab << " " << sumdab << " " << dab/sumdab << "\n";
			}
		}
		double dab2 = 1.0*dab/sumdab;
//		dab = distbin;
		dab = 1.0/pip;   //BAD THING, DONT!
		double pml = backback[2];
		double pmb = backback[3];
		double rva = RV;
        double vvh = (dab*pml*cos(glr)*facr)- (sin(gbr)*dab*pmb*sin(glr)*facr)  + cos(gbr)*rva*sin(glr) + 248.0;
        double uuh = -dab*pml*sin(glr)*facr - sin(gbr)*dab*pmb*cos(glr)*facr + cos(gbr)*rva*cos(glr) + 13.5;
        double wwh = sin(gbr)*rva + dab*pmb*facr*cos(gbr) + 7.0;
		double xx = dab*cos(glr)*cos(gbr) + 8.3;
		double yy = dab*sin(glr)*cos(gbr);
		double alpha = atan2(yy,xx);              //atan(yy/fabs(xx));
		double uur = uuh*cos(alpha) - vvh*sin(alpha);
		double vvr = vvh*cos(alpha) + uuh*sin(alpha);
/*

		int vvi = ((int)((vvh+500.0)/10.0+0.5));
		if (vvi < 0) {vvi = 0;}
		if (vvi >= 200) {vvi = 199;}
		vvs[vvi][0] += 1.0;
		if (met_K > -0.5) {
			vvs[vvi][1] += 1.0;
		}
		if (met_K < -1.0) {
			vvs[vvi][2] += 1.0;
		}
		if (fabs(wwh) < 60.0) {
			vvs[vvi][3] += 1.0;
		}
		vvs[vvi][5] += (uuh);
		vvs[vvi][6] += (wwh);
		vvs[vvi][7] += uuh*uuh;
		vvs[vvi][8] += wwh*wwh;
		vvs[vvi][9] += la.feh[jjh];;
*/
		ostr << countthem << " ";
		cout << countthem << " " << j << "\n";
		for (int jii = 0; jii < la.identifierlength[jjh]; jii++) {
		   ostr << la.identifier[jjh][jii];
		}
		ostr << " " << sin(glr)*cos(gbr)*sin(gbr) << " " << vvh << " " << wwh << " " << uuh << " " << glr << " " << gbr << " " << " " << pip << " " << pie << " " << dab << " " << dab2 << " " << la.feh[jjh] <<  " " << la.feherr[jjh] << " " << pip << " " << pie << " " << pmraerr_Tgas << " " << pmdecerr_Tgas <<  " " << eRV << " " << RV << " " << maia.excessnoisesig << " " << maia.excessnoise << " " << maia.deltaq << " " << j << " " << linei<< " " << jjh << " " << maia.phot[0] << " " << la.mag[1][jjh] << " " << la.sn[0][jjh] << " " << la.sn[1][jjh] << " " << la.sn[2][jjh] << " " << la.sn[3][jjh] << " " << la.sn[4][jjh] << " " << dra << " " << ddec << " " << dtot << " " << maia.radec[0] << " " << maia.radec[1] << " " << pml << " " << pmb << "\n";
//		cout << countthem << " " << nn << " " << dab << " " << RV << " "  << RAdeg << " " << DEdeg  << " " << uuh << " " << vvh << " " << wwh << " " << pmra_Tgas << " " << pmdec_Tgas << " " << glon << " " << glat << " " << backback[0] << " " << backback[1] << " " << backback[2] << " " << backback[3] << " | " << Teff_K << " " << Teff_IR << " " << met_K << " " << pi_Tgas << " " << pierr_Tgas << "\n";
		countthem++;


		}  //end if
//		cout << "gaia " << nn << " " << maia.phot[0] << " " << maia.lb[0] << " " << maia.lb[1] << "\n";

	    } //end if linei
	    linei++;

	} //end while  fgets
	fclose(in);

    } //end for j

   string stri = ostr.str();
   char *ch = &stri[0];
   fputs(ch, fout);
   fflush(fout);
    fclose(fout);
}




double ran3(long &idum);
void getnorm(double *Phi);


void readgaiar(rave la, int lan) {
  FILE *gaiout = fopen("TGASdist.txt", "w");
  ostringstream gaiostr;
  gaiostr << "#  " << "solution_id  " << "sourcid  " << "par" << "  " << "parerr" << "  " << "l" << "  " << "b" << "  " << "dist" << "  " << "disterr" << "\n";


#ifdef MOCK
  srand(time(NULL));
//  double Phi[2000001];
//  getnorm(Phi);

	std::random_device rd;     // only used once to initialise (seed) engine

	std::mt19937 rng(rd());    // random-number engine used (Mersenne-Twister in this case)

	std::normal_distribution<double> N(0.0,1.0); // (mean, sigma) guaranteed unbiased
#endif
  double ***veld = new double**[150];
  for (int nn = 0; nn < 150; nn++) {
    veld[nn] = new double*[60];
    for (int ii = 0; ii < 60; ii++) {
      veld[nn][ii] = new double[4];
      for (int jj = 0; jj < 4; jj++) {
	veld[nn][ii][jj] = 0.0;
      }
    }
  }
  double ***statdeg = new double**[50];
  for (int nn = 0; nn < 50; nn++) {
    statdeg[nn] = new double*[25];
    for (int mm = 0; mm < 25; mm++) {
      statdeg[nn][mm] = new double[4];
      for (int ii = 0; ii < 4; ii++) {
	statdeg[nn][mm][ii] = 0.0;
      }
    }
  }


#ifdef TESTVELBSTAT

  double velbstat[300][20];
  for (int ii = 0; ii < 300; ii++) {
	for (int jj = 0; jj < 20; jj++) {
		velbstat[ii][jj] = 0.0;
	}
  }
#endif

#ifdef PRINTSTARS
  ostringstream namostr;
  namostr << OUTDIR << "allstar" << SUFFIX;
  string namstri = namostr.str();
  char *namch = &namstri[0];
  FILE *fpi = fopen(namch, "w");
  ostringstream allstostr;
#ifdef MOCK
 ostringstream mnamostr;
  mnamostr << OUTDIR << "mockallstar" << SUFFIX;
  string mnamstri = mnamostr.str();
  char *mnamch = &mnamstri[0];
  FILE *mfpi = fopen(mnamch, "w");
  ostringstream mallstostr;
#endif

#endif

#ifdef TESTDISTPRIOR
  double dfall[5][PRIORLENGTH];
  for (int ii = 0; ii < 5; ii++) {
    for (int si = 0; si < PRIORLENGTH; si++) {
      dfall[ii][si] = 0.0;
    }
  }
#endif

  double Tcorr[3][3][3][3];
  double Tmean[3][3];
  double Tsum = 0.0;
  double Tcorrb = 0.0; double Tcorrbm1 = 0.0; double Tcorrbm2 = 0.0;
  for (int i = 0; i < 3; i++) {
    for (int ii = 0; ii < 3; ii++) {
      for (int j = 0; j < 3; j++) {
	for (int jj = 0; jj < 3; jj++) {
	  Tcorr[i][ii][j][jj] = 0.0;
	}
      }
      Tmean[i][ii] = 0.0;
    }
  }


        double facr = (parsec)*(PI)*0.001/(365.2425*24.0*3600.0*3600.0*180.0);
  double expos[200];
  double sume = 0.0;
  for (int ii = 0; ii < 200; ii++) {
	double dii = ((double)(ii-100))/25.0;
        expos[ii] = exp(-dii*dii/2.0);
	sume += expos[ii];
  }
  sume = 1.0/sume;
  for (int ii = 0; ii < 200; ii++) {
	expos[ii] *= sume;
  }
	int countthem = 0;
	for (int j = 0; j <= 15; j++) {
	FILE *in;
	if (j == 0) {in = fopen("tgas/TgasSource_000-000-000.csv", "r");}
	if (j == 1) {in = fopen("tgas/TgasSource_000-000-001.csv", "r");}
	if (j == 2) {in = fopen("tgas/TgasSource_000-000-002.csv", "r");}
	if (j == 3) {in = fopen("tgas/TgasSource_000-000-003.csv", "r");}
	if (j == 4) {in = fopen("tgas/TgasSource_000-000-004.csv", "r");}
	if (j == 5) {in = fopen("tgas/TgasSource_000-000-005.csv", "r");}
	if (j == 6) {in = fopen("tgas/TgasSource_000-000-006.csv", "r");}
	if (j == 7) {in = fopen("tgas/TgasSource_000-000-007.csv", "r");}
	if (j == 8) {in = fopen("tgas/TgasSource_000-000-008.csv", "r");}
	if (j == 9) {in = fopen("tgas/TgasSource_000-000-009.csv", "r");}
	if (j == 10) {in = fopen("tgas/TgasSource_000-000-010.csv", "r");}
	if (j == 11) {in = fopen("tgas/TgasSource_000-000-011.csv", "r");}
	if (j == 12) {in = fopen("tgas/TgasSource_000-000-012.csv", "r");}
	if (j == 13) {in = fopen("tgas/TgasSource_000-000-013.csv", "r");}
	if (j == 14) {in = fopen("tgas/TgasSource_000-000-014.csv", "r");}
	if (j == 15) {in = fopen("tgas/TgasSource_000-000-015.csv", "r");}

	gaia maia;
	char buffer[4096];
	char **pointit = NULL;
	int linei = -1;
    cout << "great new j " << j << "\n";
	while (fgets(buffer, 4095, in)) {

	    if (linei > -1) {
		int nn = maia.readbuffer(buffer, pointit);
	{
	  double par = maia.par;
	  double parerr = maia.parerr;
	  if (par/parerr > PARQUALLIMIT) {
	    double l = maia.lb[0];
	    double b = maia.lb[1];
	    if (l > 180.0) {
	      l = (l - 360.0);
		      }
	    double laccd = 360.0;
	    double baccd = 360.0;
	    double l1d = fabs(l);
	    double l2d = fabs(180.0 - fabs(l));
	    double l3d = fabs(360.0 -l);
	    if (l1d < l2d) {laccd = l1d;}
	    else {laccd = l2d;}
	    if (l3d < laccd) { laccd = l3d;}
	    baccd = fabs(b);
	    if ((laccd < ANGACCEPTl) && (baccd < ANGACCEPTb) && (b > BLIM) && (l > LLIM) && (l < LULIM)) {
	      double glr = l*PI/180.0;
	      double gbr = b*PI/180.0;
	      double RA_Tgas = maia.radec[0];
	      double DEC_Tgas = maia.radec[1];
	      double pmra_Tgas = maia.radec[2];
	      double pmdec_Tgas = maia.radec[3];
	      double pmraerr_Tgas = maia.radecerr[2];
	      double pmdecerr_Tgas = maia.radecerr[3];
	      double backback[4];
	      backback[0] = 0.0; backback[1] = 0.0; backback[2] = 0.0; backback[3] = 0.0;
	      propmotranslation(RA_Tgas, DEC_Tgas, pmra_Tgas, pmdec_Tgas, backback, 2000);
	      double pml = backback[2]; double pmb = backback[3];
	      double s = 1.0/par;
	      double dab = s;
	      double dabsq = 0.0;
	      dab = 0.0;
	      double sumdab = 0.0;
		int bacci = (int)(baccd);
		int lacci = (int)(laccd);
		int angi = bacci;
		if (lacci > bacci) {angi = lacci;}

	      if (par > parerr*PARQUALLIMITC) {
		double parerrs = sqrt(parerr*parerr + 0.00*0.00);
		double s0 = 1.0/par;
		double parsq = 1.0/(parerrs*parerrs*2.0);
		double deltas1 = s0*0.002*parerrs/par;
		double deltas2 = s0*0.002*parerrs/par;
		double s1 = s0;
		double s2 = s0;
		int flag = 0;
		int flag1 = 0;
		int flag2 = 0;
		double xfac = -cos(glr)*cos(gbr);
		double yfac = sin(glr)*cos(gbr);
		double zfac = fabs(sin(gbr));
		int counter = 0;
		int counter0 = 0;
		while (flag ==0) {
		  if (flag1 == 0) {
		    s1 -= deltas1;
		    if (s1 > 0.0) {
		      double xx = RSUN + xfac*s1;
		      double yy = yfac*s1;
		      double zz = zfac*s1;
		      double R = sqrt(xx*xx + yy*yy);
		      double rr = sqrt(R*R + zz*zz);
		      double par1 = 1.0/s1;
		      double rho = exp((RSUN-R)/2.5)*(1.0/1.12)*(exp(-zz/0.3) + 0.12*exp(-zz/0.9)) + 0.001*(pow(rr, -2.5)/(pow(RSUN, -2.5)));
		      double expo = -(par1-par)*(par1 - par)*parsq;
		      double exph = 0.0;
		      if (expo > -300.0) {
			exph = exp(expo);
		    }
		      double dfac = s1*s1*rho*exph*selfunction(s1)*deltas1;           // s1*s1*(exp(-s1/0.1) + 0.01*exp(-s1))*rho*exph*deltas1;
		      sumdab += dfac;
		      dab += dfac*s1;
		      dabsq += dfac*s1*s1;
		      if (dfac < (2.0e-5)*sumdab) {deltas1 *= 1.4;}
		      if (deltas1 > 0.005) {
			deltas1 = 0.005;
		    }
		      counter++;
		      counter0++;
		      }
		    else {
		      flag1 = 1;
		      if (flag2 == 1) {
			flag = 1;
		    }
		    }
		  }
		  if (flag2 == 0) {
		    s2 += deltas2;
		    if (s2 < 1.0e+6) {
		      double xx = RSUN + xfac*s2;
		      double yy = yfac*s2;
		      double zz = zfac*s2;
		      double R = sqrt(xx*xx + yy*yy);
		      double rr = sqrt(R*R + zz*zz);
		      double par2 = 1.0/s2;
		      double rho = exp((RSUN-R)/2.5)*(1.0/1.12)*(exp(-zz/0.3) + 0.12*exp(-zz/0.9)) + 0.001*(pow(rr, -2.5)/(pow(RSUN, -2.5)));
		      double expo = -(par2-par)*(par2 - par)*parsq;
		      double exph = 0.0;
		      if (expo > -300.0) {
			exph = exp(expo);
		      }
		      double dfac = s2*s2*rho*exph*deltas2*selfunction(s2);               //s2*s2*(exp(-s2/0.1) + 0.01*exp(-s2))*rho*exph*deltas2;
		      sumdab += dfac;
		      dab += dfac*s2;
		      dabsq += dfac*s2*s2;
		      if (dfac < (2.0e-5)*sumdab) {deltas2 *= 1.6;}
		      counter++;
		    }
		    else {
		      flag2 = 1;
		      if (flag1 == 1) {
			flag = 1;
		      }
		    }
		  }
		  //              cout << " " << par/parerr << " " << par << " " << parerr << " " << s1 << " " << s2 << " " << deltas1 << " " << deltas2 << " " << dab << " " << sumdab << " " << dab/sumdab << "\n";
		}
		//      cout << " " << counter << " " << counter0 << " " << par/parerr << " " << par << " " << parerr << " " << dab/sumdab << " " << dab << " " << sumdab << " " << deltas1 << " " << deltas2 << "\n";

		  }
	      double dab2 = 1.0*dab/sumdab;
	      //      dab = distbin;
	      dab = 1.0/par;   //BAD THING, DONT!
	      double dabp = dab;
	      if (DIAL == 0) {
                dab = dabp;
	      }
	      if (DIAL == 1) {
                dab = dab2;
		dabsq = (sqrt((dabsq/sumdab) - dab2*dab2));
		}
	      for (int ii = 0; ii <= maia.solidl; ii++) {
		gaiostr << maia.solid[ii];
	      }
		gaiostr << "  ";
		for (int ii = 0; ii <= maia.sourcidl; ii++) {
			gaiostr << maia.sourcid[ii];
		}
		gaiostr << "  " << par << "  " << parerr << "  " << l << "  " << b << "  " << dab << "  " << dabsq << "\n";

#ifdef TESTDISTPRIOR
	      if ((par/parerr > PARQUALLIMIT) && (VVCONDg)) {
		double xfac = -cos(glr)*cos(gbr);
		double yfac = sin(glr)*cos(gbr);
		double zfac = fabs(sin(gbr));
		double ddf[4][PRIORLENGTH];
		for (int si = 0; si < PRIORLENGTH; si++) {
		  ddf[0][si] = 0.0; ddf[1][si] = 0.0; ddf[2][si] = 0.0; ddf[3][si] = 0.0;
		  double s1 = 0.01*((double)(si) + 0.5);
		  double xx = RSUN + xfac*s1;
		  double yy = yfac*s1;
		  double zz = zfac*s1;
		  double R = sqrt(xx*xx + yy*yy);
		  double rr = sqrt(R*R + zz*zz);
		  double par1 = 1.0/s1;
		  double rho = exp((RSUN-R)/2.5)*(1.0/1.12)*(exp(-zz/0.3) + 0.12*exp(-zz/0.9)) + 0.001*(pow(rr, -2.5)/(pow(RSUN, -2.5)));
		  double dfac = s1*s1*rho;
		  double dfac2 = dfac*selfunction(s1);
		  double dfac3 = dfac*(exp(-s1/0.1))*s1;
		  double dfac4 = dfac*(exp(-s1/0.15));
		  ddf[0][si] += dfac;
		  ddf[1][si] += dfac2;
		  ddf[2][si] += dfac3;
		  ddf[3][si] += dfac4;
                    }
		double ddfsum[4];
		ddfsum[0] = 0.0; ddfsum[1] = 0.0; ddfsum[2] = 0.0; ddfsum[3] = 0.0;
		for (int si = 0; si < PRIORLENGTH; si++) {
		  ddfsum[0] += ddf[0][si];
		  ddfsum[1] += ddf[1][si];
		  ddfsum[2] += ddf[2][si];
		  ddfsum[3] += ddf[3][si];
                }
		ddfsum[0] = 1.0/ddfsum[0];
		ddfsum[1] = 1.0/ddfsum[1]; ddfsum[2] = 1.0/ddfsum[2]; ddfsum[3] = 1.0/ddfsum[3];
		for (int si = 0; si < PRIORLENGTH; si++) {
		  ddf[0][si] *= ddfsum[0];
		  ddf[1][si] *= ddfsum[1];
		  ddf[2][si] *= ddfsum[2]; ddf[3][si] *= ddfsum[3];
		  }
		for (int si = 0; si < PRIORLENGTH; si++) {
		  dfall[0][si] += ddf[0][si];
		  dfall[1][si] += ddf[1][si]; dfall[2][si] += ddf[2][si]; dfall[3][si] += ddf[3][si];
		  }
		int disti = (int)(dab*100.0 + 0.5);
		if (disti < 0) {disti = 0;}
		if (disti >= PRIORLENGTH) {disti = PRIORLENGTH-1;}
		dfall[4][disti] += 1.0;
	      }
#endif
	      double cb = cos(gbr); double sb = sin(gbr); double sl = sin(glr); double cl = cos(glr);
	      double cb2 = cb*cb; double sb2 = sb*sb; double sl2 = sl*sl; double cl2 = cl*cl;
	      double xx = -dab*cl*cb + 8.27;
	      double yy = dab*sl*cb;
	      double zz = dab*sb;
	      double RR = sqrt(xx*xx + yy*yy);
	      double alpha = atan2(yy,xx);                 //atan(yy/fabs(xx));
	      double beta = atan(zz/RR);
	      double calpha = cos(alpha); double salpha = sin(alpha);
	      double vvstar = VROT*cos(alpha); double uustar = VROT*sin(alpha);
//	      double rva = -WLSR*sb + (uustar - ULSR)*cb*cl + (vvstar - VSUN)*cb*sl;
		double rva = RVATERM;
		pml += PMLTERM;
		pmb += PMBTERM;

#ifdef MOCK

                long idum2 = -rand();
		double dice3 = ran3(idum2);
                double sigu = 50.0;
                double sigv = 40.0;
                double sigw = 30.0;
                double vvrot = 225.0;
                if (dice3 < 0.0001) {
                    sigu = 140.0;
                    sigv = 90.0;
                    sigw = 90.0;
                    vvrot = 0.0;
                }
//                idum2 = -rand();
//                dice3 = ran3(idum2);
//                int dicei3 = ((int)(2000000.0*dice3));
                double gauss = N(rng);   //Phi[dicei3];
                double vrr = gauss*sigu;
 //               idum2 = -rand();
 //               dicei3 = ((int)(2000000.0*dice3));
//                gauss = Phi[dicei3];
		gauss = N(rng);
               double vff = gauss*sigv + vvrot;
//                idum2 = -rand();
		gauss = N(rng);
//                gauss = Phi[dicei3];
                double vtt = gauss*sigw;
                double vuu = vrr*calpha + vff*salpha;
                double vvv = vff*calpha - vrr*salpha;
                double vww = vtt;
                double vuuh = vrr - ULSR;
                double vvvh = vff - VSUN;
                double vwwh = vtt - WLSR;
                double mrva = vuuh*cb*cl + vvvh*cb*sl + vwwh*sb;
                double mpml = (1.0/(dab*facr))*(-vuuh*sl + vvvh*cl);
                double mpmb = (1.0/(dab*facr))*(vwwh*cb - (vuuh*sb*cl + vvvh*sl*sb));
#endif


	      /*
		double wwh = (s*pmb*facr);
		if (wwh > (300.0 - 7.25)) {
		wwh = 300.0;
		}
		if (wwh < (-300.0 - 7.25)) {
		wwh = -300.0;
		}
	      */
	      double T[3][3];
	      T[0][0] = 1.0 - cb2*cl2; T[1][1] = 1.0 - cb2*sl2; T[2][2] = cb2;
	      T[0][1] = -0.5*cb2*sin(2.0*glr); T[0][2] = -0.5*sin(2.0*gbr)*cl; T[1][2] = -0.5*sin(2.0*gbr)*sl;
	      T[1][0] = T[0][1]; T[2][0] = T[0][2]; T[2][1] = T[1][2];
	      double vvh = (dab*pml*cl*facr)- (sb*dab*pmb*sl*facr);     // + cb*rva*sl + VSUN; // + cos(gbr)*rva*sin(glr) + 248.0;
	      double uuh = -dab*pml*sl*facr - sb*dab*pmb*cl*facr;     // + cb*rva*cl + ULSR;  // + 13.5;
	      double wwh = dab*pmb*facr*cb;          //  + sb*rva + WLSR;  // + 7.0;
	      int vvi = (int)((vvh + 400.0)*0.05 + 0.5);
	      if (vvi < 0) {vvi = 0;}
	      if (vvi >= 150) {vvi = 149;}
#ifdef INVERTT
	      double vvhf = 1.0/T[1][1];                           //1.0/(sqrt(T[1][1]));
	      double uuhf = T[0][0];                            //1.0 - cos(gbr)*cos(gbr)*cos(glr)*cos(glr);
	      if (uuhf < 1.0e-30) {
		uuhf = 0.0;
                }
		else {
		uuhf = 1.0/(uuhf);
			}
	      double wwhf = 1.0/T[2][2];
		double det = T[1][1]*T[2][2] - T[1][2]*T[1][2];
		double vvhcorr = (vvh*T[2][2] - wwh*T[1][2])/det;
		double wwhcorr = (wwh*T[1][1] - vvh*T[1][2])/det;
		double uuhcorr = uuh*uuhf;
	      if (uuhcorr > UUHLIM) {
		uuhcorr = UUHLIM;
		}
	      else {
		if (uuhcorr < -UUHLIM) {
			uuhcorr = -UUHLIM;
		}
	      }
	      wwhcorr += (sb*rva + WCOR);
	      vvhcorr += (cb*rva*sl + VCOR);
	      uuhcorr += (cb*rva*cl + UCOR);
#else
	      double vvhf = 1.0/T[1][1];                           //1.0/(sqrt(T[1][1]));
	      double vvhcorr = vvh*vvhf;           //        /(sqrt(1.0 - cos(gbr)*cos(gbr)*sin(glr)*sin(glr)));
	      double uuhf = T[0][0];                            //1.0 - cos(gbr)*cos(gbr)*cos(glr)*cos(glr);
	      if (uuhf < 1.0e-30) {
		uuhf = 0.0;
	      }
	      else {
		uuhf = 1.0/(uuhf);
	      }
	      double uuhcorr = uuh*uuhf;
	      double wwhf = 1.0/T[2][2];
	      double wwhcorr = wwh*wwhf;
	      //				        double vvh = (dab*pml*cl*facr)- (sb*dab*pmb*sl*facr) + cb*rva*sl + VSUN; // + cos(gbr)*rva*sin(glr) + 248.0;
	      //				        double uuh = -dab*pml*sl*facr - sb*dab*pmb*cl*facr + cb*rva*cl + ULSR;  // + 13.5;
	      //				        double wwh = dab*pmb*facr*cb  + sb*rva + WLSR;  // + 7.0;
	      vvhcorr += (cb*rva*sl + VCOR);
	      uuhcorr += (cb*rva*cl + UCOR);
	      if (uuhcorr > UUHLIM) {
		uuhcorr = UUHLIM;
	      }
	      else {
		if (uuhcorr < -UUHLIM) {
			uuhcorr = -UUHLIM;
		}
	      }
	      wwhcorr += (sb*rva + WCOR);
#endif
	      double uur = uuh*calpha - vvh*salpha;
	      double vvr = vvh*calpha + uuh*salpha;
	      double uurcorr = uuhcorr*calpha - vvhcorr*salpha;
	      double vvrcorr = vvhcorr*calpha + uuhcorr*salpha;
	      double Lzh = vvhcorr*RR;
	      double Lzr = vvrcorr*RR;
	      double Lzhv = Lzh/RSUN;
	      double Lzrv = Lzh/RSUN;

	      int vvibih = (int)((Lzhv*0.05 + 0.5));                               //(int)((vvrcorr + 150.0)*(RR/8.27)*0.05 + 0.5);
	      if (vvibih < 0) {vvibih = 0;}
	      if (vvibih >= 150) {vvibih = 149;}
	      int vvibir = (int)((Lzrv*0.05 + 0.5));                               //(int)((vvrcorr + 150.0)*(RR/8.27)*0.05 + 0.5);
	      if (vvibir < 0) {vvibir = 0;}
	      if (vvibir >= 150) {vvibir = 149;}
	      double uuhfsq = 1.0/(uuhf*uuhf); double vvhfsq = 1.0/(vvhf*vvhf);
	      double uurcorrf = 1.0/(uuhfsq*calpha*calpha + vvhfsq*salpha*salpha);
	      double vvrcorrf = 1.0/(uuhfsq*salpha*salpha + vvhfsq*calpha*calpha);

#ifdef MOCK
	      mrva = 0.0;
	      mrva += RVATERM;
	      mpml += PMLTERM;
	      mpmb += PMBTERM;
	      double mvvh = (dab*mpml*cl*facr)- (sb*dab*mpmb*sl*facr);      // + cb*rva*sl + VSUN; // + cos(gbr)*rva*sin(glr) + 248.0;
	      double muuh = -dab*mpml*sl*facr - sb*dab*mpmb*cl*facr;     // + cb*rva*cl + ULSR;  // + 13.5;
	      double mwwh = dab*mpmb*facr*cb;          //  + sb*rva + WLSR;  // + 7.0;
	      int mvvi = (int)((mvvh + 400.0)*0.05 + 0.5);
	      if (mvvi < 0) {mvvi = 0;}
	      if (mvvi >= 150) {mvvi = 149;}
#ifdef INVERTT
	      double mvvhf = 1.0/T[1][1];                           //1.0/(sqrt(T[1][1]));
	      double muuhf = T[0][0];                            //1.0 - cos(gbr)*cos(gbr)*cos(glr)*cos(glr);
	      if (muuhf < 1.0e-30) {
		muuhf = 0.0;
	      }
	      else {
		muuhf = 1.0/(muuhf);
	      }
	      double mwwhf = 1.0/T[2][2];
		double mdet = T[1][1]*T[2][2] - T[1][2]*T[1][2];
		double mvvhcorr = (mvvh*T[2][2] - mwwh*T[1][2])/mdet;
		double mwwhcorr = (mwwh*T[1][1] - mvvh*T[1][2])/mdet;
		double muuhcorr = muuh*muuhf;
	      if (muuhcorr > UUHLIM) {
		muuhcorr = UUHLIM;
	      }
	      else {
		if (muuhcorr < -UUHLIM) {
			muuhcorr = -UUHLIM;
		}
	      }
	      mvvhcorr += (cb*mrva*sl + VCOR);
	      muuhcorr += (cb*mrva*cl + UCOR);
	      mwwhcorr += (sb*mrva + WCOR);

#else
	      double mvvhf = 1.0/T[1][1];                           //1.0/(sqrt(T[1][1]));
	      double mvvhcorr = mvvh*mvvhf;           //        /(sqrt(1.0 - cos(gbr)*cos(gbr)*sin(glr)*sin(glr)));
	      double muuhf = T[0][0];                            //1.0 - cos(gbr)*cos(gbr)*cos(glr)*cos(glr);
	      if (muuhf < 1.0e-30) {
		muuhf = 0.0;
	      }
	      else {
		muuhf = 1.0/(muuhf);
	      }
	      double muuhcorr = muuh*muuhf;
	      double mwwhf = 1.0/T[2][2];
	      double mwwhcorr = mwwh*mwwhf;
	      //				        double vvh = (dab*pml*cl*facr)- (sb*dab*pmb*sl*facr) + cb*rva*sl + VSUN; // + cos(gbr)*rva*sin(glr) + 248.0;
	      //				        double uuh = -dab*pml*sl*facr - sb*dab*pmb*cl*facr + cb*rva*cl + ULSR;  // + 13.5;
	      //				        double wwh = dab*pmb*facr*cb  + sb*rva + WLSR;  // + 7.0;
	      mvvhcorr += (cb*mrva*sl + VCOR);
	      muuhcorr += (cb*mrva*cl + UCOR);
	      mwwhcorr += (sb*mrva + WCOR);
	      if (muuhcorr > UUHLIM) {
		muuhcorr = UUHLIM;
	      }
	      else {
		if (muuhcorr < -UUHLIM) {
			muuhcorr = -UUHLIM;
		}
	      }
#endif
	      double muur = muuh*calpha - mvvh*salpha;
	      double mvvr = mvvh*calpha + muuh*salpha;
	      double muurcorr = muuhcorr*calpha - mvvhcorr*salpha;
	      double mvvrcorr = mvvhcorr*calpha + muuhcorr*salpha;
	      double mLzh = mvvhcorr*RR;
	      double mLzr = mvvrcorr*RR;
	      double mLzhv = mLzh/RSUN;
	      double mLzrv = mLzh/RSUN;

	      int mvvibih = (int)((mLzhv*0.05 + 0.5));                               //(int)((vvrcorr + 150.0)*(RR/8.27)*0.05 + 0.5);
	      if (mvvibih < 0) {mvvibih = 0;}
	      if (mvvibih >= 150) {mvvibih = 149;}
	      int mvvibir = (int)((mLzrv*0.05 + 0.5));                               //(int)((vvrcorr + 150.0)*(RR/8.27)*0.05 + 0.5);
	      if (mvvibir < 0) {mvvibir = 0;}
	      if (mvvibir >= 150) {mvvibir = 149;}
	      double muuhfsq = 1.0/(muuhf*muuhf); double mvvhfsq = 1.0/(mvvhf*mvvhf);
	      double muurcorrf = 1.0/(muuhfsq*calpha*calpha + mvvhfsq*salpha*salpha);
	      double mvvrcorrf = 1.0/(muuhfsq*salpha*salpha + mvvhfsq*calpha*calpha);
	      if ((mwwh + WLSR) > WMAX) {
		mwwh = WMAX - WLSR;
	      }
	      if ((mwwh + WLSR) < -WMAX) {
		mwwh = (-WMAX) - WLSR;
	      }
	      if (mwwhcorr > WMAX) {
		mwwhcorr = WMAX;
	      }
	      if (mwwhcorr < -WMAX) {
		mwwhcorr = -WMAX;
	      }

#endif



	      if ((wwh + WLSR) > WMAX) {
		wwh = WMAX - WLSR;
	      }
	      if ((wwh + WLSR) < -WMAX) {
		wwh = (-WMAX) - WLSR;
	      }
	      if (wwhcorr > WMAX) {
		wwhcorr = WMAX;
	      }
	      if (wwhcorr < -WMAX) {
		wwhcorr = -WMAX;
	      }


	      if ((par/parerr > PARQUALLIMIT) && (vvrcorr > VVRLIMIT) && (VVCONDg)) {
#ifdef PRINTSTARS
		for (int ti = 0; ti < 3; ti++) {
		  for (int tii = 0; tii < 3; tii++) {
		    for (int tj = 0; tj < 3; tj++) {
		      for (int tjj = 0; tjj < 3; tjj++) {
			Tcorr[ti][tii][tj][tjj] += T[ti][tii]*T[tj][tjj];
		      }
		    }
		    Tmean[ti][tii] += T[ti][tii];
		  }
		}
		Tsum += 1.0;
		Tcorrb += (T[0][1]/T[1][1])*(T[0][2]/T[2][2]);
		Tcorrbm1 += (T[0][1]/T[1][1]);
		Tcorrbm2 += (T[0][2]/T[2][2]);

		for (int nii = angi; nii < 50; nii++) {
		  for (int iii = 0; iii < 20; iii++) {
		    statdeg[nii][iii][0] += 1.0;
		  }
		  statdeg[nii][0][1] += uurcorr;
		  statdeg[nii][0][2] += uurcorr*uurcorr;
		  statdeg[nii][1][1] += vvrcorr;
		  statdeg[nii][1][2] += vvrcorr*vvrcorr;
		  statdeg[nii][2][1] += wwhcorr;
		  statdeg[nii][2][2] += wwhcorr*wwhcorr;
		  statdeg[nii][3][1] += uuh;
		  statdeg[nii][3][2] += uuh*uuh;
		  statdeg[nii][4][1] += vvh;
		  statdeg[nii][4][2] += vvh*vvh;
		  statdeg[nii][5][1] += wwh;
		  statdeg[nii][5][2] += wwh*wwh;
		  statdeg[nii][6][1] += uuhcorr;
		  statdeg[nii][6][2] += uuhcorr*uuhcorr;
		  statdeg[nii][7][1] += vvhcorr;
		  statdeg[nii][7][2] += vvhcorr*vvhcorr;
		  statdeg[nii][8][1] += wwhcorr*vvrcorr;
		  statdeg[nii][8][2] += wwhcorr*vvhcorr;
		  statdeg[nii][9][1] += wwh*vvh;
		  statdeg[nii][9][2] += wwh*uuh;
		  statdeg[nii][1][3] += vvr;
		  statdeg[nii][2][3] += vvr*vvr;
		  statdeg[nii][3][3] += vvr*wwh;
		  statdeg[nii][4][3] += vvr*wwhcorr;
		  statdeg[nii][10][1] += T[0][0];
		  statdeg[nii][10][2] += T[0][0]*T[0][0];
		  statdeg[nii][11][1] += T[0][1];
		  statdeg[nii][11][2] += T[0][1]*T[0][1];
		  statdeg[nii][12][1] += T[0][2];
		  statdeg[nii][12][2] += T[0][2]*T[0][2];
		  statdeg[nii][13][1] += T[1][1];
		  statdeg[nii][13][2] += T[1][1]*T[1][1];
		  statdeg[nii][14][1] += T[1][2];
		  statdeg[nii][14][2] += T[1][2]*T[1][2];
		  statdeg[nii][15][1] += T[2][2];
		  statdeg[nii][15][2] += T[2][2]*T[2][2];
		  statdeg[nii][16][1] += (T[0][1]/T[1][1])*(T[0][2]/T[2][2]);
		  statdeg[nii][16][2] += (T[0][1]/T[1][1]);
		  statdeg[nii][16][3] += (T[0][2]/T[2][2]);
		  statdeg[nii][17][1] += (T[0][1]*T[0][2]);
		  statdeg[nii][17][2] += (T[1][2]/T[2][2]);
		  statdeg[nii][17][3] += T[1][1]*T[1][2];
		  statdeg[nii][18][1] += (T[1][2]/T[1][1]);
		  statdeg[nii][18][2] += T[1][2]*T[2][2];
		  statdeg[nii][19][1] += (T[0][1]/T[1][1])*beta;
		  statdeg[nii][19][2] += (T[0][1]/T[1][1]);
		  statdeg[nii][19][3] += beta;
		  statdeg[nii][20][1] += (T[0][2]/T[2][2])*(T[1][2]/T[1][1])*beta;
		  statdeg[nii][20][2] += (T[0][2]/T[2][2]);
		  statdeg[nii][20][3] += (T[1][2]/T[1][1])*beta;
		  statdeg[nii][21][1] += Lzh;
		  statdeg[nii][21][2] += Lzh*Lzh;
		  statdeg[nii][22][1] += Lzr;
		  statdeg[nii][22][2] += Lzr*Lzr;
		  statdeg[nii][23][1] += Lzh*wwh;
		  statdeg[nii][23][2] += Lzh*wwhcorr;
		  statdeg[nii][24][1] += Lzr*wwh;
		  statdeg[nii][24][2] += Lzr*wwhcorr;
		}
		  allstostr << RR << " " << zz << " " << l << " " << b << " " << glr << " " << gbr << " " << pml << " " << pmb << " " << uuh << " " << vvh << " " << wwh << " " << s << " " << par << " " << parerr << " " << uuhcorr << " " << vvhcorr << " " << wwhcorr << " " << uurcorr << " " << vvrcorr << " " << uuhf << " " << vvhf << " " << wwhf << " " << uurcorrf << " " << vvrcorrf << " " << rva << " tstat " << Tsum << " ";
#ifdef MOCK
                  mallstostr << RR << " " << zz << " " << l << " " << b << " " << glr << " " << gbr << " " << mpml << " " << mpmb << " " << muuh << " " << mvvh << " " << mwwh << " " << s << " " << par << " " << parerr << " " << muuhcorr << " " << mvvhcorr << " " << mwwhcorr << " " << muurcorr << " " << mvvrcorr << " " << muuhf << " " << mvvhf << " " << mwwhf << " " << muurcorrf << " " << mvvrcorrf << " " << mrva << " tstat " << Tsum << " " << vuu << " " << vvv << " " << vww << " " << vuuh << " " << vvvh << " " << vwwh << "\n";
#endif
		  for (int ti = 0; ti < 3; ti++) {
		    for (int tj = 0; tj < 3; tj++) {
		      allstostr << Tmean[ti][tj] << " ";
		    }
		  }
		  allstostr << Tcorr[0][1][0][2] << " " << Tcorr[1][1][1][2] << " " << Tcorr[1][2][2][2] << " " << Tcorrb << " " << Tcorrbm1 << " " << Tcorrbm2 << "\n";

#ifdef TESTVELBSTAT

	int vvri = (int)((vvrcorr + 100.0)/4.0 + 0.5);
	if (vvri < 0) {
		vvri = 0;
	}
	if (vvri >= 300) {
		vvri = 300;
	}
	velbstat[vvri][0] += 1.0;
	velbstat[vvri][1] += wwhcorr;
	velbstat[vvri][2] += wwhcorr*wwhcorr;
	velbstat[vvri][3] += uurcorr;
	velbstat[vvri][4] += uurcorr*uurcorr;
	velbstat[vvri][5] += zz;
	velbstat[vvri][6] += zz*zz;
	velbstat[vvri][7] += dab;
	velbstat[vvri][8] += dab*dab;
	velbstat[vvri][9] += uuhcorr;
	velbstat[vvri][10] += uuhcorr*uuhcorr;
	velbstat[vvri][11] += parerr/par;
	velbstat[vvri][12] += (parerr*parerr)/(par*par);
	velbstat[vvri][13] += par;
	velbstat[vvri][14] += par*par;
#endif

#endif

		  veld[vvibih][6][0] += 1.0;
		  veld[vvibih][6][1] += RR;
		  veld[vvibih][6][2] += wwhcorr;
		  veld[vvibih][6][3] += wwhcorr*wwhcorr;
		  veld[vvibih][7][0] += 1.0;
		  veld[vvibih][7][1] += RR;
		  veld[vvibih][7][2] += wwh;
		  veld[vvibih][7][3] += wwh*wwh;
		  veld[vvibir][12][0] += 1.0;
		  veld[vvibir][12][1] += RR;
		  veld[vvibir][12][2] += wwhcorr;
		  veld[vvibir][12][3] += wwhcorr*wwhcorr;
		  veld[vvibir][13][0] += 1.0;
		  veld[vvibir][13][1] += RR;
		  veld[vvibir][13][2] += wwh;
		  veld[vvibir][13][3] += wwh*wwh;

		  if ((fabs(b) < 10.0) && ((fabs(l) < 10.0) || (fabs(180.0 - fabs(l)) < 10.0))) {
		    veld[vvi][39][0] += 1.0;
		    veld[vvi][39][1] += 1;
		    veld[vvi][39][2] += wwh;
		    veld[vvi][39][3] += wwh*wwh;
		  }
		  if ((fabs(b) < 10.0) && ((fabs(l) < 10.0) || (fabs(180.0 - fabs(l)) < 10.0))) {
		    if (RR > 8.3) {
		      veld[vvi][0][0] += 1.0;
		      //								cout << "veld " << veld[vvi][0][0] << "\n";
		      veld[vvi][0][1] += s;
		      veld[vvi][0][2] += wwh;
		      veld[vvi][0][3] += wwh*wwh;
		    }
		    else {
		      veld[vvi][1][0] += 1.0;
		      veld[vvi][1][1] += s;
		      veld[vvi][1][2] += wwh;
		      veld[vvi][1][3] += wwh*wwh;
		    }
		  }
		  if ((fabs(b) < 20.0) && ((fabs(l) < 20.0) || (fabs(180.0 - fabs(l)) < 20.0))) {
		  veld[vvibih][8][0] += 1.0;
		  veld[vvibih][8][1] += RR;
		  veld[vvibih][8][2] += wwhcorr;
		  veld[vvibih][8][3] += wwhcorr*wwhcorr;
		  veld[vvibih][9][0] += 1.0;
		  veld[vvibih][9][1] += RR;
		  veld[vvibih][9][2] += wwh;
		  veld[vvibih][9][3] += wwh*wwh;
		  veld[vvibir][10][0] += 1.0;
		  veld[vvibir][10][1] += RR;
		  veld[vvibir][10][2] += wwhcorr;
		  veld[vvibir][10][3] += wwhcorr*wwhcorr;
		  veld[vvibir][11][0] += 1.0;
		  veld[vvibir][11][1] += RR;
		  veld[vvibir][11][2] += wwh;
		  veld[vvibir][11][3] += wwh*wwh;
		    veld[vvi][38][0] += 1.0;
		    veld[vvi][38][1] += 1;
		    veld[vvi][38][2] += wwh;
		    veld[vvi][38][3] += wwh*wwh;
		  }
		  if ((fabs(b) < 15.0) && ((fabs(l) < 15.0) || (fabs(180.0 - fabs(l)) < 15.0))) {
		    if (RR > 8.3) {
		      veld[vvi][2][0] += 1.0;
		      veld[vvi][2][1] += s;
		      veld[vvi][2][2] += wwh;
		      veld[vvi][2][3] += wwh*wwh;
		    }
		    else {
		      veld[vvi][3][0] += 1.0;
		      veld[vvi][3][1] += s;
		      veld[vvi][3][2] += wwh;
		      veld[vvi][3][3] += wwh*wwh;
		    }
		  }
		  if (fabs(b) < 15.0) {
		    if (RR > 8.3) {
		      veld[vvi][4][0] += 1.0;
		      veld[vvi][4][1] += s;
		      veld[vvi][4][2] += wwh;
		      veld[vvi][4][3] += wwh*wwh;
		    }
		    else {
		      veld[vvi][5][0] += 1.0;
		      veld[vvi][5][1] += s;
		      veld[vvi][5][2] += wwh;
		      veld[vvi][5][3] += wwh*wwh;
		    }
		  }
		}
		int RRi = 0;
		if (RR > 7.9) {
		  RRi = 1;
		  if (RR > 8.1) {
		    RRi = 2;
		    if (RR > 8.3) {
		      RRi = 3;
		      if (RR > 8.5) {
			RRi = 4;
			if (RR > 8.7) {
			  RRi = 5;
			}
		      }
		    }
		  }
		}
		int zzi = 0;
		if (zz > -0.06) {
		  zzi = 1;
		  if (zz > 0.06) {
		    zzi = 2;
		  }
		}
		int zri = zzi*6 + RRi + 20;
		veld[vvi][zri][0] += 1.0;
		veld[vvi][zri][1] += s;
		veld[vvi][zri][2] += wwh;
		veld[vvi][zri][3] += wwh*wwh;
	      }
	    }
	  }
	  int mflag = 0;
	  int jjh = 0;
	  double dold = 1.0;
	  double dra = 1.0;
	  double ddec = 1.0;
	  double dtot = 1.0;
	  /*
	    for (int jj = 0; jj < lan; jj++) {
	    dra = (maia.radec[0] - la.coords[jj][0]) - 0.00000416666667*maia.radec[2]; // to ravetest2
	    ddec = (maia.radec[1] - la.coords[jj][1]) - 0.00000416666667*maia.radec[3];
	    //			dra = (maia.radec[0] - la.TGAScoords[jj][0]); // - 0.00000416666667*maia.radec[2]; // to ravetest3 and 3p
	    //			ddec = (maia.radec[1] - la.TGAScoords[jj][1]); // - 0.00000416666667*maia.radec[3];
			dtot = sqrt(dra*dra + ddec*ddec);
			if (dtot < 0.0012) {
				if (dtot <= dold) {
//				cout << countthem << " " << linei << " " << jj << " " << maia.phot[0] << " " << la.gmag[jj] << " " << maia.radec[0]-la.coords[jj][0] << " " << maia.radec[1] - la.coords[jj][1] <<  " " << (fabs((maia.radec[0] - la.coords[jj][0]) - 0.000004166666667*maia.radec[2])) << " " << (fabs((maia.radec[1] - la.coords[jj][1]) - 0.00000416666666667*maia.radec[3])) << "\n";
                        mflag = 1;
                        jjh = jj;
                        dold = dtot;
		}
			}  // end if
		}  // end for
*/
		if (mflag == 1) {
			dra = (maia.radec[0] - la.coords[jjh][0]) - 0.00000416666667*maia.radec[2];
			ddec = (maia.radec[1] - la.coords[jjh][1]) - 0.00000416666667*maia.radec[3];
			dtot = sqrt(dra*dra + ddec*ddec);

		double backback[4];
		double pie = maia.parerr;
		double pip = maia.par; // + erfinv(fsel)*pie;
		double RA_Tgas = maia.radec[0];
		double DEC_Tgas = maia.radec[1];
		double pmra_Tgas = maia.radec[2];
		double pmdec_Tgas = maia.radec[3];
		double pmraerr_Tgas = maia.radecerr[2];
		double pmdecerr_Tgas = maia.radecerr[3];
		double glr = maia.lb[0]*PI/180.0;
		double gbr = maia.lb[1]*PI/180.0;
		double RV = la.RV[jjh];
		double RVerr = la.RVerr[jjh];
		double eRV = RVerr;
		for (int ii = 0; ii < 4; ii++) {
			backback[ii] = 0.0;
		}
		propmotranslation(RA_Tgas, DEC_Tgas, pmra_Tgas, pmdec_Tgas, backback, 2000);
//		cout << backback[0] << " " << maia.lb[0]*PI/180.0 << " "  << backback[1] << " " << maia.lb[1]*PI/180.0 << "\n";
		double dab = 0.0;
		double sumdab = 0.0;


		  /*
		for (int i = 0; i < 200; i++) {
			double di = ((double)(i-100))/25.0;
			double pip0 = pip + pie*di;
			if (pip0 > 1.0e-5) {
				double dab0 = 1.0/pip0;
				double dfac = dab0*dab0*exp(-dab0/1.0)*expos[i];
				sumdab += dfac;
				dab += dfac*dab0;
//				cout << i << " "  << dab0 << " " << pip << " " << pie << " " << pip0 << " " << dab << " " << sumdab << " " << dab/sumdab << "\n";
			}
		}
		double dab2 = 1.0*dab/sumdab;
//		dab = distbin;
		dab = 1.0/pip;   //BAD THING, DONT!
		  */

		double pml = backback[2];
		double pmb = backback[3];
		double rva = RV;
        double vvh = (dab*pml*cos(glr)*facr)- (sin(gbr)*dab*pmb*sin(glr)*facr)  + cos(gbr)*rva*sin(glr) + 248.0;
        double uuh = -dab*pml*sin(glr)*facr - sin(gbr)*dab*pmb*cos(glr)*facr + cos(gbr)*rva*cos(glr) + 13.5;
        double wwh = sin(gbr)*rva + dab*pmb*facr*cos(gbr) + 7.0;
		double xx = dab*cos(glr)*cos(gbr) + 8.3;
		double yy = dab*sin(glr)*cos(gbr);
		double RR = sqrt(xx*xx + yy*yy);
		double alpha = atan2(yy,xx);         //atan(yy/fabs(RR));
		double uur = uuh*cos(alpha) - vvh*sin(alpha);
		double vvr = vvh*cos(alpha) + uuh*sin(alpha);
	//		cout << countthem << " " << sin(glr)*cos(gbr)*sin(gbr) << " " << vvh << " " << wwh << " " << uuh << " " << glr << " " << gbr << " " << " " << pip << " " << pie << " " << dab << " " << dab2 << " " << la.met[jjh] <<  " " << la.meterr[jjh] << " " << pip << " " << pie << " " << pmraerr_Tgas << " " << pmdecerr_Tgas <<  " " << eRV << " " << RV << " " << maia.excessnoisesig << " " << maia.excessnoise << " " << maia.deltaq << " " << j << " " << linei<< " " << jjh << " " << maia.phot[0] << " " << la.gmag[jjh] << " " << la.SN[jjh] << " " << la.flagc0[jjh] << " " << la.flagc1[jjh] << " " << la.clustflag[jjh] << " " << la.chisq[jjh] << " " << dra << " " << ddec << " " << dtot << " " << maia.radec[0] << " " << maia.radec[1] << " " << pml << " " << pmb << "\n";
//		cout << countthem << " " << nn << " " << dab << " " << RV << " "  << RAdeg << " " << DEdeg  << " " << uuh << " " << vvh << " " << wwh << " " << pmra_Tgas << " " << pmdec_Tgas << " " << glon << " " << glat << " " << backback[0] << " " << backback[1] << " " << backback[2] << " " << backback[3] << " | " << Teff_K << " " << Teff_IR << " " << met_K << " " << pi_Tgas << " " << pierr_Tgas << "\n";
		countthem++;

		}  //end if
//		cout << "gaia " << nn << " " << maia.phot[0] << " " << maia.lb[0] << " " << maia.lb[1] << "\n";

	    } //end if linei
	    linei++;

	} //end while  fgets
	fclose(in);

    } //end for j

#ifdef PRINTSTARS
    string allstri = allstostr.str();
    char *allch = &allstri[0];
    fputs(allch, fpi);
    fflush(fpi);
    fclose(fpi);

    ostringstream fpostr2;
    fpostr2 << OUTDIR << "statdeg" << SUFFIX;
    string fpostri2 = fpostr2.str();
    char *fpoch2 = &fpostri2[0];
    FILE *statdegout = fopen(fpoch2, "w");
    ostringstream stostr;
    for (int ii = 0; ii < 50; ii++) {
      stostr << ii;
      for (int nn = 0; nn < 25; nn++) {
	stostr << "  " << nn;
	for (int mm = 0; mm < 4; mm++) {
	  stostr << "  " << statdeg[ii][nn][mm];
	}
      }
      stostr << "\n";
    }
    string stostri = stostr.str();
    char *stoch = &stostri[0];
    fputs(stoch, statdegout);
    fflush(statdegout);
    fclose(statdegout);
#ifdef MOCK
   string mallstri = mallstostr.str();
    char *mallch = &mallstri[0];
    fputs(mallch, mfpi);
    fflush(mfpi);
    fclose(mfpi);
#endif

#ifdef TESTVELBSTAT

  ostringstream velbnamostr;
  velbnamostr << OUTDIR << "veldistb" << SUFFIX;
  string velbnamstr = velbnamostr.str();
  char *velbnamch = &velbnamstr[0];
  FILE *fvelb = fopen(velbnamch, "w");
  ostringstream velbostr;
  for (int ii = 0; ii < 300; ii++) {
	velbostr << ii << " " << ((double)(ii))*4.0 - 100.0 << " ";
	for (int jj = 0; jj < 20; jj++) {
		velbostr << velbstat[ii][jj] << " ";
	}
	velbostr << "\n";
   }
   string velbstr = velbostr.str();
   char *velbch = &velbstr[0];
   fputs(velbch, fvelb);
   fflush(fvelb);
   fclose(fvelb);


   double velbstat2[600][20];
   for (int ii = 0; ii < 600; ii++) {
	for (int jj = 0; jj < 20; jj++) {
		velbstat2[ii][jj] = 0.0;
	}
   }


   for (int ii = 0; ii < 300; ii++) {
     if (velbstat[ii][0] > 2.0) {
	double vvr = ((double)(ii))*4.0 - 100.0;
	double parerrpar = 0.15;
	double parerr = parerrpar*3.4;           //*velbstat[ii][13]/velbstat[ii][0];
	double dab0 = 0.02;
	double glr = PI;
	double gbr = 0.0;
        cout << "velbstat ii " << ii << "\n";
	while (dab0 < 3.0) {
		double weight = 0.0;     //dab0*dab0*rho
		  double rho = 1.0;          //exp((RSUN-R)/2.5)*(1.0/1.12)*(exp(-zz/0.3) + 0.12*exp(-zz/0.9)) + 0.001*(pow(rr, -2.5)/(pow(RSUN, -2.5)));
		  double dfac = dab0*dab0*rho;
		  double weight0 = velbstat[ii][0]*dfac*selfunction(dab0);
		  double par0 = 1.0/dab0;
		for (int jj = -200; jj < 200; jj++) {
			double djj = ((double)(jj))*0.035;
			double par = par0 + parerr*djj;
		      double dab = 0.0; //1.0/par;
		      double sumdab = 0.0;
		if (par > parerr*PARQUALLIMITC) {
			double expor = ((djj*djj));
			weight = exp(-expor)*weight0;


		sumdab  = 0.0;
	      dab = 0.0;

//	      if (par > parerr*PARQUALLIMITC) {
		double parerrs = sqrt(parerr*parerr + 0.00*0.00);
		double s0 = 1.0/par;
		double parsq = 1.0/(parerrs*parerrs*2.0);
		double deltas1 = s0*0.002*parerrs/par;
		double deltas2 = s0*0.002*parerrs/par;
		double s1 = s0;
		double s2 = s0;
		int flag = 0;
		int flag1 = 0;
		int flag2 = 0;
		double xfac = -cos(glr)*cos(gbr);
		double yfac = sin(glr)*cos(gbr);
		double zfac = fabs(sin(gbr));
		int counter = 0;
		int counter0 = 0;
		while (flag ==0) {
		  if (flag1 == 0) {
		    s1 -= deltas1;
		    if (s1 > 0.0) {
		      double xx = RSUN + xfac*s1;
		      double yy = yfac*s1;
		      double zz = zfac*s1;
		      double R = sqrt(xx*xx + yy*yy);
		      double rr = sqrt(R*R + zz*zz);
		      double par1 = 1.0/s1;
		      double rho = exp((RSUN-R)/2.5)*(1.0/1.12)*(exp(-zz/0.3) + 0.12*exp(-zz/0.9)) + 0.001*(pow(rr, -2.5)/(pow(RSUN, -2.5)));
		      double expo = -(par1-par)*(par1 - par)*parsq;
		      double exph = 0.0;
		      if (expo > -300.0) {
			exph = exp(expo);
		      }
		      double dfac = s1*s1*rho*exph*selfunction(s1)*deltas1;           // s1*s1*(exp(-s1/0.1) + 0.01*exp(-s1))*rho*exph*deltas1;
		      sumdab += dfac;
		      dab += dfac*s1;
		      if (dfac < (4.0e-5)*sumdab) {deltas1 *= 1.4;}
		      if (deltas1 > 0.01) {
			deltas1 = 0.01;
		      }
		      counter++;
		      counter0++;
		    }
		    else {
		      flag1 = 1;
		      if (flag2 == 1) {
			flag = 1;
		      }
		    }
		  }
		  if (flag2 == 0) {
		    s2 += deltas2;
		    if (s2 < 1.0e+6) {
		      double xx = RSUN + xfac*s2;
		      double yy = yfac*s2;
		      double zz = zfac*s2;
		      double R = sqrt(xx*xx + yy*yy);
		      double rr = sqrt(R*R + zz*zz);
		      double par2 = 1.0/s2;
		      double rho = exp((RSUN-R)/2.5)*(1.0/1.12)*(exp(-zz/0.3) + 0.12*exp(-zz/0.9)) + 0.001*(pow(rr, -2.5)/(pow(RSUN, -2.5)));
		      double expo = -(par2-par)*(par2 - par)*parsq;
		      double exph = 0.0;
		      if (expo > -300.0) {
			exph = exp(expo);
		      }
		      double dfac = s2*s2*rho*exph*deltas2*selfunction(s2);               //s2*s2*(exp(-s2/0.1) + 0.01*exp(-s2))*rho*exph*deltas2;
		      sumdab += dfac;
		      dab += dfac*s2;
		      if (dfac < (4.0e-5)*sumdab) {deltas2 *= 1.6;}
		      counter++;
		    }
		    else {
		      flag2 = 1;
		      if (flag1 == 1) {
			flag = 1;
		      }
		    }
		  }
		  //              cout << " " << par/parerr << " " << par << " " << parerr << " " << s1 << " " << s2 << " " << deltas1 << " " << deltas2 << " " << dab << " " << sumdab << " " << dab/sumdab << "\n";
		}
		//      cout << " " << counter << " " << counter0 << " " << par/parerr << " " << par << " " << parerr << " " << dab/sumdab << " " << dab << " " << sumdab << " " << deltas1 << " " << deltas2 << "\n";

	      double dab2 = 1.0*dab/sumdab;
	      //      dab = distbin;
	      dab = 1.0/par;   //BAD THING, DONT!
	      double dabp = dab;
	      if (DIAL == 0) {
                dab = dabp;
	      }
	      if (DIAL == 1) {
                dab = dab2;
	      }

		double vvrb = (vvr - VSUN)*dab/dab0 + VSUN;
		int vvrbi = ((int)(vvrb + 100.0));
		if (vvrbi < 0) {vvrbi = 0;}
		if (vvrbi > 600) {vvrbi = 600;}
		int vvri = ((int)(vvr + 100.0));
		if (vvri < 0) {vvri = 0;}
		if (vvri > 600) {vvri = 600;}
	        velbstat2[vvrbi][0] += weight;
		velbstat2[vvrbi][1] += dab*weight;
		velbstat2[vvrbi][2] += dab*dab*weight;
		velbstat2[vvrbi][3] += (dab/dab0)*weight;
		velbstat2[vvrbi][4] += (dab*dab)/(dab0*dab0)*weight;
		velbstat2[vvrbi][5] += (parerr/par)*weight;
		velbstat2[vvrbi][6] += weight*(parerr*parerr)/(par*par);
	        velbstat2[vvri][7] += weight;
		velbstat2[vvri][8] += dab*weight;
		velbstat2[vvri][9] += dab*dab*weight;
		velbstat2[vvri][10] += (dab/dab0)*weight;
		velbstat2[vvri][11] += (dab*dab)/(dab0*dab0)*weight;
		velbstat2[vvri][12] += (parerr/par)*weight;
		velbstat2[vvri][13] += weight*(parerr*parerr)/(par*par);

	      }





		}
		dab0 += 0.04;
	}  //end while dab0

     }
   }

   ostringstream namostr5;
   namostr5 << OUTDIR << "velbstat2" << SUFFIX;
   string namstr5 = namostr5.str();
   char *namch5 = &namstr5[0];
   FILE *outf5 = fopen(namch5, "w");
   ostringstream ostr5;
   for (int ii = 0; ii < 600; ii++) {
	ostr5 << ii << " " << ((double)(ii)) - 100.0;
	for (int jj = 0; jj < 15; jj++) {
		ostr5 << " " << velbstat2[ii][jj];
	}
	ostr5 << "\n";
   }
   string stri5 = ostr5.str();
   char *ch5 = &stri5[0];
   fputs(ch5, outf5);
   fflush(outf5);
   fclose(outf5);

#endif


#endif

    ostringstream velnostr;
    velnostr << OUTDIR << "veldout" << SUFFIX;
    string velnostri = velnostr.str();
    char *velnoch = &velnostri[0];
    FILE *veldout = fopen(velnoch, "w");
    ostringstream ostr;
    for (int nn = 0; nn < 150; nn++) {
      ostr << nn << " " << ((double)(nn))*20.0 -400.0 << " ";
      for (int ii = 0; ii < 40; ii++) {
	for (int jj = 0; jj < 4; jj++) {
	  ostr << veld[nn][ii][jj] << " ";
	}
      }
      ostr << "\n";
    }
    string stri = ostr.str();
    char *ch = &stri[0];
    fputs(ch, veldout);
    fflush(veldout);
    fclose(veldout);


#ifdef TESTDISTPRIOR
    {
      ostringstream fpostr;
      fpostr << OUTDIR << "priortest" << SUFFIX;
      string fpostri = fpostr.str();
      char *fpoch = &fpostri[0];
      FILE *priout = fopen(fpoch, "w");
      ostringstream ostrp;
      for (int si = 0; si < PRIORLENGTH; si++) {
        ostrp << si << " " << ((double)(si) + 0.5)*0.01 << " " << dfall[0][si] << " " << dfall[1][si] << " " << dfall[2][si] << " " << dfall[3][si] << " " << dfall[4][si] << "\n";
      }
    string strip = ostrp.str();
    char *chp = &strip[0];
    fputs(chp, priout);
    fflush(priout);
    fclose(priout);
    }
#endif
    string gaistri = gaiostr.str();
    char *gaich = &gaistri[0];
    fputs(gaich, gaiout);
    fflush(gaiout);
    fclose(gaiout);



    for (int nn = 0; nn < 150; nn++) {
      for (int ii = 0; ii < 40; ii++) {
	delete [] veld[nn][ii];
	}
      delete [] veld[nn];
    }
    delete [] veld;
    for (int nn = 0; nn < 50; nn++) {
      for (int mm = 0; mm < 25; mm++) {
	delete [] statdeg[nn][mm];
      }
      delete [] statdeg[nn];
    }
    delete [] statdeg;

}







void getnorm(double *Phi) {
    for (int n = 0; n < 2000001; n++) {
        Phi[n] = 0.0;
    }
//    cout << "Phi set up \n";
    int Phin = 0;
    double Phid = 0.0;
    for (int n = 0; n < 100000000; n++) {
        double delirg = 0.0000001;
        double irg = ((double)(n))*delirg;
        double erg = erf(irg/(sqrt(2.0)));
        if (erg >= 1.0) {
            n == 2000000000;
        }
        erg *= 1000000.0;
        if ((erg-Phid) >= 1.0) {
            Phid += 1.0;
            Phin++;
            Phi[Phin] = irg;
            if ((erg-Phid) >= 1.0) {
                Phin--;
                int delta = ((int)(erg-Phid)) +1;
                double deltad = (double)(delta);
                for (int nn = 1; nn <= delta; nn++) {
                    Phin++;
                    double nnd = (double)(nn);
                    Phi[Phin] = irg - delirg + (delirg/deltad)*nnd;
//                    cout << n << "  " << Phin << "  " << Phi[Phin] << "  " << nn << "\n";
               }
            }
//                      cout << n << "  " << irg << "  " << erg << "  "  << Phin << "  " << Phi[Phin] << "\n";
        }
//              cout << n << "  " << irg << "  " << erg << "\n";
    }
    for (int n = 1000000; n >= 0; n--) {
        Phi[n+1000000] =        Phi[n];
    }
    for (int n = 0; n <= 1000000; n++) {
        Phi[1000000-n] = -Phi[1000000+n];
    }

}



double ran3(long &idum)
{
  const long MBIG=1000000000;
  const long MSEED=161803398;
  const long MZ=0;
  const double FAC=1.0/MBIG;

  static int inext,inextp;
  static long ma[56];
  static int iff=0;
  long mj,mk;
  int i,ii,k;

  if (idum < 0 || iff == 0) {
    iff=1;
    mj=MSEED-(idum < 0 ? -idum : idum);   // abs(idum)
    mj=(mj < 0 ? -mj : mj);               // abs(mj)
    mj %= MBIG;
    ma[55]=mj;
    mk=1;
    for (i=1;i<=54;i++) {
      ii=(21*i) % 55;
      ma[ii]=mk;
      mk=mj-mk;
      if (mk < MZ) mk += MBIG;
      mj=ma[ii];
    }
    for (k=1;k<=4;k++)
      for (i=1;i<=55;i++) {
        ma[i] -= ma[1+(i+30) % 55];
        if (ma[i] < MZ) ma[i] += MBIG;
      }
    inext=0;
    inextp=31;
    idum=1;
  }
  if (++inext == 56) inext=1;
  if (++inextp == 56) inextp=1;
  mj=ma[inext]-ma[inextp];
  if (mj < MZ) mj += MBIG;
  ma[inext]=mj;
  if  (mj*FAC < 0 || mj*FAC > 1) cout << "<<<" << mj*FAC<<endl;
  return mj*FAC;
}




void sortstar(double **star, double **starz, int n, int starcount) {
  double min = 1.0e+40;
  double max = -1.0e-40;
  int mini = -1;
  int maxi = -1;
  int minl = 0;
  double *val = new double[starcount];
  int *vali = new int[starcount];
  for (int i = 0; i < starcount; i++) {
    val[i] = star[i][n];
    vali[i] = i;
  }
  cout << "enter quicksort" << time(0) << " " << starcount << "\n";
  quicksort(val, vali, 0, starcount-1);
  cout << "left quicksort " << time(0) << "\n";
  for (int i = 0; i < starcount; i++) {
    int ii = vali[i];
    for (int j = 0; j < 30; j++) {
      starz[i][j] = star[ii][j];
    }
    //    cout << i << "  " << starz[i][n] << "  " << starz[i][3] << "  " << starz[i][4] << "\n";
  }
  delete [] val;
  delete [] vali;
}


double lindisthelpVW(double **star, int i1, int i2, double factor, int mode, double sigbin, std::ostream& ostr) {
  double facr = (parsec)*(PI)*0.001/(365.24*24.0*3600.0*3600.0*180.0);
  double ym = 0.0;
  double xm  = 0.0;
  double xym = 0.0;
  double varxy = 0.0;
  double varx = 0.0;
  double vary = 0.0;
  double den = 0.0;
  double Teffm = 0.0;
  double loggm = 0.0;
  double parpm = 0.0;
  double varparp = 0.0;
  double varTeff = 0.0;
  double varlogg = 0.0;
  double fehm = 0.0;
  double ordermean = 0.0;
  double varfeh = 0.0;
  double varxyrva = 0.0;
  double varxyrvabase = 0.0;
  double varxypm = 0.0;
  double varxyvebase = 0.0;
  double varxyve = 0.0;
  double varTvw = 0.0;
  double suu = 0.0;
  double svv = 0.0;
  double sww = 0.0;
  double um = 0.0;
  double vm = 0.0;
  double wm = 0.0;
  double suuve = 0.0;
  double svvve = 0.0;
  double swwve = 0.0;
  double umve = 0.0;
  double vmve = 0.0;
  double wmve = 0.0;
  double denve = 0.0;
  double uuhm = 0.0;
  for (int ii = i1; ii < i2; ii++) {
    	double dist0 = star[ii][0];
	double dist1 = star[ii][1];
	double dist2 = star[ii][2];
	double distm = star[ii][3];
	double vhelio = star[ii][4];
	double rva = vhelio;
	double glr = star[ii][5];
   	double gbr = star[ii][6];
	double pml = star[ii][7];
	double pmb = star[ii][8];
	double feh = star[ii][9];
	double Teff = star[ii][10];
	double logg = star[ii][11];
	double rvaerr = star[ii][12];
	double pmlerr = star[ii][13];
	double pmberr = star[ii][14];
        double orderm = star[ii][ORDERING];
	double parp = star[ii][15]/star[ii][16];
	double dab = star[ii][DIAL]*factor;
        double vvh = (dab*pml*cos(glr)*facr)- (sin(gbr)*dab*pmb*sin(glr)*facr)  + cos(gbr)*rva*sin(glr);
        double uuh = -dab*pml*sin(glr)*facr - sin(gbr)*dab*pmb*cos(glr)*facr + cos(gbr)*rva*cos(glr);
        double wwh = sin(gbr)*rva + dab*pmb*facr*cos(gbr);
		double xx = RSUN - dab*cos(glr)*cos(gbr);
		double yy = dab*sin(glr)*cos(gbr);
		double RR = sqrt(xx*xx + yy*yy);
		double zz = ((dab*sin(gbr)) + ZOFFSET);
		double rr = sqrt(RR*RR + zz*zz);
		double galpha = atan2(yy,xx);                  //atan(yy/fabs(xx));
		double gbeta = atan(zz/RR);
		double uur = (uuh+ULSR)*cos(galpha) - (vvh+VSUN)*sin(galpha);
		double vvr = (vvh+VSUN)*cos(galpha) + (uuh+ULSR)*sin(galpha);
		double wwr = wwh + WLSR;
		double Vr = (-uur*RR + wwr*zz)/rr;
		double Vz = (uur*zz + wwr*RR)/rr;
	double Tvw = -sin(gbr)*cos(gbr)*sin(glr);
	double Tuw = 0.0;
         double varxyvepart = -0.5*Tvw*sin(galpha)*sin(2.0*gbeta);

//		vvh = vvr;

	if (wwh > WWHLIM - WLSR) {wwh = WWHLIM - WLSR;};
	if (wwh < -WWHLIM - WLSR) {wwh = -WWHLIM - WLSR;}
	if (VARXYCOND) {                             // ((vvh > -150.0) && (vvh < 100.0)) {
	double x = Tvw*vvh;
#ifdef NOVELUSE
        x = -sin(gbr)*cos(gbr)*sin(glr);
#endif
	ym += wwh;
	xm += x;
	varx += (x*x);
	vary += (wwh*wwh);
	varxy += (wwh*x);
	varxyrva -= Tvw*Tvw*rvaerr*rvaerr;
	varTvw += Tvw*Tvw;
	varxyrvabase -= Tvw*Tvw;
	double vebm = fabs(varxyvepart);
	varxyvebase += varxyvepart;                    //-0.5*Tvw*sin(galpha)*sin(2.0*gbeta);
	varxypm += Tvw*Tvw*dab*dab*(pmberr*pmberr)*facr*facr;
	suu += uur*uur; svv += vvr*vvr;
	sww += wwh*wwh;
	um += uur; vm += vvr;
	uuhm += uuh;
	wm += wwh;
#ifdef USECARTVE
	suuve += uur*uur*vebm; svvve += vvr*vvr*vebm;
	swwve += wwh*wwh*vebm;
	umve += uur*vebm; vmve += vvr*vebm;
	wmve += wwh*vebm;
#endif
#ifdef USERADVE
	suuve += Vr*Vr*vebm; svvve += vvr*vvr*vebm;
	swwve += Vz*Vz*vebm;
	umve += Vr*vebm; vmve += vvr*vebm;
	wmve += Vz*vebm;
#endif
 	denve += vebm;
	Teffm += Teff;
	varTeff += (Teff*Teff);
	loggm += logg;
	parpm += parp;
	varparp += parp*parp;
	varlogg += logg*logg;
	fehm += feh;
	varfeh += feh*feh;
        ordermean += orderm;
	den += 1.0;
	}
  }
  den = 1.0/den;
  denve = 1.0/denve;
  xm *= den; ym *= den; Teffm *= den; loggm *= den; fehm *= den; parpm *= den; ordermean *= den;
  um *= den; uuhm *= den; vm *= den; wm *= den; suu *= den; sww *= den; svv *= den;
  umve *= denve; vmve *= denve; wmve *= denve; suuve *= denve; swwve *= denve; svvve *= denve;
  varx *= den; vary *= den; varxy *= den; varxyrva *= den; varxypm *= den; varxyrvabase *= den;  varTvw *= den;
  varxyvebase *= den;
//  cout << "vars " << varx << " " << vary << " " << varxy << " " << varxyrva << " " << varxypm << " " << varxyrva/varx<< " " << varxy/varx << "\n";
  varTeff *= den; varlogg *= den; varfeh *= den; varparp *= den;
  varxy -= (xm*ym);
  varx -= (xm*xm);
  vary -= (ym*ym);
  suu -= (um*um); svv -= (vm*vm); sww -= (wm*wm);
  suuve -= (umve*umve); svvve -= (vmve*vmve); swwve -= (wmve*wmve);
  varxyve = varxyvebase*(swwve - suuve);
#ifdef CORRERR
  varxy -= (varxyrva + varxypm + varxyve+varxyrvabase*(sigbin*sigbin));
#endif
  double varxbb = varx + varTvw*sww;
  varTeff -= (Teffm*Teffm); varlogg -= (loggm*loggm); varfeh -= (fehm*fehm); varparp -= (parpm*parpm);
  double beta = varxy/(varxbb);
  double alpha = ym - xm*beta;
  double varres = 0.0;
  den = 0.0;
  for (int ii = i1; ii < i2; ii++) {
   	double dist0 = star[ii][0];
	double dist1 = star[ii][1];
	double dist2 = star[ii][2];
	double distm = star[ii][3];
	double vhelio = star[ii][4];
	double rva = vhelio;
	double glr = star[ii][5];
   	double gbr = star[ii][6];
	double pml = star[ii][7];
	double pmb = star[ii][8];
	double feh = star[ii][9];
	double Teff = star[ii][10];
	double logg = star[ii][11];
	double dab = star[ii][DIAL]*factor;
        double vvh = (dab*pml*cos(glr)*facr)- (sin(gbr)*dab*pmb*sin(glr)*facr)  + cos(gbr)*rva*sin(glr);
        double uuh = -dab*pml*sin(glr)*facr - sin(gbr)*dab*pmb*cos(glr)*facr + cos(gbr)*rva*cos(glr);
       double wwh = sin(gbr)*rva + dab*pmb*facr*cos(gbr);
	double Tvw = -sin(gbr)*cos(gbr)*sin(glr);
        double Tuw = 0.0;
		double xx = RSUN - dab*cos(glr)*cos(gbr);
		double yy = dab*sin(glr)*cos(gbr);
		double RR = sqrt(xx*xx + yy*yy);
		double zz = ((dab*sin(gbr)) + ZOFFSET);
		double galpha = atan2(yy,xx);             //atan(yy/fabs(xx));
		double gbeta = atan(zz/RR);
        double varxyvepart = -0.5*Tvw*sin(galpha)*sin(2.0*gbeta);

	if (wwh > WWHLIM) {wwh = WWHLIM;};
	if (wwh < -WWHLIM) {wwh = -WWHLIM;}
	if (VARXYCOND) {                              // ((vvh > -150.0) && (vvh < 100.0)) {
       double x = sin(gbr)*cos(gbr)*sin(glr)*vvh;
#ifdef NOVELUSE
       x = sin(gbr)*cos(gbr)*sin(glr);
#endif
	double res = wwh - (alpha + beta*x);
	varres += res*res;
	den += 1.0;
	}
  }
//  double betaerr = sqrt((1.0/(den-2.0))*varres/((varxbb)*den));
  double betaerr = sqrt(((1.0/(den-2.0))*varres/((varx)*den))  + ((SYSERRRV*varxyrva*varxyrva + SYSERR*(varxypm*varxypm + varxyve*varxyve))/(varxbb*varxbb)));
  double alphaerr = sqrt((1.0/(den-2.0))*varres/den);
  if (mode == 3) {
    ostr << i1 << " " << i2 << " " << Teffm << " " << varTeff << " " << loggm << " " << varlogg << "  "  << fehm << " " << varfeh << " " << parpm << " " << varparp << " " << alpha << " " << alphaerr << " " << beta << " " << betaerr << " " << den << " " << varx << " " << xm << " " << varres << " " << varxy << " " << varxypm << " " << varxyrva << " " << ordermean << " " << varxyrvabase << " " << varxyrva/varxbb << " " << varxyrvabase/varxbb << " " << varxyvebase << " " << varxyvebase/varxbb << " " << varxyve << " " << varxyve/varxbb << " " << varxypm/varxbb << " " << suu << " " << svv << " " << sww << " " << um << " " << vm << " " << wm << " " << uuhm << " " << varTvw << " " << varxbb << " " << suuve << " " << svvve << " " << swwve << " " << umve << " " << vmve << " " <<  wmve << "\n";
  }
  if (mode == 0) {
    return beta;
  }
  if (mode == 1) {
    return (beta + betaerr);
  }
  if (mode == 2) {
    return (beta - betaerr);
  }
}

double lindisthelpUVW(double **star, int i1, int i2, double factor, int mode, double sigbin, std::ostream& ostr) {
  double facr = (parsec)*(PI)*0.001/(365.24*24.0*3600.0*3600.0*180.0);
  double ym = 0.0;
  double xm  = 0.0;
  double xym = 0.0;
  double varxy = 0.0;
  double varx = 0.0;
  double vary = 0.0;
  double den = 0.0;
  double Teffm = 0.0;
  double loggm = 0.0;
  double parpm = 0.0;
  double varparp = 0.0;
  double varTeff = 0.0;
  double varlogg = 0.0;
  double fehm = 0.0;
  double ordermean = 0.0;
  double varfeh = 0.0;
  double varxyrva = 0.0;
  double varxyrvabase = 0.0;
  double varxypm = 0.0;
  double varxyvebase = 0.0;
  double varxyve = 0.0;
  double varTvw = 0.0;
  double varTuw = 0.0;
  double suu = 0.0;
  double svv = 0.0;
  double sww = 0.0;
  double um = 0.0;
  double vm = 0.0;
  double wm = 0.0;
  double suuve = 0.0;
  double svvve = 0.0;
  double swwve = 0.0;
  double umve = 0.0;
  double vmve = 0.0;
  double wmve = 0.0;
  double denve = 0.0;
  double uuhm = 0.0;
  for (int ii = i1; ii < i2; ii++) {
    	double dist0 = star[ii][0];
	double dist1 = star[ii][1];
	double dist2 = star[ii][2];
	double distm = star[ii][3];
	double vhelio = star[ii][4];
	double rva = vhelio;
	double glr = star[ii][5];
   	double gbr = star[ii][6];
	double pml = star[ii][7];
	double pmb = star[ii][8];
	double feh = star[ii][9];
	double Teff = star[ii][10];
	double logg = star[ii][11];
	double rvaerr = star[ii][12];
	double pmlerr = star[ii][13];
	double pmberr = star[ii][14];
        double orderm = star[ii][ORDERING];
	double parp = star[ii][15]/star[ii][16];
	double dab = star[ii][DIAL]*factor;
        double vvh = (dab*pml*cos(glr)*facr)- (sin(gbr)*dab*pmb*sin(glr)*facr)  + cos(gbr)*rva*sin(glr);
        double uuh = -dab*pml*sin(glr)*facr - sin(gbr)*dab*pmb*cos(glr)*facr + cos(gbr)*rva*cos(glr);
        double wwh = sin(gbr)*rva + dab*pmb*facr*cos(gbr);
		double xx = RSUN - dab*cos(glr)*cos(gbr);
		double yy = dab*sin(glr)*cos(gbr);
		double RR = sqrt(xx*xx + yy*yy);
		double zz = ((dab*sin(gbr)) + ZOFFSET);
		double rr = sqrt(RR*RR + zz*zz);
		double galpha = atan2(yy,xx);                  //atan(yy/fabs(xx));
		double gbeta = atan(zz/RR);
		double uur = (uuh+ULSR)*cos(galpha) - (vvh+VSUN)*sin(galpha);
		double vvr = (vvh+VSUN)*cos(galpha) + (uuh+ULSR)*sin(galpha);
		double wwr = wwh + WLSR;
		double Vr = (-uur*RR + wwr*zz)/rr;
		double Vz = (uur*zz + wwr*RR)/rr;
	double Tvw = -sin(2.0*gbr)*sin(glr)*0.5;
	double Tuw = -sin(2.0*gbr)*cos(glr)*0.5;
		double varxyvepart = 0.5*sin(2.0*gbeta)*(Tuw*cos(galpha) - (Tvw*sin(galpha)));

//		vvh = vvr;
	if (wwh > WWHLIM - WLSR) {wwh = WWHLIM - WLSR;};
	if (wwh < -WWHLIM - WLSR) {wwh = -WWHLIM - WLSR;}
	if (VARXYCOND) {            // ((vvh > -150.0) && (vvh < 100.0)) {
	double x = Tvw*vvh + Tuw*uuh - WLSR*cos(gbr)*cos(gbr);
#ifdef NOVELUSE
        x = Tvw;
#endif
	ym += wwh;
	xm += x;
	varx += (x*x);
	vary += (wwh*wwh);
	varxy += (wwh*x);
	varxyrva -= (Tvw*Tvw + Tuw*Tuw)*rvaerr*rvaerr;
	varTvw += Tvw*Tvw;
	varTuw += Tuw*Tuw;
	varxyrvabase -= (Tvw*Tvw + Tuw*Tuw);
	varxypm += (Tvw*Tvw + Tuw*Tuw)*dab*dab*(pmberr*pmberr)*facr*facr;
	double vebm = fabs(varxyvepart);
	varxyvebase += varxyvepart;                // -0.5*(Tvw*sin(galpha) + Tuw*cos(galpha))*sin(2.0*gbeta);
	suu += uur*uur; svv += vvr*vvr;
	sww += wwh*wwh;
	um += uur; vm += vvr;
	uuhm += uuh;
	wm += wwh;
#ifdef USECARTVE
	suuve += uur*uur*vebm; svvve += vvr*vvr*vebm;
	swwve += wwh*wwh*vebm;
	umve += uur*vebm; vmve += vvr*vebm;
	wmve += wwh*vebm;
#endif
#ifdef USERADVE
	suuve += Vr*Vr*vebm; svvve += vvr*vvr*vebm;
	swwve += Vz*Vz*vebm;
	umve += Vr*vebm; vmve += vvr*vebm;
	wmve += Vz*vebm;
#endif
	denve += vebm;
	Teffm += Teff;
	varTeff += (Teff*Teff);
	loggm += logg;
	parpm += parp;
	varparp += parp*parp;
	varlogg += logg*logg;
	fehm += feh;
	varfeh += feh*feh;
        ordermean += orderm;
	den += 1.0;
	}
  }
  den = 1.0/den;
  denve = 1.0/denve;
  xm *= den; ym *= den; Teffm *= den; loggm *= den; fehm *= den; parpm *= den; ordermean *= den;
  um *= den; uuhm *= den; vm *= den; wm *= den; suu *= den; sww *= den; svv *= den;
  umve *= denve; vmve *= denve; wmve *= denve; suuve *= denve; swwve *= denve; svvve *= denve;
  varx *= den; vary *= den; varxy *= den; varxyrva *= den; varxypm *= den; varxyrvabase *= den;  varTvw *= den; varTuw *= den;
  varxyvebase *= den;
//  cout << "vars " << varx << " " << vary << " " << varxy << " " << varxyrva << " " << varxypm << " " << varxyrva/varx<< " " << varxy/varx << "\n";
  varTeff *= den; varlogg *= den; varfeh *= den; varparp *= den;
  varxy -= (xm*ym);
  varx -= (xm*xm);
  vary -= (ym*ym);
  suu -= (um*um); svv -= (vm*vm); sww -= (wm*wm);
  suuve -= (umve*umve); svvve -= (vmve*vmve); swwve -= (wmve*wmve);
  varxyve = varxyvebase*(swwve - suuve);
#ifdef CORRERR
  varxy -= (varxyrva + varxypm + varxyve+varxyrvabase*(sigbin*sigbin));
#endif
  double varxbb = varx + (varTvw + varTuw)*sww;
  varTeff -= (Teffm*Teffm); varlogg -= (loggm*loggm); varfeh -= (fehm*fehm); varparp -= (parpm*parpm);
  double beta = varxy/(varxbb);
  double alpha = ym - xm*beta;
  double varres = 0.0;
  den = 0.0;
  for (int ii = i1; ii < i2; ii++) {
   	double dist0 = star[ii][0];
	double dist1 = star[ii][1];
	double dist2 = star[ii][2];
	double distm = star[ii][3];
	double vhelio = star[ii][4];
	double rva = vhelio;
	double glr = star[ii][5];
   	double gbr = star[ii][6];
	double pml = star[ii][7];
	double pmb = star[ii][8];
	double feh = star[ii][9];
	double Teff = star[ii][10];
	double logg = star[ii][11];
	double dab = star[ii][DIAL]*factor;
        double vvh = (dab*pml*cos(glr)*facr)- (sin(gbr)*dab*pmb*sin(glr)*facr)  + cos(gbr)*rva*sin(glr);
        double uuh = -dab*pml*sin(glr)*facr - sin(gbr)*dab*pmb*cos(glr)*facr + cos(gbr)*rva*cos(glr);
       double wwh = sin(gbr)*rva + dab*pmb*facr*cos(gbr);
        double Tvw = -sin(2.0*gbr)*sin(glr)*0.5;
        double Tuw = -sin(2.0*gbr)*cos(glr)*0.5;
		double xx = RSUN - dab*cos(glr)*cos(gbr);
		double yy = dab*sin(glr)*cos(gbr);
		double RR = sqrt(xx*xx + yy*yy);
		double zz = ((dab*sin(gbr)) + ZOFFSET);
		double galpha = atan2(yy,xx);   // atan(yy/fabs(xx));
		double gbeta = atan(zz/RR);
		double varxyvepart = 0.5*sin(2.0*gbeta)*(Tuw*cos(galpha) - (Tvw*sin(galpha)));

	if (wwh > WWHLIM) {wwh = WWHLIM;};
	if (wwh < -WWHLIM) {wwh = -WWHLIM;}
	if (VARXYCOND) {    // ((vvh > -150.0) && (vvh < 100.0)) {
        double x = Tvw*vvh + Tuw*uuh - WLSR*cos(gbr)*cos(gbr);
#ifdef NOVELUSE
       x = sin(gbr)*cos(gbr)*sin(glr);
#endif
	double res = wwh - (alpha + beta*x);
	varres += res*res;
	den += 1.0;
	}
  }
  double betaerr = sqrt(((1.0/(den-2.0))*varres/((varx)*den)) + ((SYSERRRV*varxyrva*varxyrva + SYSERR*(varxypm*varxypm + varxyve*varxyve))/(varxbb*varxbb)));
  double alphaerr = sqrt((1.0/(den-2.0))*varres/den);
  if (mode == 3) {
    ostr << i1 << " " << i2 << " " << Teffm << " " << varTeff << " " << loggm << " " << varlogg << "  "  << fehm << " " << varfeh << " " << parpm << " " << varparp << " " << alpha << " " << alphaerr << " " << beta << " " << betaerr << " " << den << " " << varx << " " << xm << " " << varres << " " << varxy << " " << varxypm << " " << varxyrva << " " << ordermean << " " << varxyrvabase << " " << varxyrva/varxbb << " " << varxyrvabase/varxbb << " " << varxyvebase << " " << varxyvebase/varxbb << " " << varxyve << " " << varxyve/varxbb << " " << varxypm/varxbb << " " << suu << " " << svv << " " << sww << " " << um << " " << vm << " " << wm << " " << uuhm << " " << varTvw << " " << varxbb << "  " << suuve << " " << svvve << " " << swwve << " " << umve << " "<< vmve << " "<< wmve << "\n";
  }
  if (mode == 0) {
    return beta;
  }
  if (mode == 1) {
    return (beta + betaerr);
  }
  if (mode == 2) {
    return (beta - betaerr);
  }
}

double lindisthelpUW(double **star, int i1, int i2, double factor, int mode, double sigbin, std::ostream& ostr) {
  double facr = (parsec)*(PI)*0.001/(365.24*24.0*3600.0*3600.0*180.0);
  double ym = 0.0;
  double xm  = 0.0;
  double xym = 0.0;
  double varxy = 0.0;
  double varx = 0.0;
  double vary = 0.0;
  double den = 0.0;
  double Teffm = 0.0;
  double loggm = 0.0;
  double parpm = 0.0;
  double varparp = 0.0;
  double varTeff = 0.0;
  double varlogg = 0.0;
  double fehm = 0.0;
  double ordermean = 0.0;
  double varfeh = 0.0;
  double varxyrva = 0.0;
  double varxyrvabase = 0.0;
  double varxypm = 0.0;
  double varxyvebase = 0.0;
  double varxyve = 0.0;
  double varTvw = 0.0;
  double varTuw = 0.0;
  double suu = 0.0;
  double svv = 0.0;
  double sww = 0.0;
  double um = 0.0;
  double vm = 0.0;
  double wm = 0.0;
  double suuve = 0.0;
  double svvve = 0.0;
  double swwve = 0.0;
  double umve = 0.0;
  double vmve = 0.0;
  double wmve = 0.0;
  double denve = 0.0;
  double uuhm = 0.0;
  for (int ii = i1; ii < i2; ii++) {
    	double dist0 = star[ii][0];
	double dist1 = star[ii][1];
	double dist2 = star[ii][2];
	double distm = star[ii][3];
	double vhelio = star[ii][4];
	double rva = vhelio;
	double glr = star[ii][5];
   	double gbr = star[ii][6];
	double pml = star[ii][7];
	double pmb = star[ii][8];
	double feh = star[ii][9];
	double Teff = star[ii][10];
	double logg = star[ii][11];
	double rvaerr = star[ii][12];
	double pmlerr = star[ii][13];
	double pmberr = star[ii][14];
        double orderm = star[ii][ORDERING];
	double parp = star[ii][15]/star[ii][16];
	double dab = star[ii][DIAL]*factor;
        double vvh = (dab*pml*cos(glr)*facr)- (sin(gbr)*dab*pmb*sin(glr)*facr)  + cos(gbr)*rva*sin(glr);
        double uuh = -dab*pml*sin(glr)*facr - sin(gbr)*dab*pmb*cos(glr)*facr + cos(gbr)*rva*cos(glr);
        double wwh = sin(gbr)*rva + dab*pmb*facr*cos(gbr);

		double xx = RSUN - dab*cos(glr)*cos(gbr);
		double yy = dab*sin(glr)*cos(gbr);
		double RR = sqrt(xx*xx + yy*yy);
		double zz = (dab*sin(gbr) + ZOFFSET);
		double rr = sqrt(RR*RR + zz*zz);
		double galpha = atan2(yy,xx);                  //atan(yy/fabs(xx));
		double gbeta = atan(zz/RR);
		double uur = (uuh+ULSR)*cos(galpha) - (vvh+VSUN)*sin(galpha);
		double vvr = (vvh+VSUN)*cos(galpha) + (uuh+ULSR)*sin(galpha);
		double wwr = wwh + WLSR;
		double Vr = (-uur*RR + wwr*zz)/rr;
		double Vz = (uur*zz + wwr*RR)/rr;
	double Tvw = -sin(2.0*gbr)*sin(glr)*0.5;
	double Tuw = -sin(2.0*gbr)*cos(glr)*0.5;
	Tvw = 0.0;
	double varxyvepart = 0.5*(Tuw*cos(galpha))*sin(2.0*gbeta);
//		vvh = vvr;
	if (wwh > WWHLIM - WLSR) {wwh = WWHLIM - WLSR;};
	if (wwh < -WWHLIM - WLSR) {wwh = -WWHLIM - WLSR;}
	if (VARXYCOND) {    // ((vvh > -150.0) && (vvh < 100.0)) {
	double x = Tuw*uuh;
#ifdef NOVELUSE
        x = Tvw;
#endif
	ym += wwh;
	xm += x;
	varx += (x*x);
	vary += (wwh*wwh);
	varxy += (wwh*x);
	varxyrva -= (Tuw*Tuw)*rvaerr*rvaerr;
	varTvw += Tvw*Tvw;
	varTuw += Tuw*Tuw;
	varxyrvabase -= (Tuw*Tuw);
	varxypm += (Tuw*Tuw)*dab*dab*(pmberr*pmberr)*facr*facr;
        varxyvebase += varxyvepart;     //0.5*( Tuw*cos(galpha))*sin(2.0*gbeta);
	double vebm = fabs(varxyvepart);
//	varxyvebase -= 0.5*(Tvw*sin(galpha) - Tuw*cos(galpha))*sin(2.0*gbeta);
	suu += uur*uur; svv += vvr*vvr;
	sww += wwh*wwh;
	um += uur; vm += vvr;
	uuhm += uuh;
	wm += wwh;
#ifdef USECARTVE
	suuve += uur*uur*vebm; svvve += vvr*vvr*vebm;
	swwve += wwh*wwh*vebm;
	umve += uur*vebm; vmve += vvr*vebm;
	wmve += wwh*vebm;
#endif
#ifdef USERADVE
	suuve += Vr*Vr*vebm; svvve += vvr*vvr*vebm;
	swwve += Vz*Vz*vebm;
	umve += Vr*vebm; vmve += vvr*vebm;
	wmve += Vz*vebm;
#endif
	denve += vebm;
	Teffm += Teff;
	varTeff += (Teff*Teff);
	loggm += logg;
	parpm += parp;
	varparp += parp*parp;
	varlogg += logg*logg;
	fehm += feh;
	varfeh += feh*feh;
        ordermean += orderm;
	den += 1.0;
	}
  }
  den = 1.0/den;
  denve = 1.0/denve;
  xm *= den; ym *= den; Teffm *= den; loggm *= den; fehm *= den; parpm *= den; ordermean *= den;
  um *= den; uuhm *= den; vm *= den; wm *= den; suu *= den; sww *= den; svv *= den;
  umve *= denve; vmve *= denve; wmve *= denve; suuve *= denve; swwve *= denve; svvve *= denve;
  varx *= den; vary *= den; varxy *= den; varxyrva *= den; varxypm *= den; varxyrvabase *= den;  varTvw *= den; varTuw *= den;
  varxyvebase *= den;
//  cout << "vars " << varx << " " << vary << " " << varxy << " " << varxyrva << " " << varxypm << " " << varxyrva/varx<< " " << varxy/varx << "\n";
  varTeff *= den; varlogg *= den; varfeh *= den; varparp *= den;
  varxy -= (xm*ym);
  varx -= (xm*xm);
  vary -= (ym*ym);
  suu -= (um*um); svv -= (vm*vm); sww -= (wm*wm);
  suuve -= (umve*umve); svvve -= (vmve*vmve); swwve -= (wmve*wmve);
  varxyve = varxyvebase*(swwve - suuve);
#ifdef CORRERR
  varxy -= (varxyrva + varxypm + varxyve+varxyrvabase*(sigbin*sigbin));
#endif
  double varxbb = varx + (varTuw)*sww;
  varTeff -= (Teffm*Teffm); varlogg -= (loggm*loggm); varfeh -= (fehm*fehm); varparp -= (parpm*parpm);
  double beta = varxy/(varxbb);
  double alpha = ym - xm*beta;
  double varres = 0.0;
  den = 0.0;
  for (int ii = i1; ii < i2; ii++) {
   	double dist0 = star[ii][0];
	double dist1 = star[ii][1];
	double dist2 = star[ii][2];
	double distm = star[ii][3];
	double vhelio = star[ii][4];
	double rva = vhelio;
	double glr = star[ii][5];
   	double gbr = star[ii][6];
	double pml = star[ii][7];
	double pmb = star[ii][8];
	double feh = star[ii][9];
	double Teff = star[ii][10];
	double logg = star[ii][11];
	double dab = star[ii][DIAL]*factor;
        double vvh = (dab*pml*cos(glr)*facr)- (sin(gbr)*dab*pmb*sin(glr)*facr)  + cos(gbr)*rva*sin(glr);
        double uuh = -dab*pml*sin(glr)*facr - sin(gbr)*dab*pmb*cos(glr)*facr + cos(gbr)*rva*cos(glr);
       double wwh = sin(gbr)*rva + dab*pmb*facr*cos(gbr);
        double Tvw = -sin(2.0*gbr)*sin(glr)*0.5;
        double Tuw = -sin(2.0*gbr)*cos(glr)*0.5;
	Tvw = 0.0;
		double xx = RSUN - dab*cos(glr)*cos(gbr);
		double yy = dab*sin(glr)*cos(gbr);
		double RR = sqrt(xx*xx + yy*yy);
		double zz = ((dab*sin(gbr)) + ZOFFSET);
		double galpha = atan2(yy,xx);         //atan(yy/fabs(xx));
		double gbeta = atan(zz/RR);
	double varxyvepart = 0.5*(Tuw*cos(galpha))*sin(2.0*gbeta);
	if (wwh > WWHLIM - WLSR) {wwh = WWHLIM - WLSR;};
	if (wwh < -WWHLIM - WLSR) {wwh = -WWHLIM - WLSR;}
	if (VARXYCOND)  {   // ((vvh > -150.0) && (vvh < 100.0)) {
        double x = Tuw*uuh;
#ifdef NOVELUSE
       x = sin(gbr)*cos(gbr)*sin(glr);
#endif
	double res = wwh - (alpha + beta*x);
	varres += res*res;
	den += 1.0;
	}
  }
  double betaerr = sqrt(((1.0/(den-2.0))*varres/((varx)*den)) + ((SYSERRRV*varxyrva*varxyrva + SYSERR*(varxypm*varxypm + varxyve*varxyve))/(varxbb*varxbb)));
  double alphaerr = sqrt((1.0/(den-2.0))*varres/den);
  if (mode == 3) {
    ostr << i1 << " " << i2 << " " << Teffm << " " << varTeff << " " << loggm << " " << varlogg << "  "  << fehm << " " << varfeh << " " << parpm << " " << varparp << " " << alpha << " " << alphaerr << " " << beta << " " << betaerr << " " << den << " " << varx << " " << xm << " " << varres << " " << varxy << " " << varxypm << " " << varxyrva << " " << ordermean << " " << varxyrvabase << " " << varxyrva/varxbb << " " << varxyrvabase/varxbb << " " << varxyvebase << " " << varxyvebase/varxbb << " " << varxyve << " " << varxyve/varxbb << " " << varxypm/varxbb << " " << suu << " " << svv << " " << sww << " " << um << " " << vm << " " << wm << " " << uuhm << " " << varTvw << " " << varxbb << " " << suuve << " " << svvve << " " << swwve << " " << umve << " " << vmve << " " << wmve << "\n";
  }
  if (mode == 0) {
    return beta;
  }
  if (mode == 1) {
    return (beta + betaerr);
  }
  if (mode == 2) {
    return (beta - betaerr);
  }
}




void lindist(double **star, int i1, int i2, double sigbin, std::ostream& ostr) {
  double factor = 1.0;
  double back0 = -1.0e+10;
  double delta = 0.01;
  double back = back0;
  double factorsolve = 1.0;
  cout << "lindist call start \n";
//  ostringstream ostr;
  ostr << i1 << " " << i2 << " " << sigbin << " ";
  cout << "ostr first write done \n";
  for (int j = 0; j < 3;  j++) {
    back0 = -1.0e+10;
    back = back0;
    delta = 0.02;
    int flipflag = 0;
    if (factor < 0.1) {factor = 1.0;}
    while (((back0*back0) > 1.0e-11) && (factor > 1.0e-3) && (delta*delta > 1.0e-16) && (factor < 10.0)) {
     factor += delta;
#ifdef USEUV
     back = lindisthelpUVW(star, i1, i2, factor, j, sigbin, ostr);
#endif
#ifdef USEV
     back = lindisthelpVW(star, i1, i2, factor, j, sigbin, ostr);
#endif
#ifdef USEU
     back = lindisthelpUW(star, i1, i2, factor, j, sigbin, ostr);
#endif
     if ((back0*back) < 0.0) {   // || ((back*back) > (back0*back0))) {
       delta *= -0.4231;
     }
     else {
	if ((back*back) > (back0*back0)) {
		flipflag++;
		if (flipflag > 2) {
			delta *= -0.8;
		}
	}
	else {
		flipflag = 0;
	}
    }
//    if (i2 > 2000) {
//         cout << "factor " << " " << j << " " << back0 << " " << back << " " << delta << " " << factor << "\n";
//     }
     back0 = back;
    }
    if (j == 0) {
      factorsolve = factor;
    }
    ostr << j << " " << factor << " " << delta << " " << back << "  ";
  }
#ifdef USEUV
     back = lindisthelpUVW(star, i1, i2, factor, 3, sigbin, ostr);
#endif
#ifdef USEV
     back = lindisthelpVW(star, i1, i2, factor, 3, sigbin, ostr);
#endif
#ifdef USEU
     back = lindisthelpUW(star, i1, i2, factor, 3, sigbin, ostr);
#endif

}


void testfuncre(double **star, double *priorfit, double *priorfoot, int cii, double facr, double **priornew) {
        double par = 1.0/star[cii][2];
        double parerrin = par/(star[cii][21]);
	double parerr = parerrin;
#ifdef SETERR
	parerr = ERRORPLACE;
#endif

        double glr = star[cii][5];
        double gbr = star[cii][6];
        double rva = star[cii][4];
        double pml = star[cii][7];
        double pmb = star[cii][8];
        double s = 1.0/par;
         double dab = s;
         double dabsq = 0.0;
              dab = 0.0;
              double sumdab = 0.0;
              double precis = 0.0008;
	      double preclim = 1.0e-05;
        double parerrs = sqrt(parerr*parerr + 0.00*0.00);
	double s0 = 1.0/par;
	double parsq = 1.0/(parerrs*parerrs*2.0);
                double deltas1 = s0*precis*parerrs/par;
                double deltas2 = s0*precis*parerrs/par;
	double s1 = s0;
	double s2 = s0;
	int flag = 0;
	int flag1 = 0;
        int flag2 = 0;
        double xfac = -cos(glr)*cos(gbr);
        double yfac = sin(glr)*cos(gbr);
                double zfac = (sin(gbr));
        int counter = 0;
        int counter0 = 0;
        int pri1 = 0;
	int pri2 = 0;
	const double t03 = 1.0/0.3;
	const double t09 = 1.0/0.9;
	const double norm12 = 1.0/1.12;

                                double xx = RSUN + xfac*s1;
                                double yy = yfac*s1;
	                      double zz = fabs(zfac*s1 + ZOFFSET);
                                double R = hypot(xx, yy);           //sqrt(xx*xx + yy*yy);
                                double rr = hypot(R,zz);            //sqrt(R*R + zz*zz);
                                double par1 = 1.0/s1;
				double rrfac = RSUN/rr;
				double rho = exp((RSUN-R)*0.4)*norm12*(exp(-zz*t03) + 0.12*exp(-zz*t09)) + 0.001*rrfac*rrfac*sqrt(rrfac);   //(pow(rrfac, -2.5));
				double expo = -(par1-par)*(par1 - par)*parsq;
				double exph = 0.0;
//				if (expo > -300.0) {
					exph = exp(expo);
//				}
				    double sp0 = priorfoot[pri1];
				    double logs1 = log(s1);
				    while ((logs1 < sp0) && (pri1 > 0)) {
				        pri1--;
				        sp0 = priorfoot[pri1];
				    }
				    double sp1 = priorfoot[pri1+1];
				    while ((logs1 > sp1)) {
				        pri1++;
				        sp1 = priorfoot[pri1+1];
				    }
				    double facsp = (logs1 - sp0)/(sp1-sp0);
				    if (facsp < 0.0) {facsp = 0.0;}
				    if (pri1 >= PFMAX) {
				        cout << "WARNING PRI " << s1 << "\n";
				        pri1 = PFMAX - 1;
				        facsp = 1.0;
				    }
				    double facsp0 = 1.0-facsp;
				    double seldel = facsp0*priorfit[pri1] + facsp*priorfit[pri1+1];
				    double dfac00 = s1*s1*rho*selfunction(s1);      //*deltas1;
//				double sefus1 = selfunction(s1);
					double dfac = dfac00*exph*seldel;           // s1*s1*(exp(-s1/0.1) + 0.01*exp(-s1))*rho*exph*deltas1;
//					sumdab += dfac;
	double par2 = par1;
	double logs2 = logs1;
	double dfac1old0 = dfac;
	double dfac1old1 = dfac*s1;
	double dfac1old2 = dfac*s1*s1;
	double dfac2old0 = dfac;
	double dfac2old1 = dfac1old1;
	double dfac2old2 = dfac2old2;
	double dfac100old = dfac00; double dfac200old = dfac00;
	double dfac00m = dfac00; double dfacm = dfac; double dfac1 = dfac1old1; double dfac2 = dfac1old2;
	double dfacfact = 1.0; double dfac00fact = 1.0; double dfac1n = 1.0; double dfac2n = 1.0;

        while (flag ==0) {
		if (flag1 == 0) {
			s1 -= deltas1;
                        if (s1 > 0.0) {
                                xx = RSUN + xfac*s1;
                                yy = yfac*s1;
	                        zz = fabs(zfac*s1 + ZOFFSET);
                                R = hypot(xx, yy);           //sqrt(xx*xx + yy*yy);
                                rr = hypot(R,zz);            //sqrt(R*R + zz*zz);
                                par1 = 1.0/s1;
				rrfac = RSUN/rr;
				rho = exp((RSUN-R)*0.4)*norm12*(exp(-zz*t03) + 0.12*exp(-zz*t09)) + 0.001*rrfac*rrfac*sqrt(rrfac);   //(pow(rrfac, -2.5));
				expo = -(par1-par)*(par1 - par)*parsq;
				exph = 0.0;
					exph = exp(expo);
				sp0 = priorfoot[pri1];
				logs1 = log(s1);
				while ((logs1 < sp0) && (pri1 > 0)) {
				    pri1--;
				    sp0 = priorfoot[pri1];
				}
				sp1 = priorfoot[pri1+1];
				while ((logs1 > sp1)) {
				    pri1++;
				    sp1 = priorfoot[pri1+1];
				}
				facsp = (logs1 - sp0)/(sp1-sp0);
				if (facsp < 0.0) {facsp = 0.0;}
				if (pri1 >= PFMAX) {
				    cout << "WARNING PRI " << s1 << "\n";
				    pri1 = PFMAX - 1;
				    facsp = 1.0;
				}
				facsp0 = 1.0-facsp;
				seldel = facsp0*priorfit[pri1] + facsp*priorfit[pri1+1];
				dfac00 = s1*s1*rho*selfunction(s1);        //*deltas1;
//				double sefus1 = selfunction(s1);
				dfac = dfac00*exph*seldel;           // s1*s1*(exp(-s1/0.1) + 0.01*exp(-s1))*rho*exph*deltas1;
				dfacfact = (dfac1old0 + dfac)*deltas1;
				dfac00fact = (dfac100old + dfac00)*deltas1;
				dfac1n = dfac*s1;
				dfac2n = dfac*s1*s1;
				sumdab += dfacfact;
                	                priornew[pri1][0] += facsp0*dfac00fact;
                	                priornew[pri1][1] += facsp0*dfacfact;
                	                priornew[pri1+1][0] += facsp*dfac00fact;
                	                priornew[pri1+1][1] += facsp*dfacfact;
					dab += (dfac1old1 + dfac1n)*deltas1;
					dabsq += (dfac1old2 + dfac2n)*deltas1;

				  dfac1old0 = dfac;
				  dfac100old = dfac00;
				  dfac1old1 = dfac1n;
				  dfac1old2 = dfac2n;

		                      if (dfac < (preclim)*sumdab) {deltas1 *= 1.4;}
					if (deltas1 > 0.007) {
						deltas1 = 0.007;
					}
					counter++;
					counter0++;
					}
					else {
					flag1 = 1;
					if (flag2 == 1) {
						flag = 1;
					}
				}
			}
			if (flag2 == 0) {
  			    s2 += deltas2;
	                       if (s2 < 1.0e+6) {
                                xx = RSUN + xfac*s2;
                                yy = yfac*s2;
	                        zz = fabs(zfac*s2 + ZOFFSET);
                                R = hypot(xx, yy);           //sqrt(xx*xx + yy*yy);
                                rr = hypot(R,zz);            //sqrt(R*R + zz*zz);
                                par2 = 1.0/s2;
				rrfac = RSUN/rr;
				rho = exp((RSUN-R)*0.4)*norm12*(exp(-zz*t03) + 0.12*exp(-zz*t09)) + 0.001*rrfac*rrfac*sqrt(rrfac);   //(pow(rrfac, -2.5));
				expo = -(par2-par)*(par2 - par)*parsq;
				exph = exp(expo);
				sp0 = priorfoot[pri2];
				logs2 = log(s2);
				while ((logs2 < sp0) && (pri2 > 0)) {
				        pri2--;
				        sp0 = priorfoot[pri2];
				}
				sp1 = priorfoot[pri2+1];
				while ((logs2 > sp1)) {
				    pri2++;
				    sp1 = priorfoot[pri2+1];
				}
				facsp = (logs2 - sp0)/(sp1-sp0);
				if (facsp < 0.0) {facsp = 0.0;}
				if (pri2 >= PFMAX) {
				    cout << "WARNING PRI " << s2 << "\n";
				    pri2 = PFMAX - 1;
				    facsp = 1.0;
				}
				facsp0 = 1.0-facsp;
				seldel = facsp0*priorfit[pri2] + facsp*priorfit[pri2+1];

				dfac00 = s2*s2*rho*selfunction(s2);        //*deltas1;
				dfac = dfac00*exph*seldel;           // s1*s1*(exp(-s1/0.1) + 0.01*exp(-s1))*rho*exph*deltas1;
				dfacfact = (dfac2old0 + dfac)*deltas2;
				dfac00fact = (dfac200old + dfac00)*deltas2;
				dfac1n = dfac*s2;
				dfac2n = dfac*s2*s2;
				sumdab += dfacfact;
               	                priornew[pri2][0] += facsp0*dfac00fact;
                	        priornew[pri2][1] += facsp0*dfacfact;
                	        priornew[pri2+1][0] += facsp*dfac00fact;
                	        priornew[pri2+1][1] += facsp*dfacfact;
				dab += (dfac2old1 + dfac1n)*deltas2;
				dabsq += (dfac2old2 + dfac2n)*deltas2;

				  dfac2old0 = dfac;
				  dfac200old = dfac00;
				  dfac2old1 = dfac1n;
				  dfac2old2 = dfac2n;

		                      if (dfac < (preclim)*sumdab) {deltas2 *= 1.6;}
					counter++;
					}
					else {
						flag2 = 1;
						if (flag1 == 1) {
						flag = 1;
					}
				}    //end if s2 < 1.0e+6
			}    //end if flag2
//		cout << " " << par/parerr << " " << par << " " << parerr << " " << s1 << " " << s2 << " " << deltas1 << " " << deltas2 << " " << dab << " " << sumdab << " " << dab/sumdab << "\n";
//	}
//	cout << " " << counter << " " << counter0 << " " << par/parerr << " " << par << " " << parerr << " " << dab/sumdab << " " << dab << " " << sumdab << " " << deltas1 << " " << deltas2 << "\n";

	}    //end if flag
              double dab2 = 1.0*dab/sumdab;
		dab = dab2;
		if (dab < 1.0e-10) {
			cout << "WARNING dab " << dab << " " << par << " " << parerr << "\n";
		}
#ifdef USEPAR
		dab = 1.0/par;
#endif

//		dab = dist;

        double vvh = (dab*pml*cos(glr)*facr)- (sin(gbr)*dab*pmb*sin(glr)*facr)  + cos(gbr)*rva*sin(glr);
        double uuh = -dab*pml*sin(glr)*facr - sin(gbr)*dab*pmb*cos(glr)*facr + cos(gbr)*rva*cos(glr);
       double wwh = sin(gbr)*rva + dab*pmb*facr*cos(gbr);
	xx = RSUN - dab*cos(glr)*cos(gbr);
	yy = dab*sin(glr)*cos(gbr);
	double RR = sqrt(xx*xx + yy*yy);
	zz = dab*sin(gbr);
	double galpha = atan2(yy,xx);  //atan(yy/fabs(xx));
	double gbeta = atan(zz/RR);
	double uur = (uuh+ULSR)*cos(galpha) - (vvh+VSUN)*sin(galpha);
	double vvr = (vvh+VSUN)*cos(galpha) + (uuh+ULSR)*sin(galpha);
	if (1 == 1) {
//	if ((vvr > -50.0 && vvr < 300.0) && (dab > 0.03) && fabs(uur) < 450.0 && rvaerr < 10.0 && fabs(rva) < 400.0 && (ADDCOND)) {
//&& (fabs(dab-0.32) > 0.03)) {              //(fabs(RA - 150.0) > 30.0 || fabs(DEC + 60.0) > 30.0)) {
//&& (fabs(gllmc - glr) > adistlmc || fabs(gblmc - gbr) > adistlmc) && (fabs(glsmc - glr) > adistlmc || fabs(gbsmc - gbr) > adistsmc)) {
        star[cii][0] = dab;
	star[cii][1] = dab;        //1.0/par;
//	star[cii][2] = 1.0/par;
//	double dist2 = star[ii][2];
//	double distm = star[ii][3];
	star[cii][19] = vvh;
	star[cii][22] = xx;
	star[cii][23] = yy;
	star[cii][24] = zz;
	star[cii][25] = fabs(zz);
	star[cii][26] = RR;
	star[cii][27] = uur;
	star[cii][28] = vvr;
	star[cii][29] = vvr*RR;
	star[cii][20] = wwh;

//	if (1 == 1) {            //(dab > 0.25 && dab < 0.45) {
//		cout << dab << " " << glr << " " << gbr << " " << star[cii][18] << " " << uuh << " " << vvh << " " << wwh << " " << rva << " " << flagc0 << " " << flagc1 << " " << par << " " << parerr << "\n";
//	}
	}  // end if fflag
}






void priorrefit(double **star, int starti, int fini, double facr, int index, double *priorfoot, double **priorfitin) {
    int pfmax = PFMAX;
    double smaxb = 3.0;
//    double *priorfoot = new double[pfmax];
    double *priorfit = new double[pfmax];
   double **priorfitnew = new double*[pfmax];
    for (int ii = 0; ii < pfmax; ii++) {
//        priorfoot[ii] = 0.0;
        priorfit[ii] = 1.0;
        priorfitnew[ii] = new double[2];
        priorfitnew[ii][0] = 0.0;
        priorfitnew[ii][1] = 0.0;
    }
    double smax = SMAX;
//    for (int ii = 0; ii < pfmax; ii++) {
//        priorfoot[ii] = log(smax*(((double)(ii))+0.5)/((double)(pfmax)));                 //sqrt(((double)(ii+1))/((double)(pfmax))));
//    }
//    priorfoot[pfmax - 1] = 1.0e+100;
    for (int nn = 0; nn < 2; nn++) {
        for (int ii = starti; ii < fini; ii++) {
//	  if ((ii/50000)*50000 == ii) {
//		    cout << "do testfuncre " << ii << " " << index << "\n";
//	    }
	    double **priorfitinst = new double*[pfmax];
	    for (int jj = 0; jj < pfmax; jj++) {
		priorfitinst[jj] = new double[2];
		priorfitinst[jj][0] = 0.0;
		priorfitinst[jj][1] = 0.0;
	    }
            testfuncre(star, priorfit, priorfoot, ii, facr, priorfitinst);
	    double sum0 = 0.0; double sum1 = 0.0;
	    for (int jj = 0; jj < pfmax; jj++) {
		sum0 += priorfitinst[jj][0]; sum1 += priorfitinst[jj][1];
	    }
	    sum0 = 1.0/sum0; sum1 = 1.0/sum1;
	    for (int jj = 0; jj < pfmax; jj++) {
		priorfitnew[jj][0] += priorfitinst[jj][0]*sum0;
		priorfitnew[jj][1] += priorfitinst[jj][1]*sum1;
		delete [] priorfitinst[jj];
	    }
	    delete [] priorfitinst;
        }
        for (int ii = 0; ii < pfmax; ii++) {
            if (priorfitnew[ii][0] > 0.0) {
                priorfit[ii] = priorfitnew[ii][1]/priorfitnew[ii][0];
            }
            else {
                priorfit[ii] = 1.0;
            }
        }
	for (int ii = 0; ii < pfmax; ii++) {
		if (priorfoot[ii] > log(smaxb)) {
			if (ii == 0) {
				priorfit[ii] = 1.0;
			}
			else {
				priorfit[ii] = priorfit[ii - 1];
			}
		}
	}
    }
    ostringstream nostr;
    nostr << "output/priorrefit" << SUFFIX << index;
    string nstri = nostr.str();
    char *nch = &nstri[0];
    FILE *fout = fopen(nch, "w");
    ostringstream ostr;
    for (int ii = 0; ii < pfmax; ii++) {
        ostr << ii << " " << priorfoot[ii] << " " << priorfit[ii] << " " << priorfitnew[ii][0] << " " << priorfitnew[ii][1] << "\n";
    }
    string stri = ostr.str();
    char *ch = &stri[0];
    fputs(ch, fout);
    fflush(fout);
    fclose(fout);
    for (int ii = 0; ii < pfmax; ii++) {
        delete [] priorfitnew[ii];
    }
    delete [] priorfitnew;
    delete [] priorfit;
}




void quaddistb(std::ostream& ostr) {            //double **star, int *index) {       //, int *index, double *distr, int sam) {  // std::ostream& ostr) {
   int x = 0;
}

#ifdef REDOPRIOR
void quaddist(double **star, int *index, double *distr, int sam, std::ostream& ostr, double facr, double *priorfoot, double **priorfit) {
#else
void quaddist(double **star, int *index, double *distr, int sam, std::ostream& ostr, double facr) {
#endif
    cout << "quaddist create array " << sam << "\n";
//    ostringstream ostr;
    double **starsorted = new double*[sam];
    int indic = index[0];
    for (int ii = 0; ii < sam; ii++) {
      starsorted[ii] = new double[30];
      int iori = index[ii];
      for (int jj = 0; jj < 30; jj++) {
	  starsorted[ii][jj] = star[iori][jj];
      }
//      if (sam < SAMB) {
//	      cout << "TESTARRAY " << sam << " " << ii << " " << iori << " " << distr[ii] << " " << star[iori][ORDERING] << "\n";
//	}
    }
    cout << "call lindist \n";
    double sigbin = 0.0;
#ifdef REDOPRIOR
    priorrefit(starsorted, 0, (sam-1), facr, indic, priorfoot, priorfit);
#endif
    lindist(starsorted, 0, (sam-1), sigbin, ostr);
    for (int ii = 0; ii < sam; ii++) {
	delete [] starsorted[ii];
    }
    delete [] starsorted;

}





// && (par/parerr > 8.5) && (sqrt(0.5*(pmraerr*pmraerr + pmdecerr*pmdecerr)) < 3.0) && (rvaerr < 20.0) && (fabs(rva) < 300.0) && (sn1 > 25.0)) {


#ifdef TESTDISTPRIOR
void testfunc(double **star, double par, double parerrin, double glr, double gbr, double Teff, double logg, double feh, double Tefferr, double loggerr, double feherr, double dist, double pml, double pmb, int cii, double pmraerr, double pmdecerr, double rvaerr, double rva, double snr, double facr, double ***dfall, int ncounter, double *priorfoot, double **priornew) {

#ifdef REDOFULLPRIOR
	for (int ii = 0; ii < PFMAX; ii++) {
		priornew[ii][0] = 0.0; priornew[ii][1] = 0.0;
	}
#endif

#else
void testfunc(double **star, double par, double parerrin, double glr, double gbr, double Teff, double logg, double feh, double Tefferr, double loggerr, double feherr, double dist, double pml, double pmb, int cii, double pmraerr, double pmdecerr, double rvaerr, double rva, double snr, double facr) {
#endif
 	double parerr = parerrin;
#ifdef SETERR
	parerr = sqrt(ERRFACT*parerrin*parerrin + ERRORPLACE*ERRORPLACE);                      //ERRORPLACE;
#endif

//	double dab = 0.0;
            /*
	double sumdab = 0.0;
	dab = 0.0;
//	cout << par/parerr << " " << par << " " << parerr << "\n";
//	parerr = sqrt(parerr*parerr + 0.09);
        if (par > parerr*3.9) {
        double s0 = 1.0/par;
        double parsq = 1.0/(parerr*parerr*2.0);
        double deltas1 = s0*0.002*parerr/par;
        double deltas2 = s0*0.002*parerr/par;
        double s1 = s0;
        double s2 = s0;
        int flag = 0;
        int flag1 = 0;
        int flag2 = 0;
        double xfac = -cos(glr)*cos(gbr);
        double yfac = sin(glr)*cos(gbr);
        double zfac = fabs(sin(gbr));
        int counter = 0;
        int counter0 = 0;
*/


              double s = 1.0/par;
              double dab = s;
              double dabsq = 0.0;
              dab = 0.0;
              double sumdab = 0.0;
              double precis = 0.0007;
	      double preclim = 1.0e-05;
	if (par > parerrin*PARQUALLIMITC) {
        double parerrs = sqrt(parerr*parerr + 0.00*0.00);
	double s0 = 1.0/par;
	double parsq = 1.0/(parerrs*parerrs*2.0);
                double deltas1 = s0*precis*parerrs/par;
                double deltas2 = s0*precis*parerrs/par;
	double s1 = s0;
	double s2 = s0;
	int flag = 0;
	int flag1 = 0;
        int flag2 = 0;
        double xfac = -cos(glr)*cos(gbr);
        double yfac = sin(glr)*cos(gbr);
                double zfac = (sin(gbr));
        int counter = 0;
        int counter0 = 0;


        int pri1 = 0;
	int pri2 = 0;
	const double t03 = 1.0/0.3;
	const double t09 = 1.0/0.9;
	const double norm12 = 1.0/1.12;

                                double xx = RSUN + xfac*s1;
                                double yy = yfac*s1;
	                      double zz = fabs(zfac*s1 + ZOFFSET);
                                double R = hypot(xx, yy);           //sqrt(xx*xx + yy*yy);
                                double rr = hypot(R,zz);            //sqrt(R*R + zz*zz);
                                double par1 = 1.0/s1;
				double rrfac = RSUN/rr;
				double rho = exp((RSUN-R)*0.4)*norm12*(exp(-zz*t03) + 0.12*exp(-zz*t09)) + 0.001*rrfac*rrfac*sqrt(rrfac);   //(pow(rrfac, -2.5));
				double expo = -(par1-par)*(par1 - par)*parsq;
				double exph = 0.0;
//				if (expo > -300.0) {
					exph = exp(expo);
//				}


#ifdef REDOFULLPRIOR
                                   double sp0 = priorfoot[pri1];
				    double logs1 = log(s1);
				    while ((logs1 < sp0) && (pri1 > 0)) {
				        pri1--;
				        sp0 = priorfoot[pri1];
				    }
				    double sp1 = priorfoot[pri1+1];
				    while ((logs1 > sp1)) {
				        pri1++;
				        sp1 = priorfoot[pri1+1];
				    }
				    double facsp = (logs1 - sp0)/(sp1-sp0);
				    if (facsp < 0.0) {facsp = 0.0;}
				    if (pri1 >= PFMAX) {
				        cout << "WARNING PRI " << s1 << "\n";
				        pri1 = PFMAX - 1;
				        facsp = 1.0;
				    }
				    double facsp0 = 1.0-facsp;
//				    double seldel = facsp0*priorfit[pri1] + facsp*priorfit[pri1+1];

	double logs2 = logs1;

#endif

				    double dfac00 = s1*s1*rho*selfunction(s1);      //*deltas1;
//				double sefus1 = selfunction(s1);
					double dfac = dfac00*exph;           // s1*s1*(exp(-s1/0.1) + 0.01*exp(-s1))*rho*exph*deltas1;
//					sumdab += dfac;
	double par2 = par1;
	double dfac1old0 = dfac;
	double dfac1old1 = dfac*s1;
	double dfac1old2 = dfac*s1*s1;
	double dfac2old0 = dfac;
	double dfac2old1 = dfac1old1;
	double dfac2old2 = dfac2old2;
	double dfac100old = dfac00; double dfac200old = dfac00;
	double dfac00m = dfac00; double dfacm = dfac; double dfac1 = dfac1old1; double dfac2 = dfac1old2;
	double dfacfact = 1.0; double dfac00fact = 1.0; double dfac1n = 1.0; double dfac2n = 1.0;





        while (flag ==0) {
		if (flag1 == 0) {
			s1 -= deltas1;
                        if (s1 > 0.0) {

                                 xx = RSUN + xfac*s1;
                                yy = yfac*s1;
	                        zz = fabs(zfac*s1 + ZOFFSET);
                                R = hypot(xx, yy);           //sqrt(xx*xx + yy*yy);
                                rr = hypot(R,zz);            //sqrt(R*R + zz*zz);
                                par1 = 1.0/s1;
				rrfac = RSUN/rr;
				rho = exp((RSUN-R)*0.4)*norm12*(exp(-zz*t03) + 0.12*exp(-zz*t09)) + 0.001*rrfac*rrfac*sqrt(rrfac);   //(pow(rrfac, -2.5));
				expo = -(par1-par)*(par1 - par)*parsq;
				exph = 0.0;
					exph = exp(expo);

#ifdef REDOFULLPRIOR
				sp0 = priorfoot[pri1];
				logs1 = log(s1);
				while ((logs1 < sp0) && (pri1 > 0)) {
				    pri1--;
				    sp0 = priorfoot[pri1];
				}
				sp1 = priorfoot[pri1+1];
				while ((logs1 > sp1)) {
				    pri1++;
				    sp1 = priorfoot[pri1+1];
				}
				facsp = (logs1 - sp0)/(sp1-sp0);
				if (facsp < 0.0) {facsp = 0.0;}
				if (pri1 >= PFMAX) {
				    cout << "WARNING PRI " << s1 << "\n";
				    pri1 = PFMAX - 1;
				    facsp = 1.0;
				}
				facsp0 = 1.0-facsp;
//				seldel = facsp0*priorfit[pri1] + facsp*priorfit[pri1+1];
#endif

				dfac00 = s1*s1*rho*selfunction(s1);        //*deltas1;
//				double sefus1 = selfunction(s1);
				dfac = dfac00*exph;          //*seldel;           // s1*s1*(exp(-s1/0.1) + 0.01*exp(-s1))*rho*exph*deltas1;
				dfacfact = (dfac1old0 + dfac)*deltas1;
				dfac00fact = (dfac100old + dfac00)*deltas1;
				dfac1n = dfac*s1;
				dfac2n = dfac*s1*s1;
				sumdab += dfacfact;

#ifdef REDOFULLPRIOR
                	                priornew[pri1][0] += facsp0*dfac00fact;
                	                priornew[pri1][1] += facsp0*dfacfact;
                	                priornew[pri1+1][0] += facsp*dfac00fact;
                	                priornew[pri1+1][1] += facsp*dfacfact;
#endif

					dab += (dfac1old1 + dfac1n)*deltas1;
					dabsq += (dfac1old2 + dfac2n)*deltas1;

				  dfac1old0 = dfac;
				  dfac100old = dfac00;
				  dfac1old1 = dfac1n;
				  dfac1old2 = dfac2n;

		                      if (dfac < (preclim)*sumdab) {deltas1 *= 1.4;}
					if (deltas1 > 0.007) {
						deltas1 = 0.007;
					}
					counter++;
					counter0++;
					}
					else {
					flag1 = 1;
					if (flag2 == 1) {
						flag = 1;
					}
			}
		}
		if (flag2 == 0) {
			s2 += deltas2;
                        if (s2 < 1.0e+6) {

                                xx = RSUN + xfac*s2;
                                yy = yfac*s2;
	                        zz = fabs(zfac*s2 + ZOFFSET);
                                R = hypot(xx, yy);           //sqrt(xx*xx + yy*yy);
                                rr = hypot(R,zz);            //sqrt(R*R + zz*zz);
                                par2 = 1.0/s2;
				rrfac = RSUN/rr;
				rho = exp((RSUN-R)*0.4)*norm12*(exp(-zz*t03) + 0.12*exp(-zz*t09)) + 0.001*rrfac*rrfac*sqrt(rrfac);   //(pow(rrfac, -2.5));
				expo = -(par2-par)*(par2 - par)*parsq;
				exph = exp(expo);

#ifdef REDOFULLPRIOR
				sp0 = priorfoot[pri2];
				logs2 = log(s2);
				while ((logs2 < sp0) && (pri2 > 0)) {
				        pri2--;
				        sp0 = priorfoot[pri2];
				}
				sp1 = priorfoot[pri2+1];
				while ((logs2 > sp1)) {
				    pri2++;
				    sp1 = priorfoot[pri2+1];
				}
				facsp = (logs2 - sp0)/(sp1-sp0);
				if (facsp < 0.0) {facsp = 0.0;}
				if (pri2 >= PFMAX) {
				    cout << "WARNING PRI " << s2 << "\n";
				    pri2 = PFMAX - 1;
				    facsp = 1.0;
				}
				facsp0 = 1.0-facsp;
//				seldel = facsp0*priorfit[pri2] + facsp*priorfit[pri2+1];
#endif

				dfac00 = s2*s2*rho*selfunction(s2);        //*deltas1;
				dfac = dfac00*exph;   //*seldel;           // s1*s1*(exp(-s1/0.1) + 0.01*exp(-s1))*rho*exph*deltas1;
				dfacfact = (dfac2old0 + dfac)*deltas2;
				dfac00fact = (dfac200old + dfac00)*deltas2;
				dfac1n = dfac*s2;
				dfac2n = dfac*s2*s2;
				sumdab += dfacfact;

#ifdef REDOFULLPRIOR
               	                priornew[pri2][0] += facsp0*dfac00fact;
                	        priornew[pri2][1] += facsp0*dfacfact;
                	        priornew[pri2+1][0] += facsp*dfac00fact;
                	        priornew[pri2+1][1] += facsp*dfacfact;
#endif

                                dab += (dfac2old1 + dfac1n)*deltas2;
				dabsq += (dfac2old2 + dfac2n)*deltas2;

				  dfac2old0 = dfac;
				  dfac200old = dfac00;
				  dfac2old1 = dfac1n;
				  dfac2old2 = dfac2n;

		                      if (dfac < (preclim)*sumdab) {deltas2 *= 1.6;}
					counter++;
					}
					else {
						flag2 = 1;
						if (flag1 == 1) {
						flag = 1;
					}
				}    //end if s2 < 1.0e+6
		}
//		cout << " " << par/parerr << " " << par << " " << parerr << " " << s1 << " " << s2 << " " << deltas1 << " " << deltas2 << " " << dab << " " << sumdab << " " << dab/sumdab << "\n";
	}
//	cout << " " << counter << " " << counter0 << " " << par/parerr << " " << par << " " << parerr << " " << dab/sumdab << " " << dab << " " << sumdab << " " << deltas1 << " " << deltas2 << "\n";

	}
              double dab2 = 1.0*dab/sumdab;
		dab = dab2;
		if (dab < 1.0e-10) {
			cout << "WARNING dab " << dab << " " << par << " " << parerr << "\n";
	}
#ifdef USEPAR
		dab = 1.0/par;
#endif

//		dab = dist;
#ifdef TESTDISTPRIOR
      	double xfac = -cos(glr)*cos(gbr);
	double yfac = sin(glr)*cos(gbr);
		double zfac = fabs(sin(gbr));
            double ddf[4][PRIORLENGTH];
            for (int si = 0; si < PRIORLENGTH; si++) {
                ddf[0][si] = 0.0; ddf[1][si] = 0.0; ddf[2][si] = 0.0; ddf[3][si] = 0.0;
                double s1 = 0.01*((double)(si) + 0.5);
                double xx = RSUN + xfac*s1;
		double yy = yfac*s1;
		  double zz = zfac*s1;
		double R = sqrt(xx*xx + yy*yy);
		double rr = sqrt(R*R + zz*zz);
		double par1 = 1.0/s1;
		double rho = exp((RSUN-R)/2.5)*(1.0/1.12)*(exp(-zz/0.3) + 0.12*exp(-zz/0.9)) + 0.001*(pow(rr, -2.5)/(pow(RSUN, -2.5)));
                double dfac = s1*s1*rho;
                double dfac2 = dfac*selfunction(s1);
                double dfac3 = dfac*(exp(-s1/0.1))*s1;
                double dfac4 = dfac*(exp(-s1/0.15));
                ddf[0][si] += dfac;
                ddf[1][si] += dfac2;
                ddf[2][si] += dfac3;
                ddf[3][si] += dfac4;
            }
            double ddfsum[4];
            ddfsum[0] = 0.0; ddfsum[1] = 0.0; ddfsum[2] = 0.0; ddfsum[3] = 0.0;
            for (int si = 0; si < PRIORLENGTH; si++) {
                ddfsum[0] += ddf[0][si];
                ddfsum[1] += ddf[1][si];
                ddfsum[2] += ddf[2][si];
                ddfsum[3] += ddf[3][si];
            }
            ddfsum[0] = 1.0/ddfsum[0];
            ddfsum[1] = 1.0/ddfsum[1]; ddfsum[2] = 1.0/ddfsum[2]; ddfsum[3] = 1.0/ddfsum[3];
            for (int si = 0; si < PRIORLENGTH; si++) {
                ddf[0][si] *= ddfsum[0];
                ddf[1][si] *= ddfsum[1];
		  ddf[2][si] *= ddfsum[2]; ddf[3][si] *= ddfsum[3];
            }
            for (int si = 0; si < PRIORLENGTH; si++) {
		  dfall[0][si][ncounter] += ddf[0][si];
		  dfall[1][si][ncounter] += ddf[1][si]; dfall[2][si][ncounter] += ddf[2][si]; dfall[3][si][ncounter] += ddf[3][si];
            }
            int disti = (int)(dab*100.0 + 0.5);
            if (disti < 0) {disti = 0;}
            if (disti >= PRIORLENGTH) {disti = PRIORLENGTH-1;}
		dfall[4][disti][ncounter] += 1.0;
#endif

        double vvh = (dab*pml*cos(glr)*facr)- (sin(gbr)*dab*pmb*sin(glr)*facr)  + cos(gbr)*rva*sin(glr);
        double uuh = -dab*pml*sin(glr)*facr - sin(gbr)*dab*pmb*cos(glr)*facr + cos(gbr)*rva*cos(glr);
       double wwh = sin(gbr)*rva + dab*pmb*facr*cos(gbr);
	double xx = RSUN - dab*cos(glr)*cos(gbr);
	double yy = dab*sin(glr)*cos(gbr);
	double RR = sqrt(xx*xx + yy*yy);
	double zz = dab*sin(gbr);
	double galpha = atan2(yy,xx);  //atan(yy/fabs(xx));
	double gbeta = atan(zz/RR);
	double uur = (uuh+ULSR)*cos(galpha) - (vvh+VSUN)*sin(galpha);
	double vvr = (vvh+VSUN)*cos(galpha) + (uuh+ULSR)*sin(galpha);
	double gllmc = 280.46*PI/180.0; double glsmc = 302.80*PI/180.0;
	double gblmc = -32.89*PI/180.0; double gbsmc = -44.33*PI/180.0;
	double adistlmc = 7.0*PI/180.0; double adistsmc = 4.0*PI/180.0;    //7, 3
//	cout << dab << " " << RR*vvr << " " << glr << " " << gbr << " " << uur << " " << vvr << " " << wwh << " " << rva << " " << rvaerr << " " << dist << " " <<  logg << " " << Teff << " " << feh << "\n";

	if (1 == 1) {
//	if ((vvr > -50.0 && vvr < 300.0) && (dab > 0.03) && fabs(uur) < 450.0 && rvaerr < 10.0 && fabs(rva) < 400.0 && (ADDCOND)) {
//&& (fabs(dab-0.32) > 0.03)) {              //(fabs(RA - 150.0) > 30.0 || fabs(DEC + 60.0) > 30.0)) {
//&& (fabs(gllmc - glr) > adistlmc || fabs(gblmc - gbr) > adistlmc) && (fabs(glsmc - glr) > adistlmc || fabs(gbsmc - gbr) > adistsmc)) {
        star[cii][0] = dab;
	star[cii][1] = dab;        //1.0/par;
	star[cii][2] = 1.0/par;
//	double dist2 = star[ii][2];
//	double distm = star[ii][3];
	star[cii][4] = rva;
	star[cii][5] = glr;
   	star[cii][6] = gbr;
	star[cii][7] = pml;
	star[cii][8] = pmb;
	star[cii][9] = feh;
	star[cii][10] = feherr;
	star[cii][11] = snr;
	star[cii][12] = rvaerr;
	star[cii][13] = sqrt(0.5*(pmraerr*pmraerr + pmdecerr*pmdecerr)); //pmlerr
	star[cii][14] = sqrt(0.5*(pmraerr*pmraerr + pmdecerr*pmdecerr)); //pmberr
	star[cii][15] = logg;
	star[cii][16] = Teff;                //Teff;
	star[cii][17] = loggerr;
	star[cii][18] = fabs(gbr);               //fabs(sin(glr)*sin(gbr)*cos(gbr));
	star[cii][19] = vvh;
	star[cii][21] = par/parerrin;
	star[cii][22] = xx;
	star[cii][23] = yy;
	star[cii][24] = zz;
	star[cii][25] = fabs(zz);
	star[cii][26] = RR;
	star[cii][27] = uur;
	star[cii][28] = vvr;
	star[cii][29] = vvr*RR;
	star[cii][20] = wwh;

//	if (1 == 1) {            //(dab > 0.25 && dab < 0.45) {
//		cout << dab << " " << glr << " " << gbr << " " << star[cii][18] << " " << uuh << " " << vvh << " " << wwh << " " << rva << " " << flagc0 << " " << flagc1 << " " << par << " " << parerr << "\n";
//	}
	}  // end if fflag
	}




int main() {

  cout << "start SUFFIX " << SUFFIX << "\n";

  std::thread persei[NTHREADS];
  int besetzt[NTHREADS];
  for (int n = 0; n < NTHREADS; n++) {
	besetzt[n] = 0;
        }
  int ncounter = 0;

#ifdef TESTDISTPRIOR
  double smax = SMAX;
  double ***dfall = new double**[5];
  for (int ii = 0; ii < 5; ii++) {
    dfall[ii] = new double*[PRIORLENGTH];
    for (int si = 0; si < PRIORLENGTH; si++) {
     dfall[ii][si] = new double[NTHREADS];
     for (int sis = 0; sis < NTHREADS; sis++) {
      	dfall[ii][si][sis] = 0.0;
	}
    }
  }
  double **priorfit = new double*[PFMAX];
  double ***priorfitinst = new double**[NTHREADS];
  double *priorfoot = new double[PFMAX];
  for (int jj = 0; jj < NTHREADS; jj++) {
	priorfitinst[jj] = new double*[PFMAX];
	for (int ii = 0; ii < PFMAX; ii++) {
		priorfitinst[jj][ii] = new double[2];
		priorfitinst[jj][ii][0] = 0.0; priorfitinst[jj][ii][1] = 0.0;
	}
  }
  for (int ii = 0; ii < PFMAX; ii++) {
        priorfoot[ii] = log(smax*(((double)(ii))+0.5)/((double)(PFMAX)));
//	priorfoot[ii] = 0.0;
        priorfit[ii] = new double[2];
        priorfit[ii][0] = 0.0; priorfit[ii][1] = 0.0;
  }
  priorfoot[PFMAX-1] = 1.0e+100;
#endif

  double facr = (parsec)*(PI)*0.001/(365.2425*24.0*3600.0*3600.0*180.0);
  double expos[3000];
  double sume = 0.0;
  for (int ii = 0; ii < 3000; ii++) {
	double dii = ((double)(ii-1500))/200.0;
        expos[ii] = exp(-dii*dii/2.0);
	sume += expos[ii];
        }
  sume = 1.0/sume;
  for (int ii = 0; ii < 3000; ii++) {
	expos[ii] *= sume;
        }
  double **star = new double*[STARN];
  double **starsort = new double*[STARN];
  for (int ii = 0; ii < STARN; ii++) {
    star[ii] = new double[30];
    starsort[ii] = new double[30];
    for (int jj = 0; jj < 30; jj++) {
      star[ii][jj] = 0.0;
      starsort[ii][jj] = 0.0;
	}
  }
         int cii = 0;
        for (int iif = 0; iif < 8; iif++) {               //8; iif++) {
	cout << "open the file \n";
	ostringstream ostrif;
	ostrif << "../rvset/gaiarv" << iif << ".csv";
	string strif = ostrif.str();
	char *chif = &strif[0];
	FILE *troja;
	troja = fopen(chif, "r");
	char **pointit = NULL;
	char buffer[32768];
        double dab = 0.0;

	cout << "opened the file \n";
        fgets(buffer, 32767, troja);
	cout << "first line read properly \n";
	int breaker = 0;
        while (fgets(buffer, 32767, troja)/* && (breaker == 0)*/) {
	    int nn = 0;
/*
	    if (cii > 50000) {
		breaker = 1;
	    }
*/

	    int readflag = 0;
	    while (buffer[nn] == ' ') {
		 nn++;
	    } //solution_id
	    while (buffer[nn] != ',') {nn++;}
	    nn++; // designation
	    while (buffer[nn] != ',') {nn++;}
	    nn++; // source_id
	    while (buffer[nn] != ',') {nn++;}
	    nn++; // random_index
	    while (buffer[nn] != ',') {nn++;}
	    nn++; // ref_epoch
	    while (buffer[nn] != ',') {nn++;}
	    nn++; // RA
	    if (buffer[nn] == 'N') {readflag = 1;}
	    double RA = strtod(&buffer[nn], pointit);
	    while (buffer[nn] != ',') {nn++;}
	    nn++; // RA_err
	    if (buffer[nn] == 'N') {readflag = 1;}
	    double RAerr = strtod(&buffer[nn], pointit);
	    while (buffer[nn] != ',') {nn++;}
	    nn++; // DEC_obs
	    if (buffer[nn] == 'N') {readflag = 1;}
	    double DEC = strtod(&buffer[nn], pointit);
	    while (buffer[nn] != ',') {nn++;}
	    nn++; //DEC_err
	    if (buffer[nn] == 'N') {readflag = 1;}
	    double DECerr = strtod(&buffer[nn], pointit);
	    while (buffer[nn] != ',') {nn++;}
	    nn++; // par
	    if (buffer[nn] == 'N') {readflag = 1;}
	    double par = strtod(&buffer[nn], pointit);
	    while (buffer[nn] != ',') {nn++;}
	    nn++; // par_err
	    if (buffer[nn] == 'N') {readflag = 1;}
	    double parerr = strtod(&buffer[nn], pointit);
	    while (buffer[nn] != ',') {nn++;}
	    nn++; // parallax_over_error
	    double parovererr = strtod(&buffer[nn], pointit);
	    while (buffer[nn] != ',') {nn++;}
	    nn++; // pmra
	    if (buffer[nn] == 'N') {readflag = 1;}
	    double pmra = strtod(&buffer[nn], pointit);
	    while (buffer[nn] != ',') {nn++;}
	    nn++; // pmraerr
	    if (buffer[nn] == 'N') {readflag = 1;}
	    double pmraerr = strtod(&buffer[nn], pointit);
	    while (buffer[nn] != ',') {nn++;}
	    nn++; // pmdec
	    if (buffer[nn] == 'N') {readflag = 1;}
	    double pmdec = strtod(&buffer[nn], pointit);
	    while (buffer[nn] != ',') {nn++;}
	    nn++; // pmdecerr
	    if (buffer[nn] == 'N') {readflag = 1;}
	    double pmdecerr = strtod(&buffer[nn], pointit);
	    while (buffer[nn] != ',') {nn++;}
	    nn++; // ra_dec_corr
	    double radeccorr = strtod(&buffer[nn], pointit);
	    while (buffer[nn] != ',') {nn++;}
	    nn++; // ra_par_corr
	    while (buffer[nn] != ',') {nn++;}
	    nn++; // ra_pmra_corr
	    while (buffer[nn] != ',') {nn++;}
	    nn++; // ra_pmdec_corr
	    while (buffer[nn] != ',') {nn++;}
	    nn++; // dec_par_corr
	    while (buffer[nn] != ',') {nn++;}
	    nn++; // dec_pmra_corr
	    while (buffer[nn] != ',') {nn++;}
	    nn++; // dec_pmdec_corr
	    while (buffer[nn] != ',') {nn++;}
	    nn++; // par_pmra_corr
	    while (buffer[nn] != ',') {nn++;}
	    nn++; // par_pmdec_corr
	    while (buffer[nn] != ',') {nn++;}
	    nn++; // pmra_pmdec_corr
	    while (buffer[nn] != ',') {nn++;}
	    nn++; // n_obs (astrometric)
	    double nobs = strtod(&buffer[nn], pointit);
	    while (buffer[nn] != ',') {nn++;}
	    nn++; // n_obs_ac
	    while (buffer[nn] != ',') {nn++;}
	    nn++; // n_goodobs
	    double ngoodobs = strtod(&buffer[nn], pointit);
	    while (buffer[nn] != ',') {nn++;}
	    nn++; // n_badobs
	    while (buffer[nn] != ',') {nn++;}
	    nn++; // gof_al
	    while (buffer[nn] != ',') {nn++;}
	    nn++; // chi^2
	    double astchi2 = strtod(&buffer[nn], pointit);
	    double ulennart = astchi2/(ngoodobs-5.0);
	    while (buffer[nn] != ',') {nn++;}
	    nn++; // excessnoise
	    double excessnoise = strtod(&buffer[nn], pointit);
	    while (buffer[nn] != ',') {nn++;}
	    nn++; // excessnoise_sig
	    while (buffer[nn] != ',') {nn++;}
	    nn++; // params_solved
	    while (buffer[nn] != ',') {nn++;}
	    nn++; // primaryflag
	    while (buffer[nn] != ',') {nn++;}
	    nn++; // weight
	    while (buffer[nn] != ',') {nn++;}
	    nn++; // pseudocolour
	    while (buffer[nn] != ',') {nn++;}
	    nn++; // pseudocolour_err
	    while (buffer[nn] != ',') {nn++;}
	    nn++; // varpi_factor_al
	    while (buffer[nn] != ',') {nn++;}
	    nn++; //matched_obs
	    while (buffer[nn] != ',') {nn++;}
	    nn++; // visibility periods
	    double nvis = strtod(&buffer[nn], pointit);
	    while (buffer[nn] != ',') {nn++;}
	    nn++; // sigmna5d_max
	    while (buffer[nn] != ',') {nn++;}
	    nn++; // frame_rotator_obj_type
	    while (buffer[nn] != ',') {nn++;}
	    nn++; // matched obsdate
	    while (buffer[nn] != ',') {nn++;}
	    nn++; // dup_source
	    while (buffer[nn] != ',') {nn++;}
	    nn++; // phot_g_n_obs
	    double photgnobs = strtod(&buffer[nn], pointit);
	    while (buffer[nn] != ',') {nn++;}
	    nn++; // phot_g_mean_flux
	    double photg = strtod(&buffer[nn], pointit);
	    while (buffer[nn] != ',') {nn++;}
	    nn++; // phot_g_mean_flux_error
	    double photgerr = strtod(&buffer[nn], pointit);
	    while (buffer[nn] != ',') {nn++;}
	    nn++; // phot_g_mean_flux_over_error
	    while (buffer[nn] != ',') {nn++;}
	    nn++; // photgmeanmag
	    double gmag = strtod(&buffer[nn], pointit);
	    while (buffer[nn] != ',') {nn++;}
	    nn++; // phot_b_n_obs
	    double photbnobs = strtod(&buffer[nn], pointit);
	    while (buffer[nn] != ',') {nn++;}
	    nn++; // phot_b_mean_flux
	    double photb = strtod(&buffer[nn], pointit);
	    while (buffer[nn] != ',') {nn++;}
	    nn++; // phot_b_mean_flux_error
	    double photberr = strtod(&buffer[nn], pointit);
	    while (buffer[nn] != ',') {nn++;}
	    nn++; // phot_b_mean_flux_over_error
	    while (buffer[nn] != ',') {nn++;}
	    nn++; // photbmeanmag
	    double bmag = strtod(&buffer[nn], pointit);
	    while (buffer[nn] != ',') {nn++;}
	    nn++; // phot_r_n_obs
	    double photrnobs = strtod(&buffer[nn], pointit);
	    while (buffer[nn] != ',') {nn++;}
	    nn++; // phot_r_mean_flux
	    double photr = strtod(&buffer[nn], pointit);
	    while (buffer[nn] != ',') {nn++;}
	    nn++; // phot_r_mean_flux_error
	    double photrerr = strtod(&buffer[nn], pointit);
	    while (buffer[nn] != ',') {nn++;}
	    nn++; // phot_r_mean_flux_over_error
	    while (buffer[nn] != ',') {nn++;}
	    nn++; // photrmeanmag
	    double rmag = strtod(&buffer[nn], pointit);
	    while (buffer[nn] != ',') {nn++;}
	    nn++; // bp_rp_excessfactor
	    double bprpexcessfact = strtod(&buffer[nn], pointit);
	    while (buffer[nn] != ',') {nn++;}
	    nn++; // photprocmode
	    while (buffer[nn] != ',') {nn++;}
	    nn++; // bp_rp
	    while (buffer[nn] != ',') {nn++;}
	    nn++; // bp_g
	    while (buffer[nn] != ',') {nn++;}
	    nn++; // g_r
	    while (buffer[nn] != ',') {nn++;}
	    nn++; //RV
	    if (buffer[nn] == 'N') {readflag = 1;}
	    double RV = strtod(&buffer[nn], pointit);
	    while (buffer[nn] != ',') {nn++;}
	    nn++; // RVerr
	    if (buffer[nn] == 'N') {readflag = 1;}
	    double RVerr = strtod(&buffer[nn], pointit);
	    while (buffer[nn] != ',') {nn++;}
	    nn++; // RVtransits
	    while (buffer[nn] != ',') {nn++;}
	    nn++; // RVtemplateteff
	    while (buffer[nn] != ',') {nn++;}
	    nn++; // RVtemplatelogg
	    while (buffer[nn] != ',') {nn++;}
	    nn++; // RVtemplatefeh
	    while (buffer[nn] != ',') {nn++;}
	    nn++; // photvarflag
	    while (buffer[nn] != ',') {nn++;}
	    nn++; // dr2l
	    double dr2l = strtod(&buffer[nn], pointit);
	    while (buffer[nn] != ',') {nn++;}
	    nn++; // dr2b
	    double dr2b = strtod(&buffer[nn], pointit);
           while (buffer[nn] != ',') {nn++;}
            nn++; // ecliptic long
            double eclon = strtod(&buffer[nn], pointit);
           while (buffer[nn] != ',') {nn++;}
            nn++; // ec lat
            double eclat = strtod(&buffer[nn], pointit);

//	    par += PAROFFSET;
	    double paroff = PAROFFSET;
//	    double parerrb = sqrt(parerr*parerr + PARERRP*PARERRP);
	    double parerrread = parerr;
            double a2 = 0.904347;     //         +/- 0.08623      (9.535%)
            double b2 = 0.758844;  //         +/- 0.3089       (40.71%)
            double c2 = 0.116855;   //        +/- 0.08572      (73.36%)
            double d2 = -0.12;  //      +/- 0.06819      (47.59%)
	    double eclatb = eclat;
	    if (eclatb < -75.0) {eclatb = -75.0;}
            double gg = a2 + c2*cos(b2*eclatb*PI/180.0 + d2*eclatb*eclatb*(PI*PI/(180.0*180.0)));
	    double paroffdel = (gg - 1.0)*paroff*OFFCONST;        //10.0;
	    double paroffs = paroff + paroffdel;
//	    paroffs = exp(paroffs - 1.0);
	    par += paroffs;
	    double parerrb = sqrt(parerr*parerr*(PARERRF*PARERRF) + paroffs*paroffs*OFFCONE + PARERRP*PARERRP);
	    parerr = parerrb;
	    double dist = 1.0/par;
	    double rva = RV; double rvaerr = RVerr;
	    double dab = dist;
	    double *backback = new double[4];
//	    rvaerr = 7.0;
	    rva += RVOFFSET;


	    double totpmerr = sqrt(pmdecerr*pmdecerr + pmraerr*pmraerr);
	    for (int ii = 0; ii < 4; ii++)  {
		backback[ii] = 0.0;
	}
	    if (par > 10000.0) {
		cout << "WARNING parallax LARGE " << par << "\n";
		readflag = 1;
	    }
	    double bback = propmotranslation(RA, DEC, pmra, pmdec, backback, 2000);
	    double glr = backback[0]; double gbr = backback[1];
	    double pml = backback[2]; double pmb = backback[3];
//	    cout << "comparlb " << glr << " " << gbr << " " << dr2l << " " << dr2b << " " << readflag << " " << RV << " " << RVerr << " " << nn <<  "\n";
	    delete [] backback;
	   int fflag = 0;
	    if (1 == 0) {                      // (fabs(gbr*180.0/PI) < 20.0 || (RA > 0.0 && RA < 67.0 && DEC > 42.0 && DEC < 59.0)) {
		fflag = 1;
	  }
	    fflag = 0;
	if (fflag == 0 && dab < 10.0 && (readflag == 0) && (totpmerr < PMERRLIM)) {
	    if ((par/parerr > PARQUALLIMIT) && (VVCOND)) {
	    if (besetzt[ncounter] == 1) {
		persei[ncounter].join();
		besetzt[ncounter] = 0;
#ifdef REDOFULLPRIOR
	    double suminst0 = 0.0; double suminst1;
	    for (int ii = 0; ii < PFMAX; ii++) {
		suminst0 += priorfitinst[ncounter][ii][0]; suminst1 += priorfitinst[ncounter][ii][1];
	    }
	    if (suminst0 > 0.0) {suminst0 = 1.0/suminst0;}
	    if (suminst1 > 0.0) {suminst1 = 1.0/suminst1;}
            for (int ii = 0; ii < PFMAX; ii++) {
		priorfit[ii][0] += suminst0*priorfitinst[ncounter][ii][0];
		priorfit[ii][1] += suminst1*priorfitinst[ncounter][ii][1];
	    }
#endif
	}
#ifdef TESTDISTPRIOR
	    persei[ncounter] = thread(testfunc, star, par, parerr, glr, gbr, (bmag - rmag), ((double)(nvis)), ulennart, photgerr, photberr, photrerr, dist, pml, pmb, cii, pmraerr, pmdecerr, rvaerr, rva, excessnoise, facr, dfall, ncounter, priorfoot, priorfitinst[ncounter]);
#else
	    persei[ncounter] = thread(testfunc, star, par, parerr, glr, gbr, (bmag - rmag), ((double)(nvis)), ulennart, photgerr, photberr, photrerr, dist, pml, pmb, cii, pmraerr, pmdecerr, rvaerr, rva, excessnoise, facr);
#endif
	     cii++;
	     besetzt[ncounter] = 1;
	     ncounter++;
	     if (ncounter >= NTHREADS) {
		ncounter = 0;
      }
      }
	else {
		cout << "#1 " << cii << " "  << dab << " " << glr << " " << gbr << "\n";
    }

    }
	else {
		cout << "#2 " << cii << " "  << dab << " " << glr << " " << gbr << " " << RVerr << "\n";
  }


	} //end while fgets
//   }  //perhaps needed

  fclose(troja);
//  fclose(troja2);
	} //end for iif

  for (int nicounter = 0; nicounter < NTHREADS; nicounter++) {
	if (besetzt[nicounter] == 1) {
		persei[nicounter].join();
#ifdef REDOFULLPRIOR
	    double suminst0 = 0.0; double suminst1;
	    for (int ii = 0; ii < PFMAX; ii++) {
		suminst0 += priorfitinst[nicounter][ii][0]; suminst1 += priorfitinst[nicounter][ii][1];
	    }
	    if (suminst0 > 0.0) {suminst0 = 1.0/suminst0;}
	    if (suminst1 > 0.0) {suminst1 = 1.0/suminst1;}
            for (int ii = 0; ii < PFMAX; ii++) {
		priorfit[ii][0] += suminst0*priorfitinst[nicounter][ii][0];
		priorfit[ii][1] += suminst1*priorfitinst[nicounter][ii][1];
	    }
#endif
		besetzt[nicounter] = 0;
	}
   }



  cii--; //= 10;

#ifdef PRINTSTARS
  ostringstream psostrn;
  psostrn << "print" << SUFFIX;
  string psostrin = psostrn.str();
  char *psochn = &psostrin[0];
  FILE *printout = fopen(psochn, "w");
  ostringstream psostr;
  for (int iii = 0; iii < cii; iii++) {
	psostr << cii;
	for (int jjj = 0; jjj < 30; jjj++) {
		psostr << " " << star[iii][jjj];
	}
	psostr << "\n";
  }
  string psostri = psostr.str();
  char *psoch = &psostri[0];
  fputs(psoch, printout);
  fflush(printout);
  fclose(printout);
#endif
 cout << "enter quaddist stuff \n";



//  cout << "cii" << cii << "\n";
 int samb = SAMB;  // 3000;       //cii;
 int dsam = DSAM;   //1000;     //cii;
  int iisamni = 0;
  int samni = cii/dsam + 30;
  ncounter = 0;
  ostringstream *samnostr = new ostringstream[samni];
  std::ostringstream samnstr;
  int **index = new int*[NTHREADS];
  double **distr = new double*[NTHREADS];
  for (int ii = 0; ii < NTHREADS; ii++) {
	 index[ii] = new int[samb];
	 distr[ii] = new double[samb];
   }

  cout << "sortstar \n";
  sortstar(star, starsort, ORDERING, cii);
  cout << "sorted star \n";
//  for (int ii = 0; ii < cii; ii++) {
//	cout << "ii " << ii << " " << starsort[ii][ORDERING] << " " << starsort[ii][0] << " " << starsort[ii][1] << " " << starsort[ii][5] << "\n";
//  }
//  double facr = (parsec)*(PI)*0.001/(365.24*24.0*3600.0*3600.0*180.0);
  for (int ii = 1; ii < cii; ii++) {
	if (starsort[ii][ORDERING] < starsort[ii-1][ORDERING]) {
		cout << "WARNINGORDER " << ii << " " << cii << " " << starsort[ii][ORDERING] << " " << starsort[ii-1][ORDERING] << "\n";
	}
  }
  int deli = 500;
  int imin = -3*deli;
  int imax = 0;
  ostringstream linnamstr;
  linnamstr << OUTDIR << "lindist" << SUFFIX;
  string linnamstri = linnamstr.str();
  char *linnamch = &linnamstri[0];
  FILE *linout = fopen(linnamch, "w");
  ostringstream meanostr;
	cout << "begin jii loop " << cii << " " << samb << "\n";
    for (int jii = -samb + dsam; jii < cii; jii+= dsam) {

	  double means[30][3];
	for (int ii = 0; ii < 30; ii++) {
	   for (int jj = 0; jj < 3; jj++) {
		means[ii][jj] = 0.0;
	   }
	}

	cout << "jii " << jii << " " << samb << " " << dsam << "\n";
	double distmean = 0.0;
	double dn = 0.0;
	int ii = jii;
	int sam = samb;
	if (ii < 0) {
		sam = samb + ii;
		ii = 0;
	}
/*
	double Teff = starsort[ii][10];
	double feh = starsort[ii][9];
	double logg = starsort[ii][11];
	double par = starsort[ii][15];
	double parerr = starsort[ii][16];
	double parp = starsort[ii][17];
	for (int tt = 0; tt < sam; tt++) {
	    distr[tt] = 1.0e+20;
	}
      for (int jj = 0; jj < cii; jj++) {
	  double Teffb = starsort[jj][10];
	  double fehb = starsort[jj][9];
	  double loggb = starsort[jj][11];
	  double parb = starsort[jj][15];

	  double parerrb = starsort[jj][16];
	  double parpb = starsort[jj][17];
	  double dist = (parpb-parp)*(parpb-parp);               //WTEFF*(Teff-Teffb)*(Teff-Teffb)*1.0  + WFEH*(feh - fehb)*(feh-fehb)*1000.0*1000.0 + WLOGG*(logg-loggb)*(logg-loggb)*1000.0*1000.0;
	  if (dist < distr[ncounter][sam-1]) {
	      int pointi = sam-1;
	      while ((dist < distr[pointi]) && (pointi >= 0)) {
		  pointi--;
	      }
	      pointi++;
	      for (int tt = sam-1; tt > pointi; tt--) {
		index[ncounter][tt] = index[ncounter][tt-1];
		distr[ncoutner][tt] = distr[ncounter][tt-1];
	      }
	      distr[ncounter][pointi] = dist;
	      index[ncounter][pointi] = jj;
	  }
      }
*/

    if (besetzt[ncounter] == 1) {
		persei[ncounter].join();
		besetzt[ncounter] = 0;
	}
    int ttmax = ii + sam;
    if (ttmax >= cii) {
	int diff = ttmax - cii;
	for (int tt = ((cii -1) - ii); tt < sam; tt++) {
		index[ncounter][tt] = cii-1;
	}
//	cout << "reset sam " << sam << " " << cii << " " << ii + sam << " " << ttmax << "\n";
	sam = (cii -1) - ii;
     }

	for (int tt = 0; tt < sam; tt++) {
		index[ncounter][tt] = ii+tt;
/*
		if (index[ncounter][tt] >= cii) {
			sam -= (index[ncounter][tt] - cii);
			index[ncounter][tt] = cii-1;
		}
		else {
*/
			if (ADDCOND2) {                         //(starsort[ii+tt][16] < MAXTEFF) {
			for (int jj = 0; jj < 30; jj++) {
				means[jj][0] += 1.0;
				means[jj][1] += starsort[ii+tt][jj];
				means[jj][2] += starsort[ii+tt][jj]*starsort[ii+tt][jj];
			}
			}
			distmean += starsort[ii+tt][15];
			dn += 1.0;
//		}
		distr[ncounter][tt] = 1.0;
	}
	distmean = distmean/dn;

	meanostr << "means " << cii << " " << ii << " " << sam << " " << distmean << " ";
	for (int jj = 0; jj < 30; jj++) {
		meanostr << means[jj][0] << " " << means[jj][1] << " " << means[jj][2] << " ";
	}
	meanostr << "\n";

/*
	double distbias = 0.0;
	double dbn = 0.0;
	double disen = 0.0;
	for (int jii = 0; jii < cii; jii++) {
		double disi = starsort[jii][15];
		double dise = starsort[jii][16];
		double expd = exp(-(disi-distmean)*(disi-distmean)/(2.0*dise*dise))/(sqrt(dise));
		dbn += expd;
		distbias += (disi/distmean)*expd;
		disen += expd*dise;
	}
*/
//	cout << distmean << " " << dbn << " " << disen/dbn << " " <<  distbias/dbn << " ";
    if (besetzt[ncounter] == 1) {
		persei[ncounter].join();
		besetzt[ncounter] = 0;
	}

	int sampass = sam;

//	   quaddist(starsort, index, distr, sam, samnostr[iisamni]);
//	   string stri;
//	   quaddist(starsort, index, distr, sam);
//	   ostringstream samnstr;
//	   persei[ncounter] = std::thread(quaddistb, std::ref(samnstr));  //[iisamni]);
#ifdef REDOPRIOR
	    persei[ncounter] = std::thread(quaddist, starsort, index[ncounter], distr[ncounter], sampass, std::ref(samnostr[iisamni]), facr, priorfoot, priorfit);
#else
	    persei[ncounter] = std::thread(quaddist, starsort, index[ncounter], distr[ncounter], sampass, std::ref(samnostr[iisamni]), facr);
#endif
	    iisamni++;
	     besetzt[ncounter] = 1;
	     ncounter++;
	     if (ncounter >= NTHREADS) {
		ncounter = 0;
      }

//        quaddist(starsort, index, distr, sam, linout);


   }

   for (int ni = 0; ni < NTHREADS; ni++) {
	if (besetzt[ni] == 1) {
		persei[ni].join();
		besetzt[ni] = 0;
	}
   }
   cout << "try to write lindist \n";
   for (int ii = 0; ii < iisamni; ii++) {
	   cout << "try to write on file \n";
	   string ostri = samnostr[ii].str();
	   char *och = &ostri[0];
	   fputs(och, linout);
   }

   delete [] samnostr;


  cout << " finished \n";
  fclose(linout);
  ostringstream meannamstr;
  meannamstr << OUTDIR << "means" << SUFFIX;
  string meannamstri = meannamstr.str();
  char *meannamch = &meannamstri[0];
  FILE *meanout = fopen(meannamch, "w");
  string meanstri = meanostr.str();
  char *meanch = &meanstri[0];
  fputs(meanch, meanout);
  fflush(meanout);
  fclose(meanout);

 for (int ii = 0; ii < STARN; ii++) {
    delete [] star[ii];
    delete [] starsort[ii];
  }
  delete [] star;
  delete [] starsort;
   cout << "arrays deleted \n";


#ifdef TESTDISTPRIOR
    {
      double **dfall2 = new double*[5];
      for (int ii = 0; ii < 5; ii++) {
	dfall2[ii] = new double[PRIORLENGTH];
	for (int j = 0; j < PRIORLENGTH; j++) {
		dfall2[ii][j] = 0.0;
	}
      }
      for (int ii = 0; ii < 5; ii++) {
	for (int j = 0; j < PRIORLENGTH; j++) {
	   for (int pp = 0; pp < 32; pp++) {
		dfall2[ii][j] += dfall[ii][j][pp];
	   }
        }
     }


      ostringstream fpostr;
      fpostr << OUTDIR << "priortest" << SUFFIX;
      string fpostri = fpostr.str();
      char *fpoch = &fpostri[0];
      FILE *priout = fopen(fpoch, "w");
      ostringstream ostrp;
      for (int si = 0; si < PRIORLENGTH; si++) {
        ostrp << si << " " << ((double)(si) + 0.5)*0.01 << " " << dfall2[0][si] << " " << dfall2[1][si] << " " << dfall2[2][si] << " " << dfall2[3][si] << " " << dfall2[4][si] << "\n";
      }
    string strip = ostrp.str();
    char *chp = &strip[0];
    fputs(chp, priout);
    fflush(priout);
    fclose(priout);
    for (int ii = 0; ii < 5; ii++) {
	for (int jj = 0; jj < PRIORLENGTH; jj++) {
	  delete [] dfall[ii][jj];
	}
	delete [] dfall[ii];
	delete [] dfall2[ii];
    }
    delete [] dfall;
    delete [] dfall2;
    }
   for (int jj = 0; jj < NTHREADS; jj++) {
	for (int ii = 0; ii < PFMAX; ii++) {
		delete [] priorfitinst[jj][ii];
	}
	delete [] priorfitinst[jj];
   }
   delete [] priorfitinst;
   for (int ii = 0; ii < PFMAX; ii++) {
	delete [] priorfit[ii];
   }
   delete [] priorfit;
   delete [] priorfoot;
#endif

   cout << "priorarrays deleted \n";
  for (int ii = 0; ii < NTHREADS; ii++) {
	delete [] index[ii];
	delete [] distr[ii];
  }

  delete [] index;
  delete [] distr;

  return 0;
}

