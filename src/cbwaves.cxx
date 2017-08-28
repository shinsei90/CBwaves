#include <cmath>
#include <limits>
#include <fstream>
#include <sstream>
#include <map>
#include <sys/time.h>
#include <vector>
#include "include/tvector.h"
#include "include/RK4.h"
#include "include/Spline3.h"
#include "include/fft.h"
#include "include/ConfigFile.h"
#include <boost/algorithm/string.hpp>
#include <complex>
#include <stdio.h>
#include <stdarg.h>

#define H_COUNT       9
#define H_Q          (1<<0)
#define H_P05Q       (1<<1)
#define H_PQ         (1<<2)
#define H_P15Q       (1<<3)
#define H_P15Qtail   (1<<4)
#define H_PQSO       (1<<5)
#define H_P15QSO     (1<<6)
#define H_P2Q	     (1<<7)
#define H_PQSS       (1<<8)
static const char* H_NAMES[H_COUNT] = {
    "Q", "P05Q", "PQ", "P15Q", "P15Qtail", "PQSO", "P15QSO", "P2Q", "PQSS"
    // order is important!
};

#define C_COUNT       12
#define C_PN         (1<<0)
#define C_2PN        (1<<1)
#define C_SO         (1<<2)
#define C_SS         (1<<3)
#define C_RR         (1<<4)
#define C_PNSO       (1<<5)
#define C_3PN        (1<<6)
#define C_1RR        (1<<7)
#define C_2PNSO      (1<<8)
#define C_RRSO       (1<<9)
#define C_RRSS       (1<<10)
#define C_KADM       (1<<11)
static const char* C_NAMES[C_COUNT] = {
    "PN", "2PN", "SO", "SS", "RR", "PNSO", "3PN", "1RR", "2PNSO", "RRSO", "RRSS","KADM"
    // order is important!
};

#ifndef REAL
#define REAL      double
#endif /* REAL */

#define SI_G      6.67428E-11
#define SI_c      2.99792458E8
#define SI_c2     (SI_c*SI_c)
#define SI_ly     9460730472580800.0
#define SI_pc     (3.26156*SI_ly)
#define PI        3.141592653589793
//#define SI_M_SUN  1.9891e30

REAL    erx,ery,erz,er,ermin,ermax,erorbit,ecc_r;      // quantities for calculating the eccentricity;
REAL    orbits_turn;
REAL    rm2, rm1, eccm1, ecc;
string  chkpoint("no");                                // as default there is no checkpointing
string  description;                                   // default description of the run
string  chkpointfile("CBwaves.chk");                   // the chkpoint file name
string  histfile("");                                  // the previous chkpointfile
int     chkpointinterval = 600;                        // default setting for checkpointing interval in seconds
REAL    runfrom = 0;
REAL    t;

int     stepprint  = 1;                                 // output variables after stepprint step 
int     orbitprint = 0;                                 // output variables after orbitprint orbit

#define  cbwNumberOfLogLevels 7
enum     cbwLogLevels  {cbwALL = 6, cbwDEBUG = 5, cbwINFO = 4, cbwWARN =3, cbwERR = 2, cbwCRIT = 1, cbwNONE = 0 };
int      cbwLogLevel = cbwNONE;
string   adaptive        = "no";
int      adaptive_step   = 1000;
string   cbwTrace;
char *   cbwLogLevelNames[cbwNumberOfLogLevels];


template <class T> inline T sqr(T x) { return x*x; 
}

enum {
    I_rNx = 0,
    I_rNy = 1,
    I_rNz = 2,
    I_vNx = 3,
    I_vNy = 4,
    I_vNz = 5,
    I_rPNx = 6,
    I_rPNy = 7,
    I_rPNz = 8,
    I_vPNx = 9,
    I_vPNy = 10,
    I_vPNz = 11,
    I_E_PNrad = 12, // radiated energy post-Newtonian contribution
    I_r2PNx = 13,
    I_r2PNy = 14,
    I_r2PNz = 15,
    I_v2PNx = 16,
    I_v2PNy = 17,
    I_v2PNz = 18,
    I_E_2PNrad = 19, // radiated energy 2PN correction
    I_rSOx = 20,
    I_rSOy = 21,
    I_rSOz = 22,
    I_vSOx = 23,
    I_vSOy = 24,
    I_vSOz = 25,
    I_E_SOrad = 26, // radiated energy spin-orbit correction
    I_rSSx = 27,
    I_rSSy = 28,
    I_rSSz = 29,
    I_vSSx = 30,
    I_vSSy = 31,
    I_vSSz = 32,
    I_E_SSrad = 33, // radiated energy spin-spin correction
    I_rRRx = 34,
    I_rRRy = 35,
    I_rRRz = 36,
    I_vRRx = 37,
    I_vRRy = 38,
    I_vRRz = 39,
    I_E_Nrad = 40, // radiated energy Newtonian correction
    I_Jx_Nrad = 41, // radiated angular momentum Newtonian correction
    I_Jy_Nrad = 42,
    I_Jz_Nrad = 43,
    I_s1x = 44,
    I_s1y = 45,
    I_s1z = 46,
    I_s2x = 47,
    I_s2y = 48,
    I_s2z = 49,
    I_orbits = 50,
    I_r3PNx = 51,
    I_r3PNy = 52,
    I_r3PNz = 53,
    I_v3PNx = 54,
    I_v3PNy = 55,
    I_v3PNz = 56,
    I_r1RRx = 57,
    I_r1RRy = 58,
    I_r1RRz = 59,
    I_v1RRx = 60,
    I_v1RRy = 61,
    I_v1RRz = 62,
    I_rRRSOx = 63,
    I_rRRSOy = 64,
    I_rRRSOz = 65,
    I_vRRSOx = 66,
    I_vRRSOy = 67,
    I_vRRSOz = 68,
    I_rRRSSx = 69,
    I_rRRSSy = 70,
    I_rRRSSz = 71,
    I_vRRSSx = 72,
    I_vRRSSy = 73,
    I_vRRSSz = 74,
    I_rPNSOx = 75,
    I_rPNSOy = 76,
    I_rPNSOz = 77,
    I_vPNSOx = 78,
    I_vPNSOy = 79,
    I_vPNSOz = 80,
    I_r2PNSOx = 81,
    I_r2PNSOy = 82,
    I_r2PNSOz = 83,
    I_v2PNSOx = 84,
    I_v2PNSOy = 85,
    I_v2PNSOz = 86,
    I_E_PNSOrad = 87, // radiated energy PNSO correction
    I_E_25PNrad = 88, // radiated energy 2.5PN correction
    I_rHx = 89, // ADM Hamiltonian coordnates
    I_rHy = 90,
    I_rHz = 91,
    I_pHx = 92,
    I_pHy = 93,
    I_pHz = 94,
    N_COMPONENTS = 95 
};

enum {
    aux_dt                 = 0,
    aux_orbits             = 1,
    aux_orbitcount         = 2,
    aux_stepcount          = 3,
    aux_prevorbittime      = 4,
    aux_prevorbittime2     = 5,
    aux_prevorbitnumber    = 6,
    aux_prevchkpointtime   = 7,
    aux_orbitturncount     = 8,
    aux_N_COMPONENTS       = 9
};

static const string AUX_NAMES[] = {
    "aux_dt", "aux_orbits", "aux_orbitcount", "aux_stepcount", 
    "aux_prevorbittime", "aux_prevorbittime2", "aux_prevorbitnumber","aux_prevchkpointtime"
};

using namespace std;

class CBwaveODE: public virtual ODE<REAL>
{
public:
    static const string COMPONENT_NAMES[];

public:
    int corrections;
    int hterms;
    REAL m1;
    REAL m2;
    REAL m;
    REAL mu;
    REAL eta;
    REAL rmin;
    REAL rmax;
    REAL alpha;
    REAL beta;
    REAL delta1;
    REAL delta2;
    REAL delta3;
    REAL delta4;
    REAL delta5;
    REAL delta6;

public:
    CBwaveODE(REAL m1, REAL m2):
	    ODE<REAL>(N_COMPONENTS, COMPONENT_NAMES),
	    m1(m1), m2(m2) {
	corrections = C_PN | C_2PN | C_SO | C_SS | C_RR | C_PNSO | C_3PN | C_1RR | C_2PNSO | C_RRSO | C_RRSS | C_KADM;
	hterms = H_Q | H_P05Q | H_PQ | H_PQSO
			 | H_P15Q | H_P15Qtail | H_P15QSO | H_P2Q | H_PQSS;
	init_m();
    }

    bool set(const string& name, REAL value, tvalarray<REAL>& vars);

    void eval(const REAL* f, int offset, REAL x, REAL* df);

private:
    void init_m() {
	m = m1 + m2;
	mu = m1*m2/m;
	eta = mu/m;
    }
};

struct ObserverParameters {
    REAL iota, phi;
    REAL D;
    REAL theta, varphi, psi;

    ObserverParameters() {
	iota = phi = 0;
	D = 1;
	theta = varphi = psi = 0;
    }
};

// Both distance and mass are measured in meters: [G/c^2*kg] = meters

CBwaveODE           ode(SI_G/SI_c2*2e30, SI_G/SI_c2*2e30); // m1,m2
ObserverParameters  obs;
tvalarray  <REAL>   f(ode.getNumComponents());
REAL                aux[aux_N_COMPONENTS];
tvector    <string> outcolumns;
tvector    <string> auxoutcolumns;
string              datafname, ftfname;

typedef struct {
 REAL x;
 REAL y;
 REAL z;
} myvector;

typedef struct {
 REAL alpha;
 REAL beta;
} rotangles;

void printvec(const myvector invec) {
 cout << " x: " << invec.x
 << endl << " y: " << invec.y << endl << " z: "<<invec.z << endl;
}

myvector rotate(const rotangles inangles, const myvector invector) {
 // alpha rotates around the z axis
 // beta  rotates around the y axis

 REAL alpha=inangles.alpha;
 REAL beta =inangles.beta;
 myvector ret,tmpvec;

 // rotation around the z axis
 tmpvec.x = cos(alpha)*invector.x - sin(alpha)*invector.y;
 tmpvec.y = sin(alpha)*invector.x + cos(alpha)*invector.y;
 tmpvec.z = invector.z;

 // rotation around the y axis
 ret.x = cos(beta)*tmpvec.x + sin(beta)*tmpvec.z;
 ret.y = tmpvec.y;
 ret.z = - sin(beta)*tmpvec.x + cos(beta)*tmpvec.z;

 return ret;
}

rotangles findangle (const myvector invec) {

 // this routine finds the angle (alpha, beta)

 REAL length;
 rotangles angle;
  angle.alpha=0.;
  angle.beta =0.;

 length = sqrt(invec.x*invec.x+invec.y*invec.y);
 if ( length > 0.) {
  angle.alpha = acos(invec.x/length);
  if (invec.y < 0.) {
  angle.alpha = 2*PI-angle.alpha;
  }
 }
 angle.alpha = - angle.alpha;

 length = sqrt(invec.x*invec.x+invec.y*invec.y+invec.z*invec.z);
 if (length > 0.) {
  angle.beta = - acos(invec.z/length);
 }
 return angle;
}


void cbwInitLogMessages() {
    for (int i = 0 ; i < cbwNumberOfLogLevels; i ++ ) {
        cbwLogLevelNames[i] = (char *) malloc(10);
    }

    strcpy(cbwLogLevelNames[0],"NONE");
    strcpy(cbwLogLevelNames[1],"CRIT");
    strcpy(cbwLogLevelNames[2],"ERROR");
    strcpy(cbwLogLevelNames[3],"WARN");
    strcpy(cbwLogLevelNames[4],"INFO");
    strcpy(cbwLogLevelNames[5],"DEBUG");
    strcpy(cbwLogLevelNames[6],"ALL");
}

void cbwLog(const int loglevel, const char* format, ... ) {
    time_t cbwtime  = time(NULL);
    struct tm *now = localtime(&cbwtime);

    va_list args;
    if (loglevel <= cbwLogLevel) {
      fprintf(stderr, "%d/%d %d:%d:%d: ",now->tm_mon+1,now->tm_mday,
                        now->tm_hour, now->tm_min, now->tm_sec);
      fprintf(stderr, "%s: (%s) ", cbwLogLevelNames[loglevel], cbwTrace.c_str());
      va_start( args, format );
        vfprintf(stderr, format, args );
      va_end( args );
      fprintf( stderr, "\n" );
   }
}

typedef REAL (*function_t)(REAL, const REAL*, const CBwaveODE&,
			   const ObserverParameters&);
#define FUNCDEF(name)  REAL eval_##name(REAL t, const REAL* f, const CBwaveODE& ode, const ObserverParameters& obs)

static FUNCDEF(t) { return t/SI_c; }
static FUNCDEF(rNx) { return f[I_rNx]; }
static FUNCDEF(rNy) { return f[I_rNy]; }
static FUNCDEF(rNz) { return f[I_rNz]; }
static FUNCDEF(vNx) { return f[I_vNx]; }
static FUNCDEF(vNy) { return f[I_vNy]; }
static FUNCDEF(vNz) { return f[I_vNz]; }
static FUNCDEF(rPNx) { return f[I_rPNx]; }
static FUNCDEF(rPNy) { return f[I_rPNy]; }
static FUNCDEF(rPNz) { return f[I_rPNz]; }
static FUNCDEF(vPNx) { return f[I_vPNx]; }
static FUNCDEF(vPNy) { return f[I_vPNy]; }
static FUNCDEF(vPNz) { return f[I_vPNz]; }
static FUNCDEF(E_PNrad) { return f[I_E_PNrad]; } // energy loss PN correction
static FUNCDEF(r2PNx) { return f[I_r2PNx]; }
static FUNCDEF(r2PNy) { return f[I_r2PNy]; }
static FUNCDEF(r2PNz) { return f[I_r2PNz]; }
static FUNCDEF(v2PNx) { return f[I_v2PNx]; }
static FUNCDEF(v2PNy) { return f[I_v2PNy]; }
static FUNCDEF(v2PNz) { return f[I_v2PNz]; }
static FUNCDEF(E_2PNrad) { return f[I_E_2PNrad]; } // energy loss 2PN correction
static FUNCDEF(r3PNx) { return f[I_r3PNx]; }
static FUNCDEF(r3PNy) { return f[I_r3PNy]; }
static FUNCDEF(r3PNz) { return f[I_r3PNz]; }
static FUNCDEF(v3PNx) { return f[I_v3PNx]; }
static FUNCDEF(v3PNy) { return f[I_v3PNy]; }
static FUNCDEF(v3PNz) { return f[I_v3PNz]; }
static FUNCDEF(rSOx) { return f[I_rSOx]; }
static FUNCDEF(rSOy) { return f[I_rSOy]; }
static FUNCDEF(rSOz) { return f[I_rSOz]; }
static FUNCDEF(vSOx) { return f[I_vSOx]; }
static FUNCDEF(vSOy) { return f[I_vSOy]; }
static FUNCDEF(vSOz) { return f[I_vSOz]; }
static FUNCDEF(E_SOrad) { return f[I_E_SOrad]; } // energy loss spin-orbit corr.
static FUNCDEF(rSSx) { return f[I_rSSx]; }
static FUNCDEF(rSSy) { return f[I_rSSy]; }
static FUNCDEF(rSSz) { return f[I_rSSz]; }
static FUNCDEF(vSSx) { return f[I_vSSx]; }
static FUNCDEF(vSSy) { return f[I_vSSy]; }
static FUNCDEF(vSSz) { return f[I_vSSz]; }
static FUNCDEF(E_SSrad) { return f[I_E_SSrad]; } // energy loss spin-spin corr.
static FUNCDEF(rRRx) { return f[I_rRRx]; }
static FUNCDEF(rRRy) { return f[I_rRRy]; }
static FUNCDEF(rRRz) { return f[I_rRRz]; }
static FUNCDEF(vRRx) { return f[I_vRRx]; }
static FUNCDEF(vRRy) { return f[I_vRRy]; }
static FUNCDEF(vRRz) { return f[I_vRRz]; }
static FUNCDEF(rPNSOx) { return f[I_rPNSOx]; }
static FUNCDEF(rPNSOy) { return f[I_rPNSOy]; }
static FUNCDEF(rPNSOz) { return f[I_rPNSOz]; }
static FUNCDEF(vPNSOx) { return f[I_vPNSOx]; }
static FUNCDEF(vPNSOy) { return f[I_vPNSOy]; }
static FUNCDEF(vPNSOz) { return f[I_vPNSOz]; }
static FUNCDEF(r2PNSOx) { return f[I_r2PNSOx]; }
static FUNCDEF(r2PNSOy) { return f[I_r2PNSOy]; }
static FUNCDEF(r2PNSOz) { return f[I_r2PNSOz]; }
static FUNCDEF(v2PNSOx) { return f[I_v2PNSOx]; }
static FUNCDEF(v2PNSOy) { return f[I_v2PNSOy]; }
static FUNCDEF(v2PNSOz) { return f[I_v2PNSOz]; }
static FUNCDEF(E_Nrad) { return f[I_E_Nrad]; } // energy loss Newtonian corr.
static FUNCDEF(Jx_Nrad) { return f[I_Jx_Nrad]; } // ang.mom loss Newtonian corr.
static FUNCDEF(Jy_Nrad) { return f[I_Jy_Nrad]; } 
static FUNCDEF(Jz_Nrad) { return f[I_Jz_Nrad]; }
static FUNCDEF(E_25PNrad) { return f[I_E_25PNrad]; } // energy loss 2.5PN corr.
static FUNCDEF(E_PNSOrad) { return f[I_E_PNSOrad]; } // energy loss PNSO corr.
static FUNCDEF(r1RRx) { return f[I_r1RRx]; }
static FUNCDEF(r1RRy) { return f[I_r1RRy]; }
static FUNCDEF(r1RRz) { return f[I_r1RRz]; }
static FUNCDEF(v1RRx) { return f[I_v1RRx]; }
static FUNCDEF(v1RRy) { return f[I_v1RRy]; }
static FUNCDEF(v1RRz) { return f[I_v1RRz]; }
static FUNCDEF(rRRSOx) { return f[I_rRRSOx]; }
static FUNCDEF(rRRSOy) { return f[I_rRRSOy]; }
static FUNCDEF(rRRSOz) { return f[I_rRRSOz]; }
static FUNCDEF(vRRSOx) { return f[I_vRRSOx]; }
static FUNCDEF(vRRSOy) { return f[I_vRRSOy]; }
static FUNCDEF(vRRSOz) { return f[I_vRRSOz]; }
static FUNCDEF(rRRSSx) { return f[I_rRRSSx]; }
static FUNCDEF(rRRSSy) { return f[I_rRRSSy]; }
static FUNCDEF(rRRSSz) { return f[I_rRRSSz]; }
static FUNCDEF(vRRSSx) { return f[I_vRRSSx]; }
static FUNCDEF(vRRSSy) { return f[I_vRRSSy]; }
static FUNCDEF(vRRSSz) { return f[I_vRRSSz]; }
static FUNCDEF(s1x) { return f[I_s1x]; }
static FUNCDEF(s1y) { return f[I_s1y]; }
static FUNCDEF(s1z) { return f[I_s1z]; }
static FUNCDEF(s2x) { return f[I_s2x]; }
static FUNCDEF(s2y) { return f[I_s2y]; }
static FUNCDEF(s2z) { return f[I_s2z]; }
static FUNCDEF(orbits) { return f[I_orbits]; }
static FUNCDEF(rx) {
    return f[I_rNx] + f[I_rPNx] + f[I_rSOx] + f[I_r2PNx] + f[I_rSSx] + f[I_rRRx] + f[I_r3PNx] + f[I_r1RRx] 
		+ f[I_rPNSOx] + f[I_r2PNSOx] + f[I_rRRSOx] + f[I_rRRSSx];
}
static FUNCDEF(ry) {
    return f[I_rNy] + f[I_rPNy] + f[I_rSOy] + f[I_r2PNy] + f[I_rSSy] + f[I_rRRy] + f[I_r3PNy] + f[I_r1RRy] 
		+ f[I_rPNSOy] + f[I_r2PNSOy] + f[I_rRRSOy] + f[I_rRRSSy];
}
static FUNCDEF(rz) {
    return f[I_rNz] + f[I_rPNz] + f[I_rSOz] + f[I_r2PNz] + f[I_rSSz] + f[I_rRRz] + f[I_r3PNz] + f[I_r1RRz] 
		+ f[I_rPNSOz] + f[I_r2PNSOz] + f[I_rRRSOz] + f[I_rRRSSz];
}
static FUNCDEF(vx) {
    return f[I_vNx] + f[I_vPNx] + f[I_vSOx] + f[I_v2PNx] + f[I_vSSx] + f[I_vRRx] + f[I_v3PNx] + f[I_v1RRx] 
		+ f[I_vPNSOx] + f[I_v2PNSOx] + f[I_vRRSOx] + f[I_vRRSSx];
}
static FUNCDEF(vy) {
    return f[I_vNy] + f[I_vPNy] + f[I_vSOy] + f[I_v2PNy] + f[I_vSSy] + f[I_vRRy] + f[I_v3PNy] + f[I_v1RRy] 
		+ f[I_vPNSOy] + f[I_v2PNSOy] + f[I_vRRSOy] + f[I_vRRSSy];
}
static FUNCDEF(vz) {
    return f[I_vNz] + f[I_vPNz] + f[I_vSOz] + f[I_v2PNz] + f[I_vSSz] + f[I_vRRz] + f[I_v3PNz] + f[I_v1RRz] 
		+ f[I_vPNSOz] + f[I_v2PNSOz] + f[I_vRRSOz] + f[I_vRRSSz];
}

static FUNCDEF(drAdiabatic) {
	REAL m1   = ode.m1;
	REAL m2   = ode.m2;
	REAL m    = ode.m1 + ode.m2;
	REAL eta  = ode.eta;

	REAL chi1 = sqrt(f[I_s1x]*f[I_s1x] + f[I_s1y]*f[I_s1y] + f[I_s1z]*f[I_s1z]);
	REAL chi2 = sqrt(f[I_s2x]*f[I_s2x] + f[I_s2y]*f[I_s2y] + f[I_s2z]*f[I_s2z]);

	REAL cSs1 = ode.m1 * ode.m1;
	REAL s1x  = cSs1 * f[I_s1x];
	REAL s1y  = cSs1 * f[I_s1y];
	REAL s1z  = cSs1 * f[I_s1z];

	REAL cSs2 = ode.m2 * ode.m2;
	REAL s2x  = cSs2 * f[I_s2x];
	REAL s2y  = cSs2 * f[I_s2y];
	REAL s2z  = cSs2 * f[I_s2z];

	REAL rx = eval_rx(t, f, ode, obs);
	REAL ry = eval_ry(t, f, ode, obs);
	REAL rz = eval_rz(t, f, ode, obs);
	REAL rsq = rx*rx + ry*ry + rz*rz;
	REAL r = sqrt(rsq);

	REAL vx = eval_vx(t, f, ode, obs);
        REAL vy = eval_vy(t, f, ode, obs);
        REAL vz = eval_vz(t, f, ode, obs);
	
	REAL LNx = ode.mu*(ry*vz - rz*vy);
        REAL LNy = ode.mu*(rz*vx - rx*vz);
        REAL LNz = ode.mu*(rx*vy - ry*vx);


	REAL LNs1 = (LNx * s1x + LNy * s1y + LNz * s1z );
	REAL LNs2 = (LNx * s2x + LNy * s2y + LNz * s2z );
	REAL s1s2 = (s1x * s2x + s1y * s2y + s1z * s2z );
	REAL sum1 = chi1 * LNs1 * (19 * m1 * m1 / (m * m) + 15 * eta);
	REAL sum2 = chi2 * LNs2 * (19 * m2 * m2 / (m * m) + 15 * eta);

	REAL dr = -(64./5.) * eta * pow(m/r,3.) * ( 1 - 1./336. * (1751 + 588 * eta) *
                  (m/r) - (7./12. * (sum1 + sum2)  - 4 * PI) * pow(m/r,3./2.) -
                  5./48. * eta * chi1 * chi2 * (59 * s1s2 - 173 * LNs1 * LNs2) *
                  pow(m/r,2));

	return dr;
}

static FUNCDEF(r) {
    REAL rx = eval_rx(t, f, ode, obs);
    REAL ry = eval_ry(t, f, ode, obs);
    REAL rz = eval_rz(t, f, ode, obs);
    return sqrt(rx*rx + ry*ry + rz*rz);
}

static FUNCDEF(mr) {
    REAL rx = eval_rx(t, f, ode, obs);
    REAL ry = eval_ry(t, f, ode, obs);
    REAL rz = eval_rz(t, f, ode, obs);
    REAL r = sqrt(rx*rx + ry*ry + rz*rz);
    return ode.m/r;
}

static FUNCDEF(v) {
    REAL vx = eval_vx(t, f, ode, obs);
    REAL vy = eval_vy(t, f, ode, obs);
    REAL vz = eval_vz(t, f, ode, obs);
    return sqrt(vx*vx + vy*vy + vz*vz);
}

static FUNCDEF(v2) {
    REAL vx = eval_vx(t, f, ode, obs);
    REAL vy = eval_vy(t, f, ode, obs);
    REAL vz = eval_vz(t, f, ode, obs);
    return vx*vx + vy*vy + vz*vz;
}

static FUNCDEF(rdot) {
    REAL rx = eval_rx(t, f, ode, obs);
    REAL ry = eval_ry(t, f, ode, obs);
    REAL rz = eval_rz(t, f, ode, obs);
    REAL r = sqrt(rx*rx + ry*ry + rz*rz);
    REAL vx = eval_vx(t, f, ode, obs);
    REAL vy = eval_vy(t, f, ode, obs);
    REAL vz = eval_vz(t, f, ode, obs);
    return (rx*vx + ry*vy + rz*vz)/r;
}

static FUNCDEF(orbfreq) {	// orbital freq: \omega^2 = (v^2 - rdot^2)/r^2
    REAL r = eval_r(t, f, ode, obs);
    REAL v2 = eval_v2(t, f, ode, obs);
    REAL rdot = eval_rdot(t, f, ode, obs);
    return sqrt(v2-rdot*rdot)/r;
}

static FUNCDEF(x1) {
    return ode.m2/(ode.m1 + ode.m2)*eval_rx(t, f, ode, obs);
}
static FUNCDEF(y1) {
    return ode.m2/(ode.m1 + ode.m2)*eval_ry(t, f, ode, obs);
}
static FUNCDEF(z1) {
    return ode.m2/(ode.m1 + ode.m2)*eval_rz(t, f, ode, obs);
}
static FUNCDEF(x2) {
    return -ode.m1/(ode.m1 + ode.m2)*eval_rx(t, f, ode, obs);
}
static FUNCDEF(y2) {
    return -ode.m1/(ode.m1 + ode.m2)*eval_ry(t, f, ode, obs);
}
static FUNCDEF(z2) {
    return -ode.m1/(ode.m1 + ode.m2)*eval_rz(t, f, ode, obs);
}

static void eval_hij_TT(REAL t, const REAL* f, const CBwaveODE& ode,
			const ObserverParameters& obs, REAL* h_TT)
{
    REAL h[9];
    REAL N[3] = {sin(obs.iota)*cos(obs.phi),
		 sin(obs.iota)*sin(obs.phi),
		 cos(obs.iota)};
    REAL m = ode.m;

    // relative position vector
    REAL rx = eval_rx(t, f, ode, obs);
    REAL ry = eval_ry(t, f, ode, obs);
    REAL rz = eval_rz(t, f, ode, obs);
    REAL rsq = rx*rx + ry*ry + rz*rz;
    REAL r = sqrt(rsq);

    // relative velocity vector
    REAL v[3] = {eval_vx(t, f, ode, obs),
		 eval_vy(t, f, ode, obs),
		 eval_vz(t, f, ode, obs)};
    REAL vsq = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];

    // direction vector
    REAL n[3] = {rx/r, ry/r, rz/r};

    REAL Nn = N[0]*n[0] + N[1]*n[1] + N[2]*n[2];
    REAL Nv = N[0]*v[0] + N[1]*v[1] + N[2]*v[2];
    REAL rdot = (rx*v[0] + ry*v[1] + rz*v[2])/r;
    REAL rdot2 = rdot*rdot;
    REAL dm = ode.m1 - ode.m2;
    REAL eta = ode.eta;
    REAL eta2 = eta*eta;
    REAL Gm_r = ode.m/r;

    // spins
    REAL cSs1 = ode.m1*ode.m1;
    REAL S1x = cSs1*f[I_s1x];
    REAL S1y = cSs1*f[I_s1y];
    REAL S1z = cSs1*f[I_s1z];
    REAL cSs2 = ode.m2*ode.m2;
    REAL S2x = cSs2*f[I_s2x];
    REAL S2y = cSs2*f[I_s2y];
    REAL S2z = cSs2*f[I_s2z];
    REAL Sx = S1x + S2x;
    REAL Sy = S1y + S2y;
    REAL Sz = S1z + S2z;

    REAL S1[3], S2[3];
    S1[0] = S1x;
    S1[1] = S1y;
    S1[2] = S1z;
    S2[0] = S2x;
    S2[1] = S2y;
    S2[2] = S2z;

    REAL S12 = S1x*S2x+S1y*S2y+S1z*S2z;
    REAL nS1 = n[0]*S1x + n[1]*S1y + n[2]*S1z;
    REAL nS2 = n[0]*S2x + n[1]*S2y + n[2]*S2z;

    REAL Deltax = m*(S2x/ode.m2 - S1x/ode.m1);
    REAL Deltay = m*(S2y/ode.m2 - S1y/ode.m1);
    REAL Deltaz = m*(S2z/ode.m2 - S1z/ode.m1);
    REAL S[3], Delta[3];
    S[0] = Sx;
    S[1] = Sy;
    S[2] = Sz;
    Delta[0] = Deltax;
    Delta[1] = Deltay;
    Delta[2] = Deltaz;
    REAL Delta_cross_N[3];
    Delta_cross_N[0] = Deltay*N[2] - Deltaz*N[1];
    Delta_cross_N[1] = Deltaz*N[0] - Deltax*N[2];
    Delta_cross_N[2] = Deltax*N[1] - Deltay*N[0];
    REAL n_cross_v[3];
    n_cross_v[0] = n[1]*v[2] - n[2]*v[1];
    n_cross_v[1] = n[2]*v[0] - n[0]*v[2];
    n_cross_v[2] = n[0]*v[1] - n[1]*v[0];
    REAL P15QSO_A = 0;
    for(int i = 0; i < 3; ++i) {
	P15QSO_A += n_cross_v[i]*(12*S[i] + 6*dm/m*Delta[i]);
    }
    REAL v_cross_B[3];
    v_cross_B[0] = v[1]*(9*Sz + 5*dm/m*Deltaz) - v[2]*(9*Sy + 5*dm/m*Deltay);
    v_cross_B[1] = v[2]*(9*Sx + 5*dm/m*Deltax) - v[0]*(9*Sz + 5*dm/m*Deltaz);
    v_cross_B[2] = v[0]*(9*Sy + 5*dm/m*Deltay) - v[1]*(9*Sx + 5*dm/m*Deltax);
    REAL C_cross_N[3];
    C_cross_N[0] = (Sy + dm/m*Deltay)*N[2] - (Sz + dm/m*Deltaz)*N[1];
    C_cross_N[1] = (Sz + dm/m*Deltaz)*N[0] - (Sx + dm/m*Deltax)*N[2];
    C_cross_N[2] = (Sx + dm/m*Deltax)*N[1] - (Sy + dm/m*Deltay)*N[0];
    REAL n_cross_D[3];
    n_cross_D[0] = n[1]*(2*Sz + 2*dm/m*Deltaz) - n[2]*(2*Sy + 2*dm/m*Deltay);
    n_cross_D[1] = n[2]*(2*Sx + 2*dm/m*Deltax) - n[0]*(2*Sz + 2*dm/m*Deltaz);
    n_cross_D[2] = n[0]*(2*Sy + 2*dm/m*Deltay) - n[1]*(2*Sx + 2*dm/m*Deltax);
    REAL n_cross_E[3];
    n_cross_E[0] = n[1]*(12*Sz + 6*dm/m*Deltaz) - n[2]*(12*Sy + 6*dm/m*Deltay);
    n_cross_E[1] = n[2]*(12*Sx + 6*dm/m*Deltax) - n[0]*(12*Sz + 6*dm/m*Deltaz);
    n_cross_E[2] = n[0]*(12*Sy + 6*dm/m*Deltay) - n[1]*(12*Sx + 6*dm/m*Deltax);
    for(int i = 0; i < 3; ++i) {
	for(int j = 0; j < 3; ++j) {
	    REAL sum = 0;
	    if((ode.hterms & H_Q) != 0) {
		REAL Q = 2*(v[i]*v[j] - Gm_r*n[i]*n[j]);
		sum += Q;
	    }
	    if((ode.hterms & H_P05Q) != 0) {
		REAL P05Q = dm/m*(3*Gm_r*(n[i]*v[j] + n[j]*v[i]
					  - rdot*n[i]*n[j])*Nn
			+ (Gm_r*n[i]*n[j] - 2*v[i]*v[j])*Nv);
		sum += P05Q;
	    }
	    if((ode.hterms & H_PQ) != 0) {
		REAL PQ = (1-3*eta)/3.0
		    *(4*m/r*(3*rdot*n[i]*n[j] - 4*(n[i]*v[j] + n[j]*v[i]))*Nn*Nv
			    + 2*(3*v[i]*v[j] - m/r*n[i]*n[j])*Nv*Nv
			    + m/r*((3*vsq - 15*rdot*rdot + 7*m/r)*n[i]*n[j]
				+ 15*rdot*(n[i]*v[j] + n[j]*v[i])
				- 14*v[i]*v[j])*Nn*Nn)
		    + 2/3.0*m/r*rdot*(5 + 3*eta)*(n[i]*v[j] + n[j]*v[i])
		    + ((1-3*eta)*vsq - 2/3.0*(2 - 3*eta)*m/r)*v[i]*v[j]
		    + m/r*((1-3*eta)*rdot*rdot - 1/3.0*(10 + 3*eta)*vsq
			    + 29/3.0*m/r)*n[i]*n[j];
		sum += PQ;
	    }
	    if((ode.hterms & H_PQSO) != 0) {
		REAL PQSO = 1.0/rsq*(Delta_cross_N[i]*n[j]
			+ Delta_cross_N[j]*n[i]);
		sum += PQSO;
	    }
	    if((ode.hterms & H_P15Q) != 0) {
		REAL P15Q = dm/m*(1-2*eta)
		    *(0.25*m/r*((45*rdot*rdot - 9*vsq - 28*m/r)*n[i]*n[j]
				+ 58*v[i]*v[j]
				- 54*rdot*(n[i]*v[j] + n[j]*v[i]))*Nn*Nn*Nv
			    + 0.5*(m/r*n[i]*n[j] - 4*v[i]*v[j])*Nv*Nv*Nv
			    + m/r*(5/4.0*(3*vsq - 7*rdot*rdot + 6*m/r)
					*rdot*n[i]*n[j]
				- 1/12.0*(21*vsq - 105*rdot*rdot + 44*m/r)
				*(n[i]*v[j] + n[j]*v[i])
				- 17/2.0*rdot*v[i]*v[j])*Nn*Nn*Nn
			    + 1.5*m/r*(5*(n[i]*v[j] + n[j]*v[i])
				       - 3*rdot*n[i]*n[j])
			    *Nn*Nv*Nv)
		    + dm/(12*r)*Nn*(n[i]*n[j]*rdot*(rdot*rdot*(15 - 90*eta)
				- vsq*(63 - 54*eta)
				+ m/r*(242 - 24*eta))
			    - rdot*v[i]*v[j]*(186 + 24*eta)
			    + (n[i]*v[j] + n[j]*v[i])
			    *(rdot*rdot*(63 + 54*eta)
				- m/r*(128 - 36*eta) + vsq*(33 - 18*eta)))
		    + dm/m*Nv*(0.5*v[i]*v[j]*(m/r*(3 - 8*eta)
					      - 2*vsq*(1 - 5*eta))
			    - 0.5*(n[i]*v[j] + n[j]*v[i])*m/r*rdot*(7 + 4*eta)
			    - n[i]*n[j]*m/r*(3/4.0*(1 - 2*eta)*rdot*rdot
				+ 1/3.0*(26 - 3*eta)*m/r
				- 1/4.0*(7 - 2*eta)*vsq));
		sum += P15Q;
	    }
	    if((ode.hterms & H_P15Qtail) != 0) {
		REAL P15Qtail = 0; // TODO (integral es kvadrupol tenzor van benne)
		sum += P15Qtail;
	    }
	    if((ode.hterms & H_P15QSO) != 0) {
		REAL P15QSO = 1.0/rsq
		    *(2*n[i]*n[j]*P15QSO_A
			    - n[i]*v_cross_B[j] - n[j]*v_cross_B[i]
			    + (3*rdot*Nn - 2*Nv)
			    *(C_cross_N[i]*n[j] + C_cross_N[j]*n[i])
			    - (v[i]*n_cross_D[j] + v[j]*n_cross_D[i])
			    + rdot*(n[i]*n_cross_E[j] + n[j]*n_cross_E[i])
			    - 2*Nn*(C_cross_N[i]*v[j] + C_cross_N[j]*v[i]));
		sum += P15QSO;
	    }
            if((ode.hterms & H_P2Q) != 0) {			// Will, Wiseman - PRD54(96)4813, Eq. (6.11f)
                REAL P2Q = 1.0/60.0*(1 - 5*eta + 5*eta2)
                    *(24*Nv*Nv*Nv*Nv*(5*v[i]*v[j] - m/r*n[i]*n[j])
                      + m/r*Nn*Nn*Nn*Nn*( 2*(175.0*m/r - 465*rdot2 + 93*vsq)*v[i]*v[j] 
                                + 15*rdot*( -50*m/r + 63*rdot2 - 27*vsq)*(n[i]*v[j] + n[j]*v[i]) 
				+ (1155.0*m/r*rdot2 - 172.0*m*m/r/r - 945*rdot2*rdot2 -159.0*m/r*vsq + 630*rdot2*vsq - 45*vsq*vsq)*n[i]*n[j] )
                      + 24.0*m/r*Nn*Nn*Nn*Nv*( 87*rdot*v[i]*v[j]
				+ 5*rdot*(14*rdot2 - 15.0*m/r - 6*vsq)*n[i]*n[j]
                                - 8*(5.0*m/r - 10*rdot2 + 2*vsq)*(n[i]*v[j] + n[j]*v[i]) )
                      + 288.0*m/r*Nn*Nv*Nv*Nv*( rdot*n[i]*n[j] - 2*(n[i]*v[j] + n[j]*v[i]) )
                      + 24.0*m/r*Nn*Nn*Nv*Nv*( (35.0*m/r -45*rdot2 + 9*vsq)*n[i]*n[j] - 76*v[i]*v[j] + 63*rdot*(n[i]*v[j] + n[j]*v[i]) ) )
		+ 1.0/15.0*Nv*Nv*( (5*(25-78*eta+12*eta2)*m/r - (18-65*eta+45*eta2)*vsq + 9*(1-5*eta+5*eta2)*rdot2)*m/r*n[i]*n[j]
				+ 3*(5*(1-9*eta+21*eta2)*vsq - 2*(4-25*eta+45*eta2)*m/r)*v[i]*v[j]
				+ 9*(6-15*eta-10*eta2)*m/r*rdot*(n[i]*v[j] + n[j]*v[i]) )
                + 1.0/15.0*Nn*Nv*m/r*( (3*(36-145*eta+150*eta2)*vsq - 5*(127-392*eta+36*eta2)*m/r - 15*(2-15*eta+30*eta2)*rdot2)*rdot*n[i]*n[j] 
                                + 6*(98-295*eta-30*eta2)*rdot*v[i]*v[j]
                                + (5*(66-221*eta+96*eta2)*m/r - 9*(18-45*eta-40*eta2)*rdot2 - (66-265*eta+360*eta2)*vsq)*(n[i]*v[j] + n[j]*v[i]) )
                + 1.0/60.0*Nn*Nn*m/r*( (3*(33-130*eta+150*eta2)*vsq*vsq - 105*(1-10*eta+30*eta2)*rdot2*rdot2 
						+ 15*(181-572*eta+84*eta2)*m/r*rdot2 - (131-770*eta+930*eta2)*m/r*vsq
						- 60*(9-40*eta+60*eta2)*vsq*rdot2 - 8*(131-390*eta+30*eta2)*m*m/r/r)*n[i]*n[j]
                                + 4*((12+5*eta-315*eta2)*vsq - 9*(39-115*eta+35*eta2)*rdot2 + 5*(29-104*eta+84*eta2)*m/r)*v[i]*v[j]
                     + 2*(15*(18-40*eta-75*eta2)*rdot2 - 5*(197-640*eta+180*eta2)*m/r + 3*(21-130*eta+375*eta2)*vsq)*rdot*(n[i]*v[j] + n[j]*v[i]) )
                + 1.0/60.0*( ((467+780*eta-120*eta2)*m/r*vsq - 15*(61-96*eta+48*eta2)*m/r*rdot2
                                                - (144-265*eta-135*eta2)*vsq*vsq + 6*(24-95*eta+75*eta2)*vsq*rdot2
                                                - 2*(642+545*eta)*m*m/r/r - 45*(1-5*eta+5*eta2)*rdot2*rdot2)*m/r*n[i]*n[j]
                                + (4*(69+10*eta-135*eta2)*m/r*vsq - 12*(3+60*eta+25*eta2)*m/r*rdot2 
						+ 45*(1-7*eta+13*eta2)*vsq*vsq - 10*(56+165*eta-12*eta2)*m*m/r/r)*v[i]*v[j]
		     + 2*(2*(36+5*eta-75*eta2)*vsq - 6*(7-15*eta-15*eta2)*rdot2 + 5*(35+45*eta+36*eta2)*m/r)*m/r*rdot*(n[i]*v[j] + n[j]*v[i])  );
                sum += P2Q;
            }
            if((ode.hterms & H_PQSS) != 0) {			// Kidder PRD52(95)821, Eq. (3.22h)
                REAL PQSS = -6.0/ode.mu/rsq/r*( n[i]*n[j]*(S12 - 5*nS1*nS2)
                        + nS2*(n[i]*S1[j]+n[j]*S1[i]) + nS1*(n[i]*S2[j]+n[j]*S2[i]) );
                sum += PQSS;
            }
	    h[3*i + j] = 2*ode.mu/obs.D*sum;
	}
    }

    REAL P[9];
    for(int i = 0; i < 3; ++i) {
	for(int j = 0; j < 3; ++j) {
	    P[3*i + j] = -N[i]*N[j];
	}
	P[3*i + i] += 1;
    }
    for(int i = 0; i < 3; ++i) {
	for(int j = 0; j < 3; ++j) {
	    h_TT[3*i + j] = 0;
	    for(int k = 0; k < 3; ++k) {
		for(int l = 0; l < 3; ++l) {
		    REAL Lambda_ijkl = P[3*i + k]*P[3*j + l]
					- 0.5*P[3*i + j]*P[3*k + l];
		    h_TT[3*i+j] += Lambda_ijkl*h[3*k+l];
		}
	    }
	}
    }
}


#include "include/hmodes.h"

static FUNCDEF(hkp) {
    REAL hij_TT[9];
    eval_hij_TT(t, f, ode, obs, hij_TT);
    REAL hkp = 0.5*(hij_TT[3*2+2]-hij_TT[3*1+1]);
    return hkp;
}

static FUNCDEF(hp) {
    REAL hij_TT[9];
    REAL p[3] = {cos(obs.iota)*cos(obs.phi), cos(obs.iota)*sin(obs.phi),
		 -sin(obs.iota)};
    REAL q[3] = {-sin(obs.phi), cos(obs.phi), 0};
    REAL hp = 0;
    eval_hij_TT(t, f, ode, obs, hij_TT);
    for(int i = 0; i < 3; ++i) {
	for(int j = 0; j < 3; ++j) {
	    hp += 0.5*(p[i]*p[j] - q[i]*q[j])*hij_TT[3*i + j];
	}
    }
    return hp;
}

static FUNCDEF(hx) {
    REAL hij_TT[9];
    REAL p[3] = {cos(obs.iota)*cos(obs.phi), cos(obs.iota)*sin(obs.phi),
		 -sin(obs.iota)};
    REAL q[3] = {-sin(obs.phi), cos(obs.phi), 0};
    REAL hx = 0;
    eval_hij_TT(t, f, ode, obs, hij_TT);
    for(int i = 0; i < 3; ++i) {
	for(int j = 0; j < 3; ++j) {
	    hx += 0.5*(p[i]*q[j] + q[i]*p[j])*hij_TT[3*i + j];
	}
    }
    return hx;
}

static FUNCDEF(h) {
    REAL hij_TT[9];
    REAL p[3] = {cos(obs.iota)*cos(obs.phi), cos(obs.iota)*sin(obs.phi),
		 -sin(obs.iota)};
    REAL q[3] = {-sin(obs.phi), cos(obs.phi), 0};
    REAL hp = 0, hx = 0;
    eval_hij_TT(t, f, ode, obs, hij_TT);
    for(int i = 0; i < 3; ++i) {
	for(int j = 0; j < 3; ++j) {
	    hp += 0.5*(p[i]*p[j] - q[i]*q[j])*hij_TT[3*i + j];
	    hx += 0.5*(p[i]*q[j] + q[i]*p[j])*hij_TT[3*i + j];
	}
    }
    REAL ct = cos(obs.theta);
    REAL c2phi = cos(2*obs.varphi), s2phi = sin(2*obs.varphi);
    REAL c2psi = cos(2*obs.psi),    s2psi = sin(2*obs.psi);
    REAL Fp = -0.5*(1 + ct*ct)*c2phi*c2psi - ct*s2phi*s2psi;
    REAL Fx = 0.5*(1 + ct*ct)*c2phi*s2psi - ct*s2phi*c2psi;
    return Fp*hp + Fx*hx;
}

static FUNCDEF(E_N) {
    REAL r = eval_r(t, f, ode, obs);
    REAL v = eval_v(t, f, ode, obs);
    return ode.mu*(0.5*v*v - ode.m/r);
}

static FUNCDEF(E_PN) {
    if(ode.corrections & C_PN) {
	REAL m = ode.m;
	REAL eta = ode.eta;
	REAL rx = eval_rx(t, f, ode, obs);
	REAL ry = eval_ry(t, f, ode, obs);
	REAL rz = eval_rz(t, f, ode, obs);
	REAL rsq = rx*rx + ry*ry + rz*rz;
	REAL r = sqrt(rsq);
	REAL vx = eval_vx(t, f, ode, obs);
	REAL vy = eval_vy(t, f, ode, obs);
	REAL vz = eval_vz(t, f, ode, obs);
	REAL vsq = vx*vx + vy*vy + vz*vz;
	REAL rdot = (rx*vx + ry*vy + rz*vz)/r;
	return ode.mu*(3./8.0*(1 - 3*eta)*vsq*vsq + 0.5*(3 + eta)*vsq*m/r
		+ 0.5*eta*m*rdot*rdot/r + 0.5*m*m/rsq);
    } else {
	return 0;
    }
}


static FUNCDEF(E_PNtot) {
    return eval_E_PN(t, f, ode, obs) + f[I_E_PNrad];
}

static FUNCDEF(E_2PN) {
    if(ode.corrections & C_2PN) {
	REAL eta = ode.eta;
	REAL eta2 = eta*eta;
	REAL mu = ode.mu;
	REAL rx = eval_rx(t, f, ode, obs);
	REAL ry = eval_ry(t, f, ode, obs);
	REAL rz = eval_rz(t, f, ode, obs);
	REAL rsq = rx*rx + ry*ry + rz*rz;
	REAL r = sqrt(rsq);
	REAL vx = eval_vx(t, f, ode, obs);
	REAL vy = eval_vy(t, f, ode, obs);
	REAL vz = eval_vz(t, f, ode, obs);
	REAL v2 = vx*vx + vy*vy + vz*vz;
	REAL v4 = v2*v2;
	REAL rdot = (rx*vx + ry*vy + rz*vz)/r;
	REAL rdot2 = rdot*rdot;
	REAL m_r = ode.m/r;
	return mu*(5./16.*(1 - 7*eta + 13*eta2)*v2*v4
		- 3./8.*eta*(1 - 3*eta)*m_r*rdot2*rdot2
		+ 1./8.*(21 - 23*eta - 27*eta2)*m_r*v4
		+ 1./8.*(14 - 55*eta + 4*eta2)*m_r*m_r*v2
		+ 1./4.*eta*(1 - 15*eta)*m_r*v2*rdot2
		- 1./4.*(2 + 15*eta)*m_r*m_r*m_r
		+ 1./8.*(4 + 69*eta + 12*eta2)*m_r*m_r*rdot2);
    } else {
	return 0;
    }
}

static FUNCDEF(E_2PNtot) {
    return eval_E_2PN(t, f, ode, obs) + f[I_E_2PNrad];
}

static FUNCDEF(E_SO) {
    if(ode.corrections & C_SO) {
	REAL rx = eval_rx(t, f, ode, obs);
	REAL ry = eval_ry(t, f, ode, obs);
	REAL rz = eval_rz(t, f, ode, obs);
	REAL rsq = rx*rx + ry*ry + rz*rz;
	REAL r = sqrt(rsq);
	REAL vx = eval_vx(t, f, ode, obs);
	REAL vy = eval_vy(t, f, ode, obs);
	REAL vz = eval_vz(t, f, ode, obs);
	REAL cSs1 = ode.m1*ode.m1;
	REAL S1x = cSs1*f[I_s1x];
	REAL S1y = cSs1*f[I_s1y];
	REAL S1z = cSs1*f[I_s1z];
	REAL cSs2 = ode.m2*ode.m2;
	REAL S2x = cSs2*f[I_s2x];
	REAL S2y = cSs2*f[I_s2y];
	REAL S2z = cSs2*f[I_s2z];
	REAL S[3], LN[3], Delta[3];
	S[0] = S1x + S2x;
	S[1] = S1y + S2y;
	S[2] = S1z + S2z;
	LN[0] = ode.mu*(ry*vz - rz*vy);
	LN[1] = ode.mu*(rz*vx - rx*vz);
	LN[2] = ode.mu*(rx*vy - ry*vx);
	Delta[0] = ode.m*(S2x/ode.m2 - S1x/ode.m1);
	Delta[1] = ode.m*(S2y/ode.m2 - S1y/ode.m1);
	Delta[2] = ode.m*(S2z/ode.m2 - S1z/ode.m1);
	REAL dm = ode.m1 - ode.m2;
	REAL sum = 0;
	for(int i = 0; i < 3; ++i) {
	    sum += LN[i]*(S[i] + dm/ode.m*Delta[i]);
	}
	return sum/(r*rsq);
    } else {
	return 0;
    }
}

static FUNCDEF(E_SOtot) {
    return eval_E_SO(t, f, ode, obs) + f[I_E_SOrad];
}

static FUNCDEF(E_SS) {
    if(ode.corrections & C_SS) {
	REAL rx = eval_rx(t, f, ode, obs);
	REAL ry = eval_ry(t, f, ode, obs);
	REAL rz = eval_rz(t, f, ode, obs);
	REAL rsq = rx*rx + ry*ry + rz*rz;
	REAL r = sqrt(rsq);
	REAL nx = rx/r;
	REAL ny = ry/r;
	REAL nz = rz/r;
	REAL cSs1 = ode.m1*ode.m1;
	REAL S1x = cSs1*f[I_s1x];
	REAL S1y = cSs1*f[I_s1y];
	REAL S1z = cSs1*f[I_s1z];
	REAL cSs2 = ode.m2*ode.m2;
	REAL S2x = cSs2*f[I_s2x];
	REAL S2y = cSs2*f[I_s2y];
	REAL S2z = cSs2*f[I_s2z];
	return (3.*(nx*S1x + ny*S1y + nz*S1z)*(nx*S2x + ny*S2y + nz*S2z)
		- (S1x*S2x + S1y*S2y + S1z*S2z))/(r*rsq);
    } else {
	return 0;
    }
}

static FUNCDEF(E_SStot) {
    return eval_E_SS(t, f, ode, obs) + f[I_E_SSrad];
}

static FUNCDEF(E_3PN) {			// Mora, Will PRD69(04)104021, Eq. (2.11d)
    if(ode.corrections & C_3PN) {
        REAL eta = ode.eta;
        REAL eta2 = eta*eta;
        REAL mu = ode.mu;
        REAL rx = eval_rx(t, f, ode, obs);
        REAL ry = eval_ry(t, f, ode, obs);
        REAL rz = eval_rz(t, f, ode, obs);
        REAL rsq = rx*rx + ry*ry + rz*rz;
        REAL r = sqrt(rsq);
        REAL vx = eval_vx(t, f, ode, obs);
        REAL vy = eval_vy(t, f, ode, obs);
        REAL vz = eval_vz(t, f, ode, obs);
        REAL v2 = vx*vx + vy*vy + vz*vz;
        REAL v4 = v2*v2;
        REAL rdot = (rx*vx + ry*vy + rz*vz)/r;
        REAL rdot2 = rdot*rdot;
        REAL m_r = ode.m/r;
        REAL m_r2 = m_r*m_r;
        return mu*( (3.0/8.0 + 18469.0/840.0*eta)*m_r2*m_r2
                + (5.0/4.0 - (6747.0/280.0 - 41.0/64.0*PI*PI)*eta - 21.0/4.0*eta2 + 0.5*eta*eta2)*m_r2*m_r*v2
                + (3.0/2.0 + (2321.0/280.0 - 123.0/64.0*PI*PI)*eta + 51.0/4.0*eta2 + 7.0/2.0*eta*eta2)*m_r2*m_r*rdot2
                + 1.0/128.0*(35 - 413*eta + 1666*eta2 - 2261*eta*eta2)*v4*v4
                + 1.0/16.0*(135 - 194*eta + 406*eta2 - 108*eta*eta2)*m_r2*v4
                + 1.0/16.0*(12 + 248*eta - 815*eta2 - 324*eta*eta2)*m_r2*v2*rdot2
		- 1.0/48.0*eta*(731 - 492*eta - 288*eta2)*m_r2*rdot2*rdot2
		+ 1.0/16.0*(55 - 215*eta + 116*eta2 + 325*eta*eta2)*m_r*v4*v2
		+ 1.0/16.0*eta*(5 - 25*eta + 25*eta2)*m_r*rdot2*rdot2*rdot2
		- 1.0/16.0*eta*(21 + 75*eta - 375*eta2)*m_r*v4*rdot2
		- 1.0/16.0*eta*(9 - 84*eta + 165*eta2)*m_r*v2*rdot2*rdot2 );
    } else {
        return 0;
    }
}

static FUNCDEF(E_PNSO) {			// Faye,Blanchet,Buonanno PRD74(06)104033, Eq. (7.4b)
						// Boh� etal CQG30(13)075017, Eq. (3.9c)
    if(ode.corrections & C_PNSO) {
	REAL eta = ode.eta;
	REAL rx = eval_rx(t, f, ode, obs);
	REAL ry = eval_ry(t, f, ode, obs);
	REAL rz = eval_rz(t, f, ode, obs);
	REAL rsq = rx*rx + ry*ry + rz*rz;
	REAL r = sqrt(rsq);
	REAL vx = eval_vx(t, f, ode, obs);
	REAL vy = eval_vy(t, f, ode, obs);
	REAL vz = eval_vz(t, f, ode, obs);
        REAL v2 = vx*vx + vy*vy + vz*vz;
        REAL rdot = (rx*vx + ry*vy + rz*vz)/r;
	REAL cSs1 = ode.m1*ode.m1;
	REAL S1x = cSs1*f[I_s1x];
	REAL S1y = cSs1*f[I_s1y];
	REAL S1z = cSs1*f[I_s1z];
	REAL cSs2 = ode.m2*ode.m2;
	REAL S2x = cSs2*f[I_s2x];
	REAL S2y = cSs2*f[I_s2y];
	REAL S2z = cSs2*f[I_s2z];
	REAL S[3], LN[3], Delta[3];
	S[0] = S1x + S2x;
	S[1] = S1y + S2y;
	S[2] = S1z + S2z;
	LN[0] = ode.mu*(ry*vz - rz*vy);
	LN[1] = ode.mu*(rz*vx - rx*vz);
	LN[2] = ode.mu*(rx*vy - ry*vx);
	Delta[0] = ode.m*(S2x/ode.m2 - S1x/ode.m1);
	Delta[1] = ode.m*(S2y/ode.m2 - S1y/ode.m1);
	Delta[2] = ode.m*(S2z/ode.m2 - S1z/ode.m1);
	REAL dm = ode.m1 - ode.m2;
	REAL sum = 0;
	for(int i = 0; i < 3; ++i) {
	    sum += LN[i]*( S[i]*(-1.5*(1+eta)*v2 - 1.5*eta*rdot*rdot + 2*eta*ode.m/r) 
			+ dm/ode.m*Delta[i]*((0.5-2.5*eta)*v2 + 1.5*eta*ode.m/r) );
	}
	return sum/(r*rsq);
    } else {
	return 0;
    }
}
static FUNCDEF(E_2PNSO) {			// Boh� etal CQG30(13)075017, Eq. (3.9d)
    if(ode.corrections & C_2PNSO) {
	REAL eta = ode.eta;
	REAL eta2 = eta*eta;
	REAL rx = eval_rx(t, f, ode, obs);
	REAL ry = eval_ry(t, f, ode, obs);
	REAL rz = eval_rz(t, f, ode, obs);
	REAL rsq = rx*rx + ry*ry + rz*rz;
	REAL r = sqrt(rsq);
	REAL vx = eval_vx(t, f, ode, obs);
	REAL vy = eval_vy(t, f, ode, obs);
	REAL vz = eval_vz(t, f, ode, obs);
        REAL v2 = vx*vx + vy*vy + vz*vz;
        REAL rdot = (rx*vx + ry*vy + rz*vz)/r;
        REAL rdot2 = rdot*rdot;
	REAL cSs1 = ode.m1*ode.m1;
	REAL S1x = cSs1*f[I_s1x];
	REAL S1y = cSs1*f[I_s1y];
	REAL S1z = cSs1*f[I_s1z];
	REAL cSs2 = ode.m2*ode.m2;
	REAL S2x = cSs2*f[I_s2x];
	REAL S2y = cSs2*f[I_s2y];
	REAL S2z = cSs2*f[I_s2z];
	REAL S[3], LN[3], Delta[3];
	S[0] = S1x + S2x;
	S[1] = S1y + S2y;
	S[2] = S1z + S2z;
	LN[0] = ode.mu*(ry*vz - rz*vy);
	LN[1] = ode.mu*(rz*vx - rx*vz);
	LN[2] = ode.mu*(rx*vy - ry*vx);
	Delta[0] = ode.m*(S2x/ode.m2 - S1x/ode.m1);
	Delta[1] = ode.m*(S2y/ode.m2 - S1y/ode.m1);
	Delta[2] = ode.m*(S2z/ode.m2 - S1z/ode.m1);
	REAL dm = ode.m1 - ode.m2;
	REAL sum = 0;
	for(int i = 0; i < 3; ++i) {
	    sum += LN[i]*( S[i]*((15./8.*eta-45./8.*eta2)*rdot2*rdot2 + (-3./4.*eta+51./4.*eta2)*rdot2*v2 
			 + (-21./8.+31./8.*eta+55./8.*eta2)*v2*v2 
			 + ode.m/r*((-187./4.*eta-5.*eta2)*rdot2 + (-6.+143./4.*eta-10.*eta2)*v2) 
			 + ode.m*ode.m/r/r*(0.5+79./4.*eta+2.*eta2) ) 
			+ dm/ode.m*Delta[i]*(-15./8.*eta2*rdot2*rdot2 + (3.*eta+3./2.*eta2)*rdot2*v2 
			 + (3./8.-53./8.*eta+99./8.*eta2)*v2*v2 
			 + ode.m/r*((-179./8.*eta-3./2.*eta2)*rdot2 + (2.+145./8.*eta-41./4.*eta2)*v2) 
			 + ode.m*ode.m/r/r*(0.5+17./2.*eta+3./2.*eta2) ) );
	}
	return sum/(r*rsq);
    } else {
	return 0;
    }
}

/** Energy of the compact binary system (non-conserved). */
static FUNCDEF(E) {
    REAL E = eval_E_N(t, f, ode, obs);
    if(ode.corrections & C_PN) {
	E += eval_E_PN(t, f, ode, obs);
    }
    if(ode.corrections & C_SO) {
	E += eval_E_SO(t, f, ode, obs);
    }
    if(ode.corrections & C_2PN) {
	E += eval_E_2PN(t, f, ode, obs);
    }
    if(ode.corrections & C_SS) {
	E += eval_E_SS(t, f, ode, obs);
    }
    if(ode.corrections & C_3PN) {
        E += eval_E_3PN(t, f, ode, obs);
    }
    if(ode.corrections & C_PNSO) {
        E += eval_E_PNSO(t, f, ode, obs);
    }
    if(ode.corrections & C_2PNSO) {
        E += eval_E_2PNSO(t, f, ode, obs);
    }
    return E;
}

/** Radiated energy loss, sum of contributions. */
static FUNCDEF(E_rad) {
    return f[I_E_Nrad] + f[I_E_PNrad] + f[I_E_SOrad] + f[I_E_SSrad] + f[I_E_2PNrad] + f[I_E_PNSOrad] + f[I_E_25PNrad];
}

/** Total, conserved energy of the binary, with radiated energy subtracted. */
static FUNCDEF(E_tot) {
    REAL E = eval_E(t, f, ode, obs);
    REAL E_rad = eval_E_rad(t, f, ode, obs);
    return E + E_rad;
}

static FUNCDEF(Lnx) {
    REAL ry = eval_ry(t, f, ode, obs);
    REAL rz = eval_rz(t, f, ode, obs);

    REAL vy = eval_vy(t, f, ode, obs);
    REAL vz = eval_vz(t, f, ode, obs);

    return ode.mu*(ry*vz - rz*vy);
}

static FUNCDEF(Lny) {
    REAL rx = eval_rx(t, f, ode, obs);
    REAL rz = eval_rz(t, f, ode, obs);

    REAL vx = eval_vx(t, f, ode, obs);
    REAL vz = eval_vz(t, f, ode, obs);

    return ode.mu*(rz*vx - rx*vz);
}

static FUNCDEF(Lnz) {
    REAL rx = eval_rx(t, f, ode, obs);
    REAL ry = eval_ry(t, f, ode, obs);

    REAL vx = eval_vx(t, f, ode, obs);
    REAL vy = eval_vy(t, f, ode, obs);

    return ode.mu*(rx*vy - ry*vx);
}

static FUNCDEF(Ln2) {
    REAL Lnx = eval_Lnx(t, f, ode, obs);
    REAL Lny = eval_Lny(t, f, ode, obs);
    REAL Lnz = eval_Lnz(t, f, ode, obs);

    return Lnx*Lnx + Lny*Lny + Lnz*Lnz;
}

/* Orbital angular momentum, Kidder (2.9)*/
static FUNCDEF(Lx) {
    REAL m = ode.m;
    REAL mu = ode.mu;
    REAL eta = ode.eta;
    REAL eta2 = eta*eta;
    REAL m1 = ode.m1;
    REAL m2 = ode.m2;

    REAL rx = eval_rx(t, f, ode, obs);
    REAL ry = eval_ry(t, f, ode, obs);
    REAL rz = eval_rz(t, f, ode, obs);
    REAL rsq = rx*rx + ry*ry + rz*rz;
    REAL r = sqrt(rsq);
    REAL Gm_r = m/r;
    REAL nx = rx/r;
    REAL ny = ry/r;
    REAL nz = rz/r;

    REAL vx = eval_vx(t, f, ode, obs);
    REAL vy = eval_vy(t, f, ode, obs);
    REAL vz = eval_vz(t, f, ode, obs);
    REAL vsq = vx*vx + vy*vy + vz*vz;
    REAL rdot = (rx*vx + ry*vy + rz*vz)/r;

    // spins
    REAL cSs1 = ode.m1*ode.m1;
    REAL S1x = cSs1*f[I_s1x];
    REAL S1y = cSs1*f[I_s1y];
    REAL S1z = cSs1*f[I_s1z];
    REAL cSs2 = ode.m2*ode.m2;
    REAL S2x = cSs2*f[I_s2x];
    REAL S2y = cSs2*f[I_s2y];
    REAL S2z = cSs2*f[I_s2z];

    REAL S1n = S1x*nx + S1y*ny + S1z*nz;
    REAL S1v = S1x*vx + S1y*vy + S1z*vz;
    REAL S2n = S2x*nx + S2y*ny + S2z*nz;
    REAL S2v = S2x*vx + S2y*vy + S2z*vz;

    REAL Lnx = mu*(ry*vz - rz*vy);
    REAL Lx = Lnx;
    if(ode.corrections & C_PN) {
        Lx += Lnx*(0.5*vsq*(1 - 3*eta) 
                  + (3 + eta)*Gm_r);
    }
    if(ode.corrections & C_SO) {
        Lx += mu/m*(Gm_r*(nx*(S1n*(2 + m2/m1) + S2n*(2 + m1/m2))
                          - (S1x*(2 + m2/m1) + S2x*(2 + m1/m2)) )
                    - 0.5*(vx*(S1v*m2/m1 + S2v*m1/m2)
                           - (S1x*m2/m1 + S2x*m1/m2)*vsq) );
    }
    if(ode.corrections & C_2PN) {
        Lx += Lnx*(3.0/8.0*(1 - 7*eta + 13*eta*eta)*vsq*vsq
                   - 0.5*eta*(2 + 5*eta)*Gm_r*rdot*rdot
                   + 0.5*(7 - 10*eta - 9*eta*eta)*Gm_r*vsq
                   + 0.25*(14 - 41*eta + 4*eta*eta)*Gm_r*Gm_r);
    }
    if(ode.corrections & C_3PN) {			// Mora, Will PRD69(04)104021, Eq. (2.8d), (2.12d)
        Lx += Lnx*( (5.0/2.0 - (5199.0/280.0 - 41.0/32.0*PI*PI)*eta - 7*eta2 + eta*eta2)*Gm_r*Gm_r*Gm_r
	     + 1.0/16.0*(5 - 59*eta + 238*eta2 - 323*eta*eta2)*vsq*vsq*vsq
             + 1.0/12.0*(135 - 322*eta + 315*eta2 - 108*eta*eta2)*Gm_r*Gm_r*vsq
             + 1.0/24.0*(12 - 287*eta - 951*eta2 - 324*eta*eta2)*Gm_r*Gm_r*rdot*rdot
             + 1.0/8.0*(33 - 142*eta + 106*eta2 + 195*eta*eta2)*Gm_r*vsq*vsq
             - 0.25*eta*(12 - 7*eta - 75*eta2)*Gm_r*vsq*rdot*rdot
             + 3.0/8.0*eta*(2 - 2*eta - 11*eta2)*Gm_r*rdot*rdot*rdot*rdot );
    }
    return Lx;
}

static FUNCDEF(Ly) {
    REAL m = ode.m;
    REAL mu = ode.mu;
    REAL eta = ode.eta;
    REAL eta2 = eta*eta;
    REAL m1 = ode.m1;
    REAL m2 = ode.m2;

    REAL rx = eval_rx(t, f, ode, obs);
    REAL ry = eval_ry(t, f, ode, obs);
    REAL rz = eval_rz(t, f, ode, obs);
    REAL rsq = rx*rx + ry*ry + rz*rz;
    REAL r = sqrt(rsq);
    REAL Gm_r = m/r;
    REAL nx = rx/r;
    REAL ny = ry/r;
    REAL nz = rz/r;

    REAL vx = eval_vx(t, f, ode, obs);
    REAL vy = eval_vy(t, f, ode, obs);
    REAL vz = eval_vz(t, f, ode, obs);
    REAL vsq = vx*vx + vy*vy + vz*vz;
    REAL rdot = (rx*vx + ry*vy + rz*vz)/r;

    // spins
    REAL cSs1 = ode.m1*ode.m1;
    REAL S1x = cSs1*f[I_s1x];
    REAL S1y = cSs1*f[I_s1y];
    REAL S1z = cSs1*f[I_s1z];
    REAL cSs2 = ode.m2*ode.m2;
    REAL S2x = cSs2*f[I_s2x];
    REAL S2y = cSs2*f[I_s2y];
    REAL S2z = cSs2*f[I_s2z];

    REAL S1n = S1x*nx + S1y*ny + S1z*nz;
    REAL S1v = S1x*vx + S1y*vy + S1z*vz;
    REAL S2n = S2x*nx + S2y*ny + S2z*nz;
    REAL S2v = S2x*vx + S2y*vy + S2z*vz;

    REAL Lny = mu*(rz*vx - rx*vz);
    REAL Ly = Lny;
    if(ode.corrections & C_PN) {
        Ly += Lny*(0.5*vsq*(1 - 3*eta)
                  + (3 + eta)*Gm_r);
    }
    if(ode.corrections & C_SO) {
        Ly += mu/m*(Gm_r*(ny*(S1n*(2 + m2/m1) + S2n*(2 + m1/m2))
                          - (S1y*(2 + m2/m1) + S2y*(2 + m1/m2)) )
                    - 0.5*(vy*(S1v*m2/m1 + S2v*m1/m2)
                           - (S1y*m2/m1 + S2y*m1/m2)*vsq) );
    }
    if(ode.corrections & C_2PN) {
        Ly += Lny*(3.0/8.0*(1 - 7*eta + 13*eta*eta)*vsq*vsq
                   - 0.5*eta*(2 + 5*eta)*Gm_r*rdot*rdot
                   + 0.5*(7 - 10*eta - 9*eta*eta)*Gm_r*vsq
                   + 0.25*(14 - 41*eta + 4*eta*eta)*Gm_r*Gm_r);
    }
    if(ode.corrections & C_3PN) {                       // Mora, Will PRD69(04)104021, Eq. (2.8d), (2.12d)
        Ly += Lny*( (5.0/2.0 - (5199.0/280.0 - 41.0/32.0*PI*PI)*eta - 7*eta2 + eta*eta2)*Gm_r*Gm_r*Gm_r
             + 1.0/16.0*(5 - 59*eta + 238*eta2 - 323*eta*eta2)*vsq*vsq*vsq
             + 1.0/12.0*(135 - 322*eta + 315*eta2 - 108*eta*eta2)*Gm_r*Gm_r*vsq
             + 1.0/24.0*(12 - 287*eta - 951*eta2 - 324*eta*eta2)*Gm_r*Gm_r*rdot*rdot
             + 1.0/8.0*(33 - 142*eta + 106*eta2 + 195*eta*eta2)*Gm_r*vsq*vsq
             - 0.25*eta*(12 - 7*eta - 75*eta2)*Gm_r*vsq*rdot*rdot
             + 3.0/8.0*eta*(2 - 2*eta - 11*eta2)*Gm_r*rdot*rdot*rdot*rdot );
    }
    return Ly;
}

static FUNCDEF(Lz) {
    REAL m = ode.m;
    REAL mu = ode.mu;
    REAL eta = ode.eta;
    REAL eta2 = eta*eta;
    REAL m1 = ode.m1;
    REAL m2 = ode.m2;

    REAL rx = eval_rx(t, f, ode, obs);
    REAL ry = eval_ry(t, f, ode, obs);
    REAL rz = eval_rz(t, f, ode, obs);
    REAL rsq = rx*rx + ry*ry + rz*rz;
    REAL r = sqrt(rsq);
    REAL Gm_r = m/r;
    REAL nx = rx/r;
    REAL ny = ry/r;
    REAL nz = rz/r;

    REAL vx = eval_vx(t, f, ode, obs);
    REAL vy = eval_vy(t, f, ode, obs);
    REAL vz = eval_vz(t, f, ode, obs);
    REAL vsq = vx*vx + vy*vy + vz*vz;
    REAL rdot = (rx*vx + ry*vy + rz*vz)/r;

    // spins
    REAL cSs1 = ode.m1*ode.m1;
    REAL S1x = cSs1*f[I_s1x];
    REAL S1y = cSs1*f[I_s1y];
    REAL S1z = cSs1*f[I_s1z];
    REAL cSs2 = ode.m2*ode.m2;
    REAL S2x = cSs2*f[I_s2x];
    REAL S2y = cSs2*f[I_s2y];
    REAL S2z = cSs2*f[I_s2z];

    REAL S1n = S1x*nx + S1y*ny + S1z*nz;
    REAL S1v = S1x*vx + S1y*vy + S1z*vz;
    REAL S2n = S2x*nx + S2y*ny + S2z*nz;
    REAL S2v = S2x*vx + S2y*vy + S2z*vz;

    REAL Lnz =  mu*(rx*vy - ry*vx);
    REAL Lz = Lnz;
    if(ode.corrections & C_PN) {
        Lz += Lnz*(0.5*vsq*(1 - 3*eta)
                  + (3 + eta)*Gm_r);
    }
    if(ode.corrections & C_SO) {
        Lz += mu/m*(Gm_r*(nz*(S1n*(2 + m2/m1) + S2n*(2 + m1/m2))
                          - (S1z*(2 + m2/m1) + S2z*(2 + m1/m2)) )
                    - 0.5*(vz*(S1v*m2/m1 + S2v*m1/m2)
                           - (S1z*m2/m1 + S2z*m1/m2)*vsq) );
    }
    if(ode.corrections & C_2PN) {
        Lz += Lnz*(3.0/8.0*(1 - 7*eta + 13*eta*eta)*vsq*vsq
                   - 0.5*eta*(2 + 5*eta)*Gm_r*rdot*rdot
                   + 0.5*(7 - 10*eta - 9*eta*eta)*Gm_r*vsq
                   + 0.25*(14 - 41*eta + 4*eta*eta)*Gm_r*Gm_r);
    }
    if(ode.corrections & C_3PN) {                       // Mora, Will PRD69(04)104021, Eq. (2.8d), (2.12d)
        Lz += Lnz*( (5.0/2.0 - (5199.0/280.0 - 41.0/32.0*PI*PI)*eta - 7*eta2 + eta*eta2)*Gm_r*Gm_r*Gm_r
             + 1.0/16.0*(5 - 59*eta + 238*eta2 - 323*eta*eta2)*vsq*vsq*vsq
             + 1.0/12.0*(135 - 322*eta + 315*eta2 - 108*eta*eta2)*Gm_r*Gm_r*vsq
             + 1.0/24.0*(12 - 287*eta - 951*eta2 - 324*eta*eta2)*Gm_r*Gm_r*rdot*rdot
             + 1.0/8.0*(33 - 142*eta + 106*eta2 + 195*eta*eta2)*Gm_r*vsq*vsq
             - 0.25*eta*(12 - 7*eta - 75*eta2)*Gm_r*vsq*rdot*rdot
             + 3.0/8.0*eta*(2 - 2*eta - 11*eta2)*Gm_r*rdot*rdot*rdot*rdot );
    }
    return Lz;
}

static FUNCDEF(L2) {
    REAL Lx = eval_Lx(t, f, ode, obs);
    REAL Ly = eval_Ly(t, f, ode, obs);
    REAL Lz = eval_Lz(t, f, ode, obs);

    return Lx*Lx + Ly*Ly + Lz*Lz;
}

/* Total angular momentum */
static FUNCDEF(Jx) {
    REAL Lx = eval_Lx(t, f, ode, obs);
    REAL cSs1 = ode.m1*ode.m1;
    REAL S1x = cSs1*f[I_s1x];
    REAL cSs2 = ode.m2*ode.m2;
    REAL S2x = cSs2*f[I_s2x];
    return Lx + S1x + S2x;
}

static FUNCDEF(Jy) {
    REAL Ly = eval_Ly(t, f, ode, obs);
    REAL cSs1 = ode.m1*ode.m1;
    REAL S1y = cSs1*f[I_s1y];
    REAL cSs2 = ode.m2*ode.m2;
    REAL S2y = cSs2*f[I_s2y];
    return Ly + S1y + S2y;
}

static FUNCDEF(Jz) {
    REAL Lz = eval_Lz(t, f, ode, obs);
    REAL cSs1 = ode.m1*ode.m1;
    REAL S1z = cSs1*f[I_s1z];
    REAL cSs2 = ode.m2*ode.m2;
    REAL S2z = cSs2*f[I_s2z];
    return Lz + S1z + S2z;
}

static FUNCDEF(J2) {
    REAL Jx = eval_Jx(t, f, ode, obs);
    REAL Jy = eval_Jy(t, f, ode, obs);
    REAL Jz = eval_Jz(t, f, ode, obs);
    
    return Jx*Jx + Jy*Jy + Jz*Jz;
} 

/** Angular momentum of the binary, with radiated ang.mom subtracted. */
static FUNCDEF(J_tot) {
    REAL Jx = eval_Jx(t, f, ode, obs);
    REAL Jy = eval_Jy(t, f, ode, obs);
    REAL Jz = eval_Jz(t, f, ode, obs);
    REAL Jx_rad = eval_Jx_Nrad(t, f, ode, obs);
    REAL Jy_rad = eval_Jy_Nrad(t, f, ode, obs);
    REAL Jz_rad = eval_Jz_Nrad(t, f, ode, obs);

    REAL Jtot2 = (Jx + Jx_rad)*(Jx + Jx_rad) 
                +(Jy + Jy_rad)*(Jy + Jy_rad) 
                +(Jz + Jz_rad)*(Jz + Jz_rad);

    return sqrt(Jtot2);
}

/* Eccentricity defined by the Keplerian relation, Gergely etal. PRD72(05)104022 */
static FUNCDEF(ecc_N) {
    REAL m = ode.m;
    REAL mu = ode.mu;

    REAL E = eval_E(t, f, ode, obs);
    REAL L2 = eval_L2(t, f, ode, obs);
    REAL Abar = sqrt(m*m*mu*mu + 2*E*L2/mu);

    REAL ecc = Abar/m/mu;
//    if(ode.corrections & C_PN) {
//        ecc += E/4.0/m/mu/mu/Abar*((5*eta - 15)*Abar*Abar
//                  - (9 + eta)*m*m*mu*mu);
//    }

    return ecc;
//    if (ecc2 < 0) {return 0;}
//     else {return sqrt(ecc2);}
}

/* Eccentricity defined by the 2PN relation, Arun etal. PRD77(08)064035, Eq.(7.7d) */
static FUNCDEF(ecc_PN) {
    REAL m = ode.m;
    REAL mu = ode.mu;
    REAL eta = ode.eta;

    REAL E = eval_E(t, f, ode, obs);
    REAL L2 = eval_L2(t, f, ode, obs);
    REAL ep = -2*E/mu;
    REAL j = -2*E*L2/m/m/mu/mu/mu;
//    REAL Abar = sqrt(m*m*mu*mu + 2*E*L2/mu);

    REAL ecc2 = 1 - j;
    if(ode.corrections & C_PN) {
        ecc2 += ep/4.0*(24 - 4*eta + 5*j*(-3 + eta));
    }
    if(ode.corrections & C_2PN) {
        ecc2 += ep*ep/8.0*(52 + 2*eta + 2*eta*eta
                       -j*(80 - 55*eta + 4*eta*eta)
                       -8.0/j*(-17 + 11*eta));
    }

//    return ecc;
    if (ecc2 < 0) {return 0;}
     else {return sqrt(ecc2);}
}

static FUNCDEF(ecc_ell) {
    REAL rx = eval_rx(t, f, ode, obs);
    REAL ry = eval_ry(t, f, ode, obs);
    REAL rz = eval_rz(t, f, ode, obs);
    REAL rsq = rx*rx + ry*ry + rz*rz;
    REAL r = sqrt(rsq);

    REAL vx = eval_vx(t, f, ode, obs);
    REAL vy = eval_vy(t, f, ode, obs);
//    REAL vz = eval_vz(t, f, ode, obs);

    return (rx*vx + ry*vy)/r/vx;
}

/* Eccentricity calculated as (rmax-rmin)/(rmax+rmin) */
static FUNCDEF(ecc_r) {
     return ecc_r;
}

/* Orbits calculated using the turning points */
static FUNCDEF(orbits_turn) {
     return orbits_turn;
}

static FUNCDEF(ecc_RL) {  // eccentricity calculted from the length of the Laplace-Runge-Lenz vector
    REAL m = ode.m;
//    REAL mu = ode.mu;

    REAL r = eval_r(t, f, ode, obs);
    REAL v2 = eval_v2(t, f, ode, obs);
    REAL rdot = eval_rdot(t, f, ode, obs);

    return sqrt(r*r*(v2-m/r)*(v2-m/r) - r*r*rdot*rdot*(v2-2.*m/r))/m;
}

static map<string, function_t> create_function_map()
{
    map<string, function_t> m;
    m["t"] = eval_t;
    m["rNx"] = eval_rNx;
    m["rNy"] = eval_rNy;
    m["rNz"] = eval_rNz;
    m["vNx"] = eval_vNx;
    m["vNy"] = eval_vNy;
    m["vNz"] = eval_vNz;
    m["rPNx"] = eval_rPNx;
    m["rPNy"] = eval_rPNy;
    m["rPNz"] = eval_rPNz;
    m["vPNx"] = eval_vPNx;
    m["vPNy"] = eval_vPNy;
    m["vPNz"] = eval_vPNz;
    m["r2PNx"] = eval_r2PNx;
    m["r2PNy"] = eval_r2PNy;
    m["r2PNz"] = eval_r2PNz;
    m["v2PNx"] = eval_v2PNx;
    m["v2PNy"] = eval_v2PNy;
    m["v2PNz"] = eval_v2PNz;
    m["r3PNx"] = eval_r3PNx;
    m["r3PNy"] = eval_r3PNy;
    m["r3PNz"] = eval_r3PNz;
    m["v3PNx"] = eval_v3PNx;
    m["v3PNy"] = eval_v3PNy;
    m["v3PNz"] = eval_v3PNz;
    m["rSOx"] = eval_rSOx;
    m["rSOy"] = eval_rSOy;
    m["rSOz"] = eval_rSOz;
    m["vSOx"] = eval_vSOx;
    m["vSOy"] = eval_vSOy;
    m["vSOz"] = eval_vSOz;
    m["rSSx"] = eval_rSSx;
    m["rSSy"] = eval_rSSy;
    m["rSSz"] = eval_rSSz;
    m["vSSx"] = eval_vSSx;
    m["vSSy"] = eval_vSSy;
    m["vSSz"] = eval_vSSz;
    m["rRRx"] = eval_rRRx;
    m["rRRy"] = eval_rRRy;
    m["rRRz"] = eval_rRRz;
    m["vRRx"] = eval_vRRx;
    m["vRRy"] = eval_vRRy;
    m["vRRz"] = eval_vRRz;
    m["rPNSOx"] = eval_rPNSOx;
    m["rPNSOy"] = eval_rPNSOy;
    m["rPNSOz"] = eval_rPNSOz;
    m["vPNSOx"] = eval_vPNSOx;
    m["vPNSOy"] = eval_vPNSOy;
    m["vPNSOz"] = eval_vPNSOz;
    m["r2PNSOx"] = eval_r2PNSOx;
    m["r2PNSOy"] = eval_r2PNSOy;
    m["r2PNSOz"] = eval_r2PNSOz;
    m["v2PNSOx"] = eval_v2PNSOx;
    m["v2PNSOy"] = eval_v2PNSOy;
    m["v2PNSOz"] = eval_v2PNSOz;
    m["r1RRx"] = eval_r1RRx;
    m["r1RRy"] = eval_r1RRy;
    m["r1RRz"] = eval_r1RRz;
    m["v1RRx"] = eval_v1RRx;
    m["v1RRy"] = eval_v1RRy;
    m["v1RRz"] = eval_v1RRz;
    m["rRRSOx"] = eval_rRRSOx;
    m["rRRSOy"] = eval_rRRSOy;
    m["rRRSOz"] = eval_rRRSOz;
    m["vRRSOx"] = eval_vRRSOx;
    m["vRRSOy"] = eval_vRRSOy;
    m["vRRSOz"] = eval_vRRSOz;
    m["rRRSSx"] = eval_rRRSSx;
    m["rRRSSy"] = eval_rRRSSy;
    m["rRRSSz"] = eval_rRRSSz;
    m["vRRSSx"] = eval_vRRSSx;
    m["vRRSSy"] = eval_vRRSSy;
    m["vRRSSz"] = eval_vRRSSz;
    m["s1x"] = eval_s1x;
    m["s1y"] = eval_s1y;
    m["s1z"] = eval_s1z;
    m["s2x"] = eval_s2x;
    m["s2y"] = eval_s2y;
    m["s2z"] = eval_s2z;
    m["orbits"] = eval_orbits;
    m["rx"] = eval_rx;
    m["ry"] = eval_ry;
    m["rz"] = eval_rz;
    m["r"] = eval_r;
    m["mr"] = eval_mr;
    m["vx"] = eval_vx;
    m["vy"] = eval_vy;
    m["vz"] = eval_vz;
    m["v"] = eval_v;
    m["v2"] = eval_v2;
    m["h_+"] = eval_hp;
    m["h_x"] = eval_hx;
    m["h"] = eval_h;
    m["E"] = eval_E; // total non-conserved energy of the compact binary
    m["E_N"] = eval_E_N; // Newtonian energy contribution
    m["E_PN"] = eval_E_PN; // post-Newtonian energy correction
    m["E_2PN"] = eval_E_2PN; // 2nd order post-Newtonian correction
    m["E_SO"] = eval_E_SO; // energy contribution from spin-orbit coupling
    m["E_SS"] = eval_E_SS; // energy contribution from spin-spin coupling
    m["E_3PN"] = eval_E_3PN; // 3rd order post-Newtonian correction
    m["E_PNSO"] = eval_E_PNSO; // 2.5rd order post-Newtonian SO correction
    m["E_rad"] = eval_E_rad; // radiated energy loss
    m["E_PNrad"] = eval_E_PNrad; // radiated energy Newtonian contribution
    m["E_2PNrad"] = eval_E_2PNrad; // radiated energy PN correction
    m["E_SOrad"] = eval_E_SOrad; // radiated energy spin-orbit correction
    m["E_SSrad"] = eval_E_SSrad; // radiated energy spin-spin correction
    m["E_Nrad"] = eval_E_Nrad; // radiated energy Newtonian correction
    m["Jx_Nrad"] = eval_Jx_Nrad; // radiated ang.mom Newtonian correction
    m["Jy_Nrad"] = eval_Jy_Nrad;
    m["Jz_Nrad"] = eval_Jz_Nrad;
    m["E_25PNrad"] = eval_E_25PNrad; // radiated energy 2.5PN correction
    m["E_PNSOrad"] = eval_E_PNSOrad; // radiated energy PNSO correction
    m["E_PNtot"] = eval_E_PNtot; // PN correction to total conserved energy
    m["E_2PNtot"] = eval_E_2PNtot; // 2PN correction to total conserved energy
    m["E_SOtot"] = eval_E_SOtot; // spin-orbit correction to total conserved energy
    m["E_SStot"] = eval_E_SStot; // spin-spin correction to total conserved energy
    m["E_tot"] = eval_E_tot; // total conserved energy (energy loss subtracted)
    m["J_tot"] = eval_J_tot; // total conserved ang.mom (ang.mom loss subtracted)
    m["x1"] = eval_x1;
    m["y1"] = eval_y1;
    m["z1"] = eval_z1;
    m["x2"] = eval_x2;
    m["y2"] = eval_y2;
    m["z2"] = eval_z2;
    m["hkp"] = eval_hkp;
    m["hp22"] = eval_hp22; // Spherical modes
    m["hx22"] = eval_hx22;
    m["hp21"] = eval_hp21;
    m["hx21"] = eval_hx21;
    m["hp2m2"] = eval_hp2m2;
    m["hx2m2"] = eval_hx2m2;
    m["hp2m1"] = eval_hp2m1;
    m["hx2m1"] = eval_hx2m1;
    m["hp20"] = eval_hp20;
    m["hx20"] = eval_hx20;
    m["ecc_N"] = eval_ecc_N;
    m["ecc_PN"] = eval_ecc_PN;
    m["ecc_ell"] = eval_ecc_ell;
    m["ecc_r"] = eval_ecc_r;
    m["ecc_RL"] = eval_ecc_RL;
    m["orbits_turn"] = eval_orbits_turn;
    m["Lnx"] = eval_Lnx; // Newtonian orbital angular momentum
    m["Lny"] = eval_Lny;
    m["Lnz"] = eval_Lnz;
    m["Ln2"] = eval_Ln2; // square of the Newtonian orbital angular momentum
    m["Lx"] = eval_Lx; // Orbital angular momentum
    m["Ly"] = eval_Ly; 
    m["Lz"] = eval_Lz;
    m["L2"] = eval_L2; // square of the orbital angular momentum
    m["Jx"] = eval_Jx; // Total angular momentum
    m["Jy"] = eval_Jy;
    m["Jz"] = eval_Jz;
    m["J2"] = eval_J2; // magnitude of the total angular momentum
    m["drAdiab"] = eval_drAdiabatic; // Adiabatic approximations for dr
    m["rdot"] = eval_rdot; // \dot r 
    m["orbfreq"] = eval_orbfreq; // orbital frequency
    return m;
}

const string CBwaveODE::COMPONENT_NAMES[] = {
    "rNx", "rNy", "rNz", "vNx", "vNy", "vNz",
    "rPNx", "rPNy", "rPNz", "vPNx", "vPNy", "vPNz", "E_PNrad",
    "r2PNx", "r2PNy", "r2PNz", "v2PNx", "v2PNy", "v2PNz", "E_2PNrad",
    "rSOx", "rSOy", "rSOz", "vSOx", "vSOy", "vSOz", "E_SOrad",
    "rSSx", "rSSy", "rSSz", "vSSx", "vSSy", "vSSz", "E_SSrad",
    "rRRx", "rRRy", "rRRz", "vRRx", "vRRy", "vRRz", "E_Nrad", "Jx_Nrad", "Jy_Nrad", "Jz_Nrad", "E_25PNrad",
    "s1x", "s1y", "s1z", "s2x", "s2y", "s2z", "orbits",
    "rPNSOx", "rPNSOy", "rPNSOz", "vPNSOx", "vPNSOy", "vPNSOz", "E_PNSOrad",
    "r3PNx", "r3PNy", "r3PNz", "v3PNx", "v3PNy", "v3PNz",
    "r1RRx", "r1RRy", "r1RRz", "v1RRx", "v1RRy", "v1RRz",
    "rRRSOx", "rRRSOy", "rRRSOz", "vRRSOx", "vRRSOy", "vRRSOz",
    "rRRSSx", "rRRSSy", "rRRSSz", "vRRSSx", "vRRSSy", "vRRSSz",
    "r2PNSOx", "r2PNSOy", "r2PNSOz", "v2PNSOx", "v2PNSOy", "v2PNSOz"
};

static map<string, function_t> FUNCMAP = create_function_map();

bool CBwaveODE::set(const string& name, REAL value, tvalarray<REAL>& vars)
{
    if(ODE<REAL>::set(name, value, vars)) {
	return true;
    }
    if(name == "m1") {
	m1 = value;
	init_m();
	return true;
    }
    if(name == "m2") {
	m2 = value;
	init_m();
	return true;
    }
    return false;
}

void CBwaveODE::eval(const REAL* f, int offset, REAL t, REAL* df)
{
    REAL rx = eval_rx(t, f, *this, *(ObserverParameters*)0);
    REAL ry = eval_ry(t, f, *this, *(ObserverParameters*)0);
    REAL rz = eval_rz(t, f, *this, *(ObserverParameters*)0);
    REAL rsq = rx*rx + ry*ry + rz*rz;
    REAL r = sqrt(rsq);
    REAL vx = eval_vx(t, f, *this, *(ObserverParameters*)0);
    REAL vy = eval_vy(t, f, *this, *(ObserverParameters*)0);
    REAL vz = eval_vz(t, f, *this, *(ObserverParameters*)0);
    REAL nx = rx/r;
    REAL ny = ry/r;
    REAL nz = rz/r;
    REAL Gm_r2 = m/rsq;
    REAL rdot = (rx*vx + ry*vy + rz*vz)/r;
    REAL rdot2 = rdot*rdot;
    REAL vsq = vx*vx + vy*vy + vz*vz;
    REAL v = sqrt(vsq);
    df[I_vNx] = -Gm_r2*nx;
    df[I_vNy] = -Gm_r2*ny;
    df[I_vNz] = -Gm_r2*nz;
    df[I_rNx] = f[I_vNx];
    df[I_rNy] = f[I_vNy];
    df[I_rNz] = f[I_vNz];

    REAL LNx = mu*(ry*vz - rz*vy);
    REAL LNy = mu*(rz*vx - rx*vz);
    REAL LNz = mu*(rx*vy - ry*vx);
    REAL m12 = m1/m2;
    REAL m21 = m2/m1;
    REAL G_c2r3 = 1.0/(rsq*r);

    REAL cSs1 = m1*m1;
    REAL S1x = cSs1*f[I_s1x];
    REAL S1y = cSs1*f[I_s1y];
    REAL S1z = cSs1*f[I_s1z];
    REAL cSs2 = m2*m2;
    REAL S2x = cSs2*f[I_s2x];
    REAL S2y = cSs2*f[I_s2y];
    REAL S2z = cSs2*f[I_s2z];

    REAL Sx = S1x + S2x;
    REAL Sy = S1y + S2y;
    REAL Sz = S1z + S2z;
    REAL sigmax = m21*S1x + m12*S2x;
    REAL sigmay = m21*S1y + m12*S2y;
    REAL sigmaz = m21*S1z + m12*S2z;

    REAL nS1 = nx*S1x + ny*S1y + nz*S1z;
    REAL nS2 = nx*S2x + ny*S2y + nz*S2z;
    REAL vS1 = vx*S1x + vy*S1y + vz*S1z;
    REAL vS2 = vx*S2x + vy*S2y + vz*S2z;
    REAL S1S2 = S1x*S2x + S1y*S2y + S1z*S2z;

    if(corrections & C_PN) {
	REAL cPNn = -Gm_r2*((1 + 3*eta)*v*v
		- 2*(2 + eta)*m/r
		- 1.5*eta*rdot*rdot);
	REAL cPNv = Gm_r2*2*(2-eta)*rdot;
	df[I_vPNx] = cPNn*nx + cPNv*vx;
	df[I_vPNy] = cPNn*ny + cPNv*vy;
	df[I_vPNz] = cPNn*nz + cPNv*vz;
	df[I_rPNx] = f[I_vPNx];
	df[I_rPNy] = f[I_vPNy];
	df[I_rPNz] = f[I_vPNz];
	// radiated energy loss 1PN contribution, Kidder (3.25b)
        df[I_E_PNrad] = 2.0/105.0*m*m*mu*mu/(rsq*rsq)
                        *((785 - 852*eta)*vsq*vsq - 160*(17 - eta)*m/r*vsq
                          + 8*(367 - 15*eta)*m/r*rdot*rdot
                          - 2*(1487 - 1392*eta)*vsq*rdot*rdot
                          + 3*(687 - 620*eta)*rdot*rdot*rdot*rdot
                          + 16*(1 - 4*eta)*m*m/rsq);
    }
    if(corrections & C_2PN) {
	REAL Gm_r = m/r;
	REAL c2PNn = -Gm_r2*
	    ((9 + 87.0/4.0*eta)*sqr(m/r)
	     + eta*(3 - 4*eta)*v*v*v*v
	     + 15.0/8.0*eta*(1 - 3*eta)*rdot*rdot*rdot*rdot
	     - 1.5*eta*(3 - 4*eta)*sqr(v*rdot)
	     - 0.5*eta*(13 - 4*eta)*v*v*Gm_r
	     - (2 + 25*eta + 2*eta*eta)*rdot*rdot*Gm_r);
	REAL c2PNv = Gm_r2*0.5*rdot
	    *(eta*(15+4*eta)*v*v
		    - (4 + 41*eta + 8*eta*eta)*Gm_r
		    - 3*eta*(3+2*eta)*rdot*rdot);
	df[I_v2PNx] = c2PNn*nx + c2PNv*vx;
	df[I_v2PNy] = c2PNn*ny + c2PNv*vy;
	df[I_v2PNz] = c2PNn*nz + c2PNv*vz;
	df[I_r2PNx] = f[I_v2PNx];
	df[I_r2PNy] = f[I_v2PNy];
	df[I_r2PNz] = f[I_v2PNz];
	// radiated energy loss 2PN correction, Gopakumar and Iyer PRD56(97)7708, Eq. (3.5d)
	df[I_E_2PNrad] = 8.0/15.0*m*m*mu*mu/(rsq*rsq)
			*(1/42.0*(1692 - 5497*eta + 4430*eta*eta)*vsq*vsq*vsq
			  - 1/14.0*(1719 - 10278*eta + 6292*eta*eta)*vsq*vsq*rdot*rdot
			  - 1/21.0*(4446 - 5237*eta + 1393*eta*eta)*vsq*vsq*m/r
                          + 1/14.0*(2018 - 15207*eta + 7572*eta*eta)*vsq*rdot*rdot*rdot*rdot
                          + 1/7.0*(4987 - 8513*eta + 2165*eta*eta)*vsq*m/r*rdot*rdot
                          + 1/756.0*(281473 + 81828*eta + 4368*eta*eta)*vsq*m/r*m/r
                          - 1/42.0*(2501 - 20234*eta + 8404*eta*eta)*rdot*rdot*rdot*rdot*rdot*rdot
                          - 1/63.0*(33510 - 60971*eta + 14290*eta*eta)*m/r*rdot*rdot*rdot*rdot
                          - 1/252.0*(106319 + 9798*eta + 5376*eta*eta)*m/r*m/r*rdot*rdot
                          + 2/63.0*(-253 + 1026*eta - 56*eta*eta)*m/r*m/r*m/r);
    }
    if(corrections & C_SO) {
	REAL cSOn = 6*((ny*vz-vy*nz)*(Sx+sigmax)    // 6 (n x v) (S + sigma)
		+ (nz*vx-nx*vz)*(Sy+sigmay)
		+ (nx*vy-ny*vx)*(Sz+sigmaz));
	df[I_vSOx] = G_c2r3*(cSOn*nx
		- (vy*(4*Sz + 3*sigmaz) - vz*(4*Sy + 3*sigmay))
		+ 3*rdot*(ny*(2*Sz + sigmaz) - nz*(2*Sy + sigmay)));
	df[I_vSOy] = G_c2r3*(cSOn*ny
		- (vz*(4*Sx + 3*sigmax) - vx*(4*Sz + 3*sigmaz))
		+ 3*rdot*(nz*(2*Sx + sigmax) - nx*(2*Sz + sigmaz)));
	df[I_vSOz] = G_c2r3*(cSOn*nz
		- (vx*(4*Sy + 3*sigmay) - vy*(4*Sx + 3*sigmax))
		+ 3*rdot*(nx*(2*Sy + sigmay) - ny*(2*Sx + sigmax)));
	df[I_rSOx] = f[I_vSOx];
	df[I_rSOy] = f[I_vSOy];
	df[I_rSOz] = f[I_vSOz];
	// radiated energy loss spin-orbit correction, Kidder (3.25c)
	REAL Deltax = m*(S2x/m2 - S1x/m1);
	REAL Deltay = m*(S2y/m2 - S1y/m1);
	REAL Deltaz = m*(S2z/m2 - S1z/m1);
	REAL LNS = LNx*Sx + LNy*Sy + LNz*Sz;
	REAL LNDelta = LNx*Deltax + LNy*Deltay + LNz*Deltaz;
	df[I_E_SOrad] = 8.0/15.0*m*mu/(rsq*rsq*rsq)
		*(LNS*(78*rdot*rdot - 80*vsq - 8*m/r)
		  + LNDelta*(m1-m2)/m*(51*rdot*rdot - 43*vsq + 4*m/r));
    }
    if(corrections & (C_SO | C_SS)) {			// Kidder (2.4), Wang&Will PRD75(07)064017, Eq. (1.6), Blanchet,Buonanno,Faye PRD74(06)104034, Eq. (7.7)
							// Boh� etal CQG30(13)075017, Eq. (3.4)
	REAL c12 = (4.+3.*m12)/2. + m2*m2/m/m/8.*((16. + 22.*m12 + 10.*m12*m12 + m12*m12*m12)*vsq
		   - 2.*(4. + 7.*m12 + 7.*m12*m12 + 2.*m12*m12*m12)*m/r - 6.*(4. + 3.*m12)*rdot2) 
                   + m2*m2*m2*m2/m/m/m/m/16.*(15.*(4. + 3.*m12)*rdot2*rdot2 
                    + 2*(-12. - 49.*m12 - 52.*m12*m12 - 8.*m12*m12*m12 + 18.*m12*m12*m12*m12 + 7.*m12*m12*m12*m12*m12)*m*m/r/r
		    - 3*(32. + 18.*m12 + 14.*m12*m12 + 13.*m12*m12*m12)*rdot2*vsq
                    + (132. + 553.*m12 + 765.*m12*m12 + 476.*m12*m12*m12 + 120.*m12*m12*m12*m12 + 6.*m12*m12*m12*m12*m12)*vsq*m/r
                    + (-292. - 1331.*m12 - 1761.*m12*m12 - 967.*m12*m12*m12 - 185.*m12*m12*m12*m12 + 8.*m12*m12*m12*m12*m12)*rdot2*m/r
                    + (32. + 21.*m12 + 17.*m12*m12 + 22.*m12*m12*m12 + 12.*m12*m12*m12*m12 + m12*m12*m12*m12*m12)*vsq*vsq );
	REAL c21 = (4.+3.*m21)/2. + m1*m1/m/m/8.*((16. + 22.*m21 + 10.*m21*m21 + m21*m21*m21)*vsq
		   - 2.*(4. + 7.*m21 + 7.*m21*m21 + 2.*m21*m21*m21)*m/r - 6.*(4. + 3.*m21)*rdot2) 
                   + m1*m1*m1*m1/m/m/m/m/16.*(15.*(4. + 3.*m21)*rdot2*rdot2 
                    + 2*(-12. - 49.*m21 - 52.*m21*m21 - 8.*m21*m21*m21 + 18.*m21*m21*m21*m21 + 7.*m21*m21*m21*m21*m21)*m*m/r/r
		    - 3*(32. + 18.*m21 + 14.*m21*m21 + 13.*m21*m21*m21)*rdot2*vsq
                    + (132. + 553.*m21 + 765.*m21*m21 + 476.*m21*m21*m21 + 120.*m21*m21*m21*m21 + 6.*m21*m21*m21*m21*m21)*vsq*m/r
                    + (-292. - 1331.*m21 - 1761.*m21*m21 - 967.*m21*m21*m21 - 185.*m21*m21*m21*m21 + 8.*m21*m21*m21*m21*m21)*rdot2*m/r
                    + (32. + 21.*m21 + 17.*m21*m21 + 22.*m21*m21*m21 + 12.*m21*m21*m21*m21 + m21*m21*m21*m21*m21)*vsq*vsq );
	df[I_s1x] = G_c2r3*((c21*LNy - S2y + 3.*nS2*ny)*f[I_s1z]
		- (c21*LNz - S2z + 3.*nS2*nz)*f[I_s1y])
		+ m1*m2/(rsq*rsq*r)*(2./3.*vS2 + 30.*rdot*nS2)*(ny*f[I_s1z]-nz*f[I_s1y]);
	df[I_s1y] = G_c2r3*((c21*LNz - S2z + 3.*nS2*nz)*f[I_s1x]
		- (c21*LNx - S2x + 3.*nS2*nx)*f[I_s1z])
		+ m1*m2/(rsq*rsq*r)*(2./3.*vS2 + 30.*rdot*nS2)*(nz*f[I_s1x]-nx*f[I_s1z]);
	df[I_s1z] = G_c2r3*((c21*LNx - S2x + 3.*nS2*nx)*f[I_s1y]
		- (c21*LNy - S2y + 3.*nS2*ny)*f[I_s1x])
		+ m1*m2/(rsq*rsq*r)*(2./3.*vS2 + 30.*rdot*nS2)*(nx*f[I_s1y]-ny*f[I_s1x]);
	df[I_s2x] = G_c2r3*((c12*LNy - S1y + 3.*nS1*ny)*f[I_s2z]
		- (c12*LNz - S1z + 3.*nS1*nz)*f[I_s2y])
		+ m1*m2/(rsq*rsq*r)*(2./3.*vS1 + 30.*rdot*nS1)*(ny*f[I_s2z]-nz*f[I_s2y]);
	df[I_s2y] = G_c2r3*((c12*LNz - S1z + 3*nS1*nz)*f[I_s2x]
		- (c12*LNx - S1x + 3.*nS1*nx)*f[I_s2z])
		+ m1*m2/(rsq*rsq*r)*(2./3.*vS1 + 30.*rdot*nS1)*(nz*f[I_s2x]-nx*f[I_s2z]);
	df[I_s2z] = G_c2r3*((c12*LNx - S1x + 3.*nS1*nx)*f[I_s2y]
		- (c12*LNy - S1y + 3.*nS1*ny)*f[I_s2x])
		+ m1*m2/(rsq*rsq*r)*(2./3.*vS1 + 30.*rdot*nS1)*(nx*f[I_s2y]-ny*f[I_s2x]);
    }
    if(corrections & C_SS) {
	REAL cSS = 3.0/(mu*r*r*r*r);
	df[I_vSSx] = -cSS*(nx*S1S2 + S1x*nS2 + S2x*nS1 - 5*nx*nS1*nS2);
	df[I_vSSy] = -cSS*(ny*S1S2 + S1y*nS2 + S2y*nS1 - 5*ny*nS1*nS2);
	df[I_vSSz] = -cSS*(nz*S1S2 + S1z*nS2 + S2z*nS1 - 5*nz*nS1*nS2);
	df[I_rSSx] = f[I_vSSx];
	df[I_rSSy] = f[I_vSSy];
	df[I_rSSz] = f[I_vSSz];
	// radiated energy loss spin-spin correction, Kidder (3.25e)
	df[I_E_SSrad] = 4.0/15.0*m*mu/(rsq*rsq*rsq)
	    *(-3*nS1*nS2*(168*vsq - 269*rdot*rdot)
		    + 3*S1S2*(47*vsq - 55*rdot*rdot) + 71*vS1*vS2
		    - 171*rdot*(vS1*nS2 + nS1*vS2));
    }
    if(corrections & C_RR) {
	REAL cRR = 8./5.*eta*m*m/(r*r*r);
//	REAL cRR1 = cRR*rdot*(18.*v*v + 2./3.*m/r - 25.*rdot2);		// BT
//	REAL cRR2 = cRR*(6.*v*v - 2.*m/r - 15.*rdot2);
//	REAL cRR1 = cRR*rdot*(3.*v*v + 17./3.*m/r);			// DD
//	REAL cRR2 = cRR*(v*v + 3.*m/r);
// General form, Eq.(1.5) of Iyer,Will PRD52(95)6882
//	REAL alpha = 0.;		// gauge params for RR
//	REAL beta = 0.;
	REAL cRR1 = cRR*rdot*(3.*(1.+beta)*vsq + 1./3.*(23.+6.*alpha-9.*beta)*m/r - 5.*beta*rdot2);	
	REAL cRR2 = cRR*((2.+alpha)*vsq + (2.-alpha)*m/r - 3.*(1.+alpha)*rdot2);
	df[I_vRRx] = cRR1*nx - cRR2*vx;
	df[I_vRRy] = cRR1*ny - cRR2*vy;
	df[I_vRRz] = cRR1*nz - cRR2*vz;
	df[I_rRRx] = f[I_vRRx];
	df[I_rRRy] = f[I_vRRy];
	df[I_rRRz] = f[I_vRRz];
        // radiated energy loss Newtonian contribution, Kidder (3.25a)
        df[I_E_Nrad] = 8.0/15.0*m*m*mu*mu/(rsq*rsq)*(12*vsq - 11*rdot*rdot);
        // radiated angular momentum loss Newtonian contribution, Kidder (3.28a)
        df[I_Jx_Nrad] = 8.0/5.0*m*mu/(rsq*r)*LNx*(2*vsq - 3*rdot*rdot + 2*m/r);
        df[I_Jy_Nrad] = 8.0/5.0*m*mu/(rsq*r)*LNy*(2*vsq - 3*rdot*rdot + 2*m/r);
        df[I_Jz_Nrad] = 8.0/5.0*m*mu/(rsq*r)*LNz*(2*vsq - 3*rdot*rdot + 2*m/r);
        // radiated energy loss 2.5PN contribution, Arun etal. PRD77(08)064035, Eq. (5.2d)
        df[I_E_25PNrad] = 32.0/5.0*m*m*mu*mu/(rsq*rsq)*rdot*eta*( -12349.0/210.0*m/r*vsq*vsq
			+ 4524.0/35.0*m/r*vsq*rdot2
			- 2753.0/126.0*m*m/rsq*vsq
			- 985.0/14.0*m/r*rdot2*rdot2
			+ 13981.0/630.0*m*m/rsq*rdot2
			- 1.0/315.0*m*m*m/rsq/r );
    }
    if(corrections & C_PNSO) {			// Faye,Blanchet,Buonanno PRD74(06)104033, Eq. (5.7b)
						// Boh� etal CQG30(13)075017, Eq. (3.7)
	REAL cnvS   = ( (ny*vz-vy*nz)*Sx + (nz*vx-nx*vz)*Sy + (nx*vy-ny*vx)*Sz );   // (n x v) S
	REAL Sigx   = S1x*(m21-1) + S2x*(m12-1);				// (\delta m)/m * Sigma
	REAL Sigy   = S1y*(m21-1) + S2y*(m12-1);
	REAL Sigz   = S1z*(m21-1) + S2z*(m12-1);
	REAL cnvSig = ( (ny*vz-vy*nz)*Sigx + (nz*vx-nx*vz)*Sigy	+ (nx*vy-ny*vx)*Sigz );	// (n x v) (\delta m)/m * Sigma
	REAL cPNSOn = cnvS*(-30*eta*rdot2 + 24*eta*vsq - (44 + 33*eta)*m/r)
			+ cnvSig*(-15*eta*rdot2 + 12*eta*vsq - (24 + 37./2.*eta)*m/r);
	REAL cPNSOv = rdot*(10.5*cnvS*(eta-1) + cnvSig*(6.*eta-4.5));
	df[I_vPNSOx] = 1./(rsq*r)*( cPNSOn*nx + cPNSOv*vx 
			+ (ny*Sz-Sy*nz)*rdot*(-22.5*eta*rdot2 + (-1.5+22.5*eta)*vsq - (28. + 29.*eta)*m/r) 
			+ (ny*Sigz + Sigy*nz)*rdot*(-15.*eta*rdot2 + (-1.5+12.*eta)*vsq - (12. + 15.5*eta)*m/r) 
			+ (vy*Sz-Sy*vz)*((1.5+15.*eta)*rdot2 - 14.*eta*vsq + (24. + 19.*eta)*m/r) 
			+ (vy*Sigz-Sigy*vz)*((1.5+9.*eta)*rdot2 - 7.*eta*vsq + (12. + 9.5*eta)*m/r) );
	df[I_vPNSOy] = 1./(rsq*r)*( cPNSOn*ny + cPNSOv*vy 
			+ (nz*Sx-Sz*nx)*rdot*(-22.5*eta*rdot2 + (-1.5+22.5*eta)*vsq - (28. + 29.*eta)*m/r) 
			+ (nz*Sigx + Sigz*nx)*rdot*(-15.*eta*rdot2 + (-1.5+12.*eta)*vsq - (12. + 15.5*eta)*m/r) 
			+ (vz*Sx-Sz*vx)*((1.5+15.*eta)*rdot2 - 14.*eta*vsq + (24. + 19.*eta)*m/r) 
			+ (vz*Sigx-Sigz*vx)*((1.5+9.*eta)*rdot2 - 7.*eta*vsq + (12. + 9.5*eta)*m/r) );
	df[I_vPNSOz] = 1./(rsq*r)*( cPNSOn*nz + cPNSOv*vz 
			+ (nx*Sy-Sx*ny)*rdot*(-22.5*eta*rdot2 + (-1.5+22.5*eta)*vsq - (28. + 29.*eta)*m/r) 
			+ (nx*Sigy + Sigx*ny)*rdot*(-15.*eta*rdot2 + (-1.5+12.*eta)*vsq - (12. + 15.5*eta)*m/r) 
			+ (vx*Sy-Sx*vy)*((1.5+15.*eta)*rdot2 - 14.*eta*vsq + (24. + 19.*eta)*m/r) 
			+ (vx*Sigy-Sigx*vy)*((1.5+9.*eta)*rdot2 - 7.*eta*vsq + (12. + 9.5*eta)*m/r) );
	df[I_rPNSOx] = f[I_vPNSOx];
	df[I_rPNSOy] = f[I_vPNSOy];
	df[I_rPNSOz] = f[I_vPNSOz];
	// radiated energy loss PNSO correction, Blanchet,Buonanno,Faye PRD74(06)104034, Eq. (5.12)
	// formula from Erratum, PRD81(10)089901(E), Eq. (2)
	REAL Deltax = m*(S2x/m2 - S1x/m1);
	REAL Deltay = m*(S2y/m2 - S1y/m1);
	REAL Deltaz = m*(S2z/m2 - S1z/m1);
	REAL LNS = LNx*Sx + LNy*Sy + LNz*Sz;
	REAL LNDelta = LNx*Deltax + LNy*Deltay + LNz*Deltaz;
	df[I_E_PNSOrad] = 8./105.*m*mu/(rsq*rsq*rsq)
		*( LNS*( rdot2*rdot2*(-2244. + 3144.*eta) 
			+ m*m/rsq*(944. + 390.*eta) 
			+ m/r*rdot2*(-3223. + 506.*eta) 
			+ rdot2*vsq*(3519. + 5004.*eta) 
			+ m/r*vsq*(3805. - 224.*eta) 
			+ vsq*vsq*(-1207. + 1810.*eta) )
		  + LNDelta*(m1-m2)/m*( rdot2*rdot2*(-7941./4. + 2676.*eta) 
			+ m*m/rsq*(-137. + 238.*eta) 
			+ m/r*rdot2*(-7327./2. + 1199.*eta) 
			+ rdot2*vsq*(2364. - 3621.*eta) 
			+ m/r*vsq*(5387./2. - 497.*eta) 
			+ vsq*vsq*(-2603./4. + 1040.*eta) ) );
    }
    if(corrections & C_3PN) {			// Mora, Will PRD69(04)104021, Eq. (2.8d), (2.9d)
        REAL Gm_r = m/r;
        REAL c3PNn = Gm_r2*
            ( (16 + (1399.0/12.0 - 41.0/16.0*PI*PI)*eta + 71.0/2.0*eta*eta)*Gm_r*Gm_r*Gm_r
             + eta*(20827.0/840.0 + 123.0/64.0*PI*PI - eta*eta)*Gm_r*Gm_r*vsq
             - (1 + (22717.0/168.0 + 615.0/64.0*PI*PI)*eta + 11.0/8.0*eta*eta - 7*eta*eta*eta)*Gm_r*Gm_r*rdot2
             - 0.25*eta*(11 - 49*eta + 52*eta*eta)*vsq*vsq*vsq
             + 35.0/16.0*eta*(1 - 5*eta + 5*eta*eta)*rdot2*rdot2*rdot2
             - 0.25*eta*(75 + 32*eta - 40*eta*eta)*Gm_r*vsq*vsq
	     - 0.5*eta*(158 - 69*eta - 60*eta*eta)*Gm_r*rdot2*rdot2
	     + eta*(121 - 16*eta - 20*eta*eta)*Gm_r*vsq*rdot2
	     + 3.0/8.0*eta*(20 - 79*eta + 60*eta*eta)*vsq*vsq*rdot2
	     - 15.0/8.0*eta*(4 - 18*eta + 17*eta*eta)*vsq*rdot2*rdot2 );
        REAL c3PNv = Gm_r2*rdot*
            ( (4 + (5849.0/840.0 + 123.0/32.0*PI*PI)*eta - 25*eta*eta - 8*eta*eta*eta)*Gm_r*Gm_r
             - 1.0/8.0*eta*(65 - 152*eta - 48*eta*eta)*vsq*vsq
             + 15.0/8.0*eta*(3 - 8*eta - 2*eta*eta)*rdot2*rdot2
             + eta*(15 + 27*eta + 10*eta*eta)*Gm_r*vsq
             - 1.0/6.0*eta*(329 + 177*eta + 108*eta*eta)*Gm_r*rdot2
             - 3.0/4.0*eta*(16 - 37*eta - 16*eta*eta)*vsq*rdot2 );
        df[I_v3PNx] = c3PNn*nx + c3PNv*vx;
        df[I_v3PNy] = c3PNn*ny + c3PNv*vy;
        df[I_v3PNz] = c3PNn*nz + c3PNv*vz;
        df[I_r3PNx] = f[I_v3PNx];
        df[I_r3PNy] = f[I_v3PNy];
        df[I_r3PNz] = f[I_v3PNz];
    }
    if(corrections & C_1RR) {			// Iyer, Will PRD52(95)6882, Eq. (1.8c,d)
        REAL c1RR = 8.0/5.0*eta*m*m/(r*r*r);
/*        REAL c1RR1 = c1RR*rdot*( (87.0/14.0 - 48*eta)*vsq*vsq
				- (5379.0/28.0 + 136.0/3.0*eta)*vsq*m/r 
				+ 25.0/2.0*(1 + 5*eta)*vsq*rdot2 
				+ (1353.0/4.0 + 133*eta)*rdot2*m/r 
				- 35.0/2.0*(1 - eta)*rdot2*rdot2
				+ (160.0/7.0 + 55.0/3.0*eta)*m/r*m/r );
        REAL c1RR2 = c1RR*( -27.0/14.0*vsq*vsq
                                - (4861.0/84.0 + 58.0/3.0*eta)*vsq*m/r
                                + 3.0/2.0*(13 - 37*eta)*vsq*rdot2
                                + (2591.0/12.0 + 97*eta)*rdot2*m/r
                                - 25.0/2.0*(1 - 7*eta)*rdot2*rdot2
                                + 1.0/3.0*(776.0/7.0 + 55*eta)*m/r*m/r );
        REAL c1RR1 = -c1RR*rdot*( 23.0/14.0*(43 + 14*eta)*m/r*m/r	// Mora, Will PRD69(04)104021, Eq. (2.8e), (2.9e)
				+ 3.0/28.0*(61 + 70*eta)*vsq*vsq
				+ 70*rdot2*rdot2
                                + 1.0/42.0*(519 - 1276*eta)*vsq*m/r
				+ 1.0/4.0*(147 + 188*eta)*rdot2*m/r
                                - 15.0/4.0*(19 + 2*eta)*vsq*rdot2 );
        REAL c1RR2 = -c1RR*( 1.0/42.0*(1325 + 546*eta)*m/r*m/r       
                                + 1.0/28.0*(313 + 42*eta)*vsq*vsq
                                + 75*rdot2*rdot2
                                - 1.0/42.0*(205 + 777*eta)*vsq*m/r
                                + 1.0/12.0*(205 + 424*eta)*rdot2*m/r
                                - 3.0/4.0*(113 + 2*eta)*vsq*rdot2 );
*/
// General form, Eq.(2.12) of Iyer,Will PRD52(95)6882
/*        REAL alpha = 4.;		// gauge params for RR, BT
        REAL beta = 5.;
        REAL delta1 = -99./14.+27.*eta;		
	REAL delta2 = 5.*(1. - 4.*eta);
	REAL delta3 = 274./7. + 67./21.*eta;
	REAL delta4 = 5./2.*(1. - eta);
	REAL delta5 = -1./7.*(292. + 57.*eta);
	REAL delta6 = 51./28. + 71./14.*eta;	
	REAL alpha = -1.;		// DD
	REAL beta = 0.;
	REAL delta1 = 271./28. + 6.*eta;
	REAL delta2 = -77./4. - 3./2.*eta;
	REAL delta3 = 79./14. - 92./7.*eta;
	REAL delta4 = 10.;
	REAL delta5 = 5./42. + 242./21.*eta;
	REAL delta6 = -439./28. + 18./7.*eta;
	REAL alpha = 0;			// general gauge
	REAL beta = 0;
	REAL delta1 = 0.;
	REAL delta2 = 0.;
	REAL delta3 = 0.;
	REAL delta4 = 0.;
	REAL delta5 = 0.;
	REAL delta6 = 0.;
*/
	REAL c1 = 1./28.*(117. + 132.*eta) - 3./2.*beta*(1. - 3.*eta) + 3.*delta2 - 3.*delta6;
	REAL c2 = -1./42.*(297. - 310.*eta) - 3.*alpha*(1. - 4.*eta) - 3./2.*beta*(7. + 13.*eta) - 2.*delta1 - 3.*delta2 + 3.*delta5 + 3.*delta6;
	REAL c3 = 5./28.*(19. - 72.*eta) + 5./2.*beta*(1. - 3.*eta) - 5.*delta2 + 5.*delta4 + 5.*delta6;	
	REAL c4 = -1./28.*(687. - 368.*eta) - 6.*alpha*eta + 1./2.*beta*(54. + 17.*eta) - 2.*delta2 - 5.*delta4 - 6.*delta5;
	REAL c5 = -7.*delta4;
	REAL c6 = -1./21.*(1533. + 498.*eta) - alpha*(14. + 9.*eta) + 3.*beta*(7. + 4.*eta) - 2.*delta3 - 3.*delta5;
	REAL d1 = -3.*(1. - 3.*eta) - 3./2.*alpha*(1. - 3.*eta) - delta1;
	REAL d2 = -1./84.*(139. + 768.*eta) - 1./2.*alpha*(5. + 17.*eta) + delta1 - delta3;
	REAL d3 = 1./28.*(369. - 624.*eta) + 3./2.*(3.*alpha + 2.*beta)*(1. - 3.*eta) + 3.*delta1 - 3.*delta6;
	REAL d4 = 1./42.*(295. - 335.*eta) + 1./2.*alpha*(38. - 11.*eta) - 3*beta*(1. - 3.*eta) + 2.*delta1 + 4.*delta3 + 3.*delta6;
	REAL d5 = 5./28.*(19. - 72.*eta) - 5.*beta*(1. - 3.*eta) + 5.*delta6;
	REAL d6 = -1./21.*(634. - 66.*eta) + alpha*(7. + 3.*eta) + delta3;
	REAL c1RR1 = c1RR*rdot*(c1*vsq*vsq + c2*vsq*m/r + c3*vsq*rdot2 + c4*rdot2*m/r + c5*rdot2*rdot2 + c6*m*m/r/r);
	REAL c1RR2 = c1RR*(d1*vsq*vsq + d2*vsq*m/r + d3*vsq*rdot2 + d4*rdot2*m/r + d5*rdot2*rdot2 + d6*m*m/r/r);
        df[I_v1RRx] = c1RR1*nx - c1RR2*vx;
        df[I_v1RRy] = c1RR1*ny - c1RR2*vy;
        df[I_v1RRz] = c1RR1*nz - c1RR2*vz;
        df[I_r1RRx] = f[I_v1RRx];
        df[I_r1RRy] = f[I_v1RRy];
        df[I_r1RRz] = f[I_v1RRz];
    }
    if(corrections & C_RRSO) {			// Will PRD71(05)084027, Eq. (B4b)
	REAL LNS = LNx*Sx + LNy*Sy + LNz*Sz;
	REAL LNsigma = LNx*sigmax + LNy*sigmay + LNz*sigmaz;
	REAL cRRSOn = rdot/r/mu*((120*vsq + 280*rdot2 + 453*m/r)*LNS + (120*vsq + 280*rdot2 + 458*m/r)*LNsigma);
	REAL cRRSOv = 1.0/r/mu*((87*vsq - 675*rdot2 - 901.0/3.0*m/r)*LNS + 4*(18*vsq - 150*rdot2 - 66*m/r)*LNsigma);
	df[I_vRRSOx] = -mu/5.0/(rsq*rsq)*( cRRSOn*nx + cRRSOv*vx
			- 2.0/3.0*rdot*(vy*Sz-Sy*vz)*(48*v*v + 15*rdot2 + 364*m/r)
			+ 1.0/3.0*rdot*(vy*sigmaz-sigmay*vz)*(291*v*v - 705*rdot2 - 772*m/r)
			+ 0.5*(ny*Sz-Sy*nz)*(31*vsq*vsq - 260*vsq*rdot2 + 245*rdot2*rdot2 
				- 689.0/3.0*vsq*m/r + 537*rdot2*m/r + 4.0/3.0*m*m/r/r)
			+ 0.5*(ny*sigmaz-sigmay*nz)*(115*vsq*vsq - 1130*vsq*rdot2 + 1295*rdot2*rdot2 
				- 869.0/3.0*vsq*m/r + 849*rdot2*m/r + 44.0/3.0*m*m/r/r) );
	df[I_vRRSOy] = -mu/5.0/(rsq*rsq)*( cRRSOn*ny + cRRSOv*vy
			- 2.0/3.0*rdot*(vz*Sx-Sz*vx)*(48*v*v + 15*rdot2 + 364*m/r)
			+ 1.0/3.0*rdot*(vz*sigmax-sigmaz*vx)*(291*v*v - 705*rdot2 - 772*m/r)
			+ 0.5*(nz*Sx-Sz*nx)*(31*vsq*vsq - 260*vsq*rdot2 + 245*rdot2*rdot2 
				- 689.0/3.0*vsq*m/r + 537*rdot2*m/r + 4.0/3.0*m*m/r/r)
			+ 0.5*(nz*sigmax-sigmaz*nx)*(115*vsq*vsq - 1130*vsq*rdot2 + 1295*rdot2*rdot2 
				- 869.0/3.0*vsq*m/r + 849*rdot2*m/r + 44.0/3.0*m*m/r/r) );
	df[I_vRRSOz] = -mu/5.0/(rsq*rsq)*( cRRSOn*nz + cRRSOv*vz
			- 2.0/3.0*rdot*(vx*Sy-Sx*vy)*(48*v*v + 15*rdot2 + 364*m/r)
			+ 1.0/3.0*rdot*(vx*sigmay-sigmax*vy)*(291*v*v - 705*rdot2 - 772*m/r)
			+ 0.5*(nx*Sy-Sx*ny)*(31*vsq*vsq - 260*vsq*rdot2 + 245*rdot2*rdot2 
				- 689.0/3.0*vsq*m/r + 537*rdot2*m/r + 4.0/3.0*m*m/r/r)
			+ 0.5*(nx*sigmay-sigmax*ny)*(115*vsq*vsq - 1130*vsq*rdot2 + 1295*rdot2*rdot2 
				- 869.0/3.0*vsq*m/r + 849*rdot2*m/r + 44.0/3.0*m*m/r/r) );
	df[I_rRRSOx] = f[I_vRRSOx];
	df[I_rRRSOy] = f[I_vRRSOy];
	df[I_rRRSOz] = f[I_vRRSOz];
    }
    if(corrections & C_RRSS) {			// Wang, Will PRD75(07)064017, Eq. (1.3)
	REAL cRRSSn = (287*rdot2 - 99*vsq + 541.0/5.0*m/r)*rdot*S1S2 
		      - (2646*rdot2 - 714*vsq + 1961.0/5.0*m/r)*rdot*nS1*nS2
		      + (1029*rdot2 - 123*vsq + 62.9*m/r)*(nS1*vS2+nS2*vS1)
		      - 336*rdot*vS1*vS2;
	REAL cRRSSv = (171.0/5.0*vsq - 195*rdot2 - 67*m/r)*S1S2 
		      - (174*vsq - 1386*rdot2 - 1038.0/5.0*m/r)*nS1*nS2
		      - 438*rdot*(nS1*vS2+nS2*vS1)
		      + 96*vS1*vS2;
	df[I_vRRSSx] = 1.0/(rsq*rsq*r)*( cRRSSn*nx + cRRSSv*vx
			+ (2.7*vsq - 37.5*rdot2 - 50.9/3.0*m/r)*(vS2*S1x + vS1*S2x) 
			+ (7.5*vsq + 38.5*rdot2 - 19.9*m/r)*rdot*(nS2*S1x + nS1*S2x) );
	df[I_vRRSSy] = 1.0/(rsq*rsq*r)*( cRRSSn*ny + cRRSSv*vy 
			+ (2.7*vsq - 37.5*rdot2 - 50.9/3.0*m/r)*(vS2*S1y + vS1*S2y) 
			+ (7.5*vsq + 38.5*rdot2 - 19.9*m/r)*rdot*(nS2*S1y + nS1*S2y) );
	df[I_vRRSSz] = 1.0/(rsq*rsq*r)*( cRRSSn*nz + cRRSSv*vz
			+ (2.7*vsq - 37.5*rdot2 - 50.9/3.0*m/r)*(vS2*S1z + vS1*S2z) 
			+ (7.5*vsq + 38.5*rdot2 - 19.9*m/r)*rdot*(nS2*S1z + nS1*S2z) );
	df[I_rRRSSx] = f[I_vRRSSx];
	df[I_rRRSSy] = f[I_vRRSSy];
	df[I_rRRSSz] = f[I_vRRSSz];
    }
    if(corrections & C_2PNSO) {			// Boh� etal CQG30(13)075017, Eq. (3.8)
	REAL cnvS   = ( (ny*vz-vy*nz)*Sx + (nz*vx-nx*vz)*Sy + (nx*vy-ny*vx)*Sz );   // (n x v) S
	REAL Sigx   = S1x*(m21-1) + S2x*(m12-1);				// (\delta m)/m * Sigma
	REAL Sigy   = S1y*(m21-1) + S2y*(m12-1);
	REAL Sigz   = S1z*(m21-1) + S2z*(m12-1);
	REAL cnvSig = ( (ny*vz-vy*nz)*Sigx + (nz*vx-nx*vz)*Sigy	+ (nx*vy-ny*vx)*Sigz );	// (n x v) (\delta m)/m * Sigma
	REAL c2PNSOn = cnvS*((105./2. - 315./2.*eta)*eta*rdot2*rdot2 + (-60.+150.*eta)*eta*rdot2*vsq + (18.-48.*eta)*eta*vsq*vsq
			+ ((-1635./2. - 117.*eta)*eta*rdot2 + (217./4.+28.*eta)*eta*vsq)*m/r
			+ (195./2. + 749./4.*eta - 8.*eta*eta)*m*m/r/r )
		    + cnvSig*((105./4. - 315./4.*eta)*eta*rdot2*rdot2 + (-30.+75.*eta)*eta*rdot2*vsq + (9.-24.*eta)*eta*vsq*vsq
			+ ((-3147./8. - 255./4.*eta)*eta*rdot2 + (131./8.+19.*eta)*eta*vsq)*m/r
			+ (111./2. + 441./4.*eta - 5.*eta*eta)*m*m/r/r );
	REAL c2PNSOv = rdot*( cnvS*((7.5 + 195./4.*eta)*eta*rdot2 + (-3./8. - 27./8.*eta - 249./8.*eta*eta)*vsq 
				+ (777./2. + 87./2.*eta)*eta*m/r) 
			    + cnvSig*((7.5 + 105./4.*eta)*eta*rdot2 + (-3./8. - 15./4.*eta - 141./8.*eta*eta)*vsq 
				+ (381./2. + 25.*eta)*eta*m/r) );
	df[I_v2PNSOx] = 1./(rsq*r)*( c2PNSOn*nx + c2PNSOv*vx 
			+ (ny*Sz-Sy*nz)*rdot*((315./8. - 945./8.*eta)*eta*rdot2*rdot2 + (-105./2. + 585./4.*eta)*eta*rdot2*vsq
				+ (-3./8. + 165./8.*eta - 441./8.*eta*eta)*vsq*vsq 
				+ ((-1215./2. - 105.*eta)*eta*rdot2 + (1067./4.+79./2.*eta)*eta*vsq)*m/r 
				+ (121./2. + 65.*eta - 8.*eta*eta)*m*m/r/r ) 
			+ (ny*Sigz + Sigy*nz)*rdot*((105./4. - 525./8.*eta)*eta*rdot2*rdot2 + (-75./2. + 345./4.*eta)*eta*rdot2*vsq
				+ (-3./8. + 57./4.*eta - 237./8.*eta*eta)*vsq*vsq 
				+ ((-2193./8. - 279./4.*eta)*eta*rdot2 + (945./8.+23.*eta)*eta*vsq)*m/r 
				+ (57./2. + 85./4.*eta - 6.*eta*eta)*m*m/r/r ) 
			+ (vy*Sz-Sy*vz)*((-225./8. + 585./8.*eta)*eta*rdot2*rdot2 + (3./8. + 255./8.*eta - 627./8.*eta*eta)*rdot2*vsq
				+ (-21./2. + 28.*eta)*eta*vsq*vsq 
				+ ((352. + 123./2.*eta)*eta*rdot2 + (197./4.+14.*eta)*eta*vsq)*m/r 
				+ (-105./2. - 137./2.*eta)*m*m/r/r ) 
			+ (vy*Sigz-Sigy*vz)*((-15. + 315./8.*eta)*eta*rdot2*rdot2 + (3./8. + 69./4.*eta - 351./8.*eta*eta)*rdot2*vsq
				+ (-11./2. + 14.*eta)*eta*vsq*vsq 
				+ ((1325./8. + 147./4.*eta)*eta*rdot2 + (-177./8.-7.*eta)*eta*vsq)*m/r 
				+ (-57./2. - 65./2.*eta)*m*m/r/r ) );
	df[I_v2PNSOy] = 1./(rsq*r)*( c2PNSOn*ny + c2PNSOv*vy 
			+ (nz*Sx-Sz*nx)*rdot*((315./8. - 945./8.*eta)*eta*rdot2*rdot2 + (-105./2. + 585./4.*eta)*eta*rdot2*vsq
				+ (-3./8. + 165./8.*eta - 441./8.*eta*eta)*vsq*vsq 
				+ ((-1215./2. - 105.*eta)*eta*rdot2 + (1067./4.+79./2.*eta)*eta*vsq)*m/r 
				+ (121./2. + 65.*eta - 8.*eta*eta)*m*m/r/r ) 
			+ (nz*Sigx + Sigz*nx)*rdot*((105./4. - 525./8.*eta)*eta*rdot2*rdot2 + (-75./2. + 345./4.*eta)*eta*rdot2*vsq
				+ (-3./8. + 57./4.*eta - 237./8.*eta*eta)*vsq*vsq 
				+ ((-2193./8. - 279./4.*eta)*eta*rdot2 + (945./8.+23.*eta)*eta*vsq)*m/r 
				+ (57./2. + 85./4.*eta - 6.*eta*eta)*m*m/r/r )
			+ (vz*Sx-Sz*vx)*((-225./8. + 585./8.*eta)*eta*rdot2*rdot2 + (3./8. + 255./8.*eta - 627./8.*eta*eta)*rdot2*vsq
				+ (-21./2. + 28.*eta)*eta*vsq*vsq 
				+ ((352. + 123./2.*eta)*eta*rdot2 + (197./4.+14.*eta)*eta*vsq)*m/r 
				+ (-105./2. - 137./2.*eta)*m*m/r/r ) 
			+ (vz*Sigx-Sigz*vx)*((-15. + 315./8.*eta)*eta*rdot2*rdot2 + (3./8. + 69./4.*eta - 351./8.*eta*eta)*rdot2*vsq
				+ (-11./2. + 14.*eta)*eta*vsq*vsq 
				+ ((1325./8. + 147./4.*eta)*eta*rdot2 + (-177./8.-7.*eta)*eta*vsq)*m/r 
				+ (-57./2. - 65./2.*eta)*m*m/r/r ) );
	df[I_v2PNSOz] = 1./(rsq*r)*( c2PNSOn*nz + c2PNSOv*vz 
			+ (nx*Sy-Sx*ny)*rdot*((315./8. - 945./8.*eta)*eta*rdot2*rdot2 + (-105./2. + 585./4.*eta)*eta*rdot2*vsq
				+ (-3./8. + 165./8.*eta - 441./8.*eta*eta)*vsq*vsq 
				+ ((-1215./2. - 105.*eta)*eta*rdot2 + (1067./4.+79./2.*eta)*eta*vsq)*m/r 
				+ (121./2. + 65.*eta - 8.*eta*eta)*m*m/r/r ) 
			+ (nx*Sigy + Sigx*ny)*rdot*((105./4. - 525./8.*eta)*eta*rdot2*rdot2 + (-75./2. + 345./4.*eta)*eta*rdot2*vsq
				+ (-3./8. + 57./4.*eta - 237./8.*eta*eta)*vsq*vsq 
				+ ((-2193./8. - 279./4.*eta)*eta*rdot2 + (945./8.+23.*eta)*eta*vsq)*m/r 
				+ (57./2. + 85./4.*eta - 6.*eta*eta)*m*m/r/r ) 
			+ (vx*Sy-Sx*vy)*((-225./8. + 585./8.*eta)*eta*rdot2*rdot2 + (3./8. + 255./8.*eta - 627./8.*eta*eta)*rdot2*vsq
				+ (-21./2. + 28.*eta)*eta*vsq*vsq 
				+ ((352. + 123./2.*eta)*eta*rdot2 + (197./4.+14.*eta)*eta*vsq)*m/r 
				+ (-105./2. - 137./2.*eta)*m*m/r/r ) 
			+ (vx*Sigy-Sigx*vy)*((-15. + 315./8.*eta)*eta*rdot2*rdot2 + (3./8. + 69./4.*eta - 351./8.*eta*eta)*rdot2*vsq
				+ (-11./2. + 14.*eta)*eta*vsq*vsq 
				+ ((1325./8. + 147./4.*eta)*eta*rdot2 + (-177./8.-7.*eta)*eta*vsq)*m/r 
				+ (-57./2. - 65./2.*eta)*m*m/r/r ) );
	df[I_r2PNSOx] = f[I_v2PNSOx];
	df[I_r2PNSOy] = f[I_v2PNSOy];
	df[I_r2PNSOz] = f[I_v2PNSOz];
    }
    df[I_orbits] = v/(2*M_PI*r);
}

inline void print(ostream& out, REAL t, const tvalarray<REAL>& f,
		  const CBwaveODE& ode, const ObserverParameters& obs,
		  const map<string, function_t>& funcmap,
		  const tvector<string>& outcolumns)
{
    out.precision(16);
    for(unsigned i = 0; i < outcolumns.size(); ++i) {
	if(i != 0) {
	    out << ' ';
	}
	const string& name = outcolumns[i];
	map<string, function_t>::const_iterator it;
	if((it = funcmap.find(name)) != funcmap.end()) {
	    out << (*it->second)(t, f, ode, obs);
	} else {
	    out << "unknown";
	}
    }
    out << endl;
}

inline bool is_valid(REAL t, const tvalarray<REAL>& f, const CBwaveODE& ode,
		     const ObserverParameters& obs)
{
    REAL r = eval_r(t, f, ode, obs);
    REAL h = eval_h(t, f, ode, obs);
    if(r > ode.rmin && r < ode.rmax
	    && h == h) { // h cannot be NaN!
	return true;
    } else {
	return false;
    }
}

inline REAL currentTime()
{
    struct timeval time_now;
    gettimeofday(&time_now, NULL);
    long long int t = (long long int)1000000*time_now.tv_sec
		+ (long long int)time_now.tv_usec;
    return 1e-6*t;
}

static void run(CBwaveODE& ode, tvalarray<REAL>& f, REAL dt,
		REAL t, REAL tmax, REAL orbitsmax,
		const ObserverParameters& obs,
		const string& datafname, ostream& out1,
		const tvector<string>& outcolumns,
		const string& ftfname, ostream& out2, 
                ostream& out3)
{
    //
    // Integrate equations of motion
    //

    REAL runtime0             = currentTime();
    aux[aux_prevchkpointtime] = runtime0;
    int doprintoutput;

    tvector  <REAL> tArray, hArray;  // arrays for to fourier transform
//  float intpart, fractpart, intpart2, intpart3;
    RK4 <REAL> rk4(ode);

//    if (&out1 != &cout) {
//	cerr << "Integrating..." << flush;
//    }
    while(is_valid(t, f, ode, obs) && t + aux[aux_dt]/2 < tmax && f[I_orbits] < orbitsmax) {

// By default there is no output
        doprintoutput = 0;

// This unlimited push_backs are eating the memory. Comment it out for long runs, uncomment
// for template bank generation, since they are needed for the FFT

//	if( &out2) {
//	    REAL h = eval_h(t, f, ode, obs);
//	    tArray.push_back(t);
//	    hArray.push_back(h);
//	}

//	fractpart = modf (f[I_orbits], &intpart);
//	if(outcolumns.size() != 0 && modf(intpart/100.0, &intpart2) == 0 && modf( float (t/(dt*1e1)) ,&intpart3) == 0) {
//	    print(out1, t, f, ode, obs, FUNCMAP, outcolumns);
//	}

// Stepping forward with the integration
	rk4.integrate(f, t, aux[aux_dt]);

// Printing some info at the end of each X orbit
/*
        if ((f[I_orbits]-aux[aux_orbitcount] >=orbitprint) && ( orbitprint > 0)) {
            doprintoutput       = 1;
            aux[aux_orbitcount] = f[I_orbits];
        }
*/
// Printing some info at the end of ecah X orbit_turn orbit

        if ((orbits_turn - aux[aux_orbitturncount] >= orbitprint) && ( orbitprint > 0)) {
            doprintoutput       = 1;
            aux[aux_orbitturncount] = orbits_turn;
        }


// Printing some info at every X step

          if (( aux[aux_stepcount] >= stepprint) && (stepprint > 0)) {
            doprintoutput      = 1;
            aux[aux_stepcount] = 0;
          }
          aux[aux_stepcount]++;

// Printing the variables
       if (doprintoutput == 1) {
          print(out1, t, f, ode, obs, FUNCMAP, outcolumns);
       }


// Make checkpoint file if required
        REAL now = currentTime();
        if (( chkpoint == "yes" ) && (now - aux[aux_prevchkpointtime] > chkpointinterval)) {
           out3 << t/SI_c;
           for (int i = 0; i < N_COMPONENTS; i++) {
              out3 << " " << ode.COMPONENT_NAMES[i] << "=" << f[i];
           }
           for (int i = 0; i < aux_N_COMPONENTS; i++) {
              out3 << " " << AUX_NAMES[i] << "=" << aux[i];
           }
           out3 << endl;
           aux[aux_prevchkpointtime] = now;
        }



// Adjusting the time step for the actual orbiting period
        if (f[I_orbits] - 1 > aux[aux_prevorbitnumber]) {
            if (adaptive.compare("yes") == 0 ) {
              aux[aux_prevorbitnumber] = f[I_orbits];
              aux[aux_dt] = (t - aux[aux_prevorbittime2]) / adaptive_step;
              aux[aux_prevorbittime2] = t;
              aux[aux_prevorbittime]  = t;
            }
        }

// Calculating the eccentricity using (rmax-rmin)/(rmax+rmin)

        erx  = eval_rx(t, f, ode, obs);
        ery  = eval_ry(t, f, ode, obs);
        erz  = eval_rz(t, f, ode, obs);
        er = sqrt(erx*erx + ery*ery + erz*erz);

	if ((rm1 > rm2) && (rm1 > er)) {
		ermax  = rm1;
		ecc_r  = (ermax-ermin)/(ermax+ermin);
		orbits_turn++;
                ermin = 2e52;
	}
        if (er < ermin) { ermin = er; } 

//	 if ((rm1 < rm2) && (rm1 < er)) {
//		ermin  = rm1;
//	 }

	rm2 = rm1;
	rm1 = er;

// Stepping forward with the time
	t += aux[aux_dt];
    }

    REAL runtime1 = currentTime();

/*
    if(&out1 != &cout) {
	char s[256];

	sprintf(s, " %.3f orbits at t=%.6fs [%.1fs] > ",
		f[I_orbits], (t-dt)/SI_c, runtime1 - runtime0);
	cerr << s << datafname << endl;

    }
*/

    int nsteps = tArray.size();
    if(nsteps != 0) {
	//
	// Resample data. Number of points must be pow of 2 for FFT.
	//
	if(&out1 != &cout) {
	    cerr << "Resampling..." << flush;
	}
	int n;
	for(n = 1; n < nsteps; n <<= 1);
	tvalarray<REAL> dft_re_h(n), dft_im_h(n);
	spline3resample((const REAL*)tArray, (const REAL*)hArray, nsteps,
			(REAL*)dft_re_h, n);
	REAL runtime2 = currentTime();

	//
	// Fourier transform resampled data
	//
	if(&out1 != &cout) {
	    char s[256];
	    sprintf(s, " [%.1fs]", runtime2 - runtime1);
	    cerr << s << " Fourier transforming..." << flush;
	}
	for(int i = 0; i < n; ++i) {
	    dft_im_h[i] = 0;
	}
	fft((REAL*)dft_re_h, (REAL*)dft_im_h, n, FFT_FORWARD);
	REAL runtime3 = currentTime();
	if(&out1 != &cout) {
	    char s[256];
	    sprintf(s, " [%.1fs]", runtime3 - runtime2);
	    cerr << s << flush;
	}
	REAL dt = tArray[nsteps - 1]/(n - 1);
	REAL T = n*dt/SI_c;
	out2 << "# f\tRe H_k\tIm H_k\t|H_k|^2" << endl;
	for(int i = 0; i < n; ++i) {
	    REAL freq = SI_c*i/(n*dt);
	    REAL re_h = dt/SI_c/sqrt(2*M_PI)*dft_re_h[i]; // [s]
	    REAL im_h = dt/SI_c/sqrt(2*M_PI)*dft_im_h[i]; // [s]
	    REAL H2 = re_h*re_h + im_h*im_h;
	    out2 << freq << '\t' << re_h << '\t' << im_h
		 << '\t' << sqrt(H2/T) << endl;
	}
	cerr << " > " << ftfname << endl;
        out2 << EOF;
    }

}

static ostream* make_ostream(const string& fname)
{
    if(fname == "") {
	return (ostream*)0;
    } else if(fname == "-") {
	return &cout;
    } else {
	return new ofstream(fname.c_str());
    }
}


void set_checkpoint_value(string inputstring) {
   vector <string> keyvaluepair;
   boost::split(keyvaluepair,inputstring,boost::is_any_of("="));

   for (int i = 0; i < N_COMPONENTS; i++) {
    if (keyvaluepair[0] == ode.COMPONENT_NAMES[i]) {
       f[i] = atof(keyvaluepair[1].c_str());
    }
   }

   for (int i = 0; i < aux_N_COMPONENTS; i++) {
    if (keyvaluepair[0] == AUX_NAMES[i]) {
       aux[i] = atof(keyvaluepair[1].c_str());
    }
   }
   return;
}

static int parseArgs(int argc, const char* argv[],
		     CBwaveODE& ode, ObserverParameters& obs,
		     tvector<string>& outcolumns, tvalarray<REAL>& f,
		     REAL& tmax, REAL& orbitsmax,
		     REAL& dt,  REAL& freq, REAL& T,
		     REAL& epsilon, REAL& rr,
		     string& datafname, string& ftfname, string& checkpoint,
                     int * stepprint, int * orbitprint, int *loglevel, 
             string * adaptive, int * adaptive_step)
{

// Reading of the ini file variables

   if (argc < 2 ) {
     cerr << " No ini file defined! Exiting.\n" << endl;
     return -1;
   }

   // The checkpoint file is always called like the ....ini.chk
   chkpointfile = argv[1];
   chkpointfile.append(".chk");

   ConfigFile config(argv[1]);
   REAL m1   		= config.read<REAL>("m1");
   REAL m2   		= config.read<REAL>("m2");
   REAL s1x		= config.read<REAL>("s1x");
   REAL s1y		= config.read<REAL>("s1y");
   REAL s1z		= config.read<REAL>("s1z");
   REAL s2x		= config.read<REAL>("s2x");
   REAL s2y		= config.read<REAL>("s2y");
   REAL s2z		= config.read<REAL>("s2z");
   tmax      		= config.read<REAL>("tmax")*SI_c;
   orbitsmax  		= config.read<REAL>("orbitsmax");
   dt 	 		= config.read<REAL>("dt")*SI_c;
   ode.rmin   		= config.read<REAL>("rmin");
   ode.rmax   		= config.read<REAL>("rmax");
   T	   		= config.read<REAL>("T")*SI_c;
   freq 		= config.read<REAL>("f");
   epsilon 		= config.read<REAL>("epsilon");
   rr	 		= config.read<REAL>("r");
   obs.D 		= config.read<REAL>("D");
   obs.iota 		= config.read<REAL>("iota")*M_PI/180.;
   obs.phi 		= config.read<REAL>("phi")*M_PI/180.;
   obs.theta 		= config.read<REAL>("theta")*M_PI/180.;
   obs.varphi 		= config.read<REAL>("varphi")*M_PI/180.;
   obs.psi 		= config.read<REAL>("psi")*M_PI/180.;
   datafname 		= config.read<string>("outfile");
   ftfname		= config.read<string>("ftfile");
   description          = config.read<string>("description");
   string hterms	= config.read<string>("hterms");
   string corrections	= config.read<string>("corrs");
   string outvars	= config.read<string>("outvars");
   checkpoint		= config.read<string>("checkpoint");
   *orbitprint		= config.read<int>("printorbit");
   *stepprint		= config.read<int>("printstep");
   *loglevel		= config.read<int>("loglevel");
   *adaptive		= config.read<string>("adaptive");
   *adaptive_step	= config.read<int>("adaptive_step");
   ode.alpha 		= config.read<REAL>("alpha");
   ode.beta 		= config.read<REAL>("beta");
   ode.delta1 		= config.read<REAL>("delta1");
   ode.delta2 		= config.read<REAL>("delta2");
   ode.delta3 		= config.read<REAL>("delta3");
   ode.delta4 		= config.read<REAL>("delta4");
   ode.delta5 		= config.read<REAL>("delta5");
   ode.delta6 		= config.read<REAL>("delta6");


// Setting the initial values for the masses and spins

   ode.set("m1", m1, f);
   ode.set("m2", m2, f);
   ode.set("s1x", s1x, f);   
   ode.set("s1y", s1y, f);   
   ode.set("s1z", s1z, f);   
   ode.set("s2x", s2x, f);   
   ode.set("s2y", s2y, f);   
   ode.set("s2z", s2z, f);   

// Converting the 'hterms' string
   ode.hterms = 0;
   for (int i = 0; i < H_COUNT; i++) {
    string mysubstr = "'"; mysubstr.append(H_NAMES[i]); mysubstr+="'";
    if ( hterms.find(mysubstr, 0 ) != std::string::npos) {
     ode.hterms |= 1 << i;
    }
   }

// Converting the 'corrections' string
   ode.corrections = 0;
   for (int i = 0; i < C_COUNT; i++) {
    string mysubstr = "'"; mysubstr.append(C_NAMES[i]); mysubstr+="'";
    if ( corrections.find(mysubstr, 0 ) != std::string::npos) {
     ode.corrections |= 1 << i;
    }
   }

// Converting the outcolumns string
   string::size_type i = 0;
   string::size_type j = 0;
   outcolumns.resize(0);
   while((j = outvars.find(',', i + 1)) != std::string::npos) {
          outcolumns.push_back(outvars.substr(i, j - i));
          i = j + 1;
   }
   outcolumns.push_back(outvars.substr(i));

// Now, if checkpoint was required reading the checkpoint file.
   
   if (argc == 4) {
    fstream   hfile;
    string    line,field;
    vector    <string> fields;
    int       chkpointfound = 0;
    runfrom  = atof(argv[3]);
    histfile = argv[2];
    hfile.open (histfile.c_str(), fstream::in);
    while ( ! hfile.eof() ) {
      getline (hfile,line);
      boost::split(fields,line,boost::is_any_of(" "));
      if (atof(fields[0].c_str()) == runfrom) {
        t = runfrom*SI_c;
        chkpointfound = 1;
        for_each(fields.begin(),fields.end(), set_checkpoint_value); 
        hfile.close();
        return 0;
      }
    }
    if (chkpointfound == 0) {
      cerr << "The specified time hasn't been found in the checkpoint file. Exiting." << endl;
      hfile.close();
      exit(-1);
    }
   } else 
   if (argc != 2) {
      cerr << "Usage: " << argv[0] << " settings.ini [chekpointfile.chk] [timestamp]" << endl;
      exit(-1);
   }

// This line executes only when no previous checkpoint file is defined.
  return 0;
}

/******************************************************************************/
/*
 * Here is starts the main loop of CBwaves
 *
 */
/******************************************************************************/

int main(int argc, const char* argv[])
{

   string oldTrace;
   cbwInitLogMessages();
   cbwTrace = "main";


// Auxiliary variables
    int  ret;

//  Defaults values for variables

    REAL freq      = 60.0;      // [Hz]
    REAL T         = SI_c/freq; // [m] orbit time in meters
    REAL dt        = 0.;        // time step
    REAL epsilon   = 99.;       // excentricity, cannot be specified it is calculated
    REAL tmax      = 30.0*T;    // maximum evolution time
    REAL orbitsmax = 30.;       // maximum number of orbits
    REAL rr        = 253936.;   // initial separation of the two body

    REAL boostvalue = 1.0 ;     // boost value for f[I_vNy] for correct eccentricity
    REAL boostscale = 0.006;    // scaling for succ. approx. determination of boostvalue
    REAL initialspeed;          // value of f[I_vNy] without any boost
    REAL boostflag = -1.;       // aux var for succ. approx.

    ostream * devnull = make_ostream("/dev/null");
    ostream * none = NULL;

// Reading the initial parameters to determine the boostscale. The higher the
// required eccentricity is the finer the stepping set.

    ret = parseArgs(argc, argv, ode, obs, outcolumns, f, tmax, orbitsmax, dt,  freq, T,
          epsilon, rr, datafname, ftfname, chkpoint, &stepprint, &orbitprint, &cbwLogLevel,
          &adaptive, &adaptive_step);

    cbwLog(cbwINFO, "CBwaves main started. ");
    cbwLog(cbwINFO, "Configuration file has been read. ");

// Scaling for the initial velocity
    REAL m = ode.m;
    REAL eta = ode.eta;
    REAL e = epsilon;
    REAL BBB = (1.0-3.65e4/rr*(e-2.62e-1*e*e));

    REAL m1 = ode.m1;
    REAL m2 = ode.m2;
    REAL delta = m1-m2;
    REAL cSs1 = m1*m1;
    REAL S1x = cSs1*f[I_s1x];
    REAL S1y = cSs1*f[I_s1y];
    REAL S1z = cSs1*f[I_s1z];
    REAL cSs2 = m2*m2;
    REAL S2x = cSs2*f[I_s2x];
    REAL S2y = cSs2*f[I_s2y];
    REAL S2z = cSs2*f[I_s2z];


    complex<double> A(((1. + e)*m*(3. - 4.*eta)*eta)/rr,0);
    complex<double> B((-2.*pow(rr,2) + 2.*(1 + e)*m*rr*(1 + 3.*eta) +
                  (1. + e)*pow(m,2)*eta*(-13. + 4.*eta))/(2.*pow(rr,2)),0);
    complex<double>  F(-(((1. + e)*(-3.*m2*S1z*delta +
                  m1*(5.*m2*(S1z + S2z) + 3.*S2z*delta)))/(m1*m2*pow(rr,2))),0);
    complex<double> G(((1. + e)*(-24.*S1x*S2x + 12.*S1y*S2y + 12*S1z*S2z +
                       36.*pow(m,4)*eta - 16.*pow(m,3)*rr*eta +
                       4.*pow(m,2)*pow(rr,2)*eta + 87.*pow(m,4)*pow(eta,2) -
                       8.*pow(m,3)*rr*pow(eta,2)))/(4.*m*pow(rr,3)*eta),0);

    initialspeed = real(sqrt((-2.*B)/(3.*A) +
         (pow(2,0.3333333333333333)*(pow(B,2) + 12.*A*G))/
          (3.*A*pow(2.*pow(B,3) + 27.*A*pow(F,2) - 72.*A*B*G +
            sqrt(-4.*pow(pow(B,2) + 12.*A*G,3) +
             pow(2.*pow(B,3) + 27.*A*pow(F,2) - 72.*A*B*G,2)),
          0.3333333333333333)) +
         pow(2.*pow(B,3) + 27.*A*pow(F,2) - 72.*A*B*G +
          sqrt(-4.*pow(pow(B,2) + 12.*A*G,3) +
            pow(2.*pow(B,3) + 27.*A*pow(F,2) - 72.*A*B*G,2)),
         0.3333333333333333)/(3.*pow(2,0.3333333333333333)*A))/2. -
        sqrt((-4.*B)/(3.*A) - (pow(2,0.3333333333333333)*
          (pow(B,2) + 12.*A*G))/
        (3.*A*pow(2.*pow(B,3) + 27.*A*pow(F,2) - 72.*A*B*G +
            sqrt(-4.*pow(pow(B,2) + 12.*A*G,3) +
              pow(2.*pow(B,3) + 27.*A*pow(F,2) - 72.*A*B*G,2)),
           0.3333333333333333)) -
         pow(2.*pow(B,3) + 27.*A*pow(F,2) - 72.*A*B*G +
          sqrt(-4.*pow(pow(B,2) + 12.*A*G,3) +
            pow(2.*pow(B,3) + 27.*A*pow(F,2) - 72.*A*B*G,2)),
         0.3333333333333333)/(3.*pow(2,0.3333333333333333)*A) -
         (2.*F)/(A*sqrt((-2.*B)/(3.*A) +
            (pow(2,0.3333333333333333)*(pow(B,2) + 12.*A*G))/
             (3.*A*pow(2.*pow(B,3.) + 27.*A*pow(F,2) - 72.*A*B*G +
                 sqrt(-4.*pow(pow(B,2) + 12.*A*G,3) +
                   pow(2.*pow(B,3) + 27.*A*pow(F,2) - 72.*A*B*G,2)),
                0.3333333333333333)) +
            pow(2.*pow(B,3) + 27.*A*pow(F,2) - 72.*A*B*G +
               sqrt(-4.*pow(pow(B,2) + 12.*A*G,3.) +
                 pow(2.*B*B*B + 27.*A*F*F - 72.*A*B*G,2)),
              0.3333333333333333)/(3.*pow(2.,0.3333333333333333)*A))))/2.);


    cbwLog(cbwINFO, "Determining initial conditions for required eccentricity.");

    REAL ecc_prev  = 10.;
    REAL accuracy  = 1.0e-6;

    do {

        oldTrace = cbwTrace;
        cbwTrace = "eccloop";

//  Reading program parameters from .ini file
        ret = parseArgs(argc, argv, ode, obs, outcolumns, f, tmax, orbitsmax, dt,  freq, T,
             epsilon, rr, datafname, ftfname, chkpoint, &stepprint, &orbitprint, &cbwLogLevel,
             &adaptive, &adaptive_step);

        if (ret != 0) { return ret - 1;}

// Not running from a checkpoint
        if (runfrom == 0) {

        // Initial separation and velocity at periastron, Kepler orbit
        // BBB factor needs scaling for eccentricity

// Run 1.2 orbits
        f[I_rNx] = rr;
        f[I_vNy] = initialspeed * BBB * boostvalue;
        aux[aux_dt] = dt;

        auxoutcolumns.resize(1);
        auxoutcolumns[0]="ecc_r";

        ermax          = 0;         // initial values for eccentricity calculations
        ermin          = 2.e52;      //
        erorbit        = 0;
        ecc_r          = 10;

        run(ode, f, aux[aux_dt], t, tmax, 1.2, obs, datafname, *devnull, auxoutcolumns,
        ftfname, *none, *devnull);

        cbwLog(cbwDEBUG,"Boost = %g, scaler  = %g, ecc = %g", boostvalue,boostscale,ecc_r);

        if (accuracy < abs(ecc_r - epsilon)) {
         if ( abs(ecc_prev - epsilon) < abs(ecc_r - epsilon)) {
            boostflag  *=-1.0;
            boostscale /= 1.5;
         }
         ecc_prev = ecc_r;
         boostvalue *= 1.0 + boostflag*boostscale;
        }
      } // end of if runfrom == 0

      cbwTrace = oldTrace;
      // Cleaning the state for a new simulation
      for (int i = 0; i < N_COMPONENTS; i++ ) {
        f[i] = 0;
      }

      //   Zeroing the aux array
      for (int i = 0; i < aux_N_COMPONENTS; i++) {
        aux[i] = 0;
      }

    } while (accuracy < abs(ecc_r-epsilon));	// end of do-while loop

    cbwLog(cbwINFO, "Initial parameters determind.");
    // Ok, now we have the right value for the eccentricity let's start the real simulation !


    //  Reading program parameters from .ini file
    ret = parseArgs(argc, argv, ode, obs, outcolumns, f, tmax, orbitsmax, dt,  freq, T,
         epsilon, rr, datafname, ftfname, chkpoint, &stepprint, &orbitprint, &cbwLogLevel,
         &adaptive, &adaptive_step);
    if (ret != 0) { return ret - 1;}

    //  Setting the time step
    aux[aux_dt] = dt;

    // Zeroing the time before the real simulation
    t = 0;

    // Setting initial conditions
    f[I_rNx] = rr;
    f[I_vNy] = initialspeed*BBB*boostvalue;

/*
  // Initial conditions for hyperbolic orbits
    f[I_rNx] = rr;
    f[I_rNy] = 40*(ode.m1 + ode.m2);
    f[I_vNx] = - initialspeed*BBB*boostvalue - 0.1;
    f[I_vNy] = 0.0;

    REAL dist = sqrt(f[I_rNx]*f[I_rNx] + f[I_rNy]*f[I_rNy])/(ode.m1 + ode.m2);
    REAL ener = ode.mu*(0.5*f[I_vNx]*f[I_vNx] - ode.m/dist);

    cout << "E=" << ener << " d=" << dist << "M" << endl;
*/

    // Create output files

    ostream* out1 = make_ostream(datafname);
    ostream* out2 = make_ostream(ftfname);
    ostream* out3 = NULL;

    // Creat checkpoint file only if required
    if (chkpoint == "yes" ) {
       ostream* out3 = make_ostream(chkpointfile);
 	// Set the precision for the checkpoint file stream
	out3->precision(20);
    }

    // Reseting the eccentricity value for the clean run
    ermax          = 0;
    ermin          = 2e52;
    erorbit        = 0;
    ecc_r          = 10;
    orbits_turn    = 0;

    cbwLog(cbwINFO, "Running the simulation");

    // Run the simulation
    run(ode, f, aux[aux_dt], t, tmax, orbitsmax, obs, datafname, *out1, outcolumns,
        ftfname, *out2, *out3);


    cbwLog(cbwINFO, "Simulation finished, cleaning up.");

    // Free memory


    if(out1 && out1 != &cout) {
	delete out1;
    }
    if(out2 && out2 != &cout) {
	delete out2;
    }

    if(out3 && out3 != &cout) {
	delete out3;
    }

    cbwLog(cbwINFO, "Exiting.");
    return 0;

}
