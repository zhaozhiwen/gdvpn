//
//  gdvpn.cc
//
//
//  by Adam Freese, April 2012 - May 2017
//
//  This library computes the differential cross section of the reaction
//    d + gamma --> V + p + n
// for a vector meson V. Support for V=phi(1020) and V=J/Psi is present.
//
//  This program uses Cuba numerical integration routines, cf. 
//    <http://www.feynarts.de/cuba/>
// Thus -lcuba should be included as a linker flag during build.
//
//
// Standard libraries
#include <stdlib.h>
#include <unistd.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <complex>
#include <cuba.h>
#include <time.h>
// Local libraries
#include "dwf.hh"
#include "kine2.hh"
#include "nnscatter.hh"
#include "vmeson_params.hh"
#include "constants.hh"
using namespace std;

// =============================================================================
//  Kinematical Data Struct

struct kinedata {
  FourVector deuteron, gamma, vmeson, proton, neutron;
  int sd, sp, sn;
  bool pwia, singles, doubles, cauchy, vnvertex, pnvertex, soft, dsoft;
  kinedata() {
    // By default consider PWIA and pole parts of hard single rescatterings
    pwia     = true;
    singles  = true;
    doubles  = false;
    cauchy   = false;
    vnvertex = true;
    pnvertex = true;
    soft     = false;
    dsoft    = false;
    // Default initial kinematics, 8.2 GeV photon and deuteron at rest
    gamma.set(8.2,0,0,8.2);
    deuteron.set(DeuteronMass,0,0,0);
  }
};

//  A global instance that can only be seen/used internally
static kinedata gdvpn;

// =============================================================================
//  Forward declarations

// The differential cross section
double diffsig5(double,double,double,double,double);

//  Mutator and accessor for the private kinedata struct
void init_gdvpn(double,bool,bool,bool,bool,bool,bool,bool,bool);
kinedata get_gdvpn_info();

//  Feynman amplitude
static complex<double> Feynman();
//  Integrand functions for Feynman amplitude
static int ampl_intd_vn(const int*,const double[],const int*,double[]);
static int ampl_intd_pn(const int*,const double[],const int*,double[]);
static int ampl_intd_soft(const int*,const double[],const int*,double[]);
static int ampl_intd_vnpn(const int*,const double[],const int*,double[]);
static int ampl_intd_pnvn(const int*,const double[],const int*,double[]);
static int ampl_intd_dsoft(const int*,const double[],const int*,double[]);
static int ampl_intd_vnc(const int*,const double[],const int*,double[]);
static int ampl_intd_pnc(const int*,const double[],const int*,double[]);
//  Extra functions for Cauchy principal value part
static complex<double> cauchy_f_vn(double,FourVector);
static complex<double> cauchy_f_pn(double,FourVector);
static complex<double> cauchy_h_vn(double,double,FourVector);
static complex<double> cauchy_h_pn(double,double,FourVector);
//  Kinematics and phase factor
static int gdvpn_kinematics(double[],double,int);
static double phase_factor();
//  Other math routines
static void integrate(const int,const int,integrand_t,void*,double[]);
static void solvep(double,double,double,double,double&,double&,int&);
static void solvequad(double,double,double,double&,double&,int&);
//  Deuteron wave-function shortcut
static complex<double>dwf(ThreeVector);

// =============================================================================
//  Differential cross-section

double diffsig5(double t, double p_n, 
    double theta_nl, double phi_nl, double phi_v) {
  complex<double> M;
  double diffsig = 0.;
  double dsigargs[5];
  int solno, nosol;
  bool solved = false;
  //  Arguments for kinematical routine
  dsigargs[0] = t;
  dsigargs[1] = p_n;
  dsigargs[2] = theta_nl;
  dsigargs[3] = phi_nl;
  dsigargs[4] = phi_v;
  for(solno=1;solno<=2;solno++) {
    if(gdvpn_kinematics(dsigargs,vmass(),solno)==0) {
      for(gdvpn.sd=1;gdvpn.sd>=-1;gdvpn.sd-=1)
        for(gdvpn.sp=1;gdvpn.sp>=-1;gdvpn.sp-=2)
          for(gdvpn.sn=1;gdvpn.sn>=-1;gdvpn.sn-=2)
            diffsig += 0.5/kine2flux(gdvpn.deuteron,gdvpn.gamma)
              * phase_factor()*norm(Feynman());
      solved = true;
    }
  }
  diffsig *= nb_mb*mb_GeV;
  diffsig /= 3.;
  if(solved) return diffsig;
  else return -1.;
}

// =============================================================================
//  Mutator and accessor for private kinedata struct

void init_gdvpn(double E, bool pwia, bool singles, bool doubles, bool vnvertex, 
    bool pnvertex, bool soft, bool dsoft, bool cauchy) {
  gdvpn.gamma.set(E,0,0,E);
  gdvpn.deuteron.set(DeuteronMass,0,0,0);
  gdvpn.pwia     = pwia;
  gdvpn.singles  = singles;
  gdvpn.doubles  = doubles;
  gdvpn.vnvertex = vnvertex;
  gdvpn.pnvertex = pnvertex;
  gdvpn.soft     = soft;
  gdvpn.dsoft    = dsoft;
  gdvpn.cauchy   = cauchy;
  return;
}

kinedata get_gdvpn_info() {
  return gdvpn;
}

// =============================================================================
//  Integration Routine

// This is a wrapper for Vegas from the Cuba integration library
static void integrate(const int ndim, const int ncomp, integrand_t integrand, 
    void *userdata, double integral[]) {
  int neval, nregions, fail;
  double error[2], prob[2];
  Vegas(ndim,ncomp,integrand,userdata,
      1,1e-5,1e-15,0,137,0,10000,1000,500,1000,0,
      NULL,NULL,&neval,&fail,integral,error,prob);
  return;
}

// =============================================================================
//  Feynman Amplitude

static complex<double> Feynman() {
  FourVector Pp_prime;
  complex<double> ampl, wf;
  double s, t;
  double amplparts[2];
  ampl = 0;
  set_dwf(1); // initialize deuteron w/f to AV18
  if(gdvpn.pnvertex) init_nnscatter(); // need to initialize NN amplitudes
  if(gdvpn.pwia) {
    s = kine2s(gdvpn.vmeson,gdvpn.proton);
    t = kine2t(gdvpn.gamma,gdvpn.vmeson);
    Pp_prime = gdvpn.deuteron - gdvpn.neutron;
    wf = dwf(gdvpn.neutron.p);
    wf /= sqrt(Pp_prime.E/ProtonMass);
    ampl = 1.0;
    ampl *= sqrt(2*pow(2*M_PI,3)*2*gdvpn.neutron.E);
    ampl *= kine2flux(Pp_prime,gdvpn.gamma)*vnphoto(s,t);
    ampl *= wf;
    if(Pp_prime.E<=0||Pp_prime.m2()<=0) ampl = 0;
  }
  if(gdvpn.singles) {
    if(gdvpn.vnvertex) {
      integrate(2,2,(integrand_t)ampl_intd_vn,&gdvpn,amplparts);
      ampl += complex<double>(amplparts[0],amplparts[1]);
    }
    if(gdvpn.pnvertex) {
      integrate(2,2,(integrand_t)ampl_intd_pn,&gdvpn,amplparts);
      ampl += complex<double>(amplparts[0],amplparts[1]);
    }
    if(gdvpn.soft) {
      integrate(2,2,(integrand_t)ampl_intd_soft,&gdvpn,amplparts);
      ampl += complex<double>(amplparts[0],amplparts[1]);
    }
  }
  if(gdvpn.doubles) {
    integrate(4,2,(integrand_t)ampl_intd_vnpn,&gdvpn,amplparts);
    ampl += complex<double>(amplparts[0],amplparts[1]);
    integrate(4,2,(integrand_t)ampl_intd_pnvn,&gdvpn,amplparts);
    ampl += complex<double>(amplparts[0],amplparts[1]);
    if(gdvpn.dsoft) {
      integrate(4,2,(integrand_t)ampl_intd_dsoft,&gdvpn,amplparts);
      ampl += complex<double>(amplparts[0],amplparts[1]);
    }
  }
  if(gdvpn.cauchy) {
    if(gdvpn.vnvertex) {
      integrate(3,2,(integrand_t)ampl_intd_vnc,&gdvpn,amplparts);
      ampl += complex<double>(amplparts[0],amplparts[1]);
    }
    if(gdvpn.pnvertex) {
      integrate(3,2,(integrand_t)ampl_intd_pnc,&gdvpn,amplparts);
      ampl += complex<double>(amplparts[0],amplparts[1]);
    }
  }
  return ampl;
}

// =============================================================================
//  Integrand Functions for Feynman Amplitude

//  VN rescattering (pole part)
static int ampl_intd_vn(const int *ndim, const double xx[], const int *ncomp,
    double ff[]) {
  ff[0] = 0;
  ff[1] = 0;
  //  Scaled variables.
  double k = xx[0];
  double phi = 2*M_PI*xx[1];
  //  Intermediate quantities
  FourVector K, Pv_prime, Pp_prime, Pn_prime;
  double Delta, s_gn, s_vn, t_gn, t_vn;
  complex<double> ampl_intd, wf;
  //  Preliminary values for K and Pn_prime
  K.p.set(k*cos(phi),k*sin(phi),0);
  Pn_prime.p.set(gdvpn.neutron.p.x-K.p.x,gdvpn.neutron.p.y-K.p.y,0);
  Pn_prime.build(pow(ProtonMass,2),Pn_prime.p);
  K.E = gdvpn.neutron.E - Pn_prime.E;
  //  Delta & final for all four-momenta
  Delta = K.m2()/2 + K.E*gdvpn.vmeson.E - K.p*gdvpn.vmeson.p;
  Delta /= gdvpn.vmeson.p.z;
  K.p.set(K.p.x,K.p.y,Delta);
  Pn_prime.build(pow(ProtonMass,2),gdvpn.neutron.p-K.p);
  Pv_prime.build(vmass2(),gdvpn.vmeson.p+K.p);
  Pp_prime = gdvpn.deuteron-Pn_prime;
  if(Pp_prime.E<=0||Pp_prime.m2()<=0) return 0;
  //  Mandelstam variables
  s_gn = kine2s(Pv_prime,gdvpn.proton);
  s_vn = kine2s(gdvpn.vmeson,gdvpn.neutron);
  t_gn = kine2t(gdvpn.gamma,Pv_prime);
  t_vn = kine2t(Pv_prime,gdvpn.vmeson);
  //  Use w/f model
  wf = dwf(Pn_prime.p);
  wf /= sqrt(Pp_prime.E/ProtonMass);
  //  ampl_intd is unscaled integrand
  ampl_intd = complex<double>(0,0.5);
  ampl_intd /= 2*gdvpn.vmeson.p.z;
  ampl_intd *= sqrt(2*pow(2*M_PI,3)/(2*Pn_prime.E));
  ampl_intd *= kine2flux(s_gn,0,Pp_prime.m2())*vnphoto(s_gn,t_gn);
  ampl_intd *= kine2flux(s_vn,vmass2(),pow(ProtonMass,2))*vnscatter(s_vn,t_vn);
  ampl_intd *= wf;
  ampl_intd /= pow(2*M_PI,2);
  //  Scale integrand by appropriate Jacobian.
  ampl_intd *= k*(2*M_PI);
  ff[0] = real(ampl_intd);
  ff[1] = imag(ampl_intd);
  return 0;
}

//  pn rescattering (pole part)
static int ampl_intd_pn(const int *ndim, const double xx[], const int *ncomp,
    double ff[]) {
  ff[0] = 0;
  ff[1] = 0;
  //  Scaled variables.
  double k = xx[0];
  double phi = 2*M_PI*xx[1];
  //  Intermediate quantities
  FourVector K, Pp_first, Pp_middle, Pn_prime;
  double Delta, s_gn, s_pn, t_gn, t_pn;
  complex<double> ampl_intd, wf;
  //  Preliminary values for K and Pn_prime
  K.p.set(k*cos(phi),k*sin(phi),0);
  Pn_prime.p.set(gdvpn.neutron.p.x-K.p.x,gdvpn.neutron.p.y-K.p.y,0);
  Pn_prime.build(pow(ProtonMass,2),Pn_prime.p);
  K.E = gdvpn.neutron.E - Pn_prime.E;
  //  Delta & final values for all four-momenta
  Delta = K.m2()/2 + gdvpn.proton.E*K.E - K.p*gdvpn.proton.p;
  Delta /= gdvpn.proton.p.z;
  if(abs(Delta)>=1) return 0; //temporary
  K.p.set(K.p.x,K.p.y,Delta);
  Pn_prime.build(pow(ProtonMass,2),gdvpn.neutron.p-K.p);
  Pp_middle.build(pow(ProtonMass,2),gdvpn.proton.p+K.p);
  Pp_first.set(DeuteronMass-Pn_prime.E,
      Pp_middle.p+gdvpn.vmeson.p-gdvpn.gamma.p);
  if(Pp_first.E<0||Pp_first.m2()<=0) return 0;
  //  Mandelstam variables
  s_gn = kine2s(Pp_middle,gdvpn.vmeson);
  s_pn = kine2s(gdvpn.proton,gdvpn.neutron);
  t_gn = kine2t(gdvpn.gamma,gdvpn.vmeson);
  t_pn = kine2t(gdvpn.neutron,Pn_prime);
  //  Use w/f model
  wf = dwf(Pn_prime.p);
  wf /= sqrt(Pp_first.E/ProtonMass);
  //  ampl_intd is the unscaled integrand
  ampl_intd = complex<double>(0,0.5);
  ampl_intd /= 2*gdvpn.proton.p.z;
  ampl_intd *= sqrt(2*pow(2*M_PI,3)/(2*Pn_prime.E));
  ampl_intd *= kine2flux(s_gn,0,pow(ProtonMass,2))*vnphoto(s_gn,t_gn);
  ampl_intd *= kine2flux(s_pn,pow(ProtonMass,2),pow(ProtonMass,2))
    * pnscatter(s_pn,t_pn);
  ampl_intd *= wf;
  ampl_intd /= pow(2*M_PI,2);
  //  Scale integrand by appropriate Jacobian.
  ampl_intd *= k*(2*M_PI);
  ff[0] = real(ampl_intd);
  ff[1] = imag(ampl_intd);
  return 0;
}

//  Soft production (pole part)
static int ampl_intd_soft(const int *ndim, const double xx[], const int *ncomp,
    double ff[]) {
  ff[0] = 0;
  ff[1] = 0;
  //  Scaled variables.
  double k = xx[0];
  double phi = 2*M_PI*xx[1];
  //  Intermediate quantities
  FourVector K, Pv_prime, Pn_prime, Pp_prime;
  double Delta, s_gn, s_vp, t_gn, t_vp;
  complex<double> ampl_intd, wf;
  //  Preliminary values for K and Pn_prime
  K.p.set(k*cos(phi),k*sin(phi),0);
  Pn_prime.p.set(gdvpn.neutron.p.x-K.p.x,gdvpn.neutron.p.y-K.p.y,0);
  Pn_prime.build(pow(ProtonMass,2),Pn_prime.p);
  K.E = gdvpn.neutron.E - Pn_prime.E;
  //  Delta & final values for all four-momenta
  Delta = -K.m2()/2 + K.E*gdvpn.gamma.E + vmass2()/2;
  Delta /= gdvpn.gamma.p.z;
  K.p.set(K.p.x,K.p.y,Delta);
  Pv_prime.build(vmass2(),gdvpn.gamma.p-K.p);
  Pp_prime.build(pow(ProtonMass,2),gdvpn.vmeson.p+gdvpn.proton.p-Pv_prime.p);
  Pn_prime = gdvpn.deuteron-Pp_prime;
  if(Pn_prime.E<=0||Pn_prime.m2()<=0) return 0;
  //  Mandelstam variables
  s_gn = kine2s(Pv_prime,gdvpn.neutron);
  s_vp = kine2s(gdvpn.vmeson,gdvpn.proton);
  t_gn = kine2t(gdvpn.gamma,Pv_prime);
  t_vp = kine2t(Pv_prime,gdvpn.vmeson);
  //  Use w/f model
  wf = dwf(Pp_prime.p);
  wf /= sqrt(Pn_prime.E/ProtonMass);
  wf *= -1; // isospin of struck neutron
  //  ampl_intd is unscaled integrand
  ampl_intd = complex<double>(0,0.5);
  ampl_intd /= 2*gdvpn.gamma.p.z;
  ampl_intd *= sqrt(2*pow(2*M_PI,3)/(2*Pp_prime.E));
  ampl_intd *= kine2flux(s_gn,0,Pn_prime.m2())*vnphoto(s_gn,t_gn);
  ampl_intd *= kine2flux(s_vp,vmass2(),pow(ProtonMass,2))*vnscatter(s_vp,t_vp);
  ampl_intd *= wf;
  ampl_intd /= pow(2*M_PI,2);
  //  Scale integrand by appropriate Jacobian.
  ampl_intd *= k*(2*M_PI);
  ff[0] = real(ampl_intd);
  ff[1] = imag(ampl_intd);
  return 0;
}

//  Double rescattering, VN followed by pn rescattering
static int ampl_intd_vnpn(const int *ndim, const double xx[], const int *ncomp,
    double ff[]) {
  ff[0] = 0;
  ff[1] = 0;
  //  Scaled variables.
  double k_vn = xx[0];
  double k_pn = xx[1];
  double phi_vn = 2*M_PI*xx[2];
  double phi_pn = 2*M_PI*xx[3];
  //  Intermediate quantities
  FourVector K_vn, K_pn, Pp_first, Pp_middle, Pn_first, Pn_middle, Pv_prime;
  double Delta_vn, Delta_pn, s_gn, s_vn, s_pn, t_gn, t_vn, t_pn;
  complex<double> ampl_intd, wf;
  //  Preliminary values for K's and Pn's
  K_vn.p.set(k_vn*cos(phi_vn),k_vn*sin(phi_vn),0);
  K_pn.p.set(k_pn*cos(phi_pn),k_pn*sin(phi_pn),0);
  Pn_middle.p.set(gdvpn.neutron.p.x-K_pn.p.x,gdvpn.neutron.p.y-K_pn.p.y,0);
  Pn_first.p.set(Pn_middle.p.x-K_vn.p.x,Pn_middle.p.y-K_vn.p.y,0);
  Pn_middle.build(pow(ProtonMass,2),Pn_middle.p);
  Pn_first.build(pow(ProtonMass,2),Pn_first.p);
  K_vn.E = Pn_middle.E-Pn_first.E;
  K_pn.E = gdvpn.neutron.E-Pn_middle.E;
  //  Deltas & final values for all four-momenta
  Delta_vn = K_vn.m2()/2 + gdvpn.vmeson.E*K_vn.E - gdvpn.vmeson.p*K_vn.p;
  Delta_vn /= gdvpn.vmeson.p.z;
  Delta_pn = K_pn.m2()/2 + gdvpn.proton.E*K_pn.E - gdvpn.proton.p*K_pn.p;
  Delta_pn /= gdvpn.proton.p.z;
  K_vn.p.set(K_vn.p.x,K_vn.p.y,Delta_vn);
  K_pn.p.set(K_pn.p.x,K_pn.p.y,Delta_pn);
  Pn_middle.build(pow(ProtonMass,2),gdvpn.neutron.p-K_pn.p);
  Pn_first.build(pow(ProtonMass,2),Pn_middle.p-K_vn.p);
  Pp_middle.build(pow(ProtonMass,2),gdvpn.proton.p+K_pn.p);
  Pv_prime.build(vmass2(),gdvpn.vmeson.p+K_vn.p);
  Pp_first = gdvpn.deuteron - Pn_first;
  if(Pp_first.E<=0||Pp_first.m2()<=0) return 0;
  //  Mandelstam variables
  s_gn = kine2s(gdvpn.gamma,Pp_first);
  s_vn = kine2s(gdvpn.vmeson,Pn_middle);
  s_pn = kine2s(gdvpn.proton,gdvpn.neutron);
  t_gn = kine2t(gdvpn.gamma,Pv_prime);
  t_vn = kine2t(gdvpn.vmeson,Pv_prime);
  t_pn = kine2t(gdvpn.proton,Pp_middle);
  //  Use w/f model
  wf = dwf(Pn_first.p);
  wf /= sqrt(Pp_first.E/ProtonMass);
  //  ampl_intd is the unscaled integrand
  ampl_intd = -0.25;
  ampl_intd /= (2*gdvpn.vmeson.p.z)*(2*gdvpn.proton.p.z);
  ampl_intd *= (1/(2*Pn_middle.E))*sqrt(2*pow(2*M_PI,3)/(2*Pn_first.E));
  ampl_intd *= kine2flux(s_gn,0,Pp_first.m2())*vnphoto(s_gn,t_gn);
  ampl_intd *= kine2flux(s_vn,vmass2(),pow(ProtonMass,2))*vnscatter(s_vn,t_vn);
  ampl_intd *= kine2flux(s_pn,pow(ProtonMass,2),pow(ProtonMass,2))
    * pnscatter(s_pn,t_pn);
  ampl_intd *= wf;
  ampl_intd /= pow(2*M_PI,4);
  //  Scale integrand by appropriate Jacobian.
  ampl_intd *= k_vn*k_pn*(2*M_PI)*(2*M_PI);
  ff[0] = real(ampl_intd);
  ff[1] = imag(ampl_intd);
  return 0;
}

//  Double rescattering, pn followed by VN rescattering
static int ampl_intd_pnvn(const int *ndim, const double xx[], const int *ncomp,
    double ff[]) {
  ff[0] = 0;
  ff[1] = 0;
  //  Scaled variables.
  double k_pn = xx[0];
  double k_vn = xx[1];
  double phi_pn = 2*M_PI*xx[2];
  double phi_vn = 2*M_PI*xx[3];
  //  Intermediate quantities
  FourVector K_vn, K_pn, Pp_first, Pp_middle, Pn_first, Pn_middle, Pv_prime;
  double Delta_vn, Delta_pn, s_gn, s_vn, s_pn, t_gn, t_vn, t_pn;
  complex<double> ampl_intd, wf;
  //  Preliminary values for K's and Pn's
  K_pn.p.set(k_pn*cos(phi_pn),k_pn*sin(phi_pn),0);
  K_vn.p.set(k_vn*cos(phi_vn),k_vn*sin(phi_vn),0);
  Pn_middle.p.set(gdvpn.neutron.p.x-K_vn.p.x,gdvpn.neutron.p.y-K_vn.p.y,0);
  Pn_first.p.set(Pn_middle.p.x-K_pn.p.x,Pn_middle.p.y-K_pn.p.y,0);
  Pn_middle.build(pow(ProtonMass,2),Pn_middle.p);
  Pn_first.build(pow(ProtonMass,2),Pn_first.p);
  K_pn.E = Pn_middle.E - Pn_first.E;
  K_vn.E = gdvpn.neutron.E - Pn_middle.E;
  //  Deltas & final values for all four-momenta
  Delta_pn = K_pn.m2()/2 + gdvpn.proton.E*K_pn.E - gdvpn.proton.p*K_pn.p;
  Delta_pn /= gdvpn.proton.p.z;
  Delta_vn = K_vn.m2()/2 + gdvpn.vmeson.E*K_vn.E - gdvpn.vmeson.p*K_vn.p;
  Delta_vn /= gdvpn.vmeson.p.z;
  K_pn.p.set(K_pn.p.x,K_pn.p.y,Delta_pn);
  K_vn.p.set(K_vn.p.x,K_vn.p.y,Delta_vn);
  Pn_middle.build(pow(ProtonMass,2),gdvpn.neutron.p-K_vn.p);
  Pn_first.build(pow(ProtonMass,2),Pn_middle.p-K_pn.p);
  Pv_prime.build(vmass2(),gdvpn.vmeson.p+K_vn.p);
  Pp_middle.build(pow(ProtonMass,2),gdvpn.proton.p+K_pn.p);
  Pp_first = gdvpn.deuteron - Pn_first;
  if(Pp_first.E<=0||Pp_first.m2()<=0) return 0;
  //  Mandelstam variables
  s_gn = kine2s(gdvpn.gamma,Pp_first);
  s_pn = kine2s(gdvpn.proton,Pn_middle);
  s_vn = kine2s(gdvpn.neutron,gdvpn.vmeson);
  t_gn = kine2t(gdvpn.gamma,Pv_prime);
  t_pn = kine2t(gdvpn.proton,Pp_middle);
  t_vn = kine2t(gdvpn.vmeson,Pv_prime);
  //  Use w/f model
  wf = dwf(Pn_first.p);
  wf /= sqrt(Pp_first.E/ProtonMass);
  //  ampl_intd is the unscaled integrand
  ampl_intd = -0.25;
  ampl_intd /= (2*gdvpn.proton.p.z)*(2*gdvpn.vmeson.p.z);
  ampl_intd *= (1/(2*Pn_middle.E))*sqrt(2*pow(2*M_PI,3)/(2*Pn_first.E));
  ampl_intd *= kine2flux(s_gn,0,Pp_first.m2())*vnphoto(s_gn,t_gn);
  ampl_intd *= kine2flux(s_pn,pow(ProtonMass,2),pow(ProtonMass,2))
    * pnscatter(s_pn,t_pn);
  ampl_intd *= kine2flux(s_vn,vmass2(),pow(ProtonMass,2))*vnscatter(s_vn,t_vn);
  ampl_intd *= wf;
  ampl_intd /= pow(2*M_PI,4);
  //  Scale integrand by appropriate Jacobian.
  ampl_intd *= k_pn*k_vn*(2*M_PI)*(2*M_PI);
  ff[0] = real(ampl_intd);
  ff[1] = imag(ampl_intd);
  return 0;
}

//  Double rescattering, soft production followed by VN and then pn rescattering
static int ampl_intd_dsoft(const int *ndim, const double xx[], const int *ncomp,
    double ff[]) {
  ff[0] = 0;
  ff[1] = 0;
  //  Scaled variables.
  double k_vn = xx[0];
  double k_pn = xx[1];
  double phi_vn = 2*M_PI*xx[2];
  double phi_pn = 2*M_PI*xx[3];
  //  Intermediate quantities
  FourVector K_vn, K_pn, Pp_first, Pp_middle, Pn_first, Pn_middle, Pv_prime;
  double Delta_vn, Delta_pn, s_gn, s_vn, s_pn, t_gn, t_vn, t_pn;
  complex<double> ampl_intd, wf;
  //  Preliminary values for K's and Pn's
  K_vn.p.set(k_vn*cos(phi_vn),k_vn*sin(phi_vn),0);
  K_pn.p.set(k_pn*cos(phi_pn),k_pn*sin(phi_pn),0);
  Pv_prime.build(vmass2(),gdvpn.gamma.p-K_vn.p);
  Pp_middle.build(pow(ProtonMass,2),gdvpn.proton.p+K_pn.p);
  K_vn.E = gdvpn.gamma.E - Pv_prime.E;
  K_pn.E = Pp_middle.E - gdvpn.proton.E;
  //  Deltas & final values for all four-momenta
  Delta_vn = vmass2()/2 - K_vn.m2()/2 + K_vn.E*gdvpn.gamma.E;
  Delta_vn /= gdvpn.gamma.E;
  Delta_pn = K_pn.m2()/2 + gdvpn.proton.E*K_pn.E - gdvpn.proton.p*K_pn.p;
  Delta_pn /= gdvpn.proton.p.z;
  K_vn.p.set(K_vn.p.x,K_vn.p.y,Delta_vn);
  K_pn.p.set(K_pn.p.x,K_pn.p.y,Delta_pn);
  Pv_prime.build(vmass2(),gdvpn.gamma.p-K_vn.p);
  Pp_middle.build(pow(ProtonMass,2),gdvpn.proton.p+K_pn.p);
  Pn_middle.build(pow(ProtonMass,2),gdvpn.neutron.p-K_pn.p);
  Pp_first.build(pow(ProtonMass,2),K_vn.p-Pn_middle.p);
  Pn_first = gdvpn.deuteron - Pp_first;
  if(Pn_first.E<=0||Pp_first.m2()<=0) return 0;
  //  Mandelstam variables
  s_gn = kine2s(gdvpn.gamma,Pn_first);
  s_vn = kine2s(gdvpn.vmeson,Pp_middle);
  s_pn = kine2s(gdvpn.proton,gdvpn.neutron);
  t_gn = kine2t(gdvpn.gamma,Pv_prime);
  t_vn = kine2t(gdvpn.vmeson,Pv_prime);
  t_pn = kine2t(gdvpn.neutron,Pn_middle);
  //  Use w/f model
  wf = dwf(Pp_first.p);
  wf /= sqrt(Pn_first.E/ProtonMass);
  wf *= -1; // isospin of struck neutron
  //  ampl_intd is the unscaled integrand
  ampl_intd = -0.25;
  ampl_intd /= (2*gdvpn.gamma.E)*(2*gdvpn.proton.p.z);
  ampl_intd *= (1/(2*Pn_middle.E))*sqrt(2*pow(2*M_PI,3)/(2*Pp_first.E));
  ampl_intd *= kine2flux(s_gn,0,Pn_first.m2())*vnphoto(s_gn,t_gn);
  ampl_intd *= kine2flux(s_vn,vmass2(),pow(ProtonMass,2))*vnscatter(s_vn,t_vn);
  ampl_intd *= kine2flux(s_pn,pow(ProtonMass,2),pow(ProtonMass,2))
    * pnscatter(s_pn,t_pn);
  ampl_intd *= wf;
  ampl_intd /= pow(2*M_PI,4);
  //  Scale integrand by appropriate Jacobian.
  ampl_intd *= k_vn*k_pn*(2*M_PI)*(2*M_PI);
  ff[0] = real(ampl_intd);
  ff[1] = imag(ampl_intd);
  return 0;
}

//  VN rescattering, Cauchy PV part
static int ampl_intd_vnc(const int *ndim, const double xx[], const int *ncomp,
    double ff[]) {
  ff[0] = 0;
  ff[1] = 0;
  //  Scaled variables.
  double k_perp = xx[0];
  double phi = 2*M_PI*xx[1];
  double k_z = 2*(xx[2]-0.5);
  //  Intermediate quantities
  FourVector K, Pn_prime;
  complex<double> ampl_intd;
  double Delta;
  //  Figure out Delta
  K.p.set(k_perp*cos(phi),k_perp*sin(phi),0);
  Pn_prime.p.set(gdvpn.neutron.p.x-K.p.x,gdvpn.neutron.p.y-K.p.y,0);
  Pn_prime.build(pow(ProtonMass,2),Pn_prime.p);
  K.E = gdvpn.neutron.E - Pn_prime.E;
  Delta = K.m2()/2 + K.E*gdvpn.vmeson.E - K.p*gdvpn.vmeson.p;
  Delta /= gdvpn.vmeson.p.z;
  ampl_intd = cauchy_h_vn(k_z,Delta,K);
  //  Scale integrand by appropriate Jacobian.
  ampl_intd *= k_perp*(2*M_PI)*2;
  ff[0] = real(ampl_intd);
  ff[1] = imag(ampl_intd);
  return 0;
}

//  pn rescattering, Cauchy PV part
static int ampl_intd_pnc(const int *ndim, const double xx[], const int *ncomp,
    double ff[]) {
  ff[0] = 0;
  ff[1] = 0;
  //  Scaled variables.
  double k_perp = xx[0];
  double phi = 2*M_PI*xx[1];
  double k_z = 2*(xx[2]-0.5);
  //  Intermediate quantities
  FourVector K,Pn_prime;
  complex<double> ampl_intd;
  double Delta;
  //  Figure out Delta
  K.p.set(k_perp*cos(phi),k_perp*sin(phi),0);
  Pn_prime.p.set(gdvpn.neutron.p.x-K.p.x,gdvpn.neutron.p.y-K.p.y,0);
  Pn_prime.build(pow(ProtonMass,2),Pn_prime.p);
  K.E = gdvpn.neutron.E - Pn_prime.E;
  Delta = K.m2()/2 + K.E*gdvpn.vmeson.E - K.p*gdvpn.proton.p;
  Delta /= gdvpn.proton.p.z;
  ampl_intd = cauchy_h_pn(k_z,Delta,K);
  //  Scale integrand by appropriate Jacobian.
  ampl_intd *= k_perp*(2*M_PI)*2;
  ff[0] = real(ampl_intd);
  ff[1] = imag(ampl_intd);
  return 0;
}

// =============================================================================
// Extra Functions for Cauchy PV Parts

// The Cauchy PV integrand is of the form f(z)/(z-Delta).
// The trick here is that the PV integral of f(Delta)/(z-Delta) = 0,
// so we really integrate h(z)=(f(z)-f(Delta))/(z-Delta).
// The amplt_intd routines above find Delta and plug it into h,
// which (along with f) is defined here.

static complex<double> cauchy_f_vn(double k_z, FourVector K) {
  FourVector Pv_prime, Pp_prime, Pn_prime;
  double s_gn, t_gn, s_vn, t_vn, virtuality;
  complex<double> ampl_intd, wf;
  //  Kinematics
  K.p.z = k_z;
  Pn_prime.build(pow(ProtonMass,2),gdvpn.neutron.p-K.p);
  K.E = gdvpn.neutron.E-Pn_prime.E;
  if(K.m2()>=0) return 0;
  Pv_prime = gdvpn.vmeson + K;
  virtuality = abs(vmass2()-Pv_prime.m2());
  Pp_prime = gdvpn.deuteron - Pn_prime;
  if(Pp_prime.E<=0||Pp_prime.m2()<=0) return 0;
  //  Mandelstam variables
  s_gn = kine2s(Pv_prime,gdvpn.proton);
  t_gn = kine2t(gdvpn.gamma,Pv_prime);
  s_vn = kine2s(gdvpn.vmeson,gdvpn.neutron);
  t_vn = kine2t(Pv_prime,gdvpn.vmeson);
  if(t_gn>=0) return 0;
  //  Use w/f model
  wf = dwf(Pn_prime.p);
  wf /= sqrt(Pp_prime.E/ProtonMass);
  //  Get ampl_intd and return
  ampl_intd = 1;
  ampl_intd /= 2*gdvpn.vmeson.p.z;
  ampl_intd *= sqrt(2*pow(2*M_PI,3)/(2*Pn_prime.E));
  ampl_intd *= kine2flux(gdvpn.gamma,Pp_prime)*vnphoto(s_gn,t_gn);
  ampl_intd *= kine2flux(Pv_prime,Pn_prime)*vnscatter(s_vn,t_vn);
  ampl_intd *= wf;
  ampl_intd /= pow(2*M_PI,3);
  return ampl_intd;
}

static complex<double> cauchy_f_pn(double k_z, FourVector K) {
  FourVector Pp_middle, Pp_first, Pn_prime;
  double s_gn, s_pn, t_gn, t_pn, virtuality;
  complex<double> ampl_intd, wf;
  //  Kinematics
  K.p.z = k_z;
  Pn_prime.build(pow(ProtonMass,2),gdvpn.neutron.p-K.p);
  K.E = gdvpn.neutron.E - Pn_prime.E;
  if(K.m2()>=0) return 0;
  Pp_middle = gdvpn.proton + K;
  virtuality = abs(pow(ProtonMass,2)-Pp_middle.m2());
  Pp_first = gdvpn.deuteron - Pn_prime;
  if(Pp_first.E<=0||Pp_first.m2()<=0) return 0;
  //  Mandelstam variables
  s_gn = kine2s(Pp_middle,gdvpn.vmeson);
  s_pn = kine2s(gdvpn.proton,gdvpn.neutron);
  t_gn = kine2t(gdvpn.gamma,gdvpn.vmeson);
  t_pn = kine2t(Pp_middle,gdvpn.proton);
  if(t_pn>=0) return 0;
  //  Use w/f model
  wf = dwf(Pn_prime.p);
  wf /= sqrt(Pp_first.E/ProtonMass);
  //  Get ampl_intd and return
  ampl_intd = 1;
  ampl_intd /= 2*gdvpn.vmeson.p.z;
  ampl_intd *= sqrt(2*pow(2*M_PI,3)/(2*Pn_prime.E));
  ampl_intd *= kine2flux(gdvpn.gamma,Pp_first)*vnphoto(s_gn,t_gn);
  ampl_intd *= kine2flux(Pp_middle,Pn_prime)*pnscatter(s_pn,t_pn-virtuality);
  ampl_intd *= wf;
  ampl_intd /= pow(2*M_PI,3);
  return ampl_intd;
}

static complex<double> cauchy_h_vn(double k_z, double Delta, FourVector K) {
  complex<double> ampl_intd;
  static const double eps=10e-6;
  if(abs(k_z-Delta)<eps)
    ampl_intd=(cauchy_f_vn(Delta+eps,K)-cauchy_f_vn(Delta-eps,K))/(2*eps);
  else
    ampl_intd=(cauchy_f_vn(k_z,K)-cauchy_f_vn(Delta,K))/(k_z-Delta);
  return ampl_intd;
}

static complex<double> cauchy_h_pn(double k_z, double Delta, FourVector K) {
  complex<double> ampl_intd;
  static const double eps=10e-6;
  if(abs(k_z-Delta)<eps)
    ampl_intd=(cauchy_f_pn(Delta+eps,K)-cauchy_f_pn(Delta-eps,K))/(2*eps);
  else
    ampl_intd=(cauchy_f_pn(k_z,K)-cauchy_f_pn(Delta,K))/(k_z-Delta);
  return ampl_intd;
}
// =============================================================================
//  Kinematics and phase factor

// Kinematic solver
// 
// This routine takes as arguments:
//   -  the magnitude of the (spectator) neutron's momentum;
//   - the relative angle between the momentum transfer vector
//     (q=pg-pv) and the neutron momentum; 
//   -  and t. 
// The four-vectors deuteron and Pg should be set beforehand.
// solno should be input.
// This routine assumes the reaction takes place all in the xz plane 
// and that the vector meson travels in the upper half of the plane.
//
// Flag error code:
//    1 : solno went out of bounds
//    2 : neutron energy >= deuteron mass
//    3 : momentum transfer too small
//    4 : meson energy less than meson mass
//    5 : Meson angle complex
//    6 : Neutron angle complex
//
static int gdvpn_kinematics(double args[5], double mv, int solno) {
  //  args[] = {t,p_n,theta_nl,phi_nl,phi_v}
  double t = args[0], p_n = args[1], theta_nl = args[2], 
         phi_nl = args[3], phi_v = 0.;
  FourVector l;
  double a, b, c;
  double l0, lp, p1, p2;
  double p_v, theta_v;
  double theta_n, phi_n;
  int nosol;
  gdvpn.neutron.E = sqrt(pow(ProtonMass,2)+pow(p_n,2));
  if(gdvpn.neutron.E>=DeuteronMass) return 2;
  a = 2.*(DeuteronMass-gdvpn.neutron.E);
  b = -2.*p_n*cos(theta_nl);
  c = -pow(DeuteronMass,2) - t + 2.*DeuteronMass*gdvpn.neutron.E;
  solvep(a,b,c,t,p1,p2,nosol);
  if(solno>nosol) return 1; // solno went out of bounds
  if(solno==1) lp = p1;
  if(solno==2) lp = p2;
  if((-t)>=pow(lp,2)) return 3;
  l0 = sqrt(pow(lp,2)+t);
  gdvpn.vmeson.E = gdvpn.gamma.E - l0;
  if(gdvpn.vmeson.E<mv) return 4;
  p_v = sqrt(pow(gdvpn.vmeson.E,2)-pow(mv,2));
  theta_v = t + 2.*gdvpn.gamma.E*gdvpn.vmeson.E-pow(mv,2);
  theta_v /= 2.*gdvpn.gamma.E*p_v;
  if(theta_v>1.&&theta_v<1.00001) theta_v=1.;
  if(theta_v<-1.&&theta_v>-1.00001) theta_v=-1.;
  if(abs(theta_v)>1.) return 5;
  theta_v = acos(theta_v);
  gdvpn.vmeson.p.setsp(p_v,theta_v,phi_v);
  l = gdvpn.gamma - gdvpn.vmeson;
  if(l.p.theta()==0.) {
    theta_n = theta_nl;
    phi_n = phi_nl;
  }
  else {
    theta_n = cos(theta_nl)*cos(l.p.theta());
    theta_n += sin(theta_nl)*sin(l.p.theta())*cos(phi_nl-l.p.phi());
    if(theta_n>1.&&theta_n<1.00001) theta_n=1.;
    if(theta_n<-1.&&theta_n>-1.00001) theta_n=-1.;
    if(abs(theta_n)>1.) return 6;
    theta_n = acos(theta_n);
    phi_n = cos(theta_nl) - cos(theta_n)*cos(l.p.theta());
    phi_n /= sin(theta_n)*sin(l.p.theta());
    if(phi_n>1.&&phi_n<1.00001) phi_n=1.;
    if(phi_n<-1.&&phi_n>-1.00001) phi_n=-1.;
    if(abs(phi_n)>1.) return 6;
    phi_n = acos(phi_n) + l.p.phi();
  }
  gdvpn.neutron.p.setsp(p_n,theta_n,phi_n);
  gdvpn.proton = gdvpn.gamma + gdvpn.deuteron - gdvpn.vmeson - gdvpn.neutron;
  return 0;
}

//  Phase space factor from diffsig=M^2/flux*dQ.
static double phase_factor() {
  double phase;
  FourVector l;
  phase = pow(gdvpn.proton.p.r(),2) * pow(gdvpn.vmeson.p.r(),2);
  phase /= pow(2.*M_PI,5) * (8.*gdvpn.proton.E*gdvpn.vmeson.E*gdvpn.neutron.E);
  phase *= gdvpn.proton.p.r();
  phase /= abs(gdvpn.proton.p*(gdvpn.proton.beta()-gdvpn.neutron.beta()));
  return phase;
}

// =============================================================================
// Other Math Routines

// Solves aE=bp+c for p. Returns physical solutions & their number.
// Needs m^2 input.
// Uses m^2 rather than m so that it allows m^2 < 0. 
static void solvep(double a, double b, double c, double m2,
    double &p1, double &p2, int &n) {
  double A,B,C,test;
  A=a*a-b*b;
  B=-2*b*c;
  C=a*a*m2-c*c;
  solvequad(A,B,C,p1,p2,n);
  // This switch-case tests whether each p is physical.
  // p1>p2 if solvequad gave us two solutions, 
  // so if p1 fails the tests so does p2, but p2 can fail while p1 passes.
  //
  // The tests are:
  //    p must be greater than zero. (This equation is for magnitude.)
  //    (bp+c)/a must be greater than zero. (This quantity=E.)
  //    [(bp+c)/a]^2 must be greater than m^2. (This quantity=E^2.)
  //
  switch(n) {
    case 0:
      return;
      break;
    case 1:
      if(p1<0) { n=0; return; }
      test=b*p1+c;
      test/=a;
      if(test<0) { n=0; return; }
      test*=test;
      if(test<m2) { n=0; return; }
      break;
    case 2:
      if(p1<0) { n=0; return; }
      test=b*p1+c;
      test/=a;
      if(test<0) { n=0; return; }
      test*=test;
      if(test<m2) { n=0; return; }
      if(p2<0) { n=1; return; }
      test=b*p2+c;
      test/=a;
      if(test<0) { n=1; return; }
      test*=test;
      if(test<m2) { n=1; return; }
      break;
  }
}

//  Solves Ax^2+Bx+C=0 for x. Returns real solutions & their number.
static void solvequad(double A, double B, double C,
    double &x1, double &x2, int &n) {
  double disc;
  disc = pow(B,2) - 4*A*C;
  if(A==0) {
    if(B==0) { n=0; return; } // Was given C=0, which is trivial or impossible.
    n = 1;
    x1 = -C/B;
    return;
  }
  if(disc<0) {
    n = 0;
    return;
  }
  if(disc==0) {
    n = 1;
    x1 = -B/(2*A);
    return;
  }
  n = 2;
  x1 = (-B + sqrt(disc)) / (2*A);
  x2 = (-B - sqrt(disc)) / (2*A);
  return;
}

// =============================================================================
//  Deuteron Wave-Function Shortcut

static complex<double>dwf(ThreeVector p_vec) {
  return dwf(p_vec.r(),p_vec.theta(),p_vec.phi(),gdvpn.sd,gdvpn.sp,gdvpn.sn);
}

// =============================================================================
//  Routines for Screen Output

//  Outputs info for run to screen
static time_t start(double t, double p, char *filename) {
  bool stated_term = false;
  cout << "Starting gdvpn, with:" << endl;
  if(gdvpn.pwia)
    cout << "\tPWIA";
  if(gdvpn.singles) {
    cout << "\tSingle resc. ( ";
    if(gdvpn.vnvertex)
      cout << "VN ";
    if(gdvpn.pnvertex)
      cout << "pn ";
    if(gdvpn.soft)
      cout << "soft ";
    cout << ")";
  }
  if(gdvpn.cauchy)
    cout << "\tCauchy PV";
  if(gdvpn.doubles)
    cout << "\tDouble resc.";
  cout << endl;
  cout << "\tE = " << gdvpn.gamma.E << " GeV\tp = " 
    << p << " GeV\tt = " << t << " GeV**2" << endl;
  cout << "\tData being output to " << filename << endl;
  return time(NULL);
}

//  Outputs run time to screen
static void end(time_t alpha) {
  time_t omega = time(NULL);
  double runtime = (double)(omega-alpha);
  cout << "Finished. Took " << runtime << " second(s)." << endl;
}

//  Outputs progress to screen
static void progress_bar(double completed) {
  int hashes, place, percent;
  cout << "\t[";
  percent = 100*completed;
  hashes = percent/2;
  for(place=1;place<=hashes;place++)
    cout << "#";
  for(true;place<=50;place++)
    cout << " ";
  cout << "] " << percent << "%" << flush << "\r";
}
