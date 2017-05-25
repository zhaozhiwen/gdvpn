//
// vmeson_params.cc
//
//
// by Adam Freese, 2012-2017
//
// Contains methods for computing gamma N -> V N and V N -> V N cross sections.
//
#include <stdlib.h>
#include <cmath>
#include <complex>
#include "constants.hh"
using namespace std;

//  ID #'s for vector meson and parameterization
static int vmeson_id = 0; // J/Psi by default
static int param_id  = 2; // Brodsky et al.'s 2-gluon exchange by default

//  Other (real) parameters
static double vmeson_mass = JPsiMass; // J/Psi by default
static double sigma_vn    = 10.;      // J/Psi by default
static double B_vn        = -1;       // VMD by default
static double alpha_vn    = -0.2;     // J/Psi by default

//  Declarations of private functions
static complex<double> vnphoto_uncut(double,double);
static double diffsig_gnvn(double,double);
static double photo_tmin(double);
static double photo_tmax(double);
static double photo_t90(double);
static double scatter_tmax(double);
//  Meson types
static double jpsi_param(double,double);        // vmeson_id=0
static double phi_param(double,double);         // vmeson_id=1
//  Parameterizations available
//  J/Psi
static double vmp_jpsi_bchl2g(double,double);   // param_id=2
static double vmp_jpsi_bchl3g(double,double);   // param_id=3 ... gone, now v
static double vmp_jpsi_gittel(double);          // param_id=3
//  phi(1020)
static double vmp_phi_bsyp1(double);            // param_id=1
static double vmp_phi_bsyp2(double,double);     // param_id=2
static double vmp_phi_mibe(double,double);      // param_id=3

// =============================================================================
// Public accessors

//  Get meson mass, or meson mass squared
double vmass() {
  return vmeson_mass;
}

double vmass2() {
  return pow(vmeson_mass,2);
}

// =============================================================================
// Public mutator, to set meson and scattering model IDs

void vnset(int vmeson_id_arg, int param_id_arg,
    double sigma_vn_arg, double B_vn_arg) {
  switch(vmeson_id_arg) {
  case 0: // J/Psi
    vmeson_mass = JPsiMass;
    sigma_vn = 5.0;
    alpha_vn = -0.2;
    break;
  case 1: // phi(1020)
    vmeson_mass = Phi1020Mass;
    sigma_vn = 10.0;
    alpha_vn = -0.5;
    break;
  default: // not a valid meson
    exit(EXIT_FAILURE);
  }
  vmeson_id = vmeson_id_arg;
  param_id = param_id_arg;
  if(sigma_vn>0)
    sigma_vn = sigma_vn_arg;
  B_vn = B_vn_arg;
}

// =============================================================================
// Methods for dsigma/dt, which return answers in GeV**-2,
// and take arguments in GeV**2

//  Photoproduction amplitude
complex<double> vnphoto(double s, double t) {
  complex<double> ampl;
  double tmin, tmax, smin;
  smin = pow(ProtonMass+vmass(),2);
  if(s<smin) return 0;
  tmin = photo_tmin(s);
  tmax = photo_tmax(s);
  if(t>tmin||t<tmax) return 0;
  ampl = vnphoto_uncut(s,t);
  return ampl;
}

//  Elastic scattering amplitude
// B_vn < 0 means vector meson dominance is used
complex<double> vnscatter(double s, double t) {
  complex<double> ampl;
  double tmax;
  tmax = scatter_tmax(s);
  if(t<tmax) return 0;
  ampl = (sigma_vn/mb_GeV);
  if(B_vn>0) {
    ampl *= complex<double>(alpha_vn,1);
    ampl *= exp(B_vn*t/2);
  }
  else {
    ampl *= vnphoto_uncut(s,t);
    ampl /= vnphoto_uncut(s,0);
  }
  return ampl;
}

// =============================================================================
// Photoproduction amplitude wihout a cutoff

//  Photoproduction amplitude, without any t cuts
static complex<double> vnphoto_uncut(double s, double t) {
  complex<double> ampl;
  ampl = diffsig_gnvn(s,t);
  ampl *= 16.0*M_PI;
  ampl /= 1 + pow(alpha_vn,2);
  ampl = sqrt(ampl);
  ampl *= complex<double>(alpha_vn,1);
  return ampl;
}

// =============================================================================
//  Switches

//  This routine switches between mesons
static double diffsig_gnvn(double s, double t) {
  switch(vmeson_id) {
  case 0: return jpsi_param(s,t);
  case 1: return phi_param(s,t);
  default: return 0;
  }
}

//  Switches between J/Psi parameterizations
static double jpsi_param(double s, double t) {
  switch(param_id) {
/*  case 1: return vmp_jpsi_rsz(s,t);*/
  case 2: return vmp_jpsi_bchl2g(s,t);
  case 3: return vmp_jpsi_gittel(t);
  default: return 0;
  }
}

//  Switches between phi(1020) parameterizations
static double phi_param(double s, double t) {
  switch(param_id) {
  case 1: return vmp_phi_bsyp1(t);
  case 2: return vmp_phi_bsyp2(s,t);
  case 3: return vmp_phi_mibe(s,t);
  default: return 0;
  }
}

// =============================================================================
// J/Psi Parameterizations

// The following parameterizations are from Phys. Lett. B498 (2001) 23,
// by S.J. Brodsky, E. Chudakov, P. Hoyer, and J.M. Laget,
// hereafter referred to as BCHL.

// Equation (3) of BCHL, corresponding to 2-gluon exchange.
// I took the liberty of simplifying the equation by using the meanings of 
// x and nu provided in the paper, and got the normalization constant 
// by fitting to the data in Phys. Rev. Lett. 35 (1975) 483.
// The t-dependence is from Phys. Rev. D66, 031502
static double vmp_jpsi_bchl2g(double s, double t) {
  static const double N = 14.7/(nb_mb*mb_GeV);
  //  static const double B = 1.13;
  double B = 4./(1-t);
  static const double smin = pow(ProtonMass+JPsiMass,2);
  return N*pow((s-smin)/(s-pow(ProtonMass,2)),2)*exp(B*t);
}

/*
// Equation (4) of BCHL, corresponding to 3-gluon exchange.
// This is the same as the Gittelman et al. parameterization,
// after all the algebraic manipulation is done, except B=1.13.
//  I may remove this later as superfluous.
static double vmp_jpsi_bchl3g(double s, double t) {
  static const double A = 1.01/(nb_mb*mb_GeV);
  static const double B = 1.13;
  return A*exp(B*t);
}
*/

// This parameterization from Phys. Rev. Lett. 35 (1975) 1616,
// by B. Gittelman, K.M Hanson, D. Larson, E. Loh, A. Silverman, 
// and G. Theodosiou.
static double vmp_jpsi_gittel(double t) {
  static const double A = 1.01/(nb_mb*mb_GeV);
  static const double B = 1.25;
  return A*exp(B*t);
}

// =============================================================================
//  phi(1020) Parameterizations

// Two phi(1020) parameterizations from Rev. Mod. Phys. 50 (1978) 261,
// by T.H. Bauer, R.D. Spital, D.R. Yennie, and F.M. Pipkin,
// hereafter referred to as BSYP.

// Equation (3.85a) of BSYP
static double vmp_phi_bsyp1(double t) {
  static const double A = 2.59/(ub_mb*mb_GeV);
  static const double B = 5.9;
  static const double C = 1.4;
  if(t<(-0.5*B/C)) return 0;
  return A*exp(B*t+C*pow(t,2));
}

// Equations (3.85b) and (3.85c) of BSYP
static double vmp_phi_bsyp2(double s, double t) {
  static const double A = 1.34/(ub_mb*mb_GeV);
  static const double B = 4.8;
  static const double C = 1.7;
  if(t<(-0.5*B/C)) return 0;
  double alpha = 1.14 + 0.27*t;
  return pow(s,2*(alpha-1))*A*exp(B*t+C*pow(t,2));
}

// This parameterization is due to T. Mibe.
// It contains the fit for dsigma/dt(t=tmin) in Figure 3 of 
// Phys. Rev. Lett. 95 (2005) 182001, and uses B=3.38 as in the paper.
static double vmp_phi_mibe(double s, double t) {
  static const double C[6] = {-6.9937,8.2732,-3.2595,0.67068,-0.06967,0.0028878};
  static const double B = 3.38;
  double tmin = photo_tmin(s);
  double E = (s-pow(ProtonMass,2))/(2*ProtonMass);
  double A = 0;
  for(int i=0;i<6;i++)
    A += C[i]*pow(E,i);
  A *= exp(-B*tmin);
  A /= ub_mb*mb_GeV;  
  return A*exp(B*t);
}

// =============================================================================
//  Private Kinematical Routines

static double photo_tmin(double s) {
  double Ef, pf, Ei, tmin;
  Ei = (s-pow(ProtonMass,2)) / (2*sqrt(s));
  Ef = (s+vmass2()-pow(ProtonMass,2)) / (2*sqrt(s));
  pf = sqrt(pow(Ef,2)-vmass2());
  tmin = vmass2() - 2*Ei*(Ef-pf);
  return tmin;
}

static double photo_tmax(double s) {
  double Ef, pf, Ei, tmax;
  Ei = (s-pow(ProtonMass,2)) / (2*sqrt(s));
  Ef = (s+vmass2()-pow(ProtonMass,2)) / (2*sqrt(s));
  pf = sqrt(pow(Ef,2)-vmass2());
  tmax = vmass2() - 2*Ei*(Ef+pf);
  return tmax;
}

static double photo_t90(double s) {
  double Ef, Ei, t90;
  Ei = (s-pow(ProtonMass,2)) / (2*sqrt(s));
  Ef = (s+vmass2()-pow(ProtonMass,2)) / (2*sqrt(s));
  t90 = vmass2() - 2*Ei*Ef;
  return t90;
}

static double scatter_tmax(double s) {
  double tmax, E, p2;
  E = (s+vmass2()-pow(ProtonMass,2)) / (2*sqrt(s));
  p2 = pow(E,2) - vmass2();
  tmax = -4*p2;
  return tmax;
}
