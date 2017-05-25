//
//  nnscatter.cc
//  
//
//  by Adam Freese, 2012-2017, based on Fortran 77 code by Misak Sargsian.
//
//  This code computes the scattering amplitude for NN-->NN reactions, 
//  normalized to obey the optical theorem.
//  The result is in GeV**-2
//
//  init_nnscatter() *must* be called before use!
//
#include <cmath>
#include <complex>
#include <fstream>
#include <iostream>
#include <cstdlib>
using namespace std;

static const double M_PRTN=0.938272;
static const double MBPGEV=0.389379338; //mb per gev^-2

// Internal arrays
static bool B_initialized = false;
static double Bystricky_pn[10][400][91], 
              Bystricky_pp[10][400][91], 
              Bystricky_nn[10][400][91];

// Forward declarations
static double alpha(double,int);
static double B(double);
static double sigma_tot(double,int);

// =============================================================================
// Initialization - call this once before using differential cross sections!
// =============================================================================

// Call this method before using cross section!
void init_nnscatter() {
  static ifstream pn_data, pp_data, nn_data;
  if(!B_initialized) {
    B_initialized = true;
    pn_data.open("B_ampl_pn.tbl");
    pp_data.open("B_ampl_pp.tbl");
    nn_data.open("B_ampl_pp_nocoul.tbl");
    if( !pn_data.is_open() || !pp_data.is_open() || !nn_data.is_open() ) {
      cout << "NN scattering amplitude files not found, quitting." << endl;
      exit(EXIT_FAILURE);
    }
    for(int j=0;j<10;j++)
      for(int k=0;k<400;k++)
        for(int l=0;l<=90;l++) {
          pn_data >> Bystricky_pn[j][k][l];
          pp_data >> Bystricky_pp[j][k][l];
          nn_data >> Bystricky_nn[j][k][l];
        }
    pn_data.close();
    pp_data.close();
    nn_data.close();
  }
  return;
}

// =============================================================================
// Differential cross sections
//
//  d(sigma)/dt as functions of s and t, each in GeV**2
//  These are returned in units of GeV**-2
// =============================================================================

complex<double> pnscatter(double s, double t) {
  complex<double> damp;
  damp=(sigma_tot(s,1)/MBPGEV)*complex<double>(alpha(s,1),1)*exp(B(s)*t/2);
  if((-t>(s-pow(2*M_PRTN,2)))) return 0;
  return damp;
}

complex<double> ppscatter(double s, double t) {
  complex<double> damp;
  damp=(sigma_tot(s,2)/MBPGEV)*complex<double>(alpha(s,2),1)*exp(B(s)*t/2);
  if((-t>(s-pow(2*M_PRTN,2)))) return 0;
  return damp;
}

complex<double> nnscatter(double s, double t) {
  complex<double> damp;
  damp=(sigma_tot(s,3)/MBPGEV)*complex<double>(alpha(s,1),3)*exp(B(s)*t/2);
  if((-t>(s-pow(2*M_PRTN,2)))) return 0;
  return damp;
}

// =============================================================================
// Private (internal) methods
// =============================================================================

static double p_lab(double s) {
  double p;
  p=s-2*pow(M_PRTN,2);
  p/=2*M_PRTN;
  p=sqrt(pow(p,2)-pow(M_PRTN,2));
  return p;
}

static double p_com(double s) {
  double p_lb=p_lab(s);
  double p_cm=M_PRTN*p_lb/sqrt(s);
  return p_cm;
}

static void Bystricky_routine(complex<double> (&B)[5], double p_lab, 
    double a_cm, int NNtype) {
  double h_p, h_a;
  int p_mesh, a_mesh;
  double b00[10], b10[10], b01[10], b11[10];
  h_p=10;
  h_a=2;
  p_lab*=1000;
  if(a_cm>=180) a_cm=180;
  p_mesh = (int)(p_lab/h_p);
  a_mesh = (int)(a_cm/h_a);
  if(p_mesh<1) p_mesh=1;
  if(p_lab>=3995) p_mesh=399;
  for(int j=0;j<10;j++) {
    b00[j] = 0;
    b10[j] = 0;
    b01[j] = 0;
    b11[j] = 0;
    if(p_lab>=10) {
      switch(NNtype) {
        case 1: // pn
          b00[j] = Bystricky_pn[j][p_mesh-1][a_mesh];
          b10[j] = Bystricky_pn[j][p_mesh][a_mesh];
          if(a_mesh<90) {
            b01[j] = Bystricky_pn[j][p_mesh-1][a_mesh+1];
            b11[j] = Bystricky_pn[j][p_mesh][a_mesh+1];
          }
          break;
        case 2: // pp
          b00[j] = Bystricky_pp[j][p_mesh-1][a_mesh];
          b10[j] = Bystricky_pp[j][p_mesh][a_mesh];
          if(a_mesh<90) {
            b01[j] = Bystricky_pp[j][p_mesh-1][a_mesh+1];
            b11[j] = Bystricky_pp[j][p_mesh][a_mesh+1];
          }
          break;
        case 3: // nn
          b00[j] = Bystricky_nn[j][p_mesh-1][a_mesh];
          b10[j] = Bystricky_nn[j][p_mesh][a_mesh];
          if(a_mesh<90) {
            b01[j] = Bystricky_nn[j][p_mesh-1][a_mesh+1];
            b11[j] = Bystricky_nn[j][p_mesh][a_mesh+1];
          }
          break;
      }
    }
  }
  for(int j=0;j<5;j++) {
    B[j] = complex<double> (
        b00[j]*(   (h_p*(p_mesh+1)-p_lab)/h_p * (h_a*(a_mesh+1)-a_cm)/h_a )
      - b10[j]*(   (h_p* p_mesh   -p_lab)/h_p * (h_a*(a_mesh+1)-a_cm)/h_a )
      - b01[j]*(   (h_p*(p_mesh+1)-p_lab)/h_p * (h_a* a_mesh   -a_cm)/h_a )
      + b11[j]*(   (h_p* p_mesh   -p_lab)/h_p * (h_a* a_mesh   -a_cm)/h_a )
      , b00[j+5]*( (h_p*(p_mesh+1)-p_lab)/h_p * (h_a*(a_mesh+1)-a_cm)/h_a )
      - b10[j+5]*( (h_p* p_mesh   -p_lab)/h_p * (h_a*(a_mesh+1)-a_cm)/h_a )
      - b01[j+5]*( (h_p*(p_mesh+1)-p_lab)/h_p * (h_a* a_mesh   -a_cm)/h_a )
      + b11[j+5]*( (h_p* p_mesh   -p_lab)/h_p * (h_a* a_mesh   -a_cm)/h_a )
    );
  }
  return;
}

static double cubic_spline(double p, double x1, double x2, double x3,
 double y1, double y2, double y3) {
  double dyx1, dyx2, yy1, yy2, yyy2, c1, c2, gamma2;
  dyx1 = (y2-y1)/(x2-x1);
  dyx2 = (y3-y2)/(x3-x2);
  gamma2 = 6*(dyx2-dyx1);
  yyy2 = gamma2/(2*(x3-x1));
  yy1 = dyx1 - yyy2*(x2-x1)/6;
  yy2 = dyx2 - yyy2*(x3-x2)/3;
  c1 = y1 + yy1*(p-x1) + yyy2*pow(p-x1,3)/(6*(x2-x1));
  c2 = y2 + yy2*(p-x2) + yyy2*pow(p-x2,2)/2 - yyy2*pow(p-x2,3)/(6*(x3-x2));
  if(p<=x2) return c1;
  else return c2;
}

static double cubic_spline(double p, double x[3], double y[3]) {
  return cubic_spline(p,x[0],x[1],x[2],y[0],y[1],y[2]);
}

static double B(double s) {
  double p=p_lab(s);
  if(p<0.8) {
    double x[3]={0,0.545,0.8};
    double y[3]={0,1,2.232239};
    return cubic_spline(p,x,y);
  }
  if(p<1.74)
    return -0.41579 + 1.8333*p + 2.6484*pow(p,2) - 1.0031*pow(p,3);
  if(p<=3) {
    double x[3]={1.74,2.4,3};
    double y[3]={5.508093,7.202996,7.903550};
    return cubic_spline(p,x,y);
  }
  return 8.22 + 1.1*log(p/4);
}

static double alpha(double s, int NNtype) {
  double p=p_lab(s);
  if(p<1.4) {
    complex<double> B[5];
    Bystricky_routine(B,p,0,NNtype);
    return real(B[0]+B[1])/imag(B[0]+B[1]);
  }
  if(p<3) {
    double x[3]={1.4,1.9,3};
    double y[3]={-0.2845016,-0.4391529,-0.4938042};
    return cubic_spline(p,x,y);
  }
  return -0.56207 + 0.024223*p - 0.00050362*pow(p,2) + 0.0000048408*pow(p,3)
    - 0.000000017331*pow(p,4);
}

static double sigma_tot(double s, int NNtype) {
  double p=p_lab(s), p_cm=p_com(s);
  p_cm/=sqrt(MBPGEV); //Converts p_cm to units of 1/sqrt(mb)
  if(p<1.4) {
    complex<double> B[5];
    Bystricky_routine(B,p,0,NNtype);
    return 4*M_PI/p_cm * imag(B[0]+B[1])/2;
  }
  if(p<=3.1) {
    double x[3]={1.4,2.078,3.1};
    double y[3]={38.42728,41.905,43.13159};
    return cubic_spline(p,x,y);
  }
  return 47.3 - 4.27*log(p) + 0.513*pow(log(p),2);
}
