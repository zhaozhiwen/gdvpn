//
//  dwf.cc
//
//
//  Created by Adam Freese on 4/23/12.
//  Update 2/22/13 :: added AV18 and Bonn potentials.
//  Update 5/19/17 :: AV18 is default instead of Paris.
//
//  This is based on code written in FORTRAN by Misak Sargsian to calculate 
// the deuteron wave-function. 
//
//  The Paris, AV18, or Bonn potential wave function can be used. 
// Potential is chosen by calling set_dwf.
// Calling set_dwf is not necessary; if it's skipped, Paris is used by default.
//
//  The function dwf takes six arguments
//    (1) : magnitude of the relative momentum in GeV
//    (2) : colatitude (theta) of relative momentum in radians
//    (3) : azimuth (phi) of relative momentum in radians
//    (4) : deuteron spin, in units of h-bar; should be 1, 0 or -1
//    (5) : proton spin, in units of h-bar/2; should be 1 or -1
//    (6) : neutron spin, in units of h-bar/2; should be 1 or -1
// and returns a quantity in GeV^{3/2}.
//
//
#include <cmath>
#include <complex>
#include <vector>
using namespace std;

// Internal variables
static vector<double> C, D, BM;
static bool initialized = false;

// Forward declarations
static void set_Paris();
static void set_AV18();
static void set_Bonn();

static double s_wave(double);
static double d_wave(double);

// =============================================================================
// Method to initialize deuteron wave function

void set_dwf(int choice) {
  switch(choice) {
    case 1: set_AV18(); break;
    case 2: set_Paris(); break;
    case 3: set_Bonn(); break;
    default: return; break;
  }
  initialized = true;
  return;
}

// =============================================================================
//  The wave function itself

complex<double>dwf(double p, double theta, double phi, int sd, int sp, int sn) {
  static const complex<double> IMAG = complex<double>(0,1);
  static const double GEVPFM = 0.1973269718;
  if(!initialized) set_dwf(1); // Initialize to AV18 if uninitialized.
  complex<double> wf;
  p /= GEVPFM;
  switch(sd) {
    case 1:
      if(sp==1&&sn==1)
        wf = s_wave(p) + d_wave(p)*(3.0*cos(theta)*cos(theta)-1.0)/sqrt(8.0);
      else if((sp==1&&sn==-1)||(sp==-1&&sn==1))
        wf = d_wave(p)*3.0*cos(theta)*sin(theta)/sqrt(8.0)*exp(IMAG*phi);
      else if(sp==-1&&sn==-1)
        wf = d_wave(p)*3.0*sin(theta)*sin(theta)/sqrt(8.0)*exp(2.0*IMAG*phi);
      break;
    case 0:
      if(sp==1&&sn==1)
        wf = d_wave(p)*1.5*cos(theta)*sin(theta)*exp(-IMAG*phi);
      else if((sp==1&&sn==-1)||(sp==-1&&sn==1))
        wf = s_wave(p)/sqrt(2.0)-d_wave(p)*(3.0*cos(theta)*cos(theta)-1.0)/2.0;
      else if(sp==-1&&sn==-1)
        wf = -d_wave(p)*1.5*cos(theta)*sin(theta)*exp(IMAG*phi);
      break;
    case -1:
      if(sp==1&&sn==1)
        wf = d_wave(p)*3.0*sin(theta)*sin(theta)/sqrt(8.0)*exp(-2.0*IMAG*phi);
      else if((sp==1&&sn==-1)||(sp==-1&&sn==1))
        wf = -d_wave(p)*3.0*cos(theta)*sin(theta)/sqrt(8.0)*exp(-IMAG*phi);
      else if(sp==-1&&sn==-1)
        wf = s_wave(p) + d_wave(p)*(3.0*cos(theta)*cos(theta)-1.0)/sqrt(8.0);
      break;
  }
  wf *= pow(GEVPFM,-1.5);
  return wf;
}

// =============================================================================
// =============================================================================

static double s_wave(double p) {
  double wave, A;
  A = 0;
  for(unsigned int j=0;j<C.size();j++)
    A += C[j]/(pow(p,2)+pow(BM[j],2));
  wave = A/sqrt(4.0*M_PI);
  return wave;
}

static double d_wave(double p) {
  double wave, A;
  A = 0;
  for(unsigned int j=0;j<D.size();j++)
    A += D[j]/(pow(p,2)+pow(BM[j],2));
  wave = A/sqrt(4.0*M_PI);
  return wave;
}

static void set_AV18() {
  const double CC[12] = {0.706699,-0.169743,1.12368,-8.52995,19.5033,-75.7831,
    283.739,-694.734,885.257,-720.739,412.969,-103.336};
  const double DD[12] = {0.0176655,-0.124551,-1.08815,3.84848,-8.52442,20.9435,
    -49.0728,57.7382,-1.27114,-62.8361,58.1016,-17.7062};
  const double BB[12] = {0.2316,1.0,1.5,2.0,2.5,3.5,4.5,5.5,6.5,8,9.5,11.0};
  C.resize(12);
  D.resize(12);
  BM.resize(12);
  for(int j=0;j<12;j++) {
    C[j] = CC[j];
    D[j] = DD[j];
    BM[j] = BB[j];
  }
  return;
}

static void set_Paris() {
  double A, B, G;
  const double CC[12] = {0.88688076,-0.3471093,-3.050238,56.207766,-749.57334,
    5336.5279,-22706.863,60434.4690,-102920.58,112233.57,-75925.226,29059.715};
  const double DD[10] = {0.023135193,-0.85604572,5.6068193,-69.462922,416.31118,
    -1254.6621,1238.783,3373.9172,-13041.151,19512.524};
  C.resize(13);
  D.resize(13);
  BM.resize(13);
  A = 0;
  for(int j=0;j<12;j++) {
    C[j] = CC[j];
    A += C[j];
  }
  C[12] = -A;
  for(int j=0;j<13;j++)
    BM[j] = 0.23162461+j;
  A = 0;
  B = 0;
  G = 0;
  for(int j=0;j<10;j++) {
    D[j] = DD[j];
    A += D[j]/pow(BM[j],2);
    B += D[j];
    G += D[j]*pow(BM[j],2);
  }
  D[10] = pow(BM[10],2);
  D[10] /= pow(BM[12],2) - pow(BM[10],2);
  D[10] /= pow(BM[11],2) - pow(BM[10],2);
  D[10] *= -A*pow(BM[11]*BM[12],2) + B*(pow(BM[11],2)+pow(BM[12],2)) - G;
  D[11] = pow(BM[11],2);
  D[11] /= pow(BM[10],2) - pow(BM[11],2);
  D[11] /= pow(BM[12],2) - pow(BM[11],2);
  D[11] *= -A*pow(BM[12]*BM[10],2) + B*(pow(BM[12],2)+pow(BM[10],2)) - G;
  D[12] = pow(BM[12],2);
  D[12] /= pow(BM[11],2) - pow(BM[12],2);
  D[12] /= pow(BM[10],2) - pow(BM[12],2);
  D[12] *= -A*pow(BM[10]*BM[11],2) + B*(pow(BM[10],2)+pow(BM[11],2)) - G;
  //  The following is added to MS's original routine so that 
  // the same s_wave and d_wave functions can be used for all potentials
  for(int j=0;j<13;j++) {
    C[j] *= sqrt(2.0/M_PI);
    D[j] *= sqrt(2.0/M_PI);
  }
  return;
}

static void set_Bonn() {
  const double CC[10] = {0.88472985,-0.26408759,-0.044114404,-14.397512,
    85.591256,-318.76761,703.36701,-900.49586,661.45441,-259.58894};
  const double DD[8] = {0.022623762,-0.50471056,0.56278897,-16.079764,
    111.26803,-446.67490,1098.5907,-1611.4995};
  double A, TM0, TM1, TM2;
  C.resize(11);
  D.resize(11);
  BM.resize(11);
  A = 0;
  for(int j=0;j<10;j++) {
    C[j] = CC[j];
    A += C[j];
  }
  C[10] = -A;
  for(int j=0;j<11;j++)
    BM[j] = 0.2315380 + j*0.9;
  TM0 = 0;
  TM1 = 0;
  TM2 = 0;
  for(int j=1;j<8;j++) {
    D[j] = DD[j];
    TM0 += D[j];
    TM1 += D[j] / pow(BM[j],2);
    TM2 += D[j] * pow(BM[j],2);
  }
  D[8] = pow(BM[8],2);
  D[8] /= pow(BM[10],2) - pow(BM[8],2);
  D[8] /= pow(BM[9],2) - pow(BM[8],2);
  D[8] *= -TM1*pow(BM[9]*BM[10],2) + TM0*(pow(BM[9],2)+pow(BM[10],2)) - TM2;
  D[9] = pow(BM[9],2);
  D[9] /= pow(BM[8],2) - pow(BM[9],2);
  D[9] /= pow(BM[10],2) - pow(BM[9],2);
  D[9] *= -TM1*pow(BM[10]*BM[8],2) + TM0*(pow(BM[10],2)+pow(BM[8],2)) - TM2;
  D[10] = pow(BM[10],2);
  D[10] /= pow(BM[9],2) - pow(BM[10],2);
  D[10] /= pow(BM[8],2) - pow(BM[10],2);
  D[10] *= -TM1*pow(BM[8]*BM[9],2) + TM0*(pow(BM[8],2)+pow(BM[9],2)) - TM2;
  for(int j=0;j<11;j++) {
    C[j] *= sqrt(2.0/M_PI);
    D[j] *= sqrt(2.0/M_PI);
  }
  return;
}
