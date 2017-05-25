//
//  This program computes the differential cross-section of the reaction
//    d + gamma --> V + p + n
// for a vector meson V.
//
//  Created by A.J. Freese in April 2012, last updated May 2017.
//  This program uses Cuba numerical integration routines, cf.
//     <http://www.feynarts.de/cuba/>
// Thus -lcuba should be included as a linker flag during build.
//
//  The program takes user input. The following letters are boolean flags:
//      b : PWIA turned on
//      s : single resc. turned on (will do nothing without v, n or f)
//      d : double resc. turned on (does two diagrams with no additional flags)
//      v : vn vertex turned on
//      n : pn vertex turned on
//      c : Cauchy PV correction to v and n turned on
//      f : soft single resc. turned on (requires s)
//      z : soft double resc. turned on (requires d)
//  The following letters should be followed by user input:
//      E : photon energy (GeV); default is off-proton threshold
//      p : final neutron momentum (GeV); default is 0.1
//      t : -t (GeV^2); default is off-proton threshold 
//      V : integer for which vector meson (0 jpsi; 1 phi)
//      P : integer for photoproduction parameterization
//      S : total inelastic VN cross-section
//      B : slope parameter for VN elastic scattering
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
#include "constants.hh"
#include "kine2.hh"
#include "dwf.hh"
#include "nnscatter.hh"
#include "vmeson_params.hh"
#include "gdvpn.hh"
using namespace std;

// =============================================================================
// Forward declarations

// Makes datafile
void diffsig5data(double,double,double,char*);

// Screen output
time_t start(double,double,char*);
void end(time_t);
void progress_bar(double);

// =============================================================================
//  main() routine

int main(int argc, char **argv) {
//  kinedata gdvpn;
  char *datafilename="test";
  //  Kinematic variables and default values
  //  Impossible values for E and t are later corrected to threshold values
  double E=-1, p=0.1, t=1, theta=0, phi=0;
  double Ethr, sthr, tthr;
  // Whether various contributions are turned on or not
  bool pwia = false, singles = false, doubles = false, vnvertex = false,
       pnvertex = false, soft = false, dsoft = false, cauchy = false;
  // Meson properties and their defaults. Negative values will effect defaults.
  // meson ID, param ID, sigma_vn, and B_vn; see vmeson_params.cc for details.
  int Vpid=0, Vprm=2; double Vsig=-1; double Vslp=-1;
  //  Take user parameters.
  int c;
  while((c = getopt(argc,argv,"bsdvncfzE:p:t:V:P:S:B:o:")) != -1) {
    switch (c) {
      case 'b': pwia=true; break;
      case 's': singles=true; break;
      case 'd': doubles=true; break;
      case 'v': vnvertex=true; break;
      case 'n': pnvertex=true; break;
      case 'f': soft=true; break;
      case 'z': dsoft=true; break;
      case 'c': cauchy=true; break;
      case 'E': E = atof(optarg); break;
      case 'p': p = atof(optarg); break;
      case 't': t = atof(optarg); break;
      case 'V': Vpid = atoi(optarg); break;
      case 'P': Vprm = atoi(optarg); break;
      case 'S': Vsig = atof(optarg); break;
      case 'B': Vslp = atof(optarg); break;
      case 'o': datafilename = optarg; break;
      default: break;
    }
  }
  //  Prepare vector meson & get threshold kinematics
  vnset(Vpid,Vprm,Vsig,Vslp);
  kine2threshold(ProtonMass,vmass(),sthr,tthr,Ethr);
  //  If absurd E and t not overridden, use threshold values
  if(E<0.) E = Ethr;
  if(t>0.) t = tthr;
  //  Pass desired properties into gdvpn module
  init_gdvpn(E,pwia,singles,doubles,vnvertex,pnvertex,soft,dsoft,cauchy);
  diffsig5data(t,p,phi,datafilename);
  return 0;
}

// =============================================================================
// Method for making data file

void diffsig5data(double t, double p_n, double phi_nl, char *datafilename) {
  time_t chronos;
  ofstream datafile;
  int flag, maxdeg, degstep;
  const double degprad = 180./M_PI;
  const char *tab="\t";
  double diffsig, theta_nl;
  kinedata kineinfo;
  datafile.open(datafilename);
  datafile.precision(8);
  maxdeg = 180;
  degstep = 1;
  chronos = start(t,p_n,datafilename);
  progress_bar(0);
  //  Get the cross-section & write to file if kinematically allowed
  for(int itheta=0;itheta<=maxdeg;itheta+=degstep) {
    theta_nl = (double)itheta/degprad;
    diffsig = diffsig5(t,p_n,theta_nl,phi_nl,0.);
    if(diffsig>=0.) {
      kineinfo = get_gdvpn_info();
      datafile << scientific
         << theta_nl*degprad 
         << tab << kineinfo.gamma.E 
         << tab << p_n 
         << tab << (-t)
         << tab << diffsig
         << tab << kineinfo.vmeson.p.r() 
         << tab << kineinfo.vmeson.p.theta()*degprad 
         << tab << kineinfo.vmeson.p.phi()*degprad
         << tab << kineinfo.proton.p.r() 
         << tab << kineinfo.proton.p.theta()*degprad 
         << tab << kineinfo.proton.p.phi()*degprad
         << tab << kineinfo.neutron.p.r() 
         << tab << kineinfo.neutron.p.theta()*degprad 
         << tab << kineinfo.neutron.p.phi()*degprad
         << endl;
    }
    progress_bar(((double)itheta)/maxdeg);
  }
  cout << endl;
  datafile.close();
  end(chronos);
  return;
}

// =============================================================================
// Routines for Screen Output

// Outputs info for run to screen
time_t start(double t, double p, char *filename) {
  bool stated_term = false;
  kinedata kineinfo = get_gdvpn_info();
  cout << "Starting gdvpn, with:" << endl;
  if(kineinfo.pwia)
    cout << "\tPWIA";
  if(kineinfo.singles) {
    cout << "\tSingle resc. ( ";
    if(kineinfo.vnvertex)
      cout << "VN ";
    if(kineinfo.pnvertex)
      cout << "pn ";
    if(kineinfo.soft)
      cout << "soft ";
    cout << ")";
  }
  if(kineinfo.cauchy)
    cout << "\tCauchy PV";
  if(kineinfo.doubles)
    cout << "\tDouble resc.";
  cout << endl;
  cout << "\tE = " << kineinfo.gamma.E << " GeV\tp = " 
    << p << " GeV\tt = " << t << " GeV**2" << endl;
  cout << "\tData being output to " << filename << endl;
  return time(NULL);
}

//  Outputs run time to screen
void end(time_t alpha) {
  time_t omega = time(NULL);
  double runtime = (double)(omega-alpha);
  cout << "Finished. Took " << runtime << " second(s)." << endl;
}

//  Outputs progress to screen
void progress_bar(double completed) {
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
