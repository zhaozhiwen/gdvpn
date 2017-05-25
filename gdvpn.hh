
//  A struct containing all the kinematical data
struct kinedata {
  FourVector deuteron, gamma, vmeson, proton, neutron;
  int sd, sp, sn;
  bool pwia, singles, doubles, cauchy, vnvertex, pnvertex, soft, dsoft;
  kinedata() {
    pwia = false;
    singles = false;
    doubles = false;
    cauchy = false;
    vnvertex = false;
    pnvertex = false;
    soft = false;
    dsoft = false;
  }
};

// Mutator and accessor for the kinedata struct
void init_gdvpn(double,bool,bool,bool,bool,bool,bool,bool,bool);
kinedata get_gdvpn_info();

// The differential cross section
double diffsig5(double,double,double,double,double);
