//
//  Contains routines for dealing with two particle -> two particle collisions.
//  Created by A.J. Freese in April 2012.
//
#ifndef _kine2_hh
#define _kine2_hh
#include "fourvector.hh"

double kine2s(FourVector a, FourVector b) {
  FourVector c;
  double s;
  c=a+b;
  s=c*c;
  return s;
}

double kine2t(FourVector a, FourVector b) {
  double t;
  t=kine2s(a,-b);
  return t;
}

double kine2flux(FourVector Pa, FourVector Pb) {
  double temp;
  temp=(Pa*Pb)*(Pa*Pb)-Pa.m2()*Pb.m2();
  if(temp<0) temp=0;
  temp=2*sqrt(temp);
  return temp;
}

double kine2flux(double s, double ma2, double mb2) {
  double temp;
  temp=s-(ma2+mb2);
  temp*=temp;
  temp-=4*ma2*mb2;
  if(temp<0) temp=0;
  temp=sqrt(temp);
  return temp;
}

//  The following routine is given two initial four-momenta a and b, along with the Mandelstam t and the masses of the final particles, and determines what
// c and d are.
//  Kinematical flag will return a non-zero number if kinematics are bad.
//    1 : s is too small
//    2 : t is too small (too negative)
//    3 : t is too big (too positive)
//    4 : a+b is not timelike
void kine2collide(FourVector a, FourVector b, FourVector &c, FourVector &d, double mc, double md, double t, int &flag) {
  double s, pa, pc, Ea, Ec, ma2, mb2, phia, phic;
  double smin, tmin, tmax;
  ThreeVector N;
  flag=0;
  s=kine2s(a,b);
  //  First round of kinematical testing
  smin=mc+md;
  smin*=smin;
  if(s<smin) { flag=1; return; }
  //  Get magnitudes of energies and momenta for a and c in COM frame
  ma2=a.m2();
  mb2=b.m2();
  pa=s*s-2*s*(ma2+mb2)-4*ma2*mb2+(ma2+mb2)*(ma2+mb2);
  pa=sqrt(pa/(4*s));
  pc=s*s-2*s*(mc*mc+md*md)-4*mc*mc*md*md+(mc*mc+md*md)*(mc*mc+md*md);
  pc=sqrt(pc/(4*s));
  Ea=sqrt(ma2+pa*pa);
  Ec=sqrt(mc*mc+pc*pc);
  //  Second round of kinematical testing
  tmin=ma2+mc*mc-2*Ea*Ec-2*pa*pc;
  tmax=tmin+4*pa*pc;
  if(t<tmin) { flag=2; return; }
  if(t>tmax) { flag=3; return; }
  //  Next, get phic
  phic=t-ma2-mc*mc+2*Ea*Ec;
  phic=phic/(2*pa*pc);
  phic=acos(phic);
//  cout << pa << " " << pc << endl;
  //  Build c and d in COM frame with a forward along x axis
  c.E=Ec;
  c.p.setcl(pc,phic,0);
  d.build(md*md,-c.p);
  //  Rotate reference frame for a and b into one where the collision is in the xy plane
  N=a.p^b.p;
  if(N.r()>0) {
    a.p=a.p.rotate1z(N);
    b.p=b.p.rotate1z(N);
  }
  //  Boost to the COM frame (if it exists)
  FourVector P;
  ThreeVector Beta;
  P=a+b;
  Beta=P.beta();
  if(Beta.r()>=1) { flag=4; return; }
  a=a.boost(Beta);
  b=b.boost(Beta);
  phia=a.p.phi();
  //  Rotate c and d to current referene frame
  c.p=c.p.rotate1(3,-phia);
  d.p=d.p.rotate1(3,-phia);
  //  Boost back away from COM frame
  a=a.boost(-Beta);
  b=b.boost(-Beta);
  c=c.boost(-Beta);
  d=d.boost(-Beta);
  //  Rotate to original frame
  if(N.r()>0) {
    a.p=a.p.rotate1zi(N);
    b.p=b.p.rotate1zi(N);
    c.p=c.p.rotate1zi(N);
    d.p=d.p.rotate1zi(N);
  }
}

//  This routine finds the photon energy E, Mandelstam s, and Mandelstam t that will produce coherent photoproduction of C off of B.
void kine2threshold(double mb, double mc, double &s, double &t, double &E) {
  s=mb+mc;
  s*=s;
  E=s-mb*mb;
  E/=(2*mb);
  t=E/(mb+E);
  t=1-t*t;
  t=sqrt(t);
  t=1/t;
  t=2*mb*mb*(1-t);
}

//  These routines find the minimum (maximum) value of -t for which a reaction A+B->C+D can occur with a given s.
double kine2tminmax(double s, double ma2, double mb2, double mc2, double md2, bool min) {
  double pi2, pf2, ea, ec, t;
  pi2=ma2-mb2;
  pi2*=pi2;
  pi2=s*s-2*s*(ma2+mb2)+pi2;
  pi2/=4*s;
  ea=sqrt(ma2+pi2);
  pf2=mc2-md2;
  pf2*=pf2;
  pf2=s*s-2*s*(mc2+md2)+pf2;
  pf2/=4*s;
  ec=sqrt(mc2+pf2);
  t=ma2+mc2-2*ea*ec;
  if(min)
    t+=2*sqrt(pi2*pf2);
  if(!min)
    t-=2*sqrt(pi2*pf2);
  return t;
}

double kine2tmin(double s, double ma2, double mb2, double mc2, double md2) {
  return kine2tminmax(s,ma2,mb2,mc2,md2,true);
}

double kine2tmax(double s, double ma2, double mb2, double mc2, double md2) {
  return kine2tminmax(s,ma2,mb2,mc2,md2,false);
}

double kine2k2com(double s, double ma, double mb) {
  double k2;
  k2=s-(ma+mb)*(ma+mb);
  k2*=s-(ma-mb)*(ma-mb);
  k2/=4*s;
  return k2;
}

#endif













