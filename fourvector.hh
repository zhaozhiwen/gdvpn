//
//  Contains ThreeVector and FourVector classes.
//  Created by A.J. Freese in April 2012.
//
#ifndef _fourvector_hh
#define _fourvector_hh
#include <cmath>

//  Vector class.  Note that for numbered axes, 1=x, 2=y, 3=z
class ThreeVector {
public:
  //  Cartesian components
  double x,y,z;
  //  Setting routines
  void set(double,double,double);       //sets Cartesian components
  void setsp(double,double,double);     //sets spherical components
  void setcl(double,double,double);     //sets cylindrical components
  //  Routines for obtaining quantities
  double r();                 //gets modulus
  double rho();               //gets modulus of projection into xy plane
  double theta();
  double phi();
  //  ThreeVector operators
  void operator=(ThreeVector);
  ThreeVector operator+(ThreeVector);
  ThreeVector operator-(ThreeVector);
  ThreeVector operator-();
  ThreeVector operator*(double);        //scalar multiplication
  ThreeVector operator/(double);        //scalar division
  double operator*(ThreeVector);        //dot product
  double operator|(ThreeVector);        //obtians cosine of relative angle
  ThreeVector operator^(ThreeVector);     //cross product
  //  ThreeVector transformations; int denotes axis
  ThreeVector translate(int,double);
  ThreeVector rotate1(int,double);
  ThreeVector rotate3(double,double,double);  //rotates a Vector by Euler angles
  ThreeVector rotate3i(double,double,double); //inverse of Euler rotation
  ThreeVector rotate1x(ThreeVector,double); //rotates around a ThreeVector
  ThreeVector rotate1z(ThreeVector);      //rotates to frame where arg points up
  ThreeVector rotate1zi(ThreeVector);     //undoes what rotate1z does
};

//  FourVector class.  Consists of a ThreeVector and a time-like scalar component.
//  Convention is m^2 > 0 for time-like four-vectors.
class FourVector {
public:
  //  Components
  ThreeVector p;
  double E;
  //  Setting routines
  void set(double,double,double,double);    //first arg is time-like, others Cartesian
  void set(double,ThreeVector);       //first arg is time-like
  void build(double,ThreeVector);       //first arg is m^2
  //  Routines for obtaining quantities
  double m2();
  ThreeVector beta();
  //  FourVector operators
  void operator=(FourVector); 
  FourVector operator+(FourVector);
  FourVector operator-(FourVector);
  FourVector operator-();
  double operator*(FourVector);
  //  FourVector transformations
  FourVector boost(ThreeVector);
};

// ---------------------------------------------------------------------------------------
//  Routines for the Vector class.
// ---------------------------------------------------------------------------------------
//  Routines for setting
void ThreeVector::set(double x, double y, double z) {
  this->x=x;
  this->y=y;
  this->z=z;
}
void ThreeVector::setsp(double r, double theta, double phi) {
  double x,y,z;
  z=r*cos(theta);
  x=r*sin(theta)*cos(phi);
  y=r*sin(theta)*sin(phi);
  this->set(x,y,z);
}
void ThreeVector::setcl(double rho, double phi, double z) {
  double x,y;
  x=rho*cos(phi);
  y=rho*sin(phi);
  this->set(x,y,z);
}
// ---------------------------------------------------------------------------------------
//  Routines for obtaining values
double ThreeVector::r() {
  double temp;
  temp=this->x*this->x+this->y*this->y+this->z*this->z;
  temp=sqrt(temp);
  return temp;
}
double ThreeVector::rho() {
  double temp;
  temp=this->x*this->x+this->y*this->y;
  temp=sqrt(temp);
  return temp;
}
double ThreeVector::theta() {
  double temp;
  if(this->r()==0) return 0;
  temp=this->z;
  temp=temp/this->r();
  temp=acos(temp);
  return temp;
}
double ThreeVector::phi() {
  double temp;
  if ((this->y)==0) {
    if(this->x>0) return 0;
    if(this->x<0) return M_PI;
    if(this->x==0) return 0;
  }
  temp=this->y/this->x;
  temp=atan(temp);
  if ((this->x)<0) temp+=M_PI;
  if (temp<0) temp+=2*M_PI;
  return temp;
}
//----------------------------------------------------------------------------------------
//  Routines for overloading operators
void ThreeVector::operator=(ThreeVector that) {
  this->x=that.x;
  this->y=that.y;
  this->z=that.z;
}
ThreeVector ThreeVector::operator+(ThreeVector that) {
  ThreeVector temp;
  temp.x=this->x+that.x;
  temp.y=this->y+that.y;
  temp.z=this->z+that.z;
  return temp;
}
ThreeVector ThreeVector::operator-(ThreeVector that) {
  ThreeVector temp;
  temp.x=this->x-that.x;
  temp.y=this->y-that.y;
  temp.z=this->z-that.z;
  return temp;
}
ThreeVector ThreeVector::operator-() {
  ThreeVector temp;
  temp.x=-this->x;
  temp.y=-this->y;
  temp.z=-this->z;
  return temp;
}
ThreeVector ThreeVector::operator*(double s) {
  ThreeVector temp;
  temp.set(s*this->x,s*this->y,s*this->z);
  return temp;
}
ThreeVector ThreeVector::operator/(double s) {
  ThreeVector temp;
  temp.set(this->x/s,this->y/s,this->z/s);
  return temp;
}
double ThreeVector::operator*(ThreeVector that) {
  double dot;
  dot=this->x*that.x+this->y*that.y+this->z*that.z;
  return dot;
}
double ThreeVector::operator|(ThreeVector that) {
  double dot;
  ThreeVector temp;
  temp.set(this->x,this->y,this->z);
  temp=temp/this->r();
  temp=temp/that.r();
  dot=temp*that;
  return dot;
}
ThreeVector ThreeVector::operator^(ThreeVector that) {
  ThreeVector temp;
  temp.z=this->x*that.y-that.x*this->y;
  temp.x=this->y*that.z-that.y*this->z;
  temp.y=this->z*that.x-that.z*this->x;
  return temp;
}
//----------------------------------------------------------------------------------------
//  Routines for transforming ThreeVectors
//  Note that for these transformations, the object itself is not transformed.
// Each routine returns a ThreeVector. If it is desired that a ThreeVector 'vec' be
// transformed, then
//    vec=vec.trans(3,a);
// is what should be used, for instance.  Also, note that arguments are parameters
// representing the change of reference frame so translation by 'a' will actually subtract
// 'a' from a Vector component rather than add it.
ThreeVector ThreeVector::translate(int axis, double a) {
  ThreeVector temp;
  temp.set(this->x,this->y,this->z);
  switch(axis) {
  case 1:
    temp.x-=a;
    break;
  case 2:
    temp.y-=a;
    break;
  case 3:
    temp.z-=a;
    break;
  }
  return temp;
}
ThreeVector ThreeVector::rotate1(int axis, double theta) {
  ThreeVector temp;
  switch(axis) {
  case 1:
    temp.x=this->x;
    temp.y=cos(theta)*this->y+sin(theta)*this->z;
    temp.z=cos(theta)*this->z-sin(theta)*this->y;
    break;
  case 2:
    temp.y=this->y;
    temp.z=cos(theta)*this->z+sin(theta)*this->x;
    temp.x=cos(theta)*this->x-sin(theta)*this->z;
    break;
  case 3:
    temp.z=this->z;
    temp.x=cos(theta)*this->x+sin(theta)*this->y;
    temp.y=cos(theta)*this->y-sin(theta)*this->x;
    break;
  default:
    temp.set(this->x,this->y,this->z);
  }
  return temp;
}
ThreeVector ThreeVector::rotate3(double a, double b, double c) {
  ThreeVector temp;
  temp.set(this->x,this->y,this->z);
  temp=temp.rotate1(3,a);
  temp=temp.rotate1(2,b);
  temp=temp.rotate1(3,c);
  return temp;
}
ThreeVector ThreeVector::rotate3i(double a, double b, double c) {
  ThreeVector temp;
  temp.set(this->x,this->y,this->z);
  temp=temp.rotate1(3,-c);
  temp=temp.rotate1(2,-b);
  temp=temp.rotate1(3,-a);
  return temp;
}
ThreeVector ThreeVector::rotate1x(ThreeVector axis, double theta) {
  //  This is done in three steps.
  double col, azi;
  ThreeVector temp;
  temp.set(this->x,this->y,this->z);
  col=axis.theta();
  azi=axis.phi();
  //  First, we rotate to a frame where axis is pointing up.
  temp=temp.rotate3(azi,col,0);
  //  Next, we rotate by the desired angle about the current z-axis.
  temp=temp.rotate1(3,theta);
  //  Then, we undo the transformation to the frame where axis is up.
  temp=temp.rotate3i(azi,col,0);
  return temp;
}
ThreeVector ThreeVector::rotate1z(ThreeVector that) {
  ThreeVector temp;
  double col, azi;
  temp.set(this->x,this->y,this->z);
  col=that.theta();
  azi=that.phi();
  temp=temp.rotate3(azi,col,0);
  return temp;
}
ThreeVector ThreeVector::rotate1zi(ThreeVector that) {
  ThreeVector temp;
  double col, azi;
  temp.set(this->x,this->y,this->z);
  col=that.theta();
  azi=that.phi();
  temp=temp.rotate3i(azi,col,0);
  return temp;
}

// ---------------------------------------------------------------------------------------
//  Routines for the FourVector class.
// ---------------------------------------------------------------------------------------
//  Routines for setting
void FourVector::set(double t, double x, double y, double z) {
  this->p.x=x;
  this->p.y=y;
  this->p.z=z;
  this->E=t;
}
void FourVector::set(double t, ThreeVector x) {
  this->p.set(x.x,x.y,x.z);
  this->E=t;
}
void FourVector::build(double s2, ThreeVector x) {
  double t,r;
  r=x.r();
  t=s2+r*r;
  if (t<0) t=0;
  t=sqrt(t);
  this->set(t,x);
}
// ---------------------------------------------------------------------------------------
//  Routines for obtaining values
double FourVector::m2() {
  double s2;
  s2=this->E*this->E;
  s2-=this->p.r()*this->p.r();
  return s2;
}
ThreeVector FourVector::beta() {
  ThreeVector v;
  v=this->p;
  v=v/this->E;
  return v;
}
// ---------------------------------------------------------------------------------------
//  Routines for overloading operators
void FourVector::operator=(FourVector that) {
  this->p=that.p;
  this->E=that.E;
}
FourVector FourVector::operator+(FourVector that) {
  FourVector temp;
  temp.E=this->E+that.E;
  temp.p=this->p+that.p;
  return temp;
}
FourVector FourVector::operator-(FourVector that) {
  FourVector temp;
  temp.E=this->E-that.E;
  temp.p=this->p-that.p;
  return temp;
}
FourVector FourVector::operator-() {
  FourVector temp;
  temp.E=-this->E;
  temp.p=-this->p;
  return temp;
}
double FourVector::operator*(FourVector that) {
  double s;
  s=this->E*that.E;
  s-=this->p*that.p;
  return s;
}
// ---------------------------------------------------------------------------------------
//  Routines for transforming FourVectors
//  Usage notes the same as with Vectors
FourVector FourVector::boost(ThreeVector beta) {
  double b,g,pz,E;
  FourVector temp;
  temp.set(this->E,this->p);
  b=beta.r();
  if(b<0||b>=1) return temp;    //doesn't like beta out of range
  g=1-b*b;
  g=sqrt(g);
  g=1/g;
  temp.p=temp.p.rotate1z(beta);
  E=temp.E;
  pz=temp.p.z;
  temp.p.z=g*pz-b*g*E;      //boost along z axis
  temp.E=g*E-b*g*pz;
  temp.p=temp.p.rotate1zi(beta);
  return temp;
}

#endif
