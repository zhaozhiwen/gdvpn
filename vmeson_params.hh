#ifndef _vmeson_params_hh
#define _vmeson_params_hh

#include <complex>
using namespace std;

double vmass();
double vmass2();
void vnset(int,int,double,double);
complex<double> vnphoto(double,double);
complex<double> vnscatter(double,double);

#endif