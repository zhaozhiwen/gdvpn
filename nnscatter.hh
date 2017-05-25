#ifndef _nnscatter_hh
#define _nnscatter_hh

#include <complex>
using namespace std;

void init_nnscatter();

complex<double> pnscatter(double,double);
complex<double> ppscatter(double,double);
complex<double> nnscatter(double,double);

#endif
