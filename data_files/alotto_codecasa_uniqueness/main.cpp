/*This program calculates the amplitudes of the reflected and transmitted wave when 
a multi layered dielectric is illuminated by an obliquely incident transverse electric source.
*/

#include <iostream>
#include <fstream>
#include <iomanip>
#include <complex>
#include <cmath>
#include <cassert>
#include <cstring>
#include <vector>

using namespace std;
typedef complex<double> dcomplex;
#define j dcomplex(0.0,1.0)
const double pi = 3.1415926535897932;
const double C0 = 299792458.0;//speed of light in vacuum
const double mu0 = 4.0*pi*1.0e-7; 
const double eps0 = 1.0/(mu0*C0*C0);
const double planckconstantbar = 1.054571800e-34;

double mycubed(double x){
   return x*x*x;
}

int main()
{
   double kappa_0;
   double epsilon_r;
   double mu_r;
   cin >> kappa_0;
   cin >> epsilon_r;
   cin >> mu_r;
  
   double d_gamma = 4.0*kappa_0/(eps0);
   double d_chi = d_gamma;

  double C_ks = 2.0/(eps0);
  double C_kd = 1.0/mycubed(eps0);
  double C_kr = eps0;

  double C_nus = 2.0*(1.0/mu0 + 4*kappa_0*kappa_0/(eps0));
  double C_nud =(1.0/mu0)*(1.0/mu0 + kappa_0*kappa_0/(eps0))*(1.0/mu0 + 4.0*kappa_0*kappa_0/(eps0));
  double C_nur = mu0*eps0/(eps0 + 4.0*mu0*kappa_0*kappa_0);

   double d_k = -C_ks+sqrt(C_ks*C_ks+4.0*C_kd*C_kr);
   double d_nu = -C_nus+sqrt(C_nus*C_nus+4.0*C_nud*C_nur);
   double regularity_factor = 4.0*d_gamma*d_chi/(d_k*d_nu);

   cout << scientific << kappa_0 << "\t" << regularity_factor << endl;
	
   return 0;
}
