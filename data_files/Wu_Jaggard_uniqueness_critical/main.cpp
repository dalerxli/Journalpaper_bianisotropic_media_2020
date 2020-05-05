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

double calculateKu(double epsilon_r, double zeta_c);
double mycubed(double x);

int main()
{
   double zeta_c;
   double epsilon_r;
   double mu_r;
//   cin >> zeta_c;
   cin >> epsilon_r;
   cin >> mu_r;

   zeta_c = 1e-2;
   double Ku = 10.0;
  
   while( (Ku > 1) && (zeta_c > 0))
   {
      Ku = calculateKu(epsilon_r, zeta_c);
      zeta_c -= 1.0e-6;
   }

  cout << scientific << epsilon_r << "\t" << zeta_c << "\t" << Ku << endl;
	
   return 0;
}

double calculateKu(double epsilon_r, double zeta_c)
{
  double d_gamma = 2.0*zeta_c/(eps0*epsilon_r);
  double d_chi = d_gamma;

  double C_ks = 2.0/(eps0*epsilon_r);
  double C_kd = 1.0/mycubed(eps0*epsilon_r);
  double C_kr = eps0*epsilon_r;

  double C_nus = 2.0*(1.0/mu0 + zeta_c*zeta_c/(eps0*epsilon_r));
  double C_nud = mycubed(1.0/mu0 + zeta_c*zeta_c/(eps0*epsilon_r));
  double C_nur = mu0*eps0*epsilon_r/(eps0*epsilon_r + mu0*zeta_c*zeta_c);

  double d_k = -C_ks+sqrt(C_ks*C_ks+4.0*C_kd*C_kr);
  double d_nu = -C_nus+sqrt(C_nus*C_nus+4.0*C_nud*C_nur);
  double regularity_factor = 4.0*d_gamma*d_chi/(d_k*d_nu);
 
  return regularity_factor;
}

double mycubed(double x)
{
   return x*x*x;
}
