//Author: Praveen K R
//Date: 04/November/2019

#include <iostream>
#include <cmath>
#include <complex>

using namespace std;

typedef complex<double> dcomplex;
#define j dcomplex(0.0, 1.0)
const double pi = 3.1415926535897932;
const double C0 = 299792458.0;//speed of light in vacuum
const double mu0 = 4.0*pi*1.0e-7; 
const double eps0 = 1.0/(mu0*C0*C0);

double minimum(double x, double y);
double minimum(double x, double y, double z);
double maximum(double x, double y);
double maximum(double x, double y, double z);
double calculateKu(dcomplex epsilonr, double zeta0);

int main()
{
   double real_epsilon_r, imag_epsilon_r, zeta_0;
   cin >> real_epsilon_r  >> imag_epsilon_r;
   dcomplex epsilon_r = real_epsilon_r + j*imag_epsilon_r;

   zeta_0 = -1.0;
   double Ku = 1.0;
   while(zeta_0 < 0.0)
   {
      Ku = calculateKu(epsilon_r, zeta_0);
      if(Ku < 1.0)
      {
         cout << real(epsilon_r) << "\t" << imag(epsilon_r) << "\t" << abs(zeta_0) << "\t" << Ku << endl;
         break;
      }
      zeta_0 += 0.0001;
   }
  
   return 0;
}


double minimum(double x, double y)
{
   if(x < y)
   {
      return x;
   }
   else
   {
      return y;
   }
}

double minimum(double x, double y, double z)
{
   if(y < z)
   {
      return minimum(x,y);
   }
  else
  {
     return minimum(x,z);
  }
}

double maximum(double x, double y)
{
   if(x > y)
   {
      return x;
   }
   else
   {
      return y;
   }
}

double maximum(double x, double y, double z)
{
   if(y > z)
   {
      return maximum(x,y);
   }
   else 
   {
      return maximum(x,z);
   }
}


double calculateKu(dcomplex epsilonr, double zeta0)
{
   //Calculate the matrix entries
   dcomplex k11 = 1.0/(eps0*epsilonr);
   dcomplex k22 = 1.0/(eps0*epsilonr);
   dcomplex k33 = 1.0/(eps0*(epsilonr-zeta0*zeta0));

   dcomplex nu11 = 1.0/(mu0);
   dcomplex nu22 = epsilonr/(mu0*(epsilonr-zeta0*zeta0));
   dcomplex nu33 = 1.0/(mu0);
   
   dcomplex gamma23 = j*zeta0*C0/(epsilonr -zeta0*zeta0);
   dcomplex chi32 = -j*zeta0*C0/(epsilonr -zeta0*zeta0);

   //Calculate constants
   //Calculate Cgammas
   double Cgammas_inside = abs(gamma23);
   double Cgammas_outside = 0.0;
   double Cgammas = maximum(Cgammas_inside, Cgammas_outside);

   //Calculate Cchis
   double Cchis_inside = abs(chi32);
   double Cchis_outside = 0.0;
   double Cchis = maximum(Cchis_inside, Cchis_outside);

   //Calculate Cks
   double Cks_inside = abs(k11) + abs(k22) + abs(k33) - minimum(abs(k11), abs(k22), abs(k33));
   double Cks_outside = 2.0/eps0;
   double Cks = maximum(Cks_inside, Cks_outside);

   //Calculate Cnus
   double Cnus_inside = abs(nu11) + abs(nu22) + abs(nu33) - minimum(abs(nu11), abs(nu22), abs(nu33));
   double Cnus_outside = 2.0/mu0;
   double Cnus = maximum(Cnus_inside, Cnus_outside);

   //Calculate Ckd
   double Ckd_inside = abs(k11)*abs(k22)*abs(k33);
   double Ckd_outside = (1.0/eps0)*(1.0/eps0)*(1.0/eps0);
   double Ckd = minimum(Ckd_inside, Ckd_outside);

   //Calculate Cnud
   double Cnud_inside = abs(nu11)*abs(nu22)*abs(nu33);
   double Cnud_outside = (1.0/mu0)*(1.0/mu0)*(1.0/mu0);
   double Cnud = minimum(Cnud_inside, Cnud_outside);
 
   //Calculate Ckr
   double lambda_min_ks = minimum( 1.0/abs(real(k11)), 1.0/abs(real(k22)), 1.0/abs(real(k33)) ) ;
   double lambda_min_kss = minimum( 1.0/abs(imag(k11)), 1.0/abs(imag(k22)), 1.0/abs(imag(k33)) ) ;
   double Ckr_inside = sqrt(lambda_min_ks*lambda_min_ks + lambda_min_kss*lambda_min_kss);
   double Ckr_outside = eps0;
   double Ckr = minimum(Ckr_inside, Ckr_outside);

   //Calculate Cnur
   double lambda_min_nus = minimum( 1.0/abs(real(nu11)), 1.0/abs(real(nu22)), 1.0/abs(real(nu33)) ) ;
   double lambda_min_nuss = minimum( 1.0/abs(imag(nu11)), 1.0/abs(imag(nu22)), 1.0/abs(imag(nu33)) ) ;
   double Cnur_inside = sqrt(lambda_min_nus*lambda_min_nus + lambda_min_nuss*lambda_min_nuss);
   double Cnur_outside = mu0;
   double Cnur = minimum(Cnur_inside, Cnur_outside);

   //Calculate Ku
   double dk = -Cks + sqrt(Cks*Cks + 4.0*Ckd*Ckr);
   double dnu = -Cnus + sqrt(Cnus*Cnus + 4.0*Cnud*Cnur);
   double Ku = 4.0*Cchis*Cgammas/(dk*dnu);
   
   return Ku;
}
