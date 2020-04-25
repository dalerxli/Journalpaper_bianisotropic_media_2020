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
   double Cgammas = abs(gamma23);

   //Calculate Cchis
   double Cchis = abs(chi32);

   //Calculate Cks
   double Cks = abs(k11) + abs(k22) + abs(k33) - minimum(abs(k11), abs(k22), abs(k33));

   //Calculate Cnus
   double Cnus = abs(nu11) + abs(nu22) + abs(nu33) - minimum(abs(nu11), abs(nu22), abs(nu33));

   //Calculate Ckd
   double Ckd = abs(k11)*abs(k22)*abs(k33);

   //Calculate Cnud
   double Cnud = abs(nu11)*abs(nu22)*abs(nu33);

   //Calculate Ckr
   double lambda_min_ks = minimum( abs(real(1.0/k11)), abs(real(1.0/k22)), abs(real(1.0/k33)) ) ;
   double lambda_min_kss = minimum( abs(imag(1.0/k11)), abs(imag(1.0/k22)), abs(imag(1.0/k33)) ) ;
   double Ckr = sqrt(lambda_min_ks*lambda_min_ks + lambda_min_kss*lambda_min_kss);

   //Calculate Cnur
   double lambda_min_nus = minimum( abs(real(1.0/nu11)), abs(real(1.0/nu22)), abs(real(1.0/nu33)) ) ;
   double lambda_min_nuss = minimum( abs(imag(1.0/nu11)), abs(imag(1.0/nu22)), abs(imag(1.0/nu33)) ) ;
   double Cnur = sqrt(lambda_min_nus*lambda_min_nus + lambda_min_nuss*lambda_min_nuss);

   //Calculate Ku
   double dk = -Cks + sqrt(Cks*Cks + 4.0*Ckd*Ckr);
   double dnu = -Cnus + sqrt(Cnus*Cnus + 4.0*Cnud*Cnur);
   double Ku = 4.0*Cchis*Cgammas/(dk*dnu);
   
   return Ku;
}
