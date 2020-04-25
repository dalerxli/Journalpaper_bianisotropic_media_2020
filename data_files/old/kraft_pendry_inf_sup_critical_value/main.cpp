//Author: Praveen K R
//Date: 02/November/2019

#include <iostream>
#include <cmath>
#include <complex>

using namespace std;

typedef complex<double> dcomplex;
#define j dcomplex(0.0, 1.0)

int main()
{
   double real_epsilon_r, imag_epsilon_r, zeta_0;
   cin >> real_epsilon_r  >> imag_epsilon_r;
   dcomplex epsilon_r = real_epsilon_r + j*imag_epsilon_r;

   double K1, K3, eta, alpha_opt;
   zeta_0 = -1.0;
   K1 = abs(imag_epsilon_r);
   while(zeta_0 < 0.0)
   {
      K3 = abs(real_epsilon_r) + zeta_0*zeta_0;
      double b = 1.0 -(K1*K1+K3*K3);
      alpha_opt = (b + sqrt(b*b + 4.0*K3*K3))/2.0;
      eta = sqrt((1.0-alpha_opt)/2.0);
      if(abs(zeta_0) - sqrt(eta) < 0)
      {
         cout << real(epsilon_r) << "\t" << imag(epsilon_r) << "\t" << abs(zeta_0) << endl;
         break;
      }
      zeta_0 += 0.0001;
   }
  
//   cout << alpha_opt << endl;
//   cout << real(epsilon_r) << "\t" << imag(epsilon_r) << "\t" << zeta_critical << endl;

   return 0;
}
