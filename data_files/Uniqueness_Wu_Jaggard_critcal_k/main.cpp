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
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>

using namespace std;
typedef complex<double> dcomplex;
#define j dcomplex(0.0,1.0)
const double pi = 3.1415926535897932;
const double C0 = 299792458.0;//speed of light in vacuum
const double mu0 = 4.0*pi*1.0e-7; 
const double eps0 = 1.0/(mu0*C0*C0);
const double planckconstantbar = 1.054571800e-34;

double minimum2(double x1, double x2);
double mycubed(double x);

struct rparams
  {
    double epsilon_r;
    double mu_r;
  };

int rosenbrock_f (const gsl_vector * x, void *params, 
              gsl_vector * f)
{
   double epsilon_r = ((struct rparams *) params)->epsilon_r;
   double mu_r = ((struct rparams *) params)->mu_r;
   const double zeta_c = gsl_vector_get (x, 0);

   double d_gamma = 2.0*abs(zeta_c)/(eps0*epsilon_r);
   double d_chi = d_gamma;

   double C_ks = 2.0/(eps0*epsilon_r);
   double C_kd = 1.0/mycubed(eps0*epsilon_r);
   double C_kr = eps0*epsilon_r;

   double C_nus = 2.0*(1.0/mu0 + zeta_c*zeta_c/(eps0*epsilon_r));
   double C_nud = mycubed(1.0/mu0 + zeta_c*zeta_c/(eps0*epsilon_r));
   double C_nur = mu0*eps0*epsilon_r/(eps0*epsilon_r+mu0*zeta_c*zeta_c);

   double d_k = -C_ks+sqrt(C_ks*C_ks+4.0*C_kd*C_kr);
   double d_nu = -C_nus+sqrt(C_nus*C_nus+4.0*C_nud*C_nur);
   double regularity_factor = 4.0*d_gamma*d_chi/(d_k*d_nu);  



  const double y0 = regularity_factor-1;

  gsl_vector_set (f, 0, y0);

  return GSL_SUCCESS;
}

int main()
{
   double epsilon_r;
   double mu_r=1;
   
   cin >> epsilon_r;
   //cin >> mu_r;


  const gsl_multiroot_fsolver_type *T;
  gsl_multiroot_fsolver *s;

  int status;
  size_t i, iter = 0;

  const size_t n = 1;
  struct rparams p = {epsilon_r, mu_r};
  gsl_multiroot_function f = {&rosenbrock_f, n, &p};

  double x_init[1] = {0.5};
  gsl_vector *x = gsl_vector_alloc (n);
  gsl_vector_set (x, 0, x_init[0]);

  T = gsl_multiroot_fsolver_hybrids;
  s = gsl_multiroot_fsolver_alloc (T, 1);
  gsl_multiroot_fsolver_set (s, &f, x);

  //print_state (iter, s);

  do
    {
      iter++;
      status = gsl_multiroot_fsolver_iterate (s);

      //print_state (iter, s);

      if (status)   /* check if solver is stuck */
        break;

      status = 
        gsl_multiroot_test_residual (s->f, 1e-7);
    }
  while (status == GSL_CONTINUE && iter < 20000);

  cout << scientific << epsilon_r << "\t" <<  abs(gsl_vector_get (s->x, 0)) << "\t" << gsl_vector_get (s->f, 0)<< endl;

  gsl_multiroot_fsolver_free (s);
  gsl_vector_free (x);

   return 0;
}

double minimum2(double x1, double x2)
{
   if(x1 < x2)
   {
      return x1;
   }
   else return x2;

}

double mycubed(double x)
{
   return x*x*x;
}
