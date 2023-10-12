#ifndef HMC_T_H_INCLUDED
#define HMC_T_H_INCLUDED

#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

mat svmnsim(vec beta, double mu, double phi, double sigma2, double y0);
vec mvrgaussian(int n);
vec gH_mom(mat inv_G, vec p);

// theta
// Model functions
double logpost_theta(vec theta, vec h, int T);
vec glogpost_theta(vec theta, vec h, mat inv_G, mat dG_mu, mat dG_omega, mat dG_gamma, int T);
// Metrics
// hmc updates
vec hmc_theta(vec theta_cur, vec h, int L, double eps, int T, int &acc, int &div, mat inv_M, int k, double prec);

// beta
// Model functions
double logpost_b(vec b, vec h, vec l, int T, vec y_T);
vec glogpost_b(vec b, vec h, vec l, mat inv_G, mat G_b0, mat dG_delta, mat dG_b2, int T, vec y_T);
// hmc update
vec hmc_b(vec b_cur, vec h, vec l, int L, double eps, int T, vec y_T, int &acc, int &div, mat inv_M, int k, double prec);

// v
// Model functions
double logpost_v(double e, vec l, int T, double alpha, double li, double ls);
double glogpost_v(double e, vec l, double inv_G, double dG_v, int T, double alpha, double li, double ls);
// hmc update
double hmc_v(double e_cur, vec l, int L, double eps, int T, int &acc, double alpha, double li, double ls, double inv_M, int &div, int k, double prec);

// l
vec l_gibbs(double e, vec y_T, vec h, vec b, int T, double alpha, double li, double ls);

// h
// Model functions
double logpost_h(vec h, vec theta, vec b, vec l, int T, vec y_T);
vec glogpost_h(vec h, vec theta, vec b, vec l, int T, vec y_T);
// hmc update
vec hmc_h(vec h, vec theta, vec b, vec l, int L, double eps, int T, vec y_T, int &acc);

#endif // HMC_T_H_INCLUDED
