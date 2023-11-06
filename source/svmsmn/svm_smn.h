#ifndef SVM_SMN_H_INCLUDED
#define SVM_SMN_H_INCLUDED

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
int G_theta(vec theta, mat &G, mat &inv_G, mat &dG_mu, mat &dG_omega, mat &dG_gamma, int T);
// Energy functions
vec nuH_theta(mat dG_mu, mat dG_omega, mat dG_gamma, mat inv_G, vec p);
//vec gradHmom_theta(mat inv_G, vec p);
// hmc updates
vec hmc_theta(vec theta_cur, vec h, int L, vec eps, int T, int &acc);
vec rmhmc_theta(vec theta_cur, vec h, int fixp, int L, vec eps, int T, int &acc);

// beta
// Model functions
double logpost_b(vec b, vec h, vec l, int T, vec y_T);
vec glogpost_b(vec b, vec h, vec l, mat inv_G, mat G_b0, mat dG_delta, mat dG_b2, int T, vec y_T);
// Metrics
int G_b(vec b, vec h, vec l, mat &G, mat &inv_G, mat& dG_b0, mat &dG_delta, mat &dG_b2, int T, vec y_T);
// Energy functions
vec nuH_b(mat dG_b0, mat dG_delta, mat dG_b2, mat inv_G, vec p);
//vec gradHmom_b(mat inv_G, vec p);
// hmc update
vec rmhmc_b(vec b_cur, vec h, vec l, int fixp, int L, vec eps, int T, vec y_T , int &acc);

// h
// Model functions
double logpost_h(vec h, vec theta, vec b, vec l, int T, vec y_T);
vec glogpost_h(vec h, vec theta, vec b, vec l, int T, vec y_T);
// hmc update
vec hmc_h(vec h, vec theta, vec b, vec l, int L, double eps, int T, vec y_T, int &acc);

// -------------------------- SVM - t

// v
// Model functions
double logpost_v_t(double e, vec l, int T, double alpha, double li, double ls);
double glogpost_v_t(double e, vec l, double inv_G, double dG_v, int T, double alpha, double li, double ls);
// Metrics
int G_v_t(vec e, double &G, double &inv_G, double &dG_v, int T, double alpha, double li, double ls);
// Energy functions
double nuH_v_t(double dG_v, double inv_G, double p);
double gradHmom_v_t(double inv_G, double p);
// hmc update
double rmhmc_v_t(double e_cur, vec l, int fixp, int L, double eps, int T, int &acc, double alpha, double li, double ls );

// l
vec l_gibbs_t(double e, vec y_T, vec h, vec b, int T, double alpha, double li, double ls);

// -------------------------- SVM - VG

// v
// Model functions
double logpost_v_vg(double e, vec l, int T, double alpha, double li, double ls);
double glogpost_v_vg(double e, vec l, double inv_G, double dG_v, int T, double alpha, double li, double ls);
// Metrics
int G_v_vg(vec e, double &G, double &inv_G, double &dG_v, int T, double alpha, double li, double ls);
// Energy functions
double nuH_v_vg(double dG_v, double inv_G, double p);
double gradHmom_v_vg(double inv_G, double p);
// hmc update
double rmhmc_v_vg(double e_cur, vec l, int fixp, int L, double eps, int T, int &acc, double alpha, double li, double ls );

// -------------------------- SVM - S

// v
// Model functions
double logpost_v_s(double e, vec l, int T, double alpha, double li, double ls);
double glogpost_v_s(double e, vec l, double inv_G, double dG_v, int T, double alpha, double li, double ls);
// Metrics
int G_v_s(vec e, double &G, double &inv_G, double &dG_v, int T, double alpha, double li, double ls);
// Energy functions
double nuH_v_s(double dG_v, double inv_G, double p);
double gradHmom_v_s(double inv_G, double p);
// hmc update
double rmhmc_v_s(double e_cur, vec l, int fixp, int L, double eps, int T, int &acc, double alpha, double li, double ls );


#endif // SVM_SMN_H_INCLUDED
