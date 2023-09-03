// [[Rcpp::depends( RcppArmadillo )]]

#include "svm_smn_t.h"

// ############################## set seed function
void set_seed( int seed ){
  Rcpp::Environment base_env("package:base");
  Rcpp::Function set_seed_r = base_env["set.seed"];
  set_seed_r( std::floor( std::fabs( seed ) ) );
}

// [[Rcpp::export]]
List hmc_svm_t(int N, 
               int L_theta, vec eps_theta, mat M_theta,
               vec h,
               vec y_T, 
               int seed ){
  
  wall_clock timer;
  timer.tic();
  
  if( seed != 0 ) set_seed( seed );
  
  int T = y_T.n_elem - 1, a = floor( 0.1 * N );
  
  // iniciando theta
  int acc_theta = 0, div_theta = 0;
  //int L_theta;
  //vec eps_theta = zeros<vec>(3, 1);
  vec theta_cur = zeros<vec>(3, 1);
  theta_cur[ 0 ] += 1.0;
  theta_cur[ 1 ] += 0.5 * ( log( 1 + 0.98 ) - log( 1 - 0.98 ) );
  theta_cur[ 2 ] += log( sqrt( 0.017 ) );
  mat inv_M_theta = inv( M_theta );
  
  // iniciando cadeia
  mat chain_theta = zeros<mat>( 3, N + 1 );
  chain_theta.col( 0 ) += theta_cur;
  
  // chain builting  
  for(int it = 1 ; it < N + 1 ; it ++){
    
    //L_theta   = randi( distr_param(1, 20) );
    //eps_theta = randg( 3, distr_param( 10.0, 0.01 ) );
    
    theta_cur = hmc_theta( theta_cur, h, L_theta, eps_theta, T, acc_theta, 
                           div_theta, inv_M_theta );
    
    // chain update 
    chain_theta.col( it ) += theta_cur;
  
    //Progress
    if( (it % a) == 0 ) Rcout << "Progresso em " << ceil( 100 * it / N ) <<" %"<< endl;
  }
  
  // Transformations
  chain_theta.row( 1 ) = tanh( chain_theta.row( 1 ) );
  chain_theta.row( 2 ) = exp( chain_theta.row( 2 ) );
  
  vec acc = zeros<vec>(1, 1);
  acc[ 0 ] += acc_theta;

  double time = timer.toc();
  
  Rcout << "\nTempo decorrido: " << time / 60 << " min" << endl;
  Rcout << "\nTaxa aceitação:\n" << acc / N << endl;
  Rcout << "\n" << div_theta << " Divergências.\n" << endl;
  
  return List::create( Named("chain") = chain_theta, 
                       Named("acc") = acc, 
                       Named("time") = time 
  ); 
  
}
