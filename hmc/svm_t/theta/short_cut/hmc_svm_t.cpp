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
               int L_theta, double eps_theta, mat M_theta, int k_theta, double prec_theta,
               int L_b, double eps_b, mat M_b, int k_b, double prec_b,
               int L_h, double eps_h,
               int L_v, double eps_v, double M_v, int k_v, double prec_v,
               double alpha, double li, double ls,
               vec y_T, 
               bool eps_random,
               int seed ){
  
  wall_clock timer;
  timer.tic();
  
  if( seed != 0 ) set_seed( seed );
  
  int T = y_T.n_elem - 1, a = floor( 0.1 * N );
  
  // iniciando theta
  int acc_theta = 0, div_theta = 0;
  vec theta_cur = zeros<vec>(3, 1);
  theta_cur[ 0 ] += 1.0;
  theta_cur[ 1 ] += 0.5 * ( log( 1 + 0.98 ) - log( 1 - 0.98 ) );
  theta_cur[ 2 ] += log( sqrt( 0.03 ) );
  mat inv_M_theta = inv( M_theta );
  
  // iniciando b
  int acc_b = 0, div_b = 0;
  vec b_cur = zeros<vec>(3, 1);
  b_cur[ 0 ] += 0.1;
  b_cur[ 1 ] += 0.5 * ( log( 1 + 0.01 ) - log( 1 - 0.01 ) );
  b_cur[ 2 ] += -0.05;
  mat inv_M_b = inv( M_b );
  
  // iniciando h
  int acc_h = 0;
  vec h_cur = zeros<vec>(T, 1);
  h_cur[ 0 ] += 1.0 + sqrt( 0.03 ) / (1 - 0.98 * 0.98 ) * randn();
  for( int kt = 1 ; kt < T ; kt++ ){
    h_cur[ kt ] += 1.0 + 0.98 * ( h_cur[ kt - 1 ] - 1.0 ) + sqrt( 0.03 ) * randn();
  }
  
  // iniciando v
  int acc_v = 0, div_v = 0;
  double v_cur = 20;
  //double v_cur = randg( distr_param(2.0, 1 / 0.1) );
  
  // iniciando l
  vec l_cur = zeros<vec>(T, 1);
  for( int k = 0 ; k < T ; k++ ){
    //randg( distr_param(a,b) )
    l_cur[ k ] = randg( distr_param( 0.5 * v_cur, 2.0 / v_cur ) );
  }
  
  v_cur = ( 2 / alpha ) * atanh( (2 * v_cur - ls - li) / (ls - li) );
  
  // iniciando cadeia
  mat chain_theta = zeros<mat>( 3, N + 1 );
  chain_theta.col( 0 ) += theta_cur;
  
  mat chain_b = zeros<mat>( 3, N + 1 );
  chain_b.col( 0 ) += b_cur;
  
  mat chain_h = zeros<mat>( T, N + 1 );
  chain_h.col( 0 ) += h_cur;
  
  mat chain_l = zeros<mat>( T, N + 1 );
  chain_l.col( 0 ) += l_cur;
  
  mat chain_v = zeros<mat>( 1, N + 1 );
  chain_v.col( 0 ) += v_cur;
  
  double reps_theta = 0.0;
  double reps_b = 0.0;
  double reps_v = 0.0;
  
  // chain builting  
  for(int it = 1 ; it < N + 1 ; it ++){
    
    if( eps_random ){
      reps_theta = randg( distr_param( 1.0, 2e-2 ) );
      reps_b     = randg( distr_param( 100.0, 5e-5 ) );
      reps_v     = randg( distr_param( 20.0, 0.005 ) );
      
      theta_cur = hmc_theta( theta_cur, h_cur, L_theta, reps_theta, T, acc_theta, div_theta, inv_M_theta, k_theta, prec_theta );
      b_cur = hmc_b( b_cur, h_cur, l_cur, L_b, reps_b, T, y_T , acc_b, div_b, inv_M_b, k_b, prec_b );
      h_cur = hmc_h( h_cur, theta_cur, b_cur, l_cur, L_h, eps_h, T, y_T, acc_h );
      v_cur = hmc_v(v_cur, l_cur, L_v, reps_v, T, acc_v, alpha, li, ls, 1 / M_v, div_v, k_v, prec_v);
      l_cur = l_gibbs( v_cur, y_T, h_cur, b_cur, T, alpha, li, ls );
    }else{
      theta_cur = hmc_theta( theta_cur, h_cur, L_theta, eps_theta, T, acc_theta, div_theta, inv_M_theta, k_theta, prec_theta );
      b_cur = hmc_b( b_cur, h_cur, l_cur, L_b, eps_b, T, y_T , acc_b, div_b, inv_M_b, k_b, prec_b );
      h_cur = hmc_h( h_cur, theta_cur, b_cur, l_cur, L_h, eps_h, T, y_T, acc_h );
      v_cur = hmc_v(v_cur, l_cur, L_v, eps_v, T, acc_v, alpha, li, ls, 1 / M_v, div_v, k_v, prec_v);
      l_cur = l_gibbs( v_cur, y_T, h_cur, b_cur, T, alpha, li, ls );
    }
    
    // chain update 
    chain_theta.col( it ) += theta_cur;
    chain_b.col( it ) += b_cur;
    chain_h.col( it ) += h_cur;
    chain_v.col( it ) += v_cur;
    chain_l.col( it ) += l_cur;
    
    //Progress
    if( (it % a) == 0 ) Rcout << "Progresso em " << ceil( 100 * it / N ) <<" %"<< endl;
  }
  
  // Transformations
  chain_theta.row( 1 ) = tanh( chain_theta.row( 1 ) );
  chain_theta.row( 2 ) = exp( chain_theta.row( 2 ) );
  chain_b.row( 1 )     = tanh( chain_b.row( 1 ) );
  // v = 0.5 * ( (ls - li) * tanh( 0.5 * alpha * e ) + (ls + li) )
  chain_v.row( 0 ) =  0.5 * ( (ls - li) * tanh( 0.5 * alpha * chain_v.row( 0 ) ) + (ls + li) );
  
  List chain = List::create( Named("chain_theta") = chain_theta,
                             Named("chain_b") = chain_b,
                             Named("chain_h") = chain_h,
                             Named("chain_v") = chain_v,
                             Named("chain_l") = chain_l
                            );
  
  vec acc = zeros<vec>(4, 1);
  acc[ 0 ] += acc_theta;
  acc[ 1 ] += acc_b;
  acc[ 2 ] += acc_h;
  acc[ 3 ] += acc_v;
  
  vec div = zeros<vec>(3, 1);
  div[ 0 ] += div_theta;
  div[ 1 ] += div_b;
  div[ 2 ] += div_v;
  
  double time = timer.toc();
  
  Rcout << "\nTempo decorrido: " << time / 60 << " min" << endl;
  Rcout << "\nTaxa aceitação:\n" << acc / N << endl;
  Rcout << "\nDivergências:\n" << div << endl;
  
  return List::create( Named("chain") = chain, 
                       Named("acc") = acc, 
                       Named("time") = time 
  ); 
  
}
