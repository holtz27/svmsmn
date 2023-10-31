// [[Rcpp::depends( RcppArmadillo )]]
// [[Rcpp::depends( RcppGSL )]]

#include "svm_smn.h"
//#include <gsl/gsl_cdf.h>

//#############################################################################
//########################## rtgamma

double right_tgamma(double max, double shape, double rate){

  // G <- get(paste("p", spec, sep = ""), mode = "function")
  // G is the cdf function 
  // G = R::pgamma( p, shape, scale, lower, log )
  // Gin <- get(paste("q", spec, sep = ""), mode = "function")
  // Gin is the inverse cdf function 
  // Gin = R::qgamma( q, shape, scale, lower, log )
  // tt = Gin(G(a, ...) + p*(G(b, ...) - G(a, ...)), ...)
  
  double scale = 1 / rate;

  if( max <= 0 )
    stop( "argumnto b é menor ou igual a 0!" );
  
  double rtg, q, log_Gb;
  
  if( max != R_PosInf )
    log_Gb = R::pgamma( max, shape, scale, true, true ); 
    if( std::isinf( log_Gb ) )
      stop( "Trunction interval is not inside the domain of the quantile function" );
 
  double u = R::runif( 0, 1 );
  double x = R::qgamma( log( u ) + log_Gb, shape, scale, true, true );

  return x;
  
}

double left_tgamma( double min, double shape, double rate ){
  
  double scale = 1 / rate;
  
  if( min <= 0 )
    stop( "argumnto min tem que ser positivo!" );
  
  double ltg, q;
  
  double u = R::runif( 0, 1 );
  q = R::pgamma( min, shape, scale, true, false );
  q += u * R::pgamma( min, shape, scale, false, false );
  
  ltg = R::qgamma( q, shape, scale, true, false );
  
  return ltg;
}

double qtrunc( double p, double a, double b, double shape, double scale, int msn_erro ){
  
  // G <- get(paste("p", spec, sep = ""), mode = "function")
  // G is the cdf function 
  // G = R::pgamma( p, shape, scale, lower, log )
  // Gin <- get(paste("q", spec, sep = ""), mode = "function")
  // Gin is the inverse cdf function 
  // Gin = R::qgamma( q, shape, scale, lower, log )
  // tt = Gin(G(a, ...) + p*(G(b, ...) - G(a, ...)), ...)
  
  if( a >= b )
    stop( "argument a is greater than or equal to b" );
  
  double tt, q, Ga = 0.0, Gb = 1.0, log_Ga, log_Gb;
  
  if( a != 0 ) 
    Ga = R::pgamma( a, shape, scale, true, false );
    log_Ga = R::pgamma( a, shape, scale, true, true );
    //Ga = gsl_cdf_gamma_P( a, shape, scale );

  if( b != R_PosInf ) 
    Gb = R::pgamma( b, shape, scale, true, false );
    log_Gb = R::pgamma( b, shape, scale, true, true );
    //Gb = gsl_cdf_gamma_P( b, shape, scale );

  if( log_Ga == log_Gb ){ 
    cout << Gb << " " << log_Gb << endl;
    stop( "Trunction interval is not inside the domain of the quantile function" );  
  }
  
  if( msn_erro == 1 ) 
    tt = R::qgamma( log( p ) + log_Gb, shape, scale, true, true );
  else
    q = Ga + p * ( Gb - Ga );
    tt = R::qgamma( q, shape, scale, true, false );

  return tt;
} 
double rtgamma( double min, double max, double shape, double rate, int msn_erro ){
  
  double scale = 1 / rate;
  double u = R::runif( 0, 1 );
  double x = qtrunc( u, min, max, shape, scale, msn_erro );
  
  return x;
  
}
/*
//#############################################################################
//########################## rtgamma
double rtgamma( double li, double ls, double shape, double rate, int msn_erro ){
  
  double tg;
  double scale = 1.0 / rate;
  
  tg = R::pgamma( ls, shape, scale, true, false );
  tg -= R::pgamma( li, shape, scale, true, false );
  tg *= R::runif( 0, 1 );
  tg += R::pgamma( li, shape, scale, true, false );
  tg = R::qgamma( tg, shape, scale, true, false );
  
  if( std::isnan( tg ) ){
    cout << "shape: " << shape << " " << "rate: " << rate << endl;
    if( msn_erro == 0 ) stop( "Erro em v_cur!" );
      else stop( "Erro em l_cur!" );
  }else return tg;
}
*/
//########################## l
vec l_gibbs( double v, vec y_T, vec h, vec b, int T ){
  
  double b0 = b[0];
  double b1 = tanh( b[1] ); 
  double b2 = b[2];
  
  vec l_out = zeros<vec>(T, 1);
  
  vec aux = y_T.subvec( 1, T ) - b0 - b1 * y_T.subvec( 0, T - 1 ) - b2 * exp( h );
  vec u = 0.5 * exp( - h ) % aux % aux;

  // scale = 1 / rate
  for( int i = 0 ; i < T ; i++ ) 
    //l_out[ i ] = rtgamma( 0.0, 1.0, v + 0.5, u[ i ], 1 );
    l_out[ i ] = right_tgamma( 1.0, v + 0.5, u[ i ] );

  return l_out;
  
}
// ############################## set seed function
void set_seed( int seed ){
  Rcpp::Environment base_env("package:base");
  Rcpp::Function set_seed_r = base_env["set.seed"];
  set_seed_r( std::floor( std::fabs( seed ) ) );
}
// [[Rcpp::export]]
List svm_s(int N, 
           int L_theta, vec eps_theta, 
           int L_b, vec eps_b, 
           int L_h, double eps_h,
           vec y_T, 
           int seed ){
  
  wall_clock timer;
  timer.tic();
  
  if( seed != 0 ) set_seed( seed );
  
  int T = y_T.n_elem - 1, a = floor( 0.1 * N );
  
  // iniciando theta
  int acc_theta = 0;
  vec theta_cur = zeros<vec>(3, 1);
  theta_cur[ 0 ] += 0.1;
  theta_cur[ 1 ] += 0.5 * ( log( 1 + 0.98 ) - log( 1 - 0.98 ) );
  theta_cur[ 2 ] += log( 0.15 );
  
  // iniciando h
  int acc_b = 0;
  vec h_cur = zeros<vec>(T, 1);
  h_cur[ 0 ] += 0.1 + 0.15 / sqrt(1 - 0.98 * 0.98 ) * randn();
  for( int kt = 1 ; kt < T ; kt++ ){
    h_cur[ kt ] += 0.1 + 0.98 * ( h_cur[ kt - 1 ] - 0.1 ) + 0.15 * randn();
  }
  
  // iniciando b
  int acc_h = 0;
  vec b_cur = zeros<vec>(3, 1);
  b_cur[ 0 ] += 0.1;
  b_cur[ 1 ] += 0.5 * ( log( 1 + 0.03 ) - log( 1 - 0.03 ) );
  b_cur[ 2 ] += -0.10;
  
  // iniciando v
  double v_cur = 2.0;
  
  // iniciando l
  vec l_cur = zeros<vec>(T, 1);
  for( int k = 0 ; k < T ; k++ ){
    l_cur[ k ] = R::rbeta( v_cur, 1.0 );
  }
    
  // iniciando cadeia
  mat chain_theta = zeros<mat>( 3, N + 1 );
  chain_theta.col( 0 ) += theta_cur;
  
  mat chain_b = zeros<mat>( 3, N + 1 );
  chain_b.col( 0 ) += b_cur;
  
  mat chain_h = zeros<mat>( T, N + 1 );
  chain_h.col( 0 ) += h_cur;
  
  mat chain_v = zeros<mat>( 1, N + 1 );
  chain_v.col( 0 ) += v_cur;
  
  mat chain_l = zeros<mat>( T, N + 1 );
  chain_l.col( 0 ) += l_cur;
  
  // chain builting  
  for(int it = 1 ; it < N + 1 ; it ++){
    
    theta_cur = rmhmc_theta( theta_cur, h_cur, 5, L_theta, eps_theta, T, acc_theta );
    b_cur = rmhmc_b( b_cur, h_cur, l_cur, 5, L_b, eps_b, T, y_T , acc_b );
    h_cur = hmc_h( h_cur, theta_cur, b_cur, l_cur, L_h, eps_h, T, y_T, acc_h );
    //v_cur = rtgamma( 1.0, R_PosInf, T + 0.08, (0.04 - sum( log(l_cur) )), 0 );
    v_cur = left_tgamma( 1.0, T + 0.08, 0.04 - sum( log(l_cur) ) );
    l_cur = l_gibbs( v_cur, y_T, h_cur, b_cur, T );
    
    // chain update 
    chain_theta.col( it ) += theta_cur;
    chain_b.col( it ) += b_cur;
    chain_h.col( it ) += h_cur;
    chain_v.col( it ) += v_cur;
    chain_l.col( it ) += l_cur;
    
    //Progress
    if( (it % a) == 0 ) cout << "Progresso em " << ceil( 100 * it / N ) <<" %"<< endl;
  }
  
  // Transformations
  chain_theta.row( 1 ) = tanh( chain_theta.row( 1 ) );
  chain_theta.row( 2 ) = exp( chain_theta.row( 2 ) );
  chain_b.row( 1 )     = tanh( chain_b.row( 1 ) );
    
  List chain = List::create( Named("chain_theta") = chain_theta,
                             Named("chain_b") = chain_b,
                             Named("chain_h") = chain_h,
                             Named("chain_v") = chain_v,
                             Named("chain_l") = chain_l);
  
  vec acc = zeros<vec>(3, 1);
  acc[ 0 ] += acc_theta;
  acc[ 1 ] += acc_b;
  acc[ 2 ] += acc_h;
  
  double time = timer.toc();
  
  Rcout << "\nTempo decorrido: " << time / 60 << " min" << endl;
  Rcout << "\nTaxas aceitação:\n" << acc / N << endl;
  
  return List::create( Named("chain") = chain, 
                       Named("acc") = acc, 
                       Named("time") = time 
  ); 
  
}
