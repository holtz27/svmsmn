// [[Rcpp::depends( RcppArmadillo )]]

#include "svm_smn_t.h"

//#############################################################################
//########################## rtgamma
/*
double rtgamma( double min, double max, double shape, double scale ){
  
  Environment pkg = Environment::namespace_env("truncdist");
  Function f = pkg["rtrunc"];
  
  NumericVector x;
  x = f( 1, 
         Named("spec") = "gamma", 
         _["a"] = min, 
         _["b"] = max, 
         _["shape"] = shape, 
         _["scale"] = scale);
  
  //Rcout << shape << " " << scale << endl;
  
  return x[ 0 ];
}
 */
double rtgamma( double li, double ls, double shape, double scala ){
  
  double tg;
  
  tg = R::pgamma( ls, shape, scala, true, false );
  tg -= R::pgamma( li, shape, scala, true, false );
  tg *= R::runif( 0, 1 );
  tg += R::pgamma( li, shape, scala, true, false );
  tg = R::qgamma( tg, shape, scala, true, false );
  
  return tg;
}
//########################## l
vec l_gibbs(double v, vec y_T, vec h, vec b, int T){
  
  double b0 = b[0];
  double b1 = tanh( b[1] ); 
  double b2 = b[2];
  
  vec l_out = zeros<vec>(T, 1);
  
  vec aux = y_T.subvec( 1, T ) - b0 - b1 * y_T.subvec( 0, T - 1 ) - b2 * exp( h );
  vec u = 0.5 * exp( - h ) % aux % aux;
  // scale = 1 / rate
  
  for( int i = 0 ; i < T ; i++ ){

   l_out[ i ] = rtgamma( 0.0, 1.0, v + 0.5, 1.0 / u[ i ] );
    
  }
  
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
  theta_cur[ 0 ] += 0.005;
  theta_cur[ 1 ] += 0.5 * ( log( 1 + 0.98 ) - log( 1 - 0.98 ) );
  theta_cur[ 2 ] += log( sqrt( 0.017 ) );
  
  // iniciando h
  int acc_b = 0;
  vec h_cur = zeros<vec>(T, 1);
  h_cur[ 0 ] += 0.005 + sqrt( 0.03 ) / (1 - 0.95 * 0.95 ) * randn();
  for( int kt = 1 ; kt < T ; kt++ ){
    h_cur[ kt ] += 0.005 + 0.95 * ( h_cur[ kt - 1 ] -0.005 ) + sqrt( 0.03 ) * randn();
  }
  
  // iniciando b
  int acc_h = 0;
  vec b_cur = zeros<vec>(3, 1);
  b_cur[ 0 ] += 0.3;
  b_cur[ 1 ] += 0.5 * ( log( 1 + 0.03 ) - log( 1 - 0.03 ) );
  b_cur[ 2 ] += -0.025;
  
  // iniciando v
  double v_cur = 5.0;

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
    v_cur = rtgamma( 1.0, R_PosInf, T + 0.08, 1.0 / (0.04 - sum( log(l_cur) )) );
    //v_cur = rtgamma( 1.0, 40.0, T + 0.08, 1.0 / (0.04 - sum( log(l_cur) )) );
    l_cur = l_gibbs(v_cur, y_T, h_cur, b_cur, T);
   
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
  
  return List::create( Named("chain") = chain, 
                       Named("acc") = acc, 
                       Named("time") = time 
  ); 
  
}
