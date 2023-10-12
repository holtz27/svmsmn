#include "hmc_t.h"

// parameters priori
double mu_0 = 0.0, s_0 = 3.2;
double a_phi = 20.0, b_phi = 1.5;
double a_s = 2.5, b_s = 0.025;
double mu_b0 = 0.0, s_b0 = 3.2;
double a_b1 = 5.0, b_b1 = 1.5;
double mu_b2 = 0.0, s_b2 = 3.2;
double mu_e = 5;
double s_e = 10.0;
/***********************************************************************************/
/***********************************************************************************/
/***********************************************************************************/
vec mvrgaussian(int n){
    vec xr = zeros<vec>(n, 1);
    for( int j = 0 ; j < n ; j++ ){
        xr[j] = randn( );
    }
    return xr;
}
vec gH_mom(mat inv_G, vec p){
    return inv_G * p;
}
//########################## theta = (mu, phi, sigma) #########################
//######## Transformação: T(theta) = theta'                         ###########  
//######## theta' = (mu, arctanh(phi), log(sigma)) = (mu, w, gama) ############
double logpost_theta(vec theta, vec h, int T){
  //theta = (mu, w, gama)
  double L = 0.0; 
  double mu = theta[0];
  double phi = tanh( theta[1] ); 
  double sigma = exp( theta[2] );
  
  L += 0.5 * log(1 - phi * phi ) - T * log( sigma ); 
  L += - 0.5 * (1 - phi * phi ) * (h[0] - mu) * (h[0] - mu) / ( sigma * sigma ) ; 
  
  for( int j = 1 ; j < T ; j++ ){
    L += - ( 0.5 / (sigma * sigma) ) * ( h[j] - mu - phi * (h[j-1] - mu) ) * ( h[j] - mu - phi * (h[j-1] - mu) );
  }
    
  //# priori mu
  L += - 0.5 * ( mu - mu_0 ) * ( mu - mu_0 ) / ( s_0 * s_0 );
  //# priori phi
  L += (a_phi - 1) * log(1 + phi) + (b_phi - 1) * log(1 - phi);
  //# priori sigma
  L += - 2 * (a_s + 1) * log(sigma) - b_s / ( sigma * sigma );
  //# jacobiano de T
  L += log( 1 - phi * phi ) + log( sigma );
  
  return L; 
}
vec glogpost_theta(vec theta, vec h, mat inv_G, mat dG_mu, mat dG_omega, mat dG_gamma, int T){
  //theta = (mu, w, gama)
  double mu = theta[0];
  double phi = tanh( theta[1] ); 
  double sigma = exp( theta[2] );
   
  vec grad = zeros<vec>(3, 1); 
  
  // gradiente mu
  grad[0] += (1 - phi * phi ) * (h[0] - mu) / (sigma * sigma); 
  for( int j = 1 ; j < T ; j++ ){
    grad[0] += (1 - phi) / (sigma * sigma) * ( h[j] - mu - phi * (h[j-1] - mu) );
  }
  
  //priori
  grad[0] += - (mu - mu_0) / (s_0 * s_0);
  
  //gradiente w
  grad[1] += - phi + phi * ( 1 - phi * phi ) * (h[0] - mu) * (h[0] - mu) / ( sigma * sigma );
  for( int j = 1 ; j < T ; j++ ){
    grad[1] += (1 - phi * phi) / (sigma * sigma) * ( h[j] - mu - phi * (h[j-1] - mu) ) * ( h[j-1] - mu );
  }
  //priori
  grad[1] += (a_phi - 1) * (1 - phi);
  grad[1] += - (b_phi - 1) * (1 + phi);
  //jacobiano
  grad[1] += - 2 * phi;
  
  // gradiente gama
  grad[2] += - T + (1 - phi * phi ) * (h[0] - mu) * (h[0] - mu) / ( sigma * sigma ); 
  for( int j = 1 ; j < T ; j++ ){
    grad[2] += ( h[j] - mu - phi * (h[j-1] - mu) ) * ( h[j] - mu - phi * (h[j-1] - mu) ) / ( sigma * sigma );
  }
  // priori
  grad[2] += - 2 * (a_s + 1) + 2 * b_s / ( sigma * sigma );
  // jacobiano
  grad[2] += 1;
  
  // Partial Energy function
  grad[0] += ( -0.5 * trace( inv_G * dG_mu ) );
  grad[1] += ( -0.5 * trace( inv_G * dG_omega ) );
  grad[2] += ( -0.5 * trace( inv_G * dG_gamma ) );
  
  return - grad;
}
vec hmc_theta(vec theta_cur, vec h, int L, double eps, int T, int &acc, int &div, mat inv_M, int k, double prec){

    vec vH = zeros<vec>(3 , 1);
    vec pcur = mvrgaussian( 3 );
    vec grad = zeros<vec>(3, 1);
    
    mat inv_G = ones<mat>(3, 3);
    mat dG_mu = zeros<mat>(3, 3);
    mat dG_omega = zeros<mat>(3, 3);
    mat dG_gamma = zeros<mat>(3, 3);
   
    double Hm1 = 0.0;
    double Hmf = 0.0;
    double uMH = 0.0;
    double aMH = 0.0;
    
    vec auxp = pcur.t() * inv_M * pcur;
    Hm1 = - logpost_theta( theta_cur, h, T ) + 0.5 * auxp[ 0 ]; // H( theta, p )
    vec pa = pcur;
    vec pb = pcur;
    vec pc = pcur;
    vec vHa = theta_cur;
    vec vHb = theta_cur;
    vec vHc = theta_cur;

    // k numero de grupos
    vec H_star = zeros<vec>(k, 1);
    int cont = 0;

    for ( int jlp = 1 ; jlp < L + 1 ; jlp++ ){

      grad = glogpost_theta( vHa, h, inv_G, dG_mu, dG_omega, dG_gamma, T );
      if( grad.has_nan() ){
        div++;
        break;
      } 
      
      pb = pa + 0.5 * eps * grad;
      vHb = vHa + eps * ( inv_M * pb ); // eps % inv.M * pb, onde inv_M = Id

      grad = glogpost_theta( vHb, h, inv_G, dG_mu, dG_omega, dG_gamma, T );
      if( grad.has_nan() ){
        div++;
        break;
      } 
      
      pb = pb + 0.5 * eps * grad;

      vHa = vHb;
      vHc = vHb;
    	pa = pb;
    	pc = pb;

      auxp = pb.t() * inv_M * pb;
      H_star( (jlp - 1) % k ) += - logpost_theta( vHb, h, T ) + 0.5 * auxp[ 0 ];
      
      if( jlp % k == 0 ){
        if( stddev( H_star ) > prec ) break;
      } 
      
      cont++;
    }
    
    if( cont == 0 ) return theta_cur;

    vHb = vHc;
    pb = pc;

    auxp = pb.t() * inv_M * pb;
    Hmf = - logpost_theta( vHb, h, T ) + 0.5 * auxp[ 0 ];

    aMH = std::min( 1.0, exp( - Hmf + Hm1 ) );
    uMH = randu( distr_param( 0, 1 ) );
    
    if( uMH < aMH ){
    	vH = vHb;
    	acc++;
    }else{
        vH = theta_cur;
    }
    
    return vH;
}
//#############################################################################
//########################## b = (b0, b1, b2)           #######################
//######## Transformação: T(theta) = b'                 #######################
//######## b' = (b0, arctanh(b1), b2) = (b0, delta, b2) #######################
double logpost_b(vec b, vec h, vec l, int T, vec y_T){
  //b = (b0, delta, b2)
  //y = (y0, y1, ..., yT)
  
  double L = 0.0;
  double b0 = b[0];
  double b1 = tanh( b[1] ); 
  double b2 = b[2];
     
  //vec z = y_T.subvec(1, T) - b0 - b1 * y_T.subvec(0, T-1) - b2 * exp(h);
  //vec u = l % exp( -h ) % z;
  //L -= 0.5 * dot(z, u);

  for( int j = 0 ; j < T ; j++ ){
    L += - 0.5 * l[j] * exp( -h[j] ) * ( y_T[j+1] - b0 - b1 * y_T[ j ] - b2 * exp( h[j] ) ) * ( y_T[j+1] - b0 - b1 * y_T[ j ] - b2 * exp( h[j] ) );
  }

  //# priori b0  
  L += - 0.5 * (b0 - mu_b0) * (b0 - mu_b0) / ( s_b0 * s_b0 );
  //# priori b1
  L += (a_b1 - 1) * log(1 + b1) + (b_b1 - 1) * log(1 - b1);
  //# priori b2
  L += - 0.5 * (b2 - mu_b2) * (b2 - mu_b2) / ( s_b2 * s_b2 );
  //# jacobiano
  L += log( 1 - b1 * b1 );
  
  return L;
}
vec glogpost_b(vec b, vec h, vec l, mat inv_G, mat dG_b0, mat dG_delta, mat dG_b2, int T, vec y_T){
  //b = (b0, delta, b2)
  // y = (y0, y1, ..., yT)
  
  double b0 = b[0];
  double b1 = tanh( b[1] ); 
  double b2 = b[2];
  
  vec grad = zeros<vec>(3, 1);
  
  // construindo os vetores u, v e z
  //vec z = y_T.subvec(1, T) - b0 - b1 * y_T.subvec(0, T-1) - b2 * exp( h );
  //vec v = l % exp( -h );
  //vec u = l % exp( -h ) % y_T.subvec(0, T-1);
  //grad[0] += dot(v, z) - (b0 - mu_b0) / (s_b0 * s_b0);
  //grad[1] += (1 - b1 * b1 ) * ( dot(u, z) ) + (a_b1 - 1) * (1 - b1) - (b_b1 - 1) * (1 + b1);
  //grad[2] += sum( z ) - (b2 - mu_b2) / (s_b2 * s_b2);
  //grad[2] += dot(l, z) - (b2 - mu_b2) / (s_b2 * s_b2);

  for( int j = 0 ; j < T ; j++ ){
    grad[0] += l[j] * exp( -h[j] ) * ( y_T[j+1] - b0 - b1 * y_T[ j ] - b2 * exp( h[j] ) );
    grad[1] += (1 - b1 * b1 ) * l[j] * exp( -h[j] ) * ( y_T[j+1] - b0 - b1 * y_T[ j ] - b2 * exp( h[j] ) ) * y_T[ j ];
    grad[2] += l[j] * ( y_T[j+1] - b0 - b1 * y_T[ j ] - b2 * exp( h[j] ) );
  }
  // priori b0
  grad[0] += - (b0 - mu_b0) / (s_b0 * s_b0);
  // priori b1
  grad[1] += (a_b1 - 1) * (1 - b1) - (b_b1 - 1) * (1 + b1); 
  //# jacobiano b1
  grad[1] += - 2 * b1;
  // priori b0
  grad[2] += - (b2 - mu_b2) / (s_b2 * s_b2);
  
  // Partial Energy function
  grad[0] += ( -0.5 * trace( inv_G * dG_b0 ) );
  grad[1] += ( -0.5 * trace( inv_G * dG_delta ) );
  grad[2] += ( -0.5 * trace( inv_G * dG_b2 ) );
  
  return - grad;
}
vec hmc_b(vec b_cur, vec h, vec l, int L, double eps, int T, vec y_T, int &acc, int &div, mat inv_M, int k, double prec){

    vec vH = zeros<vec>(3 , 1);
    vec pcur = mvrgaussian( 3 );
    vec grad = zeros<vec>(3, 1);
    
    mat inv_G = ones<mat>(3, 3);
    mat dG_mu = zeros<mat>(3, 3);
    mat dG_omega = zeros<mat>(3, 3);
    mat dG_gamma = zeros<mat>(3, 3);
    
    double Hm1 = 0.0;
    double Hmf = 0.0;
    double uMH = 0.0;
    double aMH = 0.0;
    
    vec auxp = pcur.t() * inv_M * pcur;
    Hm1 = - logpost_b( b_cur, h, l, T, y_T ) + 0.5 * auxp[0];
    vec pa = pcur;
    vec pb = pcur;
    vec pc = pcur;
    vec vHa = b_cur;
    vec vHb = b_cur;
    vec vHc = b_cur;
    
    // k numero de grupos
    vec H_star = zeros<vec>(k, 1);
    int cont = 0;
    for (int jlp = 1 ; jlp < L + 1 ; jlp++ ){

      grad = glogpost_b( vHa, h, l, inv_G, dG_mu, dG_omega, dG_gamma, T, y_T );
      if( grad.has_nan() ){
        div++;
        break;
      }

    	pb = pa + 0.5 * eps * grad;
      
    	vHb = vHa + eps * ( inv_M * pb ); // eps % inv.M * pb, onde inv_M = Id

      grad = glogpost_b( vHb, h, l, inv_G, dG_mu, dG_omega, dG_gamma, T, y_T );
    	if( grad.has_nan() ){
        div++;
        break;
      }

      pb = pb + 0.5 * eps * grad;
      
    	vHa = vHb;
      vHc = vHb;
    	pa = pb;
    	pc = pb;
    	
      auxp = pb.t() * inv_M * pb;
      H_star( (jlp - 1) % k ) += - logpost_b( vHb, h, l, T, y_T ) + 0.5 * auxp[ 0 ];
      
      if( jlp % k == 0 ){
        if( stddev( H_star ) > prec ) break;
      } 
      cont++;
    }
    if( cont == 0 ) return b_cur;
    vHb = vHc;
    pb = pc;

    auxp = pb.t() * inv_M * pb;
    Hmf = - logpost_b( vHb, h, l, T, y_T ) + 0.5 * auxp[0];

    aMH = std::min( 1.0, exp( - Hmf + Hm1 ) );
    uMH = randu( distr_param( 0, 1 ) );
    
    if( uMH < aMH ){
    	vH = vHb;
    	acc++;
    }else{
        vH = b_cur;
    }
    
    return vH;
}
//#############################################################################
//#############################################################################
//########################## v.1
//######## Transformação: T(v) = e                         ############
//######## e = (2 / alpha) * arctanh((2 * v - ls - li) / (ls - li)), alpha != 0 e li < ls
// Const. definition
//double alpha = 2.0;
//double li = 2.0;
//double ls = 40.0;
//#############################################################################
double logpost_v( double e, vec l, int T, double alpha, double li, double ls ){
  
  double L = 0.0;
  double v = 0.5 * ( (ls - li) * tanh( 0.5 * alpha * e ) + (ls + li) );
  double mu = (2 / alpha) * atanh( (2 * mu_e - ls - li) / (ls - li) );

  //# log p( l|v )
  L += 0.5 * T * v * log( 0.5 * v  ) - T * log( tgamma( 0.5 * v )  );
  L += 0.5 * v * sum( log( l ) - l );
  //# log priori
  L += - 0.5 * ( (e - mu) / s_e ) * ( (e - mu) / s_e ); 
  
  return L;
}
double glogpost_v( double e, vec l, int T, double alpha, double li, double ls ){
  
  double v = 0.5 * ( (ls - li) * tanh( 0.5 * alpha * e ) + (ls + li) );
  double dv = 0.25 * alpha * (ls - li) / (cosh(0.5 * alpha * e) * cosh(0.5 * alpha * e));
  double grad = 0.0;
  double mu = (2 / alpha) * atanh( (2 * mu_e - ls - li) / (ls - li) );

  //# log p( l|v )
  grad += 0.5 * T * log(0.5 * v) + 0.5 * T - 0.5 * T * R::psigamma(0.5 * v, 0);
  grad += 0.5 * sum(log( l ) - l);
  // chain rule factor
  grad *= dv;
  // log priori 
  grad += - (e - mu) / (s_e * s_e);
  //cout << "e: " << e << endl; 
  //cout << "dv: " << dv( e, alpha, li, ls ) << endl;
  
  return - grad;
}
double hmc_v(double e_cur, vec l, int L, double eps, int T, int &acc, double alpha, double li, double ls, double inv_M, int &div, int k, double prec){

    double vH = 0.0;
    double pcur = mvrgaussian( 1 )[ 0 ];
    double grad;
    
    double Hm1 = 0.0;
    double Hmf = 0.0;
    double uMH = 0.0;
    double aMH = 0.0;
    
    double auxp = pcur * inv_M * pcur;
    Hm1 = - logpost_v( e_cur, l, T, alpha, li, ls ) + 0.5 * auxp;
    double pa = pcur;
    double pb = pcur;
    double pc = pcur;
    double vHa = e_cur;
    double vHb = e_cur;
    double vHc = e_cur;
    
    // k numero de grupos
    vec H_star = zeros<vec>(k, 1);
    int cont = 0;
    for (int jlp = 1 ; jlp < L + 1 ; jlp++ ){

      grad = glogpost_v( vHa, l, T, alpha, li, ls );
      if( std::isnan( grad ) ){
        div++;
        break;
      }

    	pb = pa + 0.5 * eps * grad;
      
    	vHb = vHa + eps * ( inv_M * pb ); // eps % inv.M * pb, onde inv_M = Id

      grad = glogpost_v( vHb, l, T, alpha, li, ls );
    	if( std::isnan( grad ) ){
        div++;
        break;
      }
      
      pb = pb + 0.5 * eps * grad;
      
    	vHa = vHb;
      vHc = vHb;
    	pa = pb;
    	pc = pb;

      auxp = pb * inv_M * pb;
      H_star( (jlp - 1) % k ) += - logpost_v( vHb, l, T, alpha, li, ls ) + 0.5 * auxp;
      
      if( jlp % k == 0 ){
        if( stddev( H_star ) > prec ) break;
      }
      cont++;
    }
    if( cont == 0 ) return e_cur;
    vHb = vHc;
    pb = pc;

    auxp = pb * inv_M * pb;
    Hmf = - logpost_v( vHb, l, T, alpha, li, ls ) + 0.5 * auxp;

    aMH = std::min( 1.0, exp( - Hmf + Hm1 ) );
    uMH = randu( distr_param( 0, 1 ) );
    
    if( uMH < aMH ){
    	vH = vHb;
    	acc++;
    }else{
        vH = e_cur;
    }
    
    return vH;
}
//#############################################################################
//########################## h
double logpost_h(vec h, vec theta, vec b, vec l, int T, vec y_T){
   
  double L = 0.0; 
  double mu = theta[0];
  double phi = tanh( theta[1] );
  double sigma = exp( theta[2] );
  double b0 = b[0];
  double b1 = tanh( b[1] );
  double b2 = b[2];
  
 
  L +=  - 0.5 * sum( h );
 
  for( int j = 1 ; j < T ; j++ ){
    L += - 0.5 / ( sigma * sigma ) * ( h[j] - mu - phi * (h[j-1] - mu) ) * ( h[j] - mu - phi * (h[j-1] - mu) );
    L += - 0.5 * l[j] * exp( -h[j] ) * ( y_T[j+1] - b0 - b1 * y_T[j] - b2 * exp( h[j]) ) * ( y_T[j+1] - b0 - b1 * y_T[j] - b2 * exp(h[j]) );
  }
  
  L += - 0.5 * l[0] * exp( -h[0] ) * ( y_T[1] - b0 - b1 * y_T[0] - b2 * exp( h[0]) ) * ( y_T[1] - b0 - b1 * y_T[0] - b2 * exp(h[0]) );
  L += - 0.5 * (1 - phi * phi ) * ((h[0] - mu)/sigma) * ((h[0] - mu)/sigma);
  
  return L;
}
vec glogpost_h(vec h, vec theta, vec b, vec l, int T, vec y_T){
  //h = (h1, ..., hT)
    
  double mu = theta[0];
  double phi = tanh( theta[1] );
  double sigma = exp(theta[2]);
  double b0 = b[0];
  double b1 = tanh( b[1] );
  double b2 = b[2];
  
  //# construindo o vetor s
  vec r  = zeros<vec>(T, 1);
  vec mu_t = y_T.subvec(1, T) - b0 - b1 * y_T.subvec(0, T - 1) - b2 * exp( h );
  
  vec s = - 0.5 + 0.5 * l % exp( -h ) % mu_t % mu_t + b2 * l % mu_t;
  
  //# construindo o vetor r
  r[0] += ( h[0] - phi * h[1] - mu * (1 - phi) ) / (sigma * sigma);
  r[T-1] = 1 / (sigma * sigma) * (h[T-1] - mu - phi * (h[T-2] - mu) );
  
  vec u = h.subvec(1, T - 2);
  vec v = h.subvec(2, T - 1) + h.subvec(0, T - 3);
  
  r.subvec(1, T-2) += (1 + phi * phi ) * u - phi * v - mu * (1 - phi) * (1 - phi);
  r.subvec(1, T-2) /= sigma * sigma ;
  
  return s - r;
}
vec hmc_h(vec h_cur, vec theta, vec b, vec l, int L, double eps, int T, vec y_T, int &acc){

    vec vH = zeros<vec>(T , 1);
    vec pcur = mvrgaussian( T );
    //vec pcur = zeros<vec>(T, 1);
    //pcur(span(1, T), 0) = mvrgaussian( T );

    double Hm1 = 0.0;
    double Hmf = 0.0;
    double uMH = 0.0;
    double aMH = 0.0;
    
    //vec  auxp = pcur( span(1, T), 0).t() * pcur( span(1, T), 0);
    vec auxp = pcur.t() * pcur;
    Hm1 = - logpost_h( h_cur, theta, b, l, T, y_T) + 0.5 * auxp[0];
    //int jlp=1;
    vec pa = pcur;
    vec pb = pcur;
    vec vHa = h_cur;
    vec vHb = h_cur;
    
    for (int jlp = 1 ; jlp < L + 1 ; jlp++ ){
    
    	pb = pa + 0.5 * eps * glogpost_h( vHa, theta, b, l, T, y_T);
    	vHb = vHa + eps * pb;
    	pb = pb + 0.5 * eps * glogpost_h( vHb, theta, b, l, T, y_T);
    	vHa = vHb;
    	pa = pb;
    	
    }
    
    //auxp = pb( span(1, T), 0).t() * pb( span(1, T), 0);
    auxp = pb.t() * pb;
    Hmf = - logpost_h( vHb, theta, b, l, T, y_T) + 0.5 * auxp[0];

    aMH = std::min( 1.0, exp( - Hmf + Hm1 ) );
    uMH = randu( distr_param( 0, 1 ) );
    
    if( uMH < aMH ){
    	vH = vHb;
    	acc++;
    }else{
        vH = h_cur;
    }
    
    return vH;
}
//#############################################################################
//########################## l
vec l_gibbs(double e, vec y_T, vec h, vec b, int T, double alpha, double li, double ls){
  
  //double v = exp( e );
  double v = 0.5 * ( (ls - li) * tanh( 0.5 * alpha * e ) + (ls + li) );
  double b0 = b[0];
  double b1 = tanh( b[1] ); 
  double b2 = b[2];
    
  vec l_out = zeros<vec>(T, 1);
  
  vec aux = y_T.subvec( 1, T ) - b0 - b1 * y_T.subvec( 0, T - 1 ) - b2 * exp( h );
  vec u = exp( - h ) % aux % aux;
  //R::rgamma( shape, scale = 1 / rate)
  //scale = 1 / 0.5 * ( u[ i ] + v ) = 2.0 / ( u[ i ] + v )

  for( int i = 0 ; i < T ; i++ ){
    l_out[ i ] = R::rgamma( 0.5 * (v + 1), 2.0 / ( u[ i ] + v ) ); 
  }
  
  return l_out;
  
}
