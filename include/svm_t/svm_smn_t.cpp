#include "svm_smn_t.h"

// parameters priori
double mu_0 = 0.0, s_0 = 3.2;
double a_phi = 20.0, b_phi = 1.5;
double a_s = 2.5, b_s = 0.025;
double mu_b0 = 0.0, s_b0 = 3.2;
double a_b1 = 5.0, b_b1 = 1.5;
double mu_b2 = 0.0, s_b2 = 3.2;
double mu_e = 20.0;
double s_e = 10.0;
//double a_v = 2.0, b_v = 0.1;
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
int G_theta(vec theta, mat &G, mat &inv_G, mat& dG_mu, mat &dG_omega, mat &dG_gamma, int T){
    
    //theta = (mu, w, gama)
    double phi = tanh( theta[1] );
    double sigma = exp( theta[2] );
    
    G = zeros<mat>(3, 3);
    inv_G = zeros<mat>(3, 3);
    dG_mu = zeros<mat>(3, 3);
    dG_omega = zeros<mat>(3, 3);
    dG_gamma = zeros<mat>(3, 3);

    G(0, 0) = ( ( 1 - phi ) * ( 1 - phi ) * ( T - 1 ) + ( 1 - phi * phi ) ) / ( sigma * sigma ) + ( 1.0 / ( s_0 * s_0 ) );
    G(1, 1) = 2 * phi * phi + ( T - 1 ) * ( 1 - phi * phi ) + ( a_phi + b_phi ) * ( 1 - phi * phi );
    G(1, 2) = 2 * phi;
    G(2, 1) = 2 * phi;
    G(2, 2) = 2 * T + 4 * b_s / ( sigma * sigma );
    
    //if( !G.is_sympd() ) std::cout << "theta" << theta[0] << " " << theta[1] << " " << theta[2] << std::endl;
    if( !G.is_sympd() ){
      return -1;
    }

    inv_G = inv_sympd( G );

    dG_omega(0, 0) = ( -2 * ( 1 - phi * phi ) / ( sigma * sigma ) ) * ( ( T - 1 ) * ( 1 - phi ) + phi );
    dG_omega(1, 1) = 2 * phi * ( 1 - phi * phi ) * ( 2 - T + 1 - a_phi - b_phi );
    dG_omega(1, 2) =  2 * ( 1 - phi * phi );
    dG_omega(2, 1) = 2 * ( 1 - phi * phi );

    dG_gamma(0, 0) = - ( 2 / ( sigma * sigma ) ) * ( ( T - 1 ) * ( 1 - phi ) * ( 1 - phi ) + 1 - phi * phi );
    dG_gamma(2, 2) = - 8 * b_s / ( sigma * sigma );
    
    return 1;
    
}
vec nuH_theta(mat dG_mu, mat dG_omega, mat dG_gamma, mat inv_G, vec p){
    vec iGp = inv_G * p;
    vec aux1 = iGp.t() * dG_mu * iGp;
    vec aux2 = iGp.t() * dG_omega * iGp;
    vec aux3 = iGp.t() * dG_gamma * iGp;
    vec nuH = zeros<vec>(3, 1);
    nuH[0] = aux1[0];
    nuH[1] = aux2[0];
    nuH[2] = aux3[0];
    
    return nuH;
}
vec hmc_theta(vec theta_cur, vec h, int L, vec eps, int T, int &acc){

    vec vH = zeros<vec>(3 , 1);
    vec pcur = mvrgaussian( 3 );
    
    mat inv_G = ones<mat>(3, 3);
    mat dG_mu = zeros<mat>(3, 3);
    mat dG_omega = zeros<mat>(3, 3);
    mat dG_gamma = zeros<mat>(3, 3);

    double Hm1 = 0.0;
    double Hmf = 0.0;
    double uMH = 0.0;
    double aMH = 0.0;
    
    vec auxp = pcur.t() * pcur;
    Hm1 = - logpost_theta( theta_cur, h, T ) + 0.5 * auxp[0];
    vec pa = pcur;
    vec pb = pcur;
    vec vHa = theta_cur;
    vec vHb = theta_cur;
    
    for (int jlp = 1 ; jlp < L + 1 ; jlp++ ){
    
    	pb = pa + 0.5 * eps % glogpost_theta( vHa, h, inv_G, dG_mu, dG_omega, dG_gamma, T );
    	vHb = vHa + eps % pb;
    	pb = pb + 0.5 * eps % glogpost_theta( vHb, h, inv_G, dG_mu, dG_omega, dG_gamma, T );
    	vHa = vHb;
    	pa = pb;
    	
    }
    
    auxp = pb.t() * pb;
    Hmf = - logpost_theta( vHb, h, T ) + 0.5 * auxp[0];

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
vec rmhmc_theta(vec theta_cur, vec h, int fixp, int L, vec eps, int T, int &acc){
  
    vec theta1p = zeros<vec>(3, 1);
    vec theta1 = zeros<vec>(3, 1);
    mat inv_G = zeros<mat>(3, 3);
    mat dG_mu = zeros<mat>(3, 3);
    mat dG_omega = zeros<mat>(3, 3);
    mat dG_gamma = zeros<mat>(3, 3);
    mat G = zeros<mat>(3, 3);
    int gcond = 0;
    /****************************************************************/
    gcond = G_theta(theta_cur, G, inv_G, dG_mu, dG_omega, dG_gamma, T);
    if( gcond == -1 ) return theta_cur;
    vec pcur = chol( G ).t() * mvrgaussian( 3 );
    /****************************************************************/
    vec aux = pcur.t() * inv_G * pcur;
    vec pa = zeros<vec>(3, 1);
    vec pb = zeros<vec>(3, 1);
    vec pc = pcur;
    vec gdpth1( 3 );
    vec gdmom1( 3 );
    vec gdmom2( 3 );
    /****************************************************************/
    double uMH = 0.0;
    double aMH = 0.0;
    double HMtheta1 = -logpost_theta( theta_cur, h, T ) + 0.5 * log( det(G) ) + 0.5 * aux[0];
    double HMthetaf = 0.0;

    vec gradHMCth1 = zeros<mat>(3, 1);
    vec theta1c = theta_cur;
    vec theta1a;
    vec theta1b;
    theta1a = theta1c;
    pa = pc;
    /**********************************************************/
    /**********************************************************/
    for( int jgl = 1 ; jgl < L + 1 ; jgl++ ){
	
	    gdpth1 = glogpost_theta( theta1a, h, inv_G, dG_mu, dG_omega, dG_gamma, T);
    	gradHMCth1 = gdpth1 - 0.5 * nuH_theta( dG_mu, dG_omega, dG_gamma, inv_G, pa );
      //cout << gdpth1 << ' ' << gradHMCth1 << endl;
    	/****************************************/
    	/***** pn+1/2 *****************************/
	
    	for( int jfix = 1; jfix < fixp + 1 ; jfix++ ){
        	pb = pa - 0.5 * eps % gradHMCth1;
        	gradHMCth1 = gdpth1 - 0.5 * nuH_theta( dG_mu, dG_omega, dG_gamma, inv_G, pb );
          //cout << gradHMCth1 << endl;
    	}
    
    	/**********************************************************/
    	gcond = G_theta( theta1a, G, inv_G, dG_mu, dG_omega, dG_gamma, T );
    	if( gcond == -1 ) return theta_cur;
      gdmom1 = gH_mom( inv_G, pb);
    	gdmom2 = gH_mom( inv_G, pb);
    	/****************************************************/
    	/****thetan+1****************************************/

    	for( int jfix = 1; jfix < fixp + 1; jfix ++ ){
        
        	theta1b = theta1a + 0.5 * eps % gdmom1 + 0.5 * eps % gdmom2;
        	gcond = G_theta( theta1b, G, inv_G, dG_mu, dG_omega, dG_gamma, T );
          if( gcond == -1 ) return theta_cur;
        	gdmom2 = gH_mom( inv_G, pb );
        
    	}
    /*************************************************************/
    	gdpth1 = glogpost_theta( theta1b, h, inv_G, dG_mu, dG_omega, dG_gamma, T);
    	gradHMCth1 = gdpth1 - 0.5 * nuH_theta( dG_mu, dG_omega, dG_gamma, inv_G, pb );
      //cout << gdpth1 << ' ' << gradHMCth1 << endl;

    	pb = pb - 0.5 * eps % gradHMCth1;
    	pa = pb;
    	theta1a = theta1b;
    }

    aux = pb.t() * inv_G * pb;
    HMthetaf = -logpost_theta( theta1b, h, T ) + 0.5 * log( det( G ) ) + 0.5 * aux[ 0 ];
    uMH = randu( distr_param( 0, 1 ) );
    aMH = std::min( 1.0, exp( -HMthetaf + HMtheta1 ) );
    
    if ( uMH < aMH ){ 
        theta1 = theta1b;
        acc ++; 
        }
    else{ 
        theta1 = theta_cur;
        }

    /**********************************************************/
    /**********************************************************/
    return theta1;
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
int G_b(vec b, vec h, vec l, mat &G, mat &inv_G, mat& dG_b0, mat &dG_delta, mat &dG_b2, int T, vec y_T){
   
  double b1 = tanh( b[1] );
   
  G = zeros<mat>(3, 3);
  inv_G = zeros<mat>(3, 3);
  dG_b0 = zeros<mat>(3, 3);
  dG_delta = zeros<mat>(3, 3);
  dG_b2 = zeros<mat>(3, 3);
  
  G(0, 0) += sum( l % exp( - h ) ) + 1 / ( s_b0 * s_b0 );
  G(2, 2) += sum( l %  exp( h ) ) + 1 / ( s_b2 * s_b2 );
  //G(1, 2) = ( 1 - b1 * b1 ) * sum( l % y_T.subvec(0, T - 1 ) );
  //G(2, 1) = G(1, 2);
    
  for( int j = 0 ; j < T ; j++ ){
    G(0, 1) += ( 1 - b1 * b1 ) * l[ j ] * exp( -h[j] ) * y_T[j];
    G(1, 1) += ( 1 - b1 * b1 ) * ( 1 - b1 * b1 ) * l[ j ] * exp( - h[j]  ) * y_T[j] * y_T[j];
    G(1, 2) += ( 1 - b1 * b1 ) * l[j] * y_T[j];     
  }

  dG_delta(1, 1) = - 4 * b1 * G(1, 1) - 2 * b1 * ( a_b1 + b_b1 ) * ( 1 - b1 * b1 );

  G(1, 1) += ( 1 - b1 * b1 ) * ( a_b1 + b_b1 );
  G(0, 2) += sum( l );

  G(1, 0) += G(0, 1);
  G(2, 0) += G(0, 2);
  G(2, 1) = G(1, 2);

  //if( !G.is_sympd() ) std::cout << "b" << b[0] << " " << b[1] << " " << b[2] << std::endl;
  if( !G.is_sympd() ) return -1;

  inv_G = inv_sympd( G );
  
  dG_delta(0, 1) += - 2 * b1 * G(0, 1);
  dG_delta(1, 2) += - 2 * b1 * G(1, 2);

  dG_delta(1, 0) += dG_delta(0, 1);
  dG_delta(2, 1) += dG_delta(1, 2);

  return 1;
}
vec nuH_b(mat dG_b0, mat dG_delta, mat dG_b2, mat inv_G, vec p){
    vec iGp = inv_G * p;
    vec aux1 = iGp.t() * dG_b0 * iGp;
    vec aux2 = iGp.t() * dG_delta * iGp;
    vec aux3 = iGp.t() * dG_b2 * iGp;
    vec nuH = zeros<vec>(3, 1);
    nuH[0] = aux1[0];
    nuH[1] = aux2[0];
    nuH[2] = aux3[0];
    
    return nuH;
}
vec rmhmc_b(vec b_cur, vec h, vec l, int fixp, int L, vec eps, int T, vec y_T , int &acc){
    
    vec theta2p = zeros<vec>(3, 1);
    vec b = zeros<vec>(3, 1);
    //int jfix = 1;
    //int jgl = 1;
    mat inv_G = zeros<mat>(3, 3);
    mat dG_b0 = zeros<mat>(3, 3);
    mat dG_delta = zeros<mat>(3, 3);
    mat dG_b2 = zeros<mat>(3, 3);
    mat G = zeros<mat>(3, 3);
    int gcond = 0;

    gcond = G_b( b_cur, h, l, G, inv_G, dG_b0, dG_delta, dG_b2, T, y_T );
    if( gcond == -1 ) return b_cur;

    vec pcur = chol( G ).t() * mvrgaussian( 3 );
   
    vec aux = pcur.t() * inv_G * pcur;
    vec pa = zeros<vec>(3, 1);
    vec pb = zeros<vec>(3, 1);
    vec pc = pcur;
    vec gdpth2(3);
    vec gdmom1(3);
    vec gdmom2(3);
    
    double uMH = 0.0;
    double aMH = 0.0;

    double HMtheta1 = - logpost_b( b_cur, h, l, T, y_T) + 0.5 * log( det( G ) ) + 0.5 * aux[0];
    double HMthetaf = 0.0;

    vec gradHMCth2 = zeros<mat>(3, 1);
    vec theta2c = b_cur;
    vec theta2a;
    vec theta2b;
    theta2a = theta2c;
    pa = pc;
    
    for (int jgl = 1 ; jgl < L + 1 ; jgl++ ){

    	gdpth2 = glogpost_b( theta2a, h, l, inv_G, dG_b0, dG_delta, dG_b2, T, y_T );
    	gradHMCth2 = gdpth2 - 0.5 * nuH_b(dG_b0, dG_delta, dG_b2, inv_G, pa);
    
    	for(int jfix = 1 ; jfix < fixp + 1 ; jfix++ ){
    		pb = pa - 0.5 * eps % gradHMCth2;
     		gradHMCth2 = gdpth2 - 0.5 * nuH_b( dG_b0, dG_delta, dG_b2, inv_G, pb);
    	}
    
    	gcond = G_b( theta2a, h, l, G, inv_G, dG_b0, dG_delta, dG_b2, T, y_T);
    	if( gcond == -1 ) return b_cur;
      gdmom1 = gH_mom( inv_G, pb );
    	gdmom2 = gH_mom( inv_G, pb );
    
    	for(int jfix = 1 ; jfix < fixp + 1 ; jfix++ ){
    		theta2b = theta2a + 0.5 * eps % gdmom1 + 0.5 * eps % gdmom2;
    		gcond = G_b( theta2b, h, l, G, inv_G, dG_b0, dG_delta, dG_b2, T, y_T);
        if( gcond ==-1 ) return b_cur;
    		gdmom2 = gH_mom( inv_G, pb );
    	}
    
    	gdpth2 = glogpost_b( theta2b, h, l, inv_G, dG_b0, dG_delta, dG_b2, T, y_T );
    	gradHMCth2 = gdpth2 - 0.5 * nuH_b( dG_b0, dG_delta, dG_b2, inv_G, pb);
    	pb = pb - 0.5 * eps % gradHMCth2;
    	pa = pb;
    	theta2a = theta2b;
    }

    aux = pb.t() * inv_G * pb;
    HMthetaf = - logpost_b( theta2b, h, l, T, y_T ) + 0.5 * log( det( G ) ) + 0.5 * aux[0];
    uMH = randu( distr_param( 0, 1 ) );
    aMH = std::min( 1.0, exp( - HMthetaf + HMtheta1 ) );

    if( uMH < aMH ){
    	b = theta2b;
    	acc++;
    }
    else{
    	b = b_cur;
    }

    return b;
}
//#############################################################################
//#############################################################################
//########################## v.3
//######## Transformação: T(v) = e                         ############
//######## e = (2 / alpha) * arctanh((2 * v - ls - li) / (ls - li)), alpha != 0 e li < ls
// Const. definition
// e ~ N(mu_e, sigma_e)
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
double glogpost_v( double e, vec l, double inv_G, double dG_v, int T, double alpha, double li, double ls ){
  
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
  //# Partial Energy function
  grad += - 0.5 * inv_G * dG_v;
  
  return - grad;
}
int G_v( double e, double &G, double &inv_G, double &dG_v, int T, double alpha, double li, double ls ){
  
  double v = 0.5 * ( (ls - li) * tanh( 0.5 * alpha * e ) + (ls + li) );
  double dv = 0.25 * alpha * (ls - li) / (cosh(0.5 * alpha * e) * cosh(0.5 * alpha * e));
  double dv2 = - 0.25 * alpha * alpha * (ls - li) * sinh(0.5 * alpha * e);
  dv2 /= (cosh(0.5 * alpha*e) * cosh(0.5 * alpha*e) * cosh(0.5 * alpha*e));

  G = 0.0;
  inv_G = 0.0;
  dG_v = 0.0;
  
  // log verossimilhança
  G += - dv * dv * T * ( 0.5 / v - 0.25 * R::psigamma(0.5 * v, 1) );
  // lof priori
  G += 1 / (s_e * s_e);  
  
  // inv_G
  inv_G += 1 / G;
  
  // dG_v
  dG_v += - 2 * dv * dv2 * T * ( 0.5 / v - 0.25 * R::psigamma(0.5 * v, 1) );
  dG_v += ( dv * dv * dv ) * T * ( 0.5 / (v*v) + 0.125 *  R::psigamma(0.5 * v, 2) );
    
  return 1;
}
double nuH_v(double dG_v, double inv_G, double p){
  
    double iGp = inv_G * p;
    double nuH = iGp * dG_v * iGp;
    //cout << "dG_v: " << dG_v << endl;    
    return nuH;
}
double gradHmom_v(double inv_G, double p){
    return inv_G * p;
}
double rmhmc_v(double e_cur, vec l, int fixp, int L, double eps, int T, int &acc, double alpha, double li, double ls ){
    
    //double theta2p = 0.0;
    double e = 0.0;
    double inv_G = 0.0;
    double dG_v = 0.0;
    double G = 0.0;
    
    G_v( e_cur, G, inv_G, dG_v, T, alpha, li, ls );
    double pcur =  G * randn( );
   
    double aux = pcur * inv_G * pcur;
    double pa = 0.0;
    double pb = 0.0;
    double pc = pcur;
    double gdpth2;
    double gdmom1;
    double gdmom2;
    
    double uMH = 0.0;
    double aMH = 0.0;

    double HMtheta1 = - logpost_v( e_cur, l, T, alpha, li, ls ) + 0.5 * log(  G  ) + 0.5 * aux;
    double HMthetaf = 0.0;

    double gradHMCth2 = 0.0;
    //double theta2c = e_cur;
    double theta2a;
    double theta2b = 0.0;
    theta2a = e_cur;
    pa = pc;
    
    for (int jgl = 1 ; jgl < L + 1 ; jgl++ ){

    	gdpth2 = glogpost_v( theta2a, l, inv_G, dG_v, T, alpha, li, ls );
      //cout << gdpth2 << endl;
    	gradHMCth2 = gdpth2 - 0.5 * nuH_v(dG_v, inv_G, pa);
      //cout << nuH_v(dG_v, inv_G, pa) << endl;

    	for(int jfix = 1 ; jfix < fixp + 1 ; jfix++ ){
    		pb = pa - 0.5 * eps * gradHMCth2;
     		gradHMCth2 = gdpth2 - 0.5 * nuH_v( dG_v, inv_G, pb);
        //cout << nuH_v(dG_v, inv_G, pb) << endl;
    	}
    
    	G_v( theta2a, G, inv_G, dG_v, T, alpha, li, ls );
    	gdmom1 = gradHmom_v( inv_G, pb );
    	gdmom2 = gradHmom_v( inv_G, pb );
    
    	for(int jfix = 1 ; jfix < fixp + 1 ; jfix++ ){
    		theta2b = theta2a + 0.5 * eps * gdmom1 + 0.5 * eps * gdmom2;
    		G_v( theta2b, G, inv_G, dG_v, T, alpha, li, ls );
    		gdmom2 = gradHmom_v( inv_G, pb );
    	}
    
    	gdpth2 = glogpost_v( theta2b, l, inv_G, dG_v, T, alpha, li, ls );
      //cout << gdpth2 << endl;
    	gradHMCth2 = gdpth2 - 0.5 * nuH_v( dG_v, inv_G, pb);
      //cout << nuH_v( dG_v, inv_G, pb) << endl;
    	pb = pb - 0.5 * eps * gradHMCth2;
    	pa = pb;
    	theta2a = theta2b;
    }

    aux = pb * inv_G * pb;

    HMthetaf = - logpost_v( theta2b, l, T, alpha, li, ls ) + 0.5 * log(  G  ) + 0.5 * aux;
    uMH = randu( distr_param( 0, 1 ) );

    aMH = std::min( 1.0, exp( - HMthetaf + HMtheta1 ) );

    if( uMH < aMH ){
    	e += theta2b;
    	acc++;
    }
    else{
    	e += e_cur;
    }

    return e;
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
