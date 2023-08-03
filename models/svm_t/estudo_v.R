# parameters priori
a_v = 2.0
b_v = 0.1

#############################################################################
#############################################################################
########################## v.1
######## Transformação: T(v) = e                         ############
######## e = (2 / alpha) * arctanh((2 * v - ls - li) / (ls - li)), alpha != 0 e li < ls
#Const. definition
alpha = .1
li = 2
ls = 40
v_init = 20
e_init = ( 2 / alpha ) * atanh( (2 * v_init - ls - li) / (ls - li) )

logpost_v = function( e, l, T){
  
  L = 0.0
  v = 0.5 * ( (ls - li) * tanh( 0.5 * alpha * e ) + (ls + li) )
  
  L = L + 0.5 * T * v * log( 0.5 * v  ) - T * log( gamma( 0.5 * v )  )
  L = L + 0.5 * v * sum( log( l ) - l )
  # priori
  L = L + ( a_v - 1 ) * log( v ) - b_v * v
  # log Jacobiano
  L = L + log(ls- v) + log(v - li)
  
  return( L )
}
glogpost_v = function( e, l, inv_G, dG_v, T){
  
  v = 0.5 * ( (ls - li) * tanh( 0.5 * alpha * e ) + (ls + li) )
  grad  = 0.0
  
  grad = grad + 0.5 * T * log(0.5 * v) + 0.5 * T - 0.5 * T * psigamma(0.5 * v, 0)
  grad = grad + 0.5 * sum(log( l ) - l)
  # priori 
  grad = grad + (a_v - 1) * log( v ) - b_v * v
  # log jacobiano
  grad = grad + log(ls - v) + log(v - li)
  # chain rule factor
  grad = grad * alpha * (ls - v) * (v - li) / (ls - li)
  
  # Partial Energy function
  grad = grad - 0.5 * inv_G * dG_v
  
  return( -grad )
}
G_v = function( e, T ){
  
  v = 0.5 * ( (ls - li) * tanh( 0.5 * alpha * e ) + (ls + li) )
  
  G = 0.0
  inv_G = 0.0
  dG_v = 0.0
  
  w = (2 * v - ls - li) / (ls - li)
  k0 = alpha * (ls - v) * (v - li) / (ls - li)
  k1 = - alpha * w * k0
  
  G = - k0 * k0 * (0.5 * T / v - 0.25 * T * psigamma(0.5 * v, 1))
  G = G + k0 * k0 * ((a_v - 1) / ( v * v ) + 1 / ((v - ls) * (v - ls)) + 1 / ((v - li) * (v - li)))
  G = G - k1 * ((a_v - 1) / v - b_v + 1 / (v - ls) + 1 / (v - li))
  
  #G += 0.25 * T * v * v * R::psigamma(0.5 * v, 1);
  #G +=  v * b_v - 0.5 * T * v;
  
  inv_G = 1 / G
  
  k2 = - alpha * w;
  k3 = - alpha * alpha * (6 * v * (ls + li - v) - (ls - li) * (ls - li) - 2 * ls * li) / (ls - li) * (ls - li)
  
  dG_v = - 2 * k2 * (0.5 * T * v - 0.25 * T * psigamma(0.5 * v, 1))
  dG_v = dG_v + 2 * k2 * ((a_v - 1) / ( v * v ) + 1 / ((v - ls) * (v - ls)) + 1 / ((v - li) * (v - li)))
  dG_v = dG_v - k0 * k0 * (- 0.5 * T / ( v*v ) - 0.125 * psigamma(0.5 * v, 2))
  dG_v = dG_v - 2 * k0 * k0 * ( (a_v - 1) / ( v * v * v ) + 1 / ( (v - ls) * (v - ls) * (v - ls) ) + 1 / ((v - li) * (v - li) * (v - li)) )
  dG_v = dG_v - k3 * ( (a_v - 1) / v - b_v + 1 / (v - ls) + 1 / (v - li) )
  dG_v = dG_v + k1 * k1 * ( (a_v - 1)/ (v * v) - 1 / ((v - ls) * (v - ls)) - 1 / ( (v - li) * (v - li)) )
  
  #dG_v += 0.5 * T * v * v * R::psigamma(0.5 * v, 1);
  #dG_v += 0.125 * T * v * v * v * R::psigamma(0.5 * v, 2);
  #dG_v += b_v * v - 0.5 * T * v;
  
  #return 1;
  return( c(G, inv_G, dG_v) )
}
nuH_v = function( dG_v, inv_G, p){
  
  iGp = inv_G * p
  nuH = iGp * dG_v * iGp
  
  return( nuH )
}
gradHmom_v = function( inv_G, p){
  return( inv_G * p )
}
rmhmc_v = function( e_cur, l, fixp, L, eps, T, acc){
  
  #double theta2p = 0.0;
  e = 0.0
  inv_G = 0.0
  dG_v = 0.0
  G = 0.0
  
  pointer = G_v( e_cur, T )
  G = pointer[1]
  inv_G = pointer[2]
  dG_v = pointer[3]
  pcur =  G * rnorm( 1 )
  
  aux = pcur * inv_G * pcur
  pa = 0.0
  pb = 0.0
  pc = pcur
  #gdpth2;
  #gdmom1;
  #gdmom2;
  
  uMH = 0.0
  aMH = 0.0
  
  HMtheta1 = - logpost_v( e_cur, l, T ) + 0.5 * log(  G  ) + 0.5 * aux
  HMthetaf = 0.0
  
  gradHMCth2 = 0.0
  theta2c = e_cur
  #theta2a
  #theta2b
  theta2a = theta2c
  pa = pc
  
  # int jgl = 1 ; jgl < L + 1 ; jgl++ 
  for(jgl in 1:L){
    
    gdpth2 = glogpost_v( theta2a, l, inv_G, dG_v, T )
    gradHMCth2 = gdpth2 - 0.5 * nuH_v(dG_v, inv_G, pa)
    
    # int jfix = 1 ; jfix < fixp + 1 ; jfix++ 
    for(jfix in 1:fixp){
      pb = pa - 0.5 * eps * gradHMCth2
      gradHMCth2 = gdpth2 - 0.5 * nuH_v( dG_v, inv_G, pb)
    }
    
    #G_v( theta2a, G, inv_G, dG_v, T );
    pointer = G_v( theta2a, T )
    G = pointer[1]
    inv_G = pointer[2]
    dG_v = pointer[3]
    gdmom1 = gradHmom_v( inv_G, pb )
    gdmom2 = gradHmom_v( inv_G, pb )
    
    # int jfix = 1 ; jfix < fixp + 1 ; jfix++
    for(jfix in 1:jfix){
      theta2b = theta2a + 0.5 * eps * gdmom1 + 0.5 * eps * gdmom2
      #G_v( theta2b, G, inv_G, dG_v, T );
      pointer = G_v( theta2b, T )
      G = pointer[1]
      inv_G = pointer[2]
      dG_v = pointer[3]
      gdmom2 = gradHmom_v( inv_G, pb )
    }
    
    gdpth2 = glogpost_v( theta2b, l, inv_G, dG_v, T )
    gradHMCth2 = gdpth2 - 0.5 * nuH_v( dG_v, inv_G, pb)
    pb = pb - 0.5 * eps * gradHMCth2
    pa = pb
    theta2a = theta2b
  }
  
  aux = pb * inv_G * pb
  
  HMthetaf = - logpost_v( theta2b, l, T ) + 0.5 * log(  G  ) + 0.5 * aux
  uMH = runif( 1 )
  
  aMH = min( 1.0, exp( - HMthetaf + HMtheta1 ) )
  
  if( uMH < aMH ){
    e = e + theta2b
    acc = acc + 1
  }
  else{
    e = e + e_cur
  }
  
  return( e )
}
