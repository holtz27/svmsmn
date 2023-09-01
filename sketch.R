#############################################################################
#############################################################################
########################## v.1
######## Transformação: T(v) = e                         ############
######## e = (2 / alpha) * arctanh((2 * v - ls - li) / (ls - li)), 
# alpha != 0 e li < ls
jac_v = function( v ) alpha * (ls - v) * (v - li) / (ls - li)
jac2_v = function( v ){
  w = (2 * v - ls - li) / (ls - li);
  return( - alpha * w * jac_v( v ) )
}

gjac_v = function( v ){
  w = (2 * v - ls - li) / (ls - li);
  return( - alpha * w )
}
gjac2_v = function( v ){
  
  x = 6 * v * (li + ls - v) - (ls + li)**2 - 2 * li * ls
  x = - alpha * alpha * x / (ls - li)**2 
  
  return( x )
}

#logjac = log( jac_v )
glogjac_v = function( v ) - 1 / (ls - v) + 1 / (v - li);
glogjac2_v = function( v ){
  x = - 1 / ((v - ls) * (v - ls)) - 1 / ((v - li) * (v - li))
  return( x )
}
glogjac3_v = function( v ){
  
  x = - 2 / (ls - v)**3 + 2 / (v - li)**3
  
  return( x )
}

###############
logpost_v = function( e, l, T ){
  
  L = 0.0;
  v = 0.5 * ( (ls - li) * tanh( 0.5 * alpha * e ) + (ls + li) );
  
  # log p( l|v )
  L = L + 0.5 * T * v * log( 0.5 * v  ) - T * log( gamma( 0.5 * v )  );
  L = L + 0.5 * v * sum( log( l ) - l );
  # log priori
  L = L + ( a_v - 1 ) * log( v ) - b_v * v;
  # log Jacobiano: log(ls - v) + log(v - li)
  L = L + log( jac_v( v ) )
  
  return ( L );
}
glogpost_v = function( e, l, inv_G, dG_v, T ){
  
  v = 0.5 * ( (ls - li) * tanh( 0.5 * alpha * e ) + (ls + li) );
  grad = 0.0;
  
  # log p( l|v )
  grad = grad + 0.5 * T * log(0.5 * v) + 0.5 * T - 0.5 * T * psigamma(0.5 * v, 0);
  grad = grad + 0.5 * sum(log( l ) - l);
  # log priori 
  grad = grad + (a_v - 1) / v - b_v;
  # log jacobiano
  grad = grad + glogjac_v( v )
  # chain rule factor
  grad = grad * jac_v( v )
  
  # Partial Energy function
  grad = grad - 0.5 * inv_G * dG_v ;
  
  return( -grad );
}
G_v = function( e, T ){
  
  v = 0.5 * ( (ls - li) * tanh( 0.5 * alpha * e ) + (ls + li) );
  
  # G
  #w = (2 * v - ls - li) / (ls - li);
  # dv / de
  #k0 = alpha * (ls - v) * (v - li) / (ls - li);
  # dv^2 / de^2
  #k1 = - alpha * w * k0;
  
  
  x = 0.5 * T / v - 0.25 * T * psigamma(0.5 * v, 1) - (a_v - 1)/ (v * v)
  x = x + glogjac2_v( v )
  y = (a_v - 1) / v - b_v + glogjac_v( v )
  #G = - k0 * k0 * x - k1 * y
  G = - jac_v( v ) * jac_v( v ) * x - jac2_v( v ) * y   
  
  # inv_G
  inv_G = 1 / G;
  
  # dG_v
  dG_v = - 2 * jac_v( v ) * gjac_v( v ) * x
  dG_v = dG_v - ( jac_v( v )**2 ) * ( - 0.5 * T / ( v*v ) - 0.125 * T * psigamma(0.5 * v, 2) + 2 * (a_v - 1) / ( v**3 ) + glogjac3_v( v ) )
  dG_v = dG_v - gjac2_v( v ) * y
  dG_v = dG_v - jac2_v( v ) * ( - (a_v - 1)/ (v * v) + glogjac2_v( v )   )
  dG_v = jac_v( v ) * dG_v
  
  #k2 = - alpha * w;
  #k3 = - alpha * alpha * (6 * v * (ls + li - v) - (ls - li) * (ls - li) - 2 * ls * li) / (ls - li) * (ls - li); 
  #dG_v = 2 * k0 * alpha * w * x
  #dG_v = dG_v - k0 * k0 * ( - 0.5 * T / ( v*v ) - 0.125 * T * psigamma(0.5 * v, 2) )
  #dG_v = dG_v - 2 * k0 * k0 * ( (a_v - 1) / ( v**3 ) + 1 / ((v - ls)**3) + 1 / ((v - li)**3) )
  #dG_v = dG_v - k3 * y
  #dG_v = dG_v - k1 * ( - (a_v - 1)/ (v * v) - 1 / ((v - ls) * (v - ls)) - 1 / ( (v - li) * (v - li)) )
  
  return( list(G = G, inv_G = inv_G, dG_v = dG_v) );
  #return( G )
}

nuH_v = function( dG_v, inv_G, p ){
  
  iGp = inv_G * p;
  nuH = iGp * dG_v * iGp;
  
  return( nuH );
}
gradHmom_v = function( inv_G, p ){
  return( inv_G * p );
}
rmhmc_v = function( e_cur, l, fixp, L, eps, T ){
  
  #double theta2p = 0.0;
  e = 0.0;
  acc = 0.0
  inv_G = 0.0;
  dG_v = 0.0;
  G = 0.0;
  
  # G_v( e_cur, G, inv_G, dG_v, T )
  G = G_v( e_cur, T )$G;
  inv_G = G_v( e_cur, T )$inv_G
  dG_v = G_v( e_cur, T )$dG_v
  
  pcur = rnorm( 1, sd = G );
  
  aux = pcur * inv_G * pcur;
  pa = 0.0;
  pb = 0.0;
  pc = pcur;
  #gdpth2;
  #gdmom1;
  #gdmom2;
  
  uMH = 0.0;
  aMH = 0.0;
  
  HMtheta1 = - logpost_v( e_cur, l, T ) + 0.5 * log(  G  ) + 0.5 * aux;
  HMthetaf = 0.0;
  
  gradHMCth2 = 0.0;
  theta2c = e_cur;
  #theta2a;
  #theta2b;
  theta2a = theta2c;
  pa = pc;
  
  # int jgl = 1 ; jgl < L + 1 ; jgl++
  for( jgl in 1:L ){
    
    gdpth2 = glogpost_v( theta2a, l, inv_G, dG_v, T );
    gradHMCth2 = gdpth2 - 0.5 * nuH_v(dG_v, inv_G, pa);
    
    # int jfix = 1 ; jfix < fixp + 1 ; jfix++
    for( jfix in 1:fixp ){
      pb = pa - 0.5 * eps * gradHMCth2;
      gradHMCth2 = gdpth2 - 0.5 * nuH_v( dG_v, inv_G, pb);
    }
    
    # G_v( theta2a, G, inv_G, dG_v, T );
    G = G_v( theta2a, T )$G;
    inv_G = G_v( theta2a, T )$inv_G
    dG_v = G_v( theta2a, T )$dG_v
    
    gdmom1 = gradHmom_v( inv_G, pb );
    gdmom2 = gradHmom_v( inv_G, pb );
    
    # int jfix = 1 ; jfix < fixp + 1 ; jfix++
    for( jfix in 1:fixp ){
      theta2b = theta2a + 0.5 * eps * gdmom1 + 0.5 * eps * gdmom2;
      # G_v( theta2b, G, inv_G, dG_v, T );
      G = G_v( theta2b, T )$G;
      inv_G = G_v( theta2b, T )$inv_G
      dG_v = G_v( theta2b, T )$dG_v
      gdmom2 = gradHmom_v( inv_G, pb );
    }
    
    gdpth2 = glogpost_v( theta2b, l, inv_G, dG_v, T );
    gradHMCth2 = gdpth2 - 0.5 * nuH_v( dG_v, inv_G, pb);
    pb = pb - 0.5 * eps * gradHMCth2;
    pa = pb;
    theta2a = theta2b;
  }
  
  aux = pb * inv_G * pb;
  
  HMthetaf = - logpost_v( theta2b, l, T ) + 0.5 * log(  G  ) + 0.5 * aux;
  uMH = runif( 1 )
  aMH = min( 1.0, exp( - HMthetaf + HMtheta1 ) );
  
  if( uMH < aMH ){
    e = e + theta2b;
    acc = 1;
  }else{
    e = e + e_cur;
  }
  
  return( list( e = e, acc = acc ) );
}
svm_t = function(N, L_v, eps_v, l, seed ){
  
  #wall_clock timer;
  #timer.tic();
  
  if( seed != 0 ) set.seed( seed );
  
  T = length( l )
  #a = floor( 0.1 * N );
  
  # iniciando v
  # double alpha = 2.0;
  acc_v = 0;
  v_cur = 10.0;
  v_cur = (2 / alpha) * atanh( (2 * v_cur - ls - li) / (ls - li) );
  
  chain_v = matrix(nrow = 1, ncol = (N + 1) )
  chain_v[1, 1] = v_cur;
  
  # chain builting  
  for( it in 2:(N + 1) ){
    
    v_cur = rmhmc_v(v_cur, l, 5, L_v, eps_v, T )$e;
    acc_v = acc_v + rmhmc_v(v_cur, l, 5, L_v, eps_v, T )$acc
    
    # chain update 
    chain_v[ 1, it ] = v_cur;
    
    #Progress
    #if( (it % a) == 0 ) cout << "Progresso em " << ceil( 100 * it / N ) <<" %"<< endl;
  }
  
  return( list(chain = chain_v, acc = acc_v) )
}

















####################################
# Const. definition
#### librarys
source('https://raw.githubusercontent.com/holtz27/svmsmn/main/source/data/t_data.R')
a_v = 2
b_v = .1
alpha = 1
li = 2
ls = 60
T = 2e3
v = 20
data = t_data( mu = -10.0, phi = 0.95, sigma = 0.15,
               b0 = 0.1, b1 = 0.01, b2 = -0.05,
               y0 = 0,
               v = v,
               T = 2e3,
               seed = 0 )
l = data$l
e_true = (2 / alpha) * atanh( (2 * v - ls - li) / (ls - li) )
# G
e = seq(from = -1, to = 1, by = .1)
g = apply(matrix(e, nrow = 1), MARGIN = 1, G_v, T = 2e3)
sum( g < 0 )

# log e glog
y0 = apply(matrix(e, nrow = 1), MARGIN = 1, FUN = logpost_v, l = l, T = T)
y1 = apply(matrix(e, nrow = 1), MARGIN = 1, FUN = glogpost_v, 
           l = l, inv_G = 1, dG_v = 1, T = T)

par( mfrow = c(1, 3) )
plot(e, g, type = 'l')
abline( h = 0 )
abline( v = e_true )
plot(e, y0, type = 'l')
abline(v = e_true)
plot(e, y1, type = 'l')
abline(v = e_true)
abline(h = 0 )
par( mfrow = c(1, 1) )
####################################
# H
H = function( e, p, l, T ){
  x = - logpost_v( e, l, T ) + 0.5 * log( 2 * pi * G_v( e, T ) ) + 0.5 * p * p / G_v( e, T )
  return( x )
}

p = seq(from = -10, to = 10, by = .5)
z = outer(e, p, FUN = H, l, T)
contour(e, p, z, xlab = 'e', ylab = 'p', nlevels = 25 )
abline(v = e_true )

persp(e, p, z, theta = 30, phi = 0, expand = 0.5, col = "lightblue", 
      xlab = 'e', zlab = 'H(e,p)')

par( mfrow = c(1, 2) )
h0 = apply(matrix(e, nrow = 1), MARGIN = 1, FUN = H, p = 0, l = l, T = T)
plot(e, h0, type = 'l')
abline(v = e_true )
p0 = seq(from = -10, to = 10, by = .1)
h1 = apply(matrix(p0, nrow = 1), MARGIN = 1, FUN = H, e = e_true, l = l, T = T )
plot(p0, h1, type = 'l')
par( mfrow = c(1, 1) )
##############################################################################################################################

#### librarys
source('https://raw.githubusercontent.com/holtz27/svmsmn/main/source/data/t_data.R')

a_v = 2
b_v = .1
alpha = .1
li = 2
ls = 40
T = 2e3
v = 20
e_true = (2 / alpha) * atanh( (2 * v - ls - li) / (ls - li) ) 
e_true


data = t_data( mu = -1.0, phi = 0.985, sigma = 0.15,
               b0 = 0.1, b1 = 0.01, b2 = -0.05,
               y0 = 0,
               v = v,
               T = T,
               seed = 423 )
l = data$l

N = 1e3
draws = svm_t( N = N, 
               L_v = 5, eps_v = 0.25, 
               l = l, 
               seed = 0 )

draws$acc / N
chain = draws$chain
chain = 0.5 * ( (ls - li) * tanh( 0.5 * alpha * chain ) + (ls + li) );

par( mfrow = c(1, 2) )
plot( chain[1, ], type = 'l' )
abline( h = v )
plot( acf(chain[1, ], lag.max = 100, plot = FALSE)[1:100] )
par( mfrow = c(1, 1) )


