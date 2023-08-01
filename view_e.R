# v = 20, ls = 40, l1 = 2, T = 2e3
li = 2
ls = 40
T = 2e3
v = 20
alpha = .2
(2 / alpha) * atanh((2 * v - ls - li) / (ls - li))
l = rgamma(T, shape = 0.5 * v, rate = 0.5 * v)

L = function( e ){
  v = 0.5 * ( (ls - li) * tanh( 0.5 * alpha * e ) + (ls + li) )
  a_v = 2
  b_v = .1
  # log p( l|v )
  x = T * ( 0.5 *  v * log(0.5 * v) - log( gamma(0.5 * v) ) ) 
  x = x + 0.5 * v * sum(log(l) - l) - sum( log(l) )
  # log p( v )
  x = x + ( a_v - 1 ) * log( v ) - b_v * v
  # log dv/de
  x = x + log(ls - v) + log(v - li)
  return( x )
}
G = function( e ){
  v = 0.5 * ( (ls - li) * tanh( 0.5 * alpha * e ) + (ls + li) )
  a_v = 2
  b_v = .1
  
  # dv/de
  k0 = alpha * (ls - v) * (v - li) / (ls - li)
  w = (2 * v - ls - li) / (ls - li)
  # dv²/de²
  k1 = - alpha * w * k0
  
  x = 0.5 * T / v  - 0.25 * T * psigamma(0.5 * v, 1)
  x = x - (a_v - 1) / (v * v)
  x = x - 1 / ((ls - v) * (ls - v)) - 1 / ((v - li) * (v - li)) 
  
  #x =  k0 * k0 * x
  y = k1 * ((a_v - 1) / v - b_v - 1 / (ls - v) + 1 / (v - li))
  
  return( - k0 * k0 * x - k1 * y )
}
H = function( e, p ){
  x = - L( e ) + 0.5 * log( 2 * pi * G( e ) ) + 0.5 * p * p / G( e )
  return( x )
}

e = seq(from = -10, to = 10, by = .1)
p = seq(from = -10, to = 10, by = .1)
z = outer(e, p, FUN = H)

contour(e, p, z, xlab = 'e', ylab = 'p', nlevels = 25 )
abline( v = (2 / alpha) * atanh((2 * v - ls - li) / (ls - li)) )

persp(e, p, z, theta = 30, phi = 0, expand = 0.5, col = "lightblue", 
      xlab = 'e', zlab = 'H(e,p)')

par( mfrow = c(1, 2) )
e0 = seq(from = -10, to = 10, by = .1)
h0 = apply(matrix(e0, nrow = 1), MARGIN = 1, FUN = H, p = 0)
plot(e0, h0, type = 'l')
abline(v = (2 / alpha) * atanh((2 * v - ls - li) / (ls - li)))
p0 = seq(from = -10, to = 10, by = .1)
h1 = apply(matrix(p0, nrow = 1), MARGIN = 1, FUN = H, e = (2 / alpha) * atanh((2 * v - ls - li) / (ls - li)))
plot(p0, h1, type = 'l')
par( mfrow = c(1, 1) )
