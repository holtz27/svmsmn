t_data = function(mu, phi, sigma,
                   b0, b1, b2,
                   y0,
                   v, 
                   T,
                   seed = NULL){
  ############################################################################
  ################################# SVM-SMN Model ############################
  # y_t = b0 + b1 * y_{t-1} + b2 * e(h_t) + e_t/lambda_t^{1/2}
  # h_t = mu + phi * ( h_{t-1} - mu ) + sigma * n_{t}
  ############################################################################
  ###################### t-Student: lambda|v ~ Gama(v/2, v/2) 
  if( is.null( seed ) ) seed = sample( 1:1e6, 1 )
  y = h = l = rep(0, T)
  set.seed( seed )
  for(t in 1:T){
    if(t == 1){
      l[t] = rgamma(1, shape = v/2, rate = v/2)
      h[t] = rnorm(1, mean = mu, sd = sigma * 1 / sqrt( (1 - phi * phi) ) )
      y[t] = rnorm(1, b0 + b1 * y0 + b2 * exp( h[t] ), exp(h[t]/2) / sqrt( l[t] ))
    }else{
      l[t] = rgamma(1, shape = v/2, rate = v/2)
      h[t] = rnorm(1, mean = (mu + phi * ( h[t-1] - mu )), sd = sigma)
      y[t] = rnorm(1, b0 + b1 * y[t-1] + b2 * exp(h[t]), exp(h[t]/2) / sqrt( l[t] ))
    }
  }
  return( list(y = y, h = h, l = l) )
}
