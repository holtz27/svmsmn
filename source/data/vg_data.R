vg_data = function(mu, phi, sigma,
                   b0, b1, b2,
                   y0,
                   v, 
                   T,
                   seed = NULL){
  library( invgamma )
  ############################################################################
  ################################# SVM-SMN Model ############################
  # y_t = b0 + b1 * y_{t-1} + b2 * e(h_t) + e_t/lambda_t^{1/2}
  # h_t = mu + phi * ( h_{t-1} - mu ) + sigma * n_{t}
  ############################################################################
  ###################### Variance Gamma: lambda|v ~ IGamma(v/2, v/2) 
  ### Obs.: X ~ Gamma(a, rate = b) -> 1/X ~ IGamma(a, 1/rate = scale = b)
  # https://en.wikipedia.org/wiki/Inverse-gamma_distribution#Related_distributions
  if( is.null( seed ) ) seed = sample( 1:1e6, 1 )
  y = h = l = rep(0, T)
  set.seed( seed )
  for(t in 1:T){
    if(t == 1){
      #l[t] = 1 / rgamma(1, shape = v/2, rate = v/2)
      l[t] = rinvgamma(1, shape = 0.5 * v, rate = 0.5 * v)
      h[t] = rnorm(1, mean = mu, sd = sigma * 1 / sqrt( (1 - phi * phi) ) )
      y[t] = rnorm(1, b0 + b1 * y0 + b2 * exp( h[t] ), exp(h[t]/2) / sqrt( l[t] ))
    }else{
      #l[t] = 1 / rgamma(1, shape = v/2, rate = v/2)
      l[t] = rinvgamma(1, shape = 0.5 * v, rate = 0.5 * v)
      h[t] = rnorm(1, mean = (mu + phi * ( h[t-1] - mu )), sd = sigma)
      y[t] = rnorm(1, b0 + b1 * y[t-1] + b2 * exp(h[t]), exp(h[t]/2) / sqrt( l[t] ))
    }
  }
  return( list(y = y, h = h, l = l) )
}
