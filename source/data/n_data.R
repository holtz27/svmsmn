normal_data = function(mu, phi, sigma,
                       b0, b1, b2,
                       y0, 
                       T,
                       seed = NULL){
  ############################################################################
  ################################# SVM-SMN Model ############################
  # y_t = b0 + b1*y_{t-1} + b2*e(h_t) + e_t/lambda_t^{1/2}
  # h_t = mu + phi(h_{t-1} - mu) + sigma*n_{t}
  ############################################################################
  ###################### Normal: l ~ 1 
  if( is.null(seed) ) seed = sample(1:1e6, 1)
  y = h = rep(0, T)
  set.seed( seed )
  for(t in 1:T){
    if(t == 1){
      h[t] = rnorm(1, mean = mu, sd = sigma * sqrt( 1 / ( 1 - phi**2 ) ))
      y[t] = rnorm(1, b0 + b1 * y0 + b2 * exp(h[t]), exp(h[t]/2) )
    }else{
      h[t] = rnorm(1, mean = (mu + phi * ( h[t-1] - mu) ), sd = sigma)
      y[t] = rnorm(1, b0 + b1*y[t - 1] + b2 * exp(h[t]), exp(h[t]/2) )
    }
  }
  return( y = y, h = h )
}
