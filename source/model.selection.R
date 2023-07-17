############### Some auxiliar functions
svmsmn_loglik = function(param, data, y0){
  #param_hat = ( b = b_hat, h = h_hat, l = l_hat )
  param = as.vector( param )
  T  = length( data )
  b0 = param[1]
  b1 = param[2]
  b2 = param[3]
  h  = param[4:(T+3)]
  l  = param[(T + 4):(2 * T + 3)]
  
  #mu
  mu = matrix(c(b0 + b1 * c(y0, data[1:(T-1)]) + b2 * exp( h )), 
              ncol = 1)
  # Sigma
  Sigma = diag( l * exp( -h ), nrow = T )
  
  #log( p (y | b, h, l) )
  k = -0.5 * T * log( 2 * pi ) - 0.5 * sum(  h - log(l)  )
  loglik = as.numeric( -0.5 * t(data - mu) %*% Sigma %*% (data - mu) )  
  
  return( k + loglik )
  
}
#Rcpp::sourceCpp('C:/Users/8936381/Documents/source/svmsmn_loglik.cpp')
svmsmn_loglik_i = function(data, y0, param_draws){
  
  x = apply(X = param_draws, MARGIN = 2, FUN = svmsmn_loglik, data, y0)
  
  return( x )
}
lik = function(data_i, draws, data_, data_past, T){
  # data_ = ( y_{1}, y_{2}, ..., y_{T} )
  # data_past = ( y_{0}, y_{1}, ..., y_{T-1} )
  # draws (2 * T + 3) x M
  data_past = matrix(data_past, ncol = 1)
  k = which( data_ == as.numeric( data_i ) )
  N = ncol( draws )
  log_l = rep(0, N)
  
  for(col in 1:N){
    
    b0 = as.numeric( draws[1, col] )
    b1 = as.numeric( draws[2, col] )
    b2 = as.numeric( draws[3, col] )
    h  = as.numeric( draws[3 + k, col] )
    l  = as.numeric( draws[T + 3 + k, col] )
    
    log_l[ col ] = dnorm(as.numeric( data_i ), 
                         mean = b0 + b1 * as.numeric( data_past[k, 1]) + b2 * exp( h ), 
                         sd = exp( 0.5 * h ) / sqrt( l ), log = FALSE )
  }
  
  return( log_l )
}
loglik = function(data_i, draws, data_, data_past, T){
  # data_ = ( y_{1}, y_{2}, ..., y_{T} )
  # data_past = ( y_{0}, y_{1}, ..., y_{T-1} )
  # draws (2 * T + 3) x M
  data_past = matrix(data_past, ncol = 1)
  k = which( data_ == as.numeric( data_i ) )
  N = ncol( draws )
  log_l = rep(0, N)
  
  for(col in 1:N){
    
    b0 = as.numeric( draws[1, col] )
    b1 = as.numeric( draws[2, col] )
    b2 = as.numeric( draws[3, col] )
    h  = as.numeric( draws[3 + k, col] )
    l  = as.numeric( draws[T + 3 + k, col] )
    
    log_l[ col ] = dnorm(as.numeric( data_i ), 
                         mean = b0 + b1 * as.numeric( data_past[k, 1]) + b2 * exp( h ), 
                         sd = exp( 0.5 * h ) / sqrt( l ), log = TRUE )
  }
  
  return( log_l )
}
############### DIC
svmsmn_dic = function(data, y0, param_draws){
  
  param_hat = apply(X = param_draws, MARGIN = 1, mean)
  
  pD = 2 * ( svmsmn_loglik( param_hat, data, y0 ) - mean( svmsmn_loglik_i(data, y0, param_draws) ) )      
  DIC = - 2 * svmsmn_loglik( param_hat, data, y0 ) + 2 * pD  
  
  return( DIC )
}
############### waic
svmsmn_waic = function(data, y0, draws){
  T = length( data )
  X = sapply(X = data, 
             FUN = loglik,
             draws,
             data,
             c( y0, data[1:(T-1)] ),
             length( data )
  ) 
  
  waic = loo::waic( X )
  
  return( waic )
}
############### loo
svmsmn_loo = function(data, y0, draws, cores){
  
  T = length( data )
  data_past = c( y0, data[1:(T-1)] )
  
  r_eff = loo::relative_eff(lik,
                            chain_id = rep(1, ncol( draws ) ),
                            data = matrix( data , ncol = 1), 
                            draws = draws,
                            data_ = data,
                            data_past = data_past,
                            T = T)
  
  loo = loo::loo(loglik, 
                 r_eff = r_eff, 
                 data = matrix( data , ncol = 1), 
                 draws = draws,
                 data_ = data,
                 data_past = data_past,
                 T = T,
                 cores = getOption('mc.cores', cores)
  )
  return( loo )
}
