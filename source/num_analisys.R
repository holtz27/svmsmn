num_analisys = function( draws, burn = 0, lags = 1, names, digits ){
  if( !require(coda) ) install.packages("coda")
  ############### Numeric Analysis
  Draws = draws
  N = ncol( Draws ) 
  if( burn != 0 ) Draws = Draws[, -c( 1:burn )]
  jumps = seq(1, N - burn, by = lags)
  Draws = Draws[, jumps ]
  
  chain = coda::as.mcmc( t( Draws ) )
  CD = coda::geweke.diag( chain )
  N_eff = coda::effectiveSize( chain )
  IF = ncol( Draws ) / N_eff
  mc_error = apply( chain, 
                    MARGIN = 2, 
                    FUN = sd) / sqrt( N_eff )
  
  theta_hat = apply( chain, MARGIN = 2, FUN = mean )
  theta_sd  = apply( chain, MARGIN = 2, FUN = sd )
  theta_min = apply( chain, MARGIN = 2, FUN = quantile, probs = c(0.025) )
  theta_max = apply( chain, MARGIN = 2, FUN = quantile, probs = c(0.975) )
  data = matrix(
    c(theta_hat,
      theta_sd,
      theta_min,
      theta_max,
      CD$z,
      IF,
      mc_error), nrow = nrow( Draws ), byrow = FALSE
  )
  row.names( data ) = names
  colnames( data ) = c( 'media', 'sd', '2.5%', '97.5%', 'CD', 'IF', 'MC erro')
  
  return( round( data, digits ) )
}
