trace_plots = function(draws, burn = 0, lags = 1, names){
  N = ncol( draws )
  Draws = draws[, -c( 1:burn )]
  jumps = seq(1, N - burn, by = lags)
  Draws = Draws[, jumps ]
  
  n = nrow( Draws )
  par(mfrow = c(1,1))
  mat = matrix(seq(1, 3 * n), nrow = 3, byrow = FALSE)
  layout( mat )
  
  for(i in 1:n){
    plot( Draws[i, ], type = 'l', main = names[i], xlab = '', ylab = '')
    plot(acf( Draws[i, ], lag.max = 100, plot = FALSE)[1:100], main ='', xlab = '', ylab = 'ACF')
    plot(density( Draws[i, ] ), main ='', xlab = '', ylab = '')
  }
  
  par(mfrow = c(1,1))
}

abs_plots = function( draws_h, date, y ){
  library(ggplot2)
  #T = length( y )
  e.vol_hat = apply( exp( 0.5 * draws_h ) , MARGIN = 1, FUN = mean )
  e.vol_min = apply( exp( 0.5 * draws_h ) , MARGIN = 1, FUN = quantile, probs = c(0.025) )
  e.vol_max = apply( exp( 0.5 * draws_h ) , MARGIN = 1, FUN = quantile, probs = c(0.975) )
  data = matrix(c(abs( y ), e.vol_hat, e.vol_min, e.vol_max), ncol = 4)
  data = data.frame(data)
  data = cbind( date, data )
  names(data) = c('date', 'y.abs', 'e.hat', 'e.min','e.max')
  h = ggplot(data)
  h = h + geom_line(aes(date, y.abs), color = 'grey70')
  h = h + geom_ribbon(aes(x = date, ymax = e.vol_max, ymin = e.vol_min), 
                      fill = 'blue' ,alpha = 0.2)
  h = h + geom_line(aes(date, e.hat), linewidth = 0.75)
  h = h + theme_test() + xlab('') + ylab('|Retornos|')
  h = h + theme(axis.title.y = element_text(size = 18),
                axis.text.x = element_text(size = 16),
                axis.text.y = element_text(size = 18))
  h = h + xlim(as.Date(c(date[1], tail(date, 1) ) ) )
  h
}

tail_plot = function(draws, date, model_name){
  library(ggplot2)
  if( !require(latex2exp) ) install.packages("latex2exp")
  library(latex2exp)
  l_hat = apply(draws, MARGIN = 1, mean)
  df = data.frame(date = date, l = l_hat)
  h = ggplot(df) + geom_line(aes(x=date, y = l)) 
  h = h + theme_test() + xlab('')
  #h = h + scale_x_date(date_breaks = '7 year', date_labels = '%Y')
  h = h + ylab(TeX(paste0('SVM-', model_name, ': ', '$\\lambda_{t}$')))
  h = h + theme(axis.title.y = element_text(size = 18),
                axis.text.x = element_text(size = 14),
                axis.text.y = element_text(size = 14))
  #h = h + xlim(as.Date(c(date[1], tail(date, 1) ) ) )
  h
}
