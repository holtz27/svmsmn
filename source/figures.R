trace_plots = function(draws, burn = 0, lags = 1, names, lag.max = 100){
  Draws = draws
  N = ncol( Draws ) 
  if( burn != 0 ) Draws = Draws[, -c( 1:burn )]
  jumps = seq(1, N - burn, by = lags)
  Draws = Draws[, jumps ]
  
  n = nrow( Draws ) 
  par(mfrow = c(1,1))
  mat = matrix(seq(1, 3 * n), nrow = 3, byrow = FALSE)
  layout( mat )
  
  for(i in 1:n){
    plot( Draws[i, ], type = 'l', 
          main = names[ i ], 
          xlab = '', ylab = '', 
          cex.main = 2.5, 
          cex.axis = 1.5)
    plot(acf( Draws[i, ], 
              lag.max = lag.max, plot = FALSE)[1:lag.max], 
         main ='', xlab = '', ylab = 'ACF', 
         cex.lab = 1.3, 
         cex.axis = 1.5)
    hist( Draws[i, ], breaks = 20, freq = FALSE, col = 'white',
          main = '', xlab = '', ylab = '',
          cex.axis = 1.5 )
    lines(density(Draws[i, ]))
  }
  par(mfrow = c(1,1))
}


abs_plots = function( draws_h, burn = 0, lags = 1, dates = NULL, y ){
  if( !require(ggplot2) ) install.packages("ggplot2")
  library(ggplot2)
  
  Draws_h = draws_h
  N = ncol( Draws_h )
  if( (burn != 0) ) Draws_h = Draws_h[, -c( 1:burn )]
  jumps = seq(1, N - burn, by = lags)
  Draws_h = Draws_h[, jumps ]
  
  if( is.null(dates) ) dates = seq.Date(from = Sys.Date(), 
                                        length.out = length(y), 
                                        by = 'day')
  
  e.vol_hat = apply( exp( 0.5 * Draws_h ) , MARGIN = 1, FUN = mean )
  e.vol_min = apply( exp( 0.5 * Draws_h ) , MARGIN = 1, FUN = quantile, probs = c(0.025) )
  e.vol_max = apply( exp( 0.5 * Draws_h ) , MARGIN = 1, FUN = quantile, probs = c(0.975) )
  data = matrix(c(abs( y ), e.vol_hat, e.vol_min, e.vol_max), ncol = 4)
  data = data.frame(data)
  data = cbind( dates, data )
  names(data) = c('date', 'y.abs', 'e.hat', 'e.min','e.max')
  h = ggplot(data)
  h = h + geom_line(aes(date, y.abs), color = 'grey70')
  h = h + geom_ribbon(aes(x = date, ymax = e.vol_max, ymin = e.vol_min), 
                      fill = 'blue' ,alpha = 0.2)
  h = h + geom_line(aes(date, e.hat), linewidth = 0.75)
  h = h + theme_test() + xlab('') + ylab('|Retornos|')
  h = h + theme(axis.title.y = element_text(size = 26),
                axis.text.x = element_text(size = 24),
                axis.text.y = element_text(size = 26))
  h = h + xlim(as.Date(c(dates[1], tail(dates, 1) ) ) )
  h
}


tail_plot = function(draws_l, burn = 0, lags = 1, dates = NULL, model_name){
  if( !require(ggplot2) ) install.packages("ggplot2")
  library(ggplot2)
  if( !require(latex2exp) ) install.packages('latex2exp', repos='http://cran.us.r-project.org', dependencies = TRUE)
  library(latex2exp)
  
  Draws_l = draws_l[, -1]
  N = ncol( Draws_l ) 
  if( burn != 0 ) Draws_l = Draws_l[, -c( 1:burn )]
  jumps = seq(1, N - burn, by = lags)
  Draws_l = Draws_l[, jumps ]
  
  if( is.null(dates) ) dates = seq.Date(from = Sys.Date(), 
                                        length.out = nrow( Draws_l ), 
                                        by = 'day')
  
  l_hat = apply(Draws_l, MARGIN = 1, mean)
  df = data.frame(dates = dates, l = l_hat)
  
  h = ggplot( df ) + geom_line(aes( x = dates, y = l )) 
  h = h + theme_test() + xlab('')
  #h = h + scale_x_date(date_breaks = '7 year', date_labels = '%Y')
  h = h + ylab(TeX(paste0('SVM-', model_name, ': ', '$\\lambda_{t}$')))
  h = h + theme(axis.title.y = element_text(size = 18),
                axis.text.x = element_text(size = 14),
                axis.text.y = element_text(size = 14))
  h = h + xlim(as.Date(c(dates[1], tail(dates, 1) ) ) )
  h
}
