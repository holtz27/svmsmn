
cambio = read.csv('https://raw.githubusercontent.com/holtz27/svmsmn/main/aplica%C3%A7%C3%A3o/cambio/gbp_usd.csv')
cambio = cambio[, c('Date', 'Adj.Close')]
cambio[, 2] = as.numeric( cambio[, 2] ) 
cambio = na.omit( cambio )
dates = as.Date( cambio[, 1], "%Y-%m-%d" )
View(cambio)
T = nrow(cambio)
log.ret = 100 * ( log( cambio[2:T, 2] ) - log( cambio[1:(T-1), 2] ) )

# Plots
library(ggplot2)
df = data.frame( Retorno = log.ret, Tempo = dates[-1] )

g = ggplot(df) + geom_line(aes(x = Tempo, y = Retorno))
g = g + scale_x_date(date_breaks = "32 month", date_labels = "%b %Y")
g = g + theme_test() + theme(axis.title.y = element_text(size = 18),
                             axis.text.x = element_text(size = 16),
                             axis.text.y = element_text(size = 18))
g = g + xlab('')
h = ggplot( df, aes(Retorno) )
h = h + geom_histogram(aes(y = after_stat(density)), bins = 40, color = 'white')
h = h + theme_test() + ylab('')
h = h + theme_test() + theme(axis.title.x = element_text(size = 18),
                             axis.text.x = element_text(size = 18),
                             axis.text.y = element_text(size = 18))
gridExtra::grid.arrange(g, h, nrow = 1, ncol = 2) 

data_summary = matrix(c( mean( log.ret ),
                         sd( log.ret ),
                         min( log.ret ),
                         max( log.ret ),
                         moments::skewness( log.ret ),
                         moments::kurtosis( log.ret ) ), nrow = 1)
colnames( data_summary ) = c( 'mean', 'sd', 'min', 'max', 'skewness', 'kurtosis')
round( data_summary, digits = 4 )

Box.test(log.ret, lag = 12, type = 'Ljung-Box')
Box.test(log.ret**2, lag = 12, type = 'Ljung-Box')
acf(log.ret, lag.max = 5, plot = FALSE)
