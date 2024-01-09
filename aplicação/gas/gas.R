

# Gás Natural
#Banco coletado em: 
#https://br.investing.com/
gas = read.csv('gas.csv', dec = ',')
gas = gas[, c('Data', 'Último')]
#gas[, 2] = stringr::str_remove(gas[, 2], "[.]")
#gas[, 2] = stringr::str_replace(gas[, 2], "[,]", ".")
gas[, 2] = as.numeric( gas[, 2] )
gas[, 1] = stringr::str_replace_all(gas[, 1], "[.]", "-")
dates = as.Date( gas[, 1], "%d-%m-%Y" )
gas$Data = dates
gas = gas[ order(gas$Data), ]
View(gas)
T = nrow( gas )
log.ret = 100 * ( log( gas[2:T, 2] ) - log( gas[1:(T-1), 2] ) )

# Plots
library(ggplot2)
df = data.frame( Retorno = log.ret , Tempo = dates[-1]  ) #Tempo = dates[-1]

g = ggplot(df) + geom_line(aes(x = Tempo, y = Retorno))
g = g + scale_x_date(date_breaks = "60 month", date_labels = "%b %Y")
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


