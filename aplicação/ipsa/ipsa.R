#Banco coletado em: 
#https://br.investing.com/indices/
# getwd()
ipsa = read.csv('https://raw.githubusercontent.com/holtz27/svmsmn/main/aplica%C3%A7%C3%A3o/ipsa/IPSA.csv')
ipsa = ipsa[, c('Data', 'Último')]
ipsa[, 2] = stringr::str_remove(ipsa[, 2], "[.]")
ipsa[, 2] = stringr::str_replace(ipsa[, 2], "[,]", ".")
ipsa[, 2] = as.numeric( ipsa[, 2] ) 
#View(ipsa)
T = nrow( ipsa )
log.ret = 100 * ( log( ipsa[2:T, 2] ) - log( ipsa[1:(T-1), 2] ) )
ipsa[, 1] = stringr::str_replace_all(ipsa[, 1], "[.]", "-")
dates = as.Date( ipsa[, 1], "%d-%m-%Y" )

# Plots
library(ggplot2)
df = data.frame( Retorno = log.ret , Tempo = dates[-1]  ) #Tempo = dates[-1]

g = ggplot(df) + geom_line(aes(x = Tempo, y = Retorno))
g = g + scale_x_date(date_breaks = "18 month", date_labels = "%b %Y")
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
round( data_summary, digits = 2 )

Box.test(log.ret, lag = 12, type = 'Ljung-Box')
Box.test(log.ret**2, lag = 12, type = 'Ljung-Box')
acf(log.ret, lag.max = 5, plot = FALSE)


