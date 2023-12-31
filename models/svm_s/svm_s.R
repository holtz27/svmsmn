################################################################################
#### librarys
library( loo )
library( ggplot2 )
source( 'https://raw.githubusercontent.com/holtz27/svmsmn/main/source/num_analisys.R' )
source( 'https://raw.githubusercontent.com/holtz27/svmsmn/main/source/figures.R' )
source( 'https://raw.githubusercontent.com/holtz27/svmsmn/main/source/data/s_data.R' )
source( 'https://raw.githubusercontent.com/holtz27/svmsmn/main/source/model.selection.R' )

data = s_data(mu = -1, phi = 0.975, sigma = 0.15,
              b0 = 0.1, b1 = 0.01, b2 = -0.05,
              y0 = 0,
              v = 2,
              T = 2e3,
              seed = 0)
y = data$y
l = data$l
h = data$h

# getwd()
path = 'C:/Users/8936381/Documents/svm_smn/Simulacao/Estudos_Simulacao/slash/svm_s.cpp'
Rcpp::sourceCpp( path )

# Sampling
N = 5e4
samples = svm_s(N,
                L_theta = 20, eps_theta = c( 0.5, 0.5, 0.5 ), 
                L_b = 20, eps_b = c( 0.1, 0.1, 0.1 ), 
                L_h = 50, eps_h = 0.015,
                y_T = c(0, y), 
                seed = 0 )
samples$time / 60
samples$acc / N
################## Save outputs
#save(samples, file = 'C:/Users/8936381/Documents/Simulacao/Estudos_Simulacao/ts/ts_ES.RDara')
#load('C:/Users/8936381/Documents/Simulacao/Estudos_Simulacao/ts/ts_ES.RDara')
chain_theta = samples$chain$chain_theta
chain_b = samples$chain$chain_b
chain_h = samples$chain$chain_h
chain_v = samples$chain$chain_v
chain_l = samples$chain$chain_l
# draws
draws = matrix(c( chain_theta[1, ],
                  chain_theta[2, ],
                  chain_theta[3, ],
                  chain_b[1, ],
                  chain_b[2, ],
                  chain_b[3, ],
                  chain_v
), nrow = 7, byrow = TRUE)
draws = rbind( draws, chain_h )
draws = rbind( draws, chain_l )
############################### Convergence analysis
################### Trace plots
### burn
burn = 1e3
draws = draws[, -c( 1:burn )]
lags = 20
jumps = seq(1, N - burn, by = lags)
draws = draws[, jumps ]
trace_plots(draws[1:7, ], 
            names = c('mu', 'phi', 'sigma', 'b0', 'b1', 'b2', 'v') )
tail_plot(chain_l, 1:T, 'S')
################### Numeric Analysis
num_analisys(draws[1:7, ], 
             names = c( 'mu', 'phi', 'sigma', 'b0', 'b1', 'b2', 'v'),
             digits = 4 )
################################################################################
############################## Model Selection
# p( y | theta ) = p( y | b, theta_h, v, h ) = p( y | b, h )
# data = ( y0 , y )
# param_hat = ( b = b_hat, h = h_hat, l = l_hat )
################################################################################
############### dic
ts_dic = svmsmn_dic(data = y, y0 = 0, 
                    param_draws = draws[ c( 4:6, 
                                            8:(T+7),
                                            (T+8):(nrow(draws)) ), ]
)
############### waic
ts_waic = svmsmn_waic(data = y, y0 = 0,
                      draws = draws[ c( 4:6, 
                                        8:(T+7),
                                        (T+8):(nrow(draws)) ), ]
)
############### loo
ts_loo = svmsmn_loo(data = y, y0 = 0, 
                    draws = draws[ c( 4:6, 
                                      8:(T+7),
                                      (T+8):(nrow(draws)) ), ],
                    cores = 4
)
ts_criterium = matrix( c(ts_dic, 
                         ts_loo$estimates[3,1],
                         ts_waic$estimates[3,1]), nrow = 1)
row.names( ts_criterium ) = c('ts')
colnames( ts_criterium ) = c('dic', 'loo', 'waic')
ts_criterium
###############################################################################
###############################################################################
############################### l
l_hat = apply(chain_l, MARGIN = 1, FUN = mean)
l_min = apply(chain_l, MARGIN = 1, FUN = quantile, probs = c(0.025) )
l_max = apply(chain_l, MARGIN = 1, FUN = quantile, probs = c(0.975) )
data2 = matrix(c(1:length(l), l, l_hat, l_min, l_max), ncol = 5)
data2 = data.frame(data2)
names(data2) = c('obs', 'l', 'media', 'min', 'max')
# plot1
#a = sample(1:(length(l) - 101), 1)
a = 1
f = ggplot(data2[a:(a + 100), ])
#f = ggplot(data2)
f = f + geom_line(aes(obs, l), color = 'red')
f = f + geom_line(aes(obs, media), color = 'black')
f = f + geom_line(aes(obs, min), linetype = 'dashed')
f = f + geom_line(aes(obs, max), linetype = 'dashed')
f
############### Numeric Analysis
mcmcchain_l = coda::as.mcmc( t( chain_l ) ) 
####### Geweke Statistic
# |G| > 1.96 evidencia não convergencia
CD_l = coda::geweke.diag( mcmcchain_l )
# Fração de valores que est]ao no intervalo ( -1.96 , 1.96 )
geweke_l = sum( abs( CD_l$z ) < 1.96 ) / 2e3
geweke_l
####### Fator de ineficiência (IF)
# IF = N / N_eff, onde N_eff = effective Sample Size
# Se IF >> 1, indica uma má mistura da cadeia gerada
N_eff_l = coda::effectiveSize( mcmcchain_l )
IF_l = N_new / N_eff_l
mean(IF_l)
####### MCMC error
# MCerror = sd( Variavel ) / sqrt( N_eff )
mc_error_l = round( apply( chain_l, 
                           MARGIN = 1, 
                           FUN = sd) / sqrt( N_eff_l ), 5 )
mean(mc_error_l)
# plots
mat = matrix(c(1,1,1,2,3,4), nrow = 2, byrow = TRUE)
layout( mat )
plot(l_hat, type = 'l', xlab = 'Tempo', ylab = '')
plot( CD_l$z, main = 'Geweke diagnostic', xlab = '', ylab = '' )
abline(h = -1.96, lty = 2, lwd = 5, col = 'blue')
abline(h = 1.96, lty = 2, lwd = 5, col = 'blue')
plot( IF_l, main = 'Inefficiency factors', xlab = '', ylab = '' )
plot( mc_error_l, main = 'MC errors', xlab = '', ylab = '' )
###############################################################################
###############################################################################
############################### h
h_hat = apply(chain_h, MARGIN = 1, FUN = mean)
h_min = apply(chain_h, MARGIN = 1, FUN = quantile, probs = c(0.025) )
h_max = apply(chain_h, MARGIN = 1, FUN = quantile, probs = c(0.975) )
data = matrix(c(1:T, h, h_hat, h_min, h_max), ncol = 5)
#data = matrix(c(1:T, h_hat, h_min, h_max), ncol = 4)
data = data.frame(data)
names(data) = c('obs', 'vdd', 'media', 'min','max')
#names(data) = c('obs', 'media', 'min','max')
#plots
g1 = ggplot(data[ 1:1000, ])
g1 = g1 + geom_line(aes(obs, media))
g1 = g1 + geom_line(aes(obs, vdd), color = 'red')
g1 = g1 + geom_line(aes(obs, min), linetype = 'dashed')
g1 = g1 + geom_line(aes(obs, max), linetype = 'dashed')
g1 = g1 + theme_test()
g2 = ggplot(data[ 1001:2000, ])
g2 = g2 + geom_line(aes(obs, media))
g2 = g2 + geom_line(aes(obs, vdd), color = 'red')
g2 = g2 + geom_line(aes(obs, min), linetype = 'dashed')
g2 = g2 + geom_line(aes(obs, max), linetype = 'dashed')
g2 = g2 + theme_test()
g1 / g2
############### Numeric Analysis
mcmcchain_h = coda::as.mcmc( t( chain_h ) ) 
####### Geweke Statistic
# |G| > 1.96 evidencia não convergencia
CD_h = coda::geweke.diag( mcmcchain_h )
# Fração de valores que est]ao no intervalo ( -1.96 , 1.96 )
geweke_h = sum( abs( CD_h$z ) < 1.96 ) / T
geweke_h
####### Fator de ineficiência (IF)
# IF = N / N_eff, onde N_eff = effective Sample Size
# Se IF >> 1, indica uma má mistura da cadeia gerada
N_eff_h = coda::effectiveSize( mcmcchain_h )
IF_h = N_new / N_eff_h
####### MCMC error
# MCerror = sd( Variavel ) / sqrt( N_eff )
mc_error_h = round( apply( chain_h, 
                           MARGIN = 1, 
                           FUN = sd) / sqrt( N_eff_h ), 
                    5 )
# convergence plots
par( mfrow = c(1,3) )
plot( CD_h$z, main = 'Geweke diagnostic', xlab = '', ylab = '' )
abline(h = -1.96)
abline(h = 1.96)
plot( IF_h, main = 'Inefficiency factors', xlab = '', ylab = '' )
abline(h = 1)
plot( mc_error_h, main = 'MCMC errors', xlab = '', ylab = '' )
par( mfrow = c(1,1) )
