################################################################################
#### librarys
library(ggplot2)
library(patchwork)
path = 'C:/Users/8936381/Documents/svm_smn/Simulacao/Estudos_Simulacao/normal/normal_model.cpp'
Rcpp::sourceCpp( path )

N = 3e1
samples = svmn(N,
               L_theta = 20, eps_theta = 0.5,
               L_b = 20, eps_b = 0.1,
               L_h = 50, eps_h = 0.015,
               y_T = c(0, y),
               seed = 976335 ) 

samples$time / 60
samples$acc / N
################### Save outputs
#output = list( data = y, fit = samples, N = N )
#save(output, file = 'C:/Users/8936381/Documents/Simulacao/Estudos_Simulacao/normal/ES.RDara')
#load('C:/Users/8936381/Documents/Simulacao/Estudos_Simulacao/normal/ES.RDara')
#y = output$data
#samples = output$fit
#N = output$N
#T = length( y )
################### Draws
chain_theta = samples$chain$chain_theta
chain_b     = samples$chain$chain_b
chain_h     = samples$chain$chain_h
# Transformations
chain_theta[2, ] = tanh( chain_theta[2, ] )
chain_theta[3, ] = exp( chain_theta[3, ] )
chain_b[2, ]     = tanh( chain_b[2, ] )
### burn
burn = 0
chain_theta  = chain_theta[, - c( 1:burn ) ] 
chain_b     = chain_b[, - c( 1:burn ) ]
chain_h     = chain_h[, - c( 1:burn ) ]
# Jumps
lags = 1
jumps = seq(1, N - burn, by = lags)
chain_theta = chain_theta[, jumps ]
chain_b     = chain_b[, jumps ]
chain_h     = chain_h[, jumps ]
N_new = length( jumps )
draws = rbind( chain_theta, chain_b, chain_h )
draws = rbind( draws, matrix(1, nrow = T, ncol = N_new) ) 
draws 
################################################################################
############################## Trace plots
source('source/figures.R')
trace_plots(draws[1:6, ], names = c('mu', 'phi', 'sigma',
                                    'b0', 'b1', 'b2') 
            )
################################################################################
############################## Numeric Analysis
source('source/mcmc_diagnostic.R')
mcmc_diagnostic( chain = t( draws[1:6, ] ), 
                 dec = 3,
                 names = c('mu', 'phi', 'sigma', 
                           'b0', 'b1', 'b2')
               )
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
# other plots
e.vol_hat = apply( exp( chain_h ), MARGIN = 1, FUN = mean )
e.vol_min = apply( exp( chain_h ), MARGIN = 1, FUN = quantile, probs = c(0.025) )
e.vol_max = apply( exp( chain_h ), MARGIN = 1, FUN = quantile, probs = c(0.975) )
data = matrix(c(1:T, abs(y), e.vol_hat, e.vol_min, e.vol_max), ncol = 5)
data = data.frame(data)
names(data) = c('obs', 'y.abs', 'e.hat', 'e.min','e.max')
h = ggplot(data)
h = h + geom_line(aes(obs, y.abs), color = 'grey70')
h = h + geom_ribbon(aes(x = 1:T, ymax = e.vol_max, ymin = e.vol_min), 
                    fill = 'blue' ,alpha = 0.2)
h = h + geom_line(aes(obs, e.hat), linewidth = 0.75)
h = h + theme_test() + xlab('Tempo') + ylab('|retornos|')
h
# convergence plots
par( mfrow = c(1,3) )
plot( CD_h$z, main = 'Geweke diagnostic', xlab = '', ylab = '' )
abline(h = -1.96)
abline(h = 1.96)
plot( IF_h, main = 'Inefficiency factors', xlab = '', ylab = '' )
abline(h = 1)
plot( mc_error_h, main = 'MCMC errors', xlab = '', ylab = '' )
par( mfrow = c(1,1) )
###############################################################################

###############################################################################
############################## Model Selection
# p( y | theta ) = p( y | b, theta_h, v, h ) = p( y | b, h, l )
###############################################################################
############### DIC
source('source/model_selection.R')
normal_dic = svmsmn_dic(data = y, y0 = 0, 
                        param_draws = draws[4:(nrow(draws)), ])
############### waic
normal_waic = svmsmn_waic(data = y, y0 = 0,
                          draws = draws[4:(nrow(draws)), ])
############### loo
normal_loo = svmsmn_loo(data = y, y0 = 0, 
                        draws = draws[4:(nrow(draws)), ])
