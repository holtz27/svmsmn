################################################################################
#### librarys
library( loo )
source( 'https://raw.githubusercontent.com/holtz27/svmsmn/main/source/num_analisys.R' )
source( 'https://raw.githubusercontent.com/holtz27/svmsmn/main/source/figures.R' )
source( 'https://raw.githubusercontent.com/holtz27/svmsmn/main/source/data/t_data.R' )
source( 'https://raw.githubusercontent.com/holtz27/svmsmn/main/source/model.selection.R' )


# getwd()
path = 'svm_smn/Simulacao/Estudos_Simulacao/ts/svm_t.cpp'
Rcpp::sourceCpp( path )

#data
data = t_data(mu = 1.0, phi = 0.975, sigma = 0.15,
              b0 = 0.1, b1 = 0.01, b2 = -0.05,
              y0 = 0,
              v = 20,
              T = 2e3,
              seed = 8936381)
y = data$y
h = data$h

# Sampling
N = 5e4
samples = svm_t(N,
                L_theta = 20, eps_theta = c( 0.5, 0.5, 0.5 ), 
                L_b = 20, eps_b = c( 0.1, 0.1, 0.1 ), 
                L_h = 50, eps_h = 0.015,
                L_v = 10, eps_v = 0.5, alpha = 0.01, li = 2, ls = 40,
                y_T = c(0, y), 
                seed = 67092 )
samples$time / 60
samples$acc / N
################## Save outputs
#save(samples, file = '.RDara')
#load('.RDara')
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
burn = 5e3
draws = draws[, -c( 1:burn )]
lags = 20
jumps = seq(1, N - burn, by = lags)
draws = draws[, jumps ]
# plots
trace_plots(draws[1:7, ], 
            names = c('mu', 'phi', 'sigma', 'b0', 'b1', 'b2', 'v') )
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
############################### h
h_hat = apply(chain_h, MARGIN = 1, FUN = mean)
h_min = apply(chain_h, MARGIN = 1, FUN = quantile, probs = c(0.025) )
h_max = apply(chain_h, MARGIN = 1, FUN = quantile, probs = c(0.975) )
data = matrix(c(1:nrow(chain_h), h, h_hat, h_min, h_max), ncol = 5)
data = data.frame( data )
names(data) = c('obs', 'vdd', 'media', 'min','max')
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
