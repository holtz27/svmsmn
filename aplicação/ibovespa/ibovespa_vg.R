# Reading data
ibovespa = read.csv('https://raw.githubusercontent.com/holtz27/svmsmn/main/aplica%C3%A7%C3%A3o/ibovespa/ibovespa.csv')
ibovespa = ibovespa[, c('Date', 'Adj.Close')]
ibovespa[, 2] = as.numeric( ibovespa[, 2] ) 
ibovespa = na.omit(ibovespa)
#ibovespa = tail(ibovespa, n = 3992)
ibovespa = ibovespa[1:3500,]
#View(ibovespa)
T = nrow(ibovespa)
log.ret = 100 * ( log( ibovespa[2:T, 2] ) - log( ibovespa[1:(T-1), 2] ) ) 
################################################################################
#### librarys
library( loo )
source( 'https://raw.githubusercontent.com/holtz27/svmsmn/main/source/num_analisys.R' )
source( 'https://raw.githubusercontent.com/holtz27/svmsmn/main/source/figures.R' )
source( 'https://raw.githubusercontent.com/holtz27/svmsmn/main/source/model.selection.R' )

# getwd()
path = 'svm_smn/Simulacao/Estudos_Simulacao/vgamma/svm_vg.cpp'
Rcpp::sourceCpp( path )

#### Sampling
N = 1e3
samples = svm_vg(N,
                 L_theta = 20, eps_theta = c( 0.5, 0.5, 0.5 ), 
                 L_b = 20, eps_b = c( 0.1, 0.1, 0.1 ), 
                 L_h = 50, eps_h = 0.015,
                 L_v = 20, eps_v = 0.1, alpha = 0.1, li = 2, ls = 40,
                 y_T = c(0, log.ret), 
                 seed = 0 
                 )
#samples$time / 60
#samples$acc / N
################## Outputs
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
burn = 0
lags = 20
# plots
trace_plots(draws[1:7, ],
            burn = burn, lags = lags,
            names = c('mu', 'phi', 'sigma', 'b0', 'b1', 'b2', 'v') )
abs_plots( chain_h, 
           burn = burn, lags = lags, 
           date = NULL, 
           log.ret )

tail_plot( chain_l, 
           date = NULL,
           burn = burn, lags = lags,
           model_name = 'VG'
)
################### Numeric Analysis
num_analisys(draws[1:7, ],
             burn = burn, lags = lags,
             names = c( 'mu', 'phi', 'sigma', 'b0', 'b1', 'b2', 'v'),
             digits = 4 )
################################################################################
############################## Model Selection
# p( y | theta ) = p( y | b, theta_h, v, h ) = p( y | b, h )
# data = ( y0 , y )
# param_hat = ( b = b_hat, h = h_hat, l = l_hat )
################################################################################
############### dic
ts_dic = svmsmn_dic(data = log.ret, y0 = 0, 
                    param_draws = draws[ c( 4:6, 
                                            8:(T+7),
                                            (T+8):(nrow(draws)) ), ]
)
############### waic
ts_waic = svmsmn_waic(data = log.ret, y0 = 0,
                      draws = draws[ c( 4:6, 
                                        8:(T+7),
                                        (T+8):(nrow(draws)) ), ]
)
############### loo
ts_loo = svmsmn_loo(data = log.ret, y0 = 0, 
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

################## Save
#save(samples, file = '.RDara')
