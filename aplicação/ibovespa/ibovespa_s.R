# Reading data
ibovespa = read.csv('https://raw.githubusercontent.com/holtz27/svmsmn/main/aplica%C3%A7%C3%A3o/ibovespa/ibovespa.csv')
ibovespa = ibovespa[, c('Date', 'Close')]
ibovespa[, 2] = as.numeric( ibovespa[, 2] ) 
ibovespa = na.omit(ibovespa)
#View(ibovespa)
T = nrow(ibovespa)
log.ret = 100 * ( log( ibovespa[2:T, 2] ) - log( ibovespa[1:(T-1), 2] ) )
dates = as.Date( ibovespa[, 1] )
T = length( log.ret )
par( mfrow = c(1 , 2) )
plot( log.ret, type = 'l' )
hist( log.ret, breaks = 40 )
par( mfrow = c(1, 1) )
################################################################################
#### librarys
library( loo )
source( 'https://raw.githubusercontent.com/holtz27/svmsmn/main/source/num_analisys.R' )
source( 'https://raw.githubusercontent.com/holtz27/svmsmn/main/source/figures.R' )
source( 'https://raw.githubusercontent.com/holtz27/svmsmn/main/source/model.selection.R' )

# getwd()
path = 'svm_s/svm_s.cpp'
Rcpp::sourceCpp( path )

# Sampling
N = 1e5
samples = svm_s(N,
                L_theta = 20, eps_theta = c( 0.5, 0.5, 0.5 ), 
                L_b = 20, eps_b = c( 0.1, 0.1, 0.1 ), 
                L_h = 50, eps_h = 0.015,
                y_T = c(0, log.ret), 
                seed = 893
                )
samples$time / 60
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
burn = 2e4
lags = 20
# plots
trace_plots(draws[1:7, ],
            burn = burn, lags = lags,
            names = c( expression( mu ), 
                       expression( phi ), 
                       expression( sigma ), 
                       expression( b[0] ), 
                       expression( b[1] ), 
                       expression( b[2] ), 
                       expression( v )
                       )
            )
abs_plots( chain_h, burn = burn, lags = lags, date = dates[-1], log.ret )
tail_plot( chain_l, date = dates[-1], burn = burn, lags = lags, model_name = 't')
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
slash_dic = svmsmn_dic(data = log.ret, y0 = 0, 
                       param_draws = draws[ c( 4:6, 
                                               8:(T+7),
                                               (T+8):(nrow(draws)) ), ]
)
############### waic
slash_waic = svmsmn_waic(data = log.ret, y0 = 0,
                         draws = draws[ c( 4:6, 
                                           8:(T+7),
                                           (T+8):(nrow(draws)) ), ]
)
############### loo
slash_loo = svmsmn_loo(data = log.ret, y0 = 0, 
                       draws = draws[ c( 4:6, 
                                         8:(T+7),
                                         (T+8):(nrow(draws)) ), ],
                       cores = 4
)
slash_criterium = matrix( c(slash_dic, 
                            slash_loo$estimates[3,1],
                            slash_waic$estimates[3,1]), nrow = 1)
row.names( slash_criterium ) = c('slash')
colnames( slash_criterium ) = c('dic', 'loo', 'waic')
slash_criterium
################## Save
save( samples, file = 'svm_s/ibovespa_s.RData' )
#load( file = 'svm_s/ibovespa_s.RData' )
rm( list = ls() )
