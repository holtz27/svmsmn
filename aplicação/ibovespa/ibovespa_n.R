# Reading data
ibovespa = read.csv('https://raw.githubusercontent.com/holtz27/rbras/main/aplica%C3%A7%C3%B5es/%5EBVSP.csv')
ibovespa = ibovespa[, c('Date', 'Close')]
ibovespa[, 2] = as.numeric( ibovespa[, 2] ) 
ibovespa = na.omit(ibovespa)
#ibovespa = tail(ibovespa, n = 1996) #n = 3992
ibovespa = ibovespa[1:3500,]
View(ibovespa)
T = nrow(ibovespa)
log.ret = 100 * ( log( ibovespa[2:T, 2] ) - log( ibovespa[1:(T-1), 2] ) )
dates = as.Date( ibovespa[, 1] )
################################################################################
#### librarys
library( loo )
source( 'https://raw.githubusercontent.com/holtz27/svmsmn/main/source/num_analisys.R' )
source( 'https://raw.githubusercontent.com/holtz27/svmsmn/main/source/figures.R' )
source( 'https://raw.githubusercontent.com/holtz27/svmsmn/main/source/model.selection.R' )

path = 'svm_smn/Simulacao/Estudos_Simulacao/normal/svm_n.cpp'
Rcpp::sourceCpp( path )

N = 3e2
samples = svmn(N,
               L_theta = 20, eps_theta = 0.5,
               L_b = 20, eps_b = 0.1,
               L_h = 50, eps_h = 0.015,
               y_T = c(0, log.ret),
               seed = 0 
               ) 

samples$time / 60
samples$acc / N
################## Save outputs
#save(samples, file = '.RDara')
#load('.RDara')
chain_theta = samples$chain$chain_theta
chain_b = samples$chain$chain_b
chain_h = samples$chain$chain_h
# draws
draws = matrix(c( chain_theta[1, ],
                  chain_theta[2, ],
                  chain_theta[3, ],
                  chain_b[1, ],
                  chain_b[2, ],
                  chain_b[3, ]
                  ), nrow = 6, byrow = TRUE )
draws = rbind( draws, chain_h )
draws = rbind( draws, matrix(1, nrow = T, ncol = N + 1) )
############################### Convergence analysis
################### Trace plots
### burn
burn = 0
lags = 1
# plots
trace_plots(draws[1:6, ],
            burn = burn, lags = lags,
            names = c('mu', 'phi', 'sigma', 'b0', 'b1', 'b2') )
abs_plots(chain_h, dates[-1], log.ret)
################### Numeric Analysis
num_analisys(draws[1:6, ],
             burn = burn, lags = lags,
             names = c( 'mu', 'phi', 'sigma', 'b0', 'b1', 'b2'),
             digits = 4 )
###############################################################################
############################## Model Selection
# p( y | theta ) = p( y | b, theta_h, v, h ) = p( y | b, h )
# data = ( y0 , y )
# param_hat = ( b = b_hat, h = h_hat, l = l_hat )
###############################################################################
############### DIC
normal_dic = svmsmn_dic(data = log.ret, y0 = 0, 
                        param_draws = draws[4:(nrow(draws)), ])
############### waic
normal_waic = svmsmn_waic(data = log.ret, y0 = 0,
                          draws = draws[4:(nrow(draws)), ])
############### loo
normal_loo = svmsmn_loo(data = log.ret, y0 = 0, 
                        draws = draws[4:(nrow(draws)), ])
