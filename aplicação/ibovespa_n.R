
# Reading data
ibovespa = read.csv('https://raw.githubusercontent.com/holtz27/rbras/main/aplica%C3%A7%C3%B5es/%5EBVSP.csv')
ibovespa = ibovespa[, c('Date', 'Close')]
ibovespa[, 2] = as.numeric( ibovespa[, 2] ) 
ibovespa = na.omit(ibovespa)
ibovespa = tail(ibovespa, n = 3992)
#View(ibovespa)
T = nrow(ibovespa)
log.ret = 100 * ( log( ibovespa[2:T, 2] ) - log( ibovespa[1:(T-1), 2] ) ) 
################################################################################
#### librarys
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
               seed = 976335 ) 

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
############################### Convergence analysis
################### Trace plots
### burn
burn = 0
draws = draws[, -c( 1:burn )]
lags = 1
jumps = seq(1, N - burn, by = lags)
draws = draws[, jumps ]
# plots
trace_plots(draws[1:6, ], 
            names = c('mu', 'phi', 'sigma', 'b0', 'b1', 'b2') )
################### Numeric Analysis
num_analisys(draws[1:6, ], 
             names = c( 'mu', 'phi', 'sigma', 'b0', 'b1', 'b2'),
             digits = 4 )
