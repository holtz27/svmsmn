################################################################################
#### librarys
source( 'https://raw.githubusercontent.com/holtz27/svmsmn/main/source/num_analisys.R' )
source( 'https://raw.githubusercontent.com/holtz27/svmsmn/main/source/figures.R' )
source( 'https://raw.githubusercontent.com/holtz27/svmsmn/main/source/data/t_data.R' )
# getwd()
path = 'svm_smn/Simulacao/Estudos_Simulacao/ts/hmc/short_cut/hmc_svm_t.cpp'
Rcpp::sourceCpp( path )

#data
data = t_data(mu = 0.10, phi = 0.98, sigma = 0.15,
              b0 = 0.1, b1 = 0.03, b2 = -0.10,
              y0 = 0,
              v = 5,
              T = 2e3,
              seed = 0
              )
y = data$y
############################################################################## 
# Plots
library(ggplot2)
df = data.frame( Retorno = y)
g = ggplot(df) + geom_line(aes(x = 1:length( y ), y = Retorno))
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
######

# Sampling
N = 2e4
samples = hmc_svm_t(N,
                    # theta
                    L_theta = 50, eps_theta = 0.01, M_theta = diag(1, 3, 3),
                    k_theta = 2, prec_theta = 0.005,
                    # b
                    L_b = 100, eps_b = 0.005, M_b = diag(1, 3, 3),
                    k_b = 2, prec_b = 0.005,
                    # h
                    L_h = 50, eps_h = 0.015,
                    # v
                    L_v = 20, eps_v = 0.1, M_v = 1.0,
                    k_v = 2, prec_v = 0.01, 
                    alpha = 0.1, li = 2, ls = 40,
                    y_T = c(0, y), 
                    seed = 0
                    )
################## utputs
chain_theta = samples$chain$chain_theta
chain_b = samples$chain$chain_b
chain_h = samples$chain$chain_h
chain_v = samples$chain$chain_v
#chain_l = samples$chain$chain_l
# draws
draws = matrix(c( chain_theta[1, ],
                  chain_theta[2, ],
                  chain_theta[3, ],
                  chain_b[1, ],
                  chain_b[2, ],
                  chain_b[3, ],
                  chain_v[1, ]
), nrow = 7, byrow = TRUE)
#draws = rbind( draws, chain_h )
#draws = rbind( draws, chain_l )
############################### Convergence analysis
### burn
burn = 1e4
################### Numeric Analysis
num_analisys(draws[1:7, ],
             burn = burn, lags = 1,
             names = c( 'mu', 'phi', 'sigma', 'b0', 'b1', 'b2', 'v'),
             digits = 4 )
################### Plots
# Trace plot
trace_plots( draws[1:7, ], 
             burn = burn, lags = 1,
             lag.max = 400,
             names = c('mu', 'phi', 'sigma', 'b0', 'b1', 'b2', 'v') 
             )
