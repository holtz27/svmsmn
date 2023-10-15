################################################################################
#### librarys
source( 'https://raw.githubusercontent.com/holtz27/svmsmn/main/source/num_analisys.R' )
source( 'https://raw.githubusercontent.com/holtz27/svmsmn/main/source/figures.R' )
source( 'https://raw.githubusercontent.com/holtz27/svmsmn/main/source/data/vg_data.R' )
# getwd()
path = 'bruno/estudo2/hmc_svm_vg.cpp'
Rcpp::sourceCpp( path )

#data
data = vg_data(mu = 0.10, phi = 0.98, sigma = 0.15,
               b0 = 0.1, b1 = 0.03, b2 = -0.10,
               y0 = 0,
               v = 10,
               T = 2e3,
               seed = 42
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
N = 1e3
samples = hmc_svm_vg(N,
                    # theta
                    L_theta = 20, eps_theta = 0.005, M_theta = diag(1, 3, 3),
                    k_theta = 10, prec_theta = 0.1,
                    # b
                    L_b = 20, eps_b = 0.005, M_b = diag(1, 3, 3),
                    k_b = 10, prec_b = 0.05,
                    # h
                    L_h = 50, eps_h = 0.03,
                    # v
                    L_v = 20, eps_v = 0.025, M_v = 1.0,
                    k_v = 10, prec_v = 0.15, 
                    alpha = 0.1, li = 2, ls = 40,
                    y_T = c(0, y), 
                    seed = 301072
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
burn = 0
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
