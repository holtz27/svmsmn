################################################################################
#### librarys
source( 'https://raw.githubusercontent.com/holtz27/svmsmn/main/source/num_analisys.R' )
source( 'https://raw.githubusercontent.com/holtz27/svmsmn/main/source/figures.R' )
source( 'https://raw.githubusercontent.com/holtz27/svmsmn/main/source/data/t_data.R' )

# getwd()
path = 'svm_smn/Simulacao/Estudos_Simulacao/ts/hmc/theta/hmc_theta.cpp'
Rcpp::sourceCpp( path )

#data
data = t_data(mu = 1.0, phi = 0.98, sigma = 0.15,
              b0 = 0.1, b1 = 0.01, b2 = 0.1,
              y0 = 0,
              v = 20,
              T = 2e3,
              seed = 0
              )
y = data$y
h = data$h
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
M_theta = diag(1, 3, 3)

# Sampling
N = 5e4
samples = hmc_svm_t(N,
                    L_theta = 10, eps_theta = c(0.01, 0.01, 0.01), 
                    M_theta,
                    h = h,
                    y_T = c(0, y), 
                    seed = 0 )
draws = samples$chain
############################### Convergence analysis
### burn
burn = 0
lags = 1
################### Plots
# Trace plot
trace_plots( draws, 
             burn = burn, lags = lags,
             names = c('mu', 'phi', 'sigma', 'b0', 'b1', 'b2', 'v') )

# Covariance matrix
draws = draws[, -c(1:burn)]
draws = draws[, seq(1, N - burn, by = lags)]
x = data.frame( mu = draws[1, ],
                phi = draws[2,],
                sigma = draws[3, ]
)
M_theta = cov( x )
M_theta
#solve( M_theta )
