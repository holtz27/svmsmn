################################################################################
#### librarys
source( 'https://raw.githubusercontent.com/holtz27/svmsmn/main/source/data/t_data.R' )
# getwd()
path = 'svm_smn/Simulacao/Estudos_Simulacao/ts/hmc/short_cut/hmc_svm_t.cpp'
Rcpp::sourceCpp( path )
path = 'svm_smn/Simulacao/Estudos_Simulacao/ts/svm_t.cpp'
Rcpp::sourceCpp( path )

M_theta = diag(c(1, 1, 1), 3)
M_b = diag(1, 3, 3)
M_v = 1.0

n_rep = 10
# Set seed
seed = 1836398
nseeds = n_rep
set.seed( seed )
seeds = sample(1:1e6, nseeds)

hmc_times = rmhmc_times = matrix(nrow = 1, ncol = n_rep)
hmc_eff_s = hmc_eff = rmhmc_eff_s = rmhmc_eff = matrix(nrow = 7, ncol = n_rep)

N = 3e4
burn = 1e4

for(i in 1:n_rep){
  #data
  data = t_data(mu = 1.0, phi = 0.98, sigma = 0.15,
                b0 = 0.1, b1 = 0.01, b2 = -0.05,
                y0 = 0,
                v = 20,
                T = 2e3,
                seed = seeds[ i ]
  )
  y = data$y
  ##################################### hmc
  # Sampling
  samples = hmc_svm_t(N,
                      # theta
                      L_theta = 50, eps_theta = 0.01, M_theta,
                      k_theta = 2, prec_theta = 0.005,
                      # b
                      L_b = 50, eps_b = 0.005, M_b,
                      k_b = 2, prec_b = 0.005,
                      # h
                      L_h = 50, eps_h = 0.015,
                      # v
                      L_v = 20, eps_v = 0.05, M_v,
                      k_v = 2, prec_v = 0.01, alpha = 1.0, li = 2, ls = 40,
                      y_T = c(0, y), 
                      seed = 0,
                      eps_random = FALSE
  )
  chain_theta = samples$chain$chain_theta
  chain_b = samples$chain$chain_b
  chain_h = samples$chain$chain_h
  chain_v = samples$chain$chain_v
  # draws
  draws = matrix(c( chain_theta[1, ],
                    chain_theta[2, ],
                    chain_theta[3, ],
                    chain_b[1, ],
                    chain_b[2, ],
                    chain_b[3, ],
                    chain_v[1, ]
  ), nrow = 7, byrow = TRUE)
  draws = draws[, -c( 1:burn )]
  chain = coda::as.mcmc( t( draws ) )
  hmc_eff[ ,i ] = coda::effectiveSize( chain )
  hmc_times[ ,i ] = samples$time
  
  ##################################### rmhmc
  samples = svm_t(N,
                  L_theta = 20, eps_theta = c( 0.5, 0.5, 0.5 ), 
                  L_b = 20, eps_b = c( 0.1, 0.1, 0.1 ), 
                  L_h = 50, eps_h = 0.015,
                  L_v = 20, eps_v = 0.5, alpha = 1.0, li = 2, ls = 40,
                  y_T = c(0, y), 
                  seed = 0
  )
  chain_theta = samples$chain$chain_theta
  chain_b = samples$chain$chain_b
  chain_h = samples$chain$chain_h
  chain_v = samples$chain$chain_v
  # draws
  draws = matrix(c( chain_theta[1, ],
                    chain_theta[2, ],
                    chain_theta[3, ],
                    chain_b[1, ],
                    chain_b[2, ],
                    chain_b[3, ],
                    chain_v
  ), nrow = 7, byrow = TRUE)
  draws = draws[, -c( 1:burn )]
  chain = coda::as.mcmc( t( draws ) )
  rmhmc_eff[ ,i ] = coda::effectiveSize( chain )
  rmhmc_times[ ,i ] = samples$time
  
  # summary
  hmc_eff_s[, i ] = hmc_eff[, i] / hmc_times[ i ]
  rmhmc_eff_s[, i ] = rmhmc_eff[, i] / rmhmc_times[ i ]
}

data = matrix(c(apply(hmc_eff_s, 1, mean),
                apply(rmhmc_eff_s, 1, mean)),ncol = 2)
data = round(data, 3)
row.names(data) = c('mu','phi','sigma','b0','b1','b2','v')
colnames(data) = c('hmc', 'rmhmc')
data




