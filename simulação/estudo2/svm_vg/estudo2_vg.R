################################################################################
#### librarys
source( 'https://raw.githubusercontent.com/holtz27/svmsmn/main/source/data/vg_data.R' )
# getwd()

Rcpp::sourceCpp( 'bruno/estudo2/hmc_svm_vg.cpp' )
Rcpp::sourceCpp( 'bruno/estudo2/svm_vg.cpp' )

n_rep = 2
# Set seed
seed = 7836398
nseeds = n_rep
set.seed( seed )
seeds = sample(1:1e6, nseeds)

hmc_times = rmhmc_times = matrix(nrow = 1, ncol = n_rep)
hmc_eff_s = hmc_eff = rmhmc_eff_s = rmhmc_eff = matrix(nrow = 7, ncol = n_rep)

N = 2e2
burn = 1e2

for(i in 1:n_rep){
  #data
  data = vg_data(mu = 0.1, phi = 0.98, sigma = 0.15,
                 b0 = 0.1, b1 = 0.03, b2 = -0.1,
                 y0 = 0,
                 v = 10,
                 T = 2e3,
                 seed = seeds[ i ]
  )
  y = data$y
  ##################################### hmc
  # Sampling
  samples = hmc_svm_vg(N,
                      # theta
                      L_theta = 10, eps_theta = 0.01, M_theta = diag(1, 3, 3),
                      k_theta = 10, prec_theta = 0.015,
                      # b
                      L_b = 10, eps_b = 0.0025, M_b = diag(1, 3, 3),
                      k_b = 10, prec_b = 0.015,
                      # h
                      L_h = 50, eps_h = 0.015,
                      # v
                      L_v = 10, eps_v = 0.05, M_v = 1.0,
                      k_v = 10, prec_v = 0.01, alpha = 0.1, li = 2, ls = 40,
                      y_T = c(0, y), 
                      seed = seeds[ i ] + 1
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
  samples = svm_vg(N,
                  L_theta = 20, eps_theta = c( 0.5, 0.5, 0.5 ), 
                  L_b = 20, eps_b = c( 0.1, 0.1, 0.1 ), 
                  L_h = 50, eps_h = 0.015,
                  L_v = 20, eps_v = 0.5, alpha = 0.1, li = 2, ls = 40,
                  y_T = c(0, y), 
                  seed = seeds[ i ] + 2
  )
  chain_theta = samples$chain$chain_theta
  chain_b = samples$chain$chain_b
  #chain_h = samples$chain$chain_h
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

Data = list( hmc_times = hmc_times, hmc_eff = hmc_eff,
             rmhmc_times = rmhmc_times, rmhmc_eff = rmhmc_eff
)

save(Data, file = 'estudo2_t.RData')
rm(list = ls())

