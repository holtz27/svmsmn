################################################################################
#### librarys
source( 'https://raw.githubusercontent.com/holtz27/svmsmn/main/source/num_analisys.R' )
source( 'https://raw.githubusercontent.com/holtz27/svmsmn/main/source/data/vg_data.R' )
# getwd()
path = 'svm_smn/Simulacao/Estudos_Simulacao/vgamma/svm_vg.cpp'
Rcpp::sourceCpp( path )

# N° réplicas
n_rep = 3
vies = smse = matrix(nrow = 7, ncol = n_rep)
# Set semente
seed = 8936381
nseeds = n_rep
set.seed( seed )
seeds = sample(1:1e6, nseeds)

sumario = NULL
# data setting
theta = c(mu = 1.0,
          phi = 0.98, 
          sigma = 0.15,
          b0 = 0.1,
          b1 = 0.01,
          b2 = -0.05,
          v = 20
)
T = 1e3
alpha = c(0.01, 0.05, 0.1)
resultados = list( )
# mcmc setting
N = 5e3
burn = 1e3
lags = 20

# sampling
for( k in 1:length( alpha ) ){
  for( it in 1:n_rep ){
    if( it == 1 ) time = Sys.time()
    #data
    data = vgamma_data(mu = theta[1], phi = theta[2], sigma = theta[3],
                       b0 = theta[4], b1 = theta[5], b2 = theta[6],
                       y0 = 0,
                       v = theta[7],
                       T = T,
                       seed = seeds[ it ] )
    y = data$y
    
    cat( paste0('réplica: ', it, '; ', 'alpha: ', alpha[ k ],'\n' ) )
    
    # Sampling
    #N = 
    samples = svm_vg(N,
                     L_theta = 20, eps_theta = c( 0.5, 0.5, 0.5 ), 
                     L_b = 20, eps_b = c( 0.1, 0.1, 0.1 ), 
                     L_h = 50, eps_h = 0.015,
                     L_v = 10, eps_v = 0.1, alpha = alpha[ k ], li = 0, ls = 40,
                     y_T = c(0, y), 
                    seed = 0 )
    chain_theta = samples$chain$chain_theta
    chain_b = samples$chain$chain_b
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
    ############################### Convergence analysis
    ################### Trace plots
    ### burn
    draws = draws[, -c( 1:burn )]
    jumps = seq(1, N - burn, by = lags)
    draws = draws[, jumps ]
    ################### Numeric Analysis
    #resultados[[ (k - 1) * n_rep + it ]] = num_analisys(draws, 
    #                                                    names = c( 'mu', 'phi', 'sigma', 'b0', 'b1', 'b2', 'v'),
    #                                                    digits = 4 )
    #vies[ , it ] = resultados[[ (k - 1) * n_rep + it ]][ , 1] - theta
    #smse[ , it ] = vies[ , it ] ** 2
    resultados[[ (k - 1) * n_rep + it ]] = apply( draws, MARGIN = 1, FUN = mean )
    vies[ , it ] = resultados[[ (k - 1) * n_rep + it ]] - theta
    smse[ , it ] = vies[ , it ] ** 2
    
    cat( '\n' )
    
    if( it == n_rep ){
      # metrics calculation
      Vies = apply( vies, MARGIN = 1, mean )
      Smse = apply( smse, MARGIN = 1, mean )
      Smse = sqrt( Smse )
      # summary
      data = matrix( c(Vies, Smse), ncol = 2 )
      data = cbind( data , rep(alpha[ k ], 7))
      row.names( data ) = names = c( 'mu', 'phi', 'sigma', 'b0', 'b1', 'b2', 'v')
      colnames( data ) = c( 'vies', 'smse', 'alpha' )
      sumario = cbind( sumario, data )
      time = Sys.time() - time
    } 
  }
}

time 
round( sumario, 4 )
resultados
