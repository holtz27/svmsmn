################################################################################
#### librarys
source( 'https://raw.githubusercontent.com/holtz27/svmsmn/main/source/num_analisys.R' )
source( 'https://raw.githubusercontent.com/holtz27/svmsmn/main/source/data/t_data.R' )
require( coda )
# getwd()
path = 'svm_smn/Simulacao/Estudos_Simulacao/ts/svm_t.cpp'
Rcpp::sourceCpp( path )

# N° réplicas
n_rep = 200
vies = smse = matrix(nrow = 7, ncol = n_rep)
# Set semente
seed = 8936381
nseeds = n_rep
set.seed( seed )
seeds = sample(1:1e6, nseeds)

sumario = NULL
ruim = 0
# data setting
theta = c(mu = 0.1,
          phi = 0.95, 
          sigma = 0.05,
          b0 = 0.1,
          b1 = 0.03,
          b2 = -0.10,
          v = 5
)
T = 1e3
resultados = list( )
# mcmc setting
N = 2e4
burn = 0.5 * N
lags = 10

# sampling
for( it in 1:n_rep ){
  if( it == 1 ) time = Sys.time()
  
  cat( paste0('réplica: ', it, '\n' ) )
  
  repeat{
    x = 0
    
    x = tryCatch(
      {
        #data
        data = t_data(mu = theta[1], phi = theta[2], sigma = theta[3],
                      b0 = theta[4], b1 = theta[5], b2 = theta[6],
                      y0 = 0,
                      v = theta[7],
                      T = T,
                      seed = seeds[ it ] )
        y = data$y
        # Sampling
        samples = svm_t(N,
                        L_theta = 20, eps_theta = c( 0.5, 0.5, 0.5 ), 
                        L_b = 20, eps_b = c( 0.1, 0.1, 0.1 ), 
                        L_h = 50, eps_h = 0.015,
                        L_v = 10, eps_v = 0.5, alpha = 0.1, li = 2, ls = 40,
                        y_T = c(0, y), 
                        seed = seeds[ it ] )
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
        ################### Numeric Analysis
        resultados[[ it ]] = num_analisys(draws,
                                          burn = burn, lags = lags,
                                          names = c( 'mu', 'phi', 'sigma', 'b0', 'b1', 'b2', 'v'),
                                          digits = 4 )
        x = 0
      },
      error = function( e ){
        message('Ocorreu um erro :(')
        print( e )
        return( -1 )
      }
    )
    
    if( x == -1 ){
      seeds[ it ] = sample(1:1e6, 1)
      ruim = ruim + 1
    }else{
      break
    } 
  }
  
  
  vies[ , it ] = resultados[[ it ]][ , 1] - theta
  smse[ , it ] = vies[ , it ] ** 2
  
  cat( '\n' )
  
  if( it == n_rep ){
    # metrics calculation
    Vies = apply( vies, MARGIN = 1, mean )
    Smse = apply( smse, MARGIN = 1, mean )
    Smse = sqrt( Smse )
    # summary
    data = matrix( c(Vies, Smse), ncol = 2 )
    row.names( data ) = names = c( 'mu', 'phi', 'sigma', 'b0', 'b1', 'b2', 'v')
    colnames( data ) = c( 'vies', 'smse' )
    sumario = cbind( sumario, data )
    time = Sys.time() - time
  } 
}

resultados = list(time = time,
                  sumario = round( sumario, 4 ),
                  resultados = resultados,
                  ruim = ruim)
#save(resultados, file = 'estudo1_s.RData')
