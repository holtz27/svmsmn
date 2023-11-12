
load("~/Documentos/svm_smn/svm_t/estudo2/es2_t.RData")
#load("~/Documentos/svm_smn/svm_vg/estudo2/es2_vg.RData")

# Ess médio
hmc_ess = apply( Data$hmc_eff, MARGIN = 1, mean ) 
rmhmc_ess = apply( Data$rmhmc_eff, MARGIN = 1, mean ) 

# SD
hmc_sd = apply( Data$hmc_eff, MARGIN = 1, sd )  
rmhmc_sd = apply( Data$rmhmc_eff, MARGIN = 1, sd )  

# ess_min
ess_min_hmc = ess_min_rmhmc = 0
T = length( Data$hmc_times )

for( t in 1:T ){
  ess_min_hmc = ess_min_hmc + 60 * Data$hmc_eff[, t] / Data$hmc_times[ t ]
  ess_min_rmhmc = ess_min_rmhmc + 60 * Data$rmhmc_eff[, t] / Data$rmhmc_times[ t ]
}

ess_min_hmc = ess_min_hmc / T
ess_min_rmhmc = ess_min_rmhmc / T

data = NULL
data = matrix( c(hmc_ess, rmhmc_ess), ncol = 2 )
data = data.frame( data )
colnames( data ) = c( 'ess_hmc', 'ess_rmhmc' )

data = rbind( data, matrix( c(hmc_sd, rmhmc_sd), ncol = 2 ) )
data = data.frame( data )
colnames( data ) = c( 'sd_hmc', 'sd_rmhmc' )

data = round( data, 2 )
row.names( data ) = c( 'mu', 'phi', 'sigma', 'b0', 'b1', 'b2', 'v' )
data

# boxplot
par( mfrow = c(1, 7) )
# hmc
boxplot( as.numeric( Data$hmc_eff[1, ] / Data$hmc_times[ 1 ] ), 
         xlab = expression(mu), 
         ylab = 'eff/s', cex.axis = 1.3, cex.lab = 1.5 )
boxplot( as.numeric( Data$hmc_eff[2, ] / Data$hmc_times[ 2 ] ), xlab = expression(phi),
         ylab = '', cex.axis = 1.3, cex.lab = 1.5 )
boxplot( as.numeric( Data$hmc_eff[3, ] / Data$hmc_times[ 3 ] ), xlab = expression(sigma),
         ylab = '', cex.axis = 1.3, cex.lab = 1.5 )
boxplot( as.numeric( Data$hmc_eff[4, ] / Data$hmc_times[ 4 ] ), xlab = expression(b[0]),
         ylab = '', cex.axis = 1.3, cex.lab = 1.5 )
boxplot( as.numeric( Data$hmc_eff[5, ] / Data$hmc_times[ 5 ] ), xlab = expression(b[1]),
         ylab = '', cex.axis = 1.3, cex.lab = 1.5 )
boxplot( as.numeric( Data$hmc_eff[6, ] / Data$hmc_times[ 6 ] ), xlab = expression(b[2]), 
         ylab = '', cex.axis = 1.3, cex.lab = 1.5 )
boxplot( as.numeric( Data$hmc_eff[7, ] / Data$hmc_times[ 7 ] ), xlab = expression(nu),
         ylab = '', cex.axis = 1.3, cex.lab = 1.5 )
mtext( 'HMC', line = -2, outer = TRUE )
par( mfrow = c(1, 1) )

par( mfrow = c(1, 7) )
# rmhmc
boxplot( as.numeric( Data$rmhmc_eff[1, ] / Data$rmhmc_times[ 1 ] ), xlab = expression(mu), 
         ylab = 'eff/s', cex.axis = 1.3, cex.lab = 1.5 )
boxplot( as.numeric( Data$rmhmc_eff[2, ] / Data$rmhmc_times[ 2 ] ), xlab = expression(phi),
         ylab = '', cex.axis = 1.3, cex.lab = 1.5 )
boxplot( as.numeric( Data$rmhmc_eff[3, ] / Data$rmhmc_times[ 3 ] ), xlab = expression(sigma),
         ylab = '', cex.axis = 1.3, cex.lab = 1.5 )
boxplot( as.numeric( Data$rmhmc_eff[4, ] / Data$rmhmc_times[ 4 ] ), xlab = expression(b[0]),
         ylab = '', cex.axis = 1.3, cex.lab = 1.5 )
boxplot( as.numeric( Data$rmhmc_eff[5, ] / Data$rmhmc_times[ 5 ] ), xlab = expression(b[1]),
         ylab = '', cex.axis = 1.3, cex.lab = 1.5 )
boxplot( as.numeric( Data$rmhmc_eff[6, ] / Data$rmhmc_times[ 6 ] ), xlab = expression(b[2]),
         ylab = '', cex.axis = 1.3, cex.lab = 1.5 )
boxplot( as.numeric( Data$rmhmc_eff[7, ] / Data$rmhmc_times[ 7 ] ), xlab = expression(nu),
         ylab = '', cex.axis = 1.3, cex.lab = 1.5 )
mtext( 'RMHMC', line = -2, outer = TRUE )
par( mfrow = c(1, 1) )
