
load("~/Documentos/svm_smn/svm_t/estudo2/es2_t.RData")
#load("~/Documentos/svm_smn/svm_vg/estudo2/es2_vg.RData")

x1 = x2 = 0
T = length( Data$hmc_times )

for( t in 1:T ){
  x1 = x1 + Data$hmc_eff[, t] / Data$hmc_times[ t ]
  x2 = x2 + Data$rmhmc_eff[, t] / Data$rmhmc_times[ t ]
}

x1 = x1 / t
x2 = x2 / t

data = matrix( c(x1, x2), ncol = 2 )
data = data.frame( data )
row.names( data ) = c( 'mu', 'phi', 'sigma', 'b0', 'b1', 'b2', 'v' )
colnames( data ) = c( 'hmc', 'rmhmc' )
data = round( data, 3 )
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
