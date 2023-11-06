
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


par( mfrow = c(2, 7) )
# hmc
boxplot( as.numeric( Data$hmc_eff[1, ] / Data$hmc_times ) )
boxplot( as.numeric( Data$hmc_eff[2, ] / Data$hmc_times ) )
boxplot( as.numeric( Data$hmc_eff[3, ] / Data$hmc_times ) )
boxplot( as.numeric( Data$hmc_eff[4, ] / Data$hmc_times ) )
boxplot( as.numeric( Data$hmc_eff[5, ] / Data$hmc_times ) )
boxplot( as.numeric( Data$hmc_eff[6, ] / Data$hmc_times ) )
boxplot( as.numeric( Data$hmc_eff[7, ] / Data$hmc_times ) )
# rmhmc
boxplot( as.numeric( Data$rmhmc_eff[1, ] / Data$rmhmc_times ) )
boxplot( as.numeric( Data$rmhmc_eff[2, ] / Data$rmhmc_times ) )
boxplot( as.numeric( Data$rmhmc_eff[3, ] / Data$rmhmc_times ) )
boxplot( as.numeric( Data$rmhmc_eff[4, ] / Data$rmhmc_times ) )
boxplot( as.numeric( Data$rmhmc_eff[5, ] / Data$rmhmc_times ) )
boxplot( as.numeric( Data$rmhmc_eff[6, ] / Data$rmhmc_times ) )
boxplot( as.numeric( Data$rmhmc_eff[7, ] / Data$rmhmc_times ) )
par( mfrow = c(1, 1) )
