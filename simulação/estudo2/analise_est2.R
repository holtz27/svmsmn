
x1 = x2 = 0

for( t in 1:10 ){
  x1 = x1 + Data$hmc_eff[, t] / Data$hmc_times[ t ]
  x2 = x2 + Data$rmhmc_eff[, t] / Data$rmhmc_times[ t ]
}

mean( x1 )
mean( x2 )

par( mfrow = c(1, 2) )
boxplot( x1 )
boxplot( x2 )
par( mfrow = c(1, 1) )
