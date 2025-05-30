densv=function(v,a,b,c){
  z=(v-c)/a
  jac=1/(a*b*(1-z*z))
  
  e=atanh(z)/b
  dens=dnorm(e, mean=-10, sd=10)*jac
  return(dens)
}


alpha=0.08
a=0.5*(40-2)
b=0.5*alpha
c=0.5*(40+2)

v=seq(2,40,0.001)
plot(v, densv(v,a,b,c), type='l')

