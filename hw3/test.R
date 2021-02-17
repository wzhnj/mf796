#install.packages("NMOF")
library(NMOF)

S = 250.0
r = 0.02
tau = 6/12
q = 0.0/100.0
sigma = 0.2
theta = 0.1
v0 = 0.08
vT=0.08
k=0.7
rho=-0.4
X=250
callHestoncf(S, X, tau, r, q, v0, vT, rho, k, sigma, implVol = TRUE)
compute.implied.volatility(r=r,te= tau,s0= S,k=X,y=0,call.price= 9.276415,lower= 0.001,upper= 10)
price.bsm.option(r = r, te = tau, s0 = S, X, 
                 sigma=c(0.01,0.99,sigma) , y = 0)$call

