#install.packages("RND")
library(RND)
library(pracma)
#define a function that performs fft on vg process
vg_fft<-function(alpha, n, B,K, params){
  r=params[4]
  T=params[5]
  S0=params[7]
  
  N=2^n
  Eta=B/N
  Lambda_Eta = 2*pi/N
  Lambda=Lambda_Eta/Eta
  
  J = 1:N
  vj = (J-1)*Eta
  m=1:N
  Beta = log(S0)-Lambda*N/2
  km=Beta+(m-1)*Lambda
  
  ii<-complex(real=0,imaginary=1)
  #calculate values of characteristic function
  Psi_vj = rep(0,length(J))
  for (zz in 1:N){
    u<-(vj[[zz]]-(alpha+1.0)*ii)
    numer<-vg_cf(u,params)
    denom <- ((alpha+ii*vj[[zz]])*(alpha+1.0+ii*vj[[zz]]))
    
    Psi_vj[[zz]] <- (numer/denom)
  }
  
  #compute fft
  XX=(Eta/2)*exp(-ii*Beta*vj)*Psi_vj*(2 - dirac(J-1))
  ZZ = fft(XX)
  
  #calculate option prices
  Multiplier = exp(-alpha*km)/pi
  ZZ2 <- Multiplier*ZZ
  Km<-exp(km)
  
  #discard strikes that are 0 or infinity to avoid errors in interpolation
  inds <- (!is.infinite(Km)& !is.nan(Km) & (Km>1e-16) & (Km<1e16))
  px_interp <- approx(Km[inds], Re(ZZ2[inds]),method = "linear",xout=K)
  
  fft_price = Re(exp(-r*T)*px_interp$y)
  return (fft_price)
}

#define a dirac delta function
dirac<-function(n){
  y <- rep(0,length(n))
  y[n==0] = 1
  return(y)
}

#define a function that computes the characteristic function for variance gamma
vg_cf<-function (u,params){
  sigma<-params[1]
  nu<-params[2]
  theta<-params[3]
  r<-params[4]
  t<-params[5]
  q<-params[6]
  S0<-params[7]
  kappa<-params[8]
  rho<-params[9]
  
  ii <- complex(real=0,imaginary = 1)
  lambda=sqrt(sigma*sigma*(u*u+ii*u)+(kappa-ii*rho*sigma*u)^2)
  omega_numer=exp(ii*u*log(S0)+ii*u*(r-q)*t+kappa*theta*t*(kappa-ii*rho*sigma*u)/(sigma^2))
  omega_denom=(cosh(lambda*t/2)+(kappa-ii*rho*sigma*u)/lambda*sinh(lambda*t/2))^(2*kappa*theta/(sigma^2))
  omega=omega_numer/omega_denom
  
  y=omega*exp(-1*(u^2+ii*u)*nu/(lambda*coth(lambda*t/2)+kappa-ii*rho*sigma*u))
  return (y)
}

plot_change_alpha<-function(alpha, n, B,K, params){
  r=params[4]
  T=params[5]
  S0=params[7]
  
  
  
  
  vj = seq(0.01,10,0.01)
  
  
  ii<-complex(real=0,imaginary=1)
  #calculate values of characteristic function
  Psi_vj = rep(0,length(vj))
  for (zz in 1:length(vj)){
    u<-(vj[[zz]]-(alpha+1.0)*ii)
    numer<-vg_cf(u,params)*exp(-r*T)
    denom <- ((alpha+ii*vj[[zz]])*(alpha+1.0+ii*vj[[zz]]))
    
    Psi_vj[[zz]] <- (numer/denom)
  }
  
  
  return(Psi_vj*exp(-alpha*vj))
}





alpha = 1.5
n = 14
B = 150.0

S0 = 250.0
r = 0.02
expT = 6/12.0
K = 100.0
q = 0.0/100.0
sig = 0.2
theta = 0.1
nu = 0.08
kappa=0.7
rho=-0.4

params = c(sig,nu,theta,r,expT,q,S0,kappa,rho)

K = 260
call_px=vg_fft(alpha,n,B,K,params)
print(call_px)
# put_px = vg_fft(-alpha,n,B,K,params)
# print(put_px)
price.bsm.option(r = r, te = expT, s0 = S0, k = 250, sigma =c(0.01,0.2,0.99) , y = 0)$call

#First plotsK = 250
# a=plot_change_alpha(0.01,n,B,K,params)
# b=plot_change_alpha(0.25,n,B,K,params)
# c=plot_change_alpha(0.5,n,B,K,params)
# d=plot_change_alpha(1,n,B,K,params)
# e=plot_change_alpha(1.5,n,B,K,params)
# f=plot_change_alpha(5,n,B,K,params)
# g=plot_change_alpha(10,n,B,K,params)
# h=plot_change_alpha(20,n,B,K,params)
# 
# 
# attach(mtcars)
# v=seq(0.01,10,0.01)
# par(mfrow=c(4,2))
# plot(v,a, main="¦Á = 0.01")
# plot(v,b, main="¦Á = 0.25")
# plot(v,c, main="¦Á = 0.5")
# plot(v,d, main="¦Á = 1")
# plot(v,e, main="¦Á = 1.5",xlim=c(0,4))
# plot(v,f, main="¦Á = 5",xlim=c(0,4))
# plot(v,g, main="¦Á = 10",xlim=c(0,1))
# plot(v,h, main="¦Á = 20",xlim=c(0,1))


# # second plots
# K = 250
# alpha=seq(0.01,1.5,0.01)
# price = rep(0,length(alpha))
# 
# for (zz in 1:length(alpha)){
# 
#   price[[zz]]=vg_fft(alpha[[zz]],n,B,K,params)
# 
# }
# plot(alpha,price)
# price
#third plots
# K = 250
# n=seq(10,20)
# alpha=1.5
# price = rep(0,length(n))
# 
# for (zz in 1:length(n)){
# 
#   price[[zz]]=vg_fft(alpha,n[[zz]],B,K,params)
# 
# }


#lo <- loess(price~n)
#plot(n,price)
# print(price)
# print(price)
#lines(predict(lo), col='red', lwd=2)

#fourth plot
# K = 250
# n=14
# alpha=1.5
# 
# B=seq(50,500,50)
# price = rep(0,length(B))
# for (zz in 1:length(B)){
# 
#   price[[zz]]=vg_fft(alpha,n,B[[zz]],K,params)
# 
# }
# 
# 
# 
# plot(B,price)
# price


#change K=260
#first plot for K=260
# K = 260
# alpha=1.5
# price = rep(0,length(alpha))
# 
# for (zz in 1:length(alpha)){
# 
#   price[[zz]]=vg_fft(alpha[[zz]],n,B,K,params)
# 
# }
# plot(alpha,price)




#second plots for K=260
# K = 260
# n=seq(1,20)
# alpha=1.5
# price = rep(0,length(n))
# 
# for (zz in 1:length(n)){
# 
#   price[[zz]]=vg_fft(alpha,n[[zz]],B,K,params)
# 
# }
# 
# 
# lo <- loess(price~n)
# plot(n[10:20],price[10:20])
# lines(predict(lo), col='red', lwd=2)

#third plot for K=260
# K = 260
# n=14
# alpha=1.5
# 
# B=seq(100,1000,50)
# price = rep(0,length(B))
# for (zz in 1:length(B)){
# 
#   price[[zz]]=vg_fft(alpha,n,B[[zz]],K,params)
# 
# }
# 
# 
# 
# plot(B,price)


#Exploring Heston Parameters Change K
# alpha = 1.5
# n = 14
# B = 700.0
# 
# S0 = 150.0
# r = 0.025
# expT = 3/12
# K = 150.0
# q = 0.0/100.0
# sig = 0.4
# theta = 0.12
# nu = 0.09
# kappa=0.5
# rho=0.25
# 
# params = c(sig,nu,theta,r,expT,q,S0,kappa,rho)
# 
# K = seq(100,200,by=10)
# call_px=vg_fft(alpha,n,B,K,params)
# print(call_px)
# implied_vol=compute.implied.volatility(r=r,te=expT,s0= S0,k=K,y=0,call.price= call_px,lower= 0.001,upper= 0.999)
# plot(K, implied_vol)
# call_px



#change expT
# alpha = 1.5
# n = 14
# B = 700.0
# 
# S0 = 150.0
# r = 0.025
# expT = 3/12
# K = 150.0
# q = 0.0/100.0
# sig = 0.4
# theta = 0.12
# nu = 0.09
# kappa=0.5
# rho=0.25
# call_px=rep(0,24)
# implied_vol=rep(0,24)
# for (m in 1:24){
# params = c(sig,nu,theta,r,m/12,q,S0,kappa,rho)
# call_px[[m]]=vg_fft(alpha,n,B,K,params)
# implied_vol[[m]]=compute.implied.volatility(r=r,te=m/12,s0= S0,k=K,y=0,call.price= call_px[[m]],lower= 0.001,upper= 0.999)
# }
# print(call_px)
# print(implied_vol)
# plot(1:24, implied_vol)
# plot(1:24, call_px)
# 
# 
# 
# 
# 
# 
# 
# 




#change sigma
# alpha = 1.5
# n = 14
# B = 700.0
# 
# S0 = 150.0
# r = 0.025
# expT = 3/12
# K = 150.0
# q = 0.0/100.0
# sig = 0.4
# theta = 0.12
# nu = 0.09
# kappa=0.5
# rho=0.25
# call_px=rep(0,9)
# implied_vol=rep(0,9)
# for (sigma in 1:9){
#   params = c(sigma/10,nu,theta,r,expT,q,S0,kappa,rho)
#   call_px[[sigma]]=vg_fft(alpha,n,B,K,params)
#   implied_vol[[sigma]]=compute.implied.volatility(r=r,te=expT,s0= S0,k=K,y=0,call.price= call_px[[sigma]],lower= 0.001,upper= 0.999)
# }
# print(call_px)
# print(implied_vol)
# plot(seq(0.1,0.9,0.1), implied_vol, xlab = "sigma")
#plot(seq(0.1,0.9,0.1), call_px)



#change nu
# alpha = 1.5
# n = 14
# B = 700.0
# 
# S0 = 150.0
# r = 0.025
# expT = 3/12
# K = 150.0
# q = 0.0/100.0
# sig = 0.4
# theta = 0.12
# nu = 0.09
# kappa=0.5
# rho=0.25
# call_px=rep(0,20)
# implied_vol=rep(0,20)
# for (nu in 1:20){
#   params = c(sig,nu/100,theta,r,expT,q,S0,kappa,rho)
#   call_px[[nu]]=vg_fft(alpha,n,B,K,params)
#   implied_vol[[nu]]=compute.implied.volatility(r=r,te=expT,s0= S0,k=K,y=0,call.price= call_px[[nu]],lower= 0.001,upper= 0.999)
# }
# print(call_px)
# print(implied_vol)
# plot(seq(0.01,0.2,0.01), implied_vol, xlab = "nu")
#plot(seq(0.1,0.9,0.1), call_px)

#change kappa
# alpha = 1.5
# n = 14
# B = 700.0
# 
# S0 = 150.0
# r = 0.025
# expT = 3/12
# K = 150.0
# q = 0.0/100.0
# sig = 0.4
# theta = 0.12
# nu = 0.09
# kappa=0.5
# rho=0.25
# call_px=rep(0,9)
# implied_vol=rep(0,9)
# for (kappa in 1:9){
#   params = c(sig,nu,theta,r,expT,q,S0,kappa/10,rho)
#   call_px[[kappa]]=vg_fft(alpha,n,B,K,params)
#   implied_vol[[kappa]]=compute.implied.volatility(r=r,te=expT,s0= S0,k=K,y=0,call.price= call_px[[kappa]],lower= 0.001,upper= 0.999)
# }
# print(call_px)
# print(implied_vol)
# plot(seq(0.1,0.9,0.1), implied_vol, xlab = "kappa")

#change rho
# alpha = 1.5
# n = 14
# B = 700.0
# 
# S0 = 150.0
# r = 0.025
# expT = 3/12
# K = 150.0
# q = 0.0/100.0
# sig = 0.4
# theta = 0.12
# nu = 0.09
# kappa=0.5
# rho=0.25
# call_px=rep(0,99)
# implied_vol=rep(0,99)
# for (rho in 1:99){
#   params = c(sig,nu,theta,r,expT,q,S0,kappa,rho/100)
#   call_px[[rho]]=vg_fft(alpha,n,B,K,params)
#   implied_vol[[rho]]=compute.implied.volatility(r=r,te=expT,s0= S0,k=K,y=0,call.price= call_px[[rho]],lower= 0.001,upper= 0.999)
# }
# print(call_px)
# print(implied_vol)
# plot(seq(0.01,0.99,0.01), implied_vol, xlab = "rho")


#change theta
# alpha = 1.5
# n = 14
# B = 700.0
# 
# S0 = 150.0
# r = 0.025
# expT = 3/12
# K = 150.0
# q = 0.0/100.0
# sig = 0.4
# theta = 0.12
# nu = 0.09
# kappa=0.5
# rho=0.25
# call_px=rep(0,99)
# implied_vol=rep(0,99)
# for (theta in 1:99){
#   params = c(sig,nu,theta/100,r,expT,q,S0,kappa,rho)
#   call_px[[theta]]=vg_fft(alpha,n,B,K,params)
#   implied_vol[[theta]]=compute.implied.volatility(r=r,te=expT,s0= S0,k=K,y=0,call.price= call_px[[theta]],lower= 0.001,upper= 0.999)
# }
# print(call_px)
# print(implied_vol)
# plot(seq(0.01,0.99,0.01), implied_vol, xlab = "theta")

attach(mtcars)

par(mfrow=c(3,2))

alpha = 1.5
n = 14
B = 700.0

#original
S0 = 150.0
r = 0.025
expT = 3/12
K = 150.0
q = 0.0/100.0
sig = 0.4
theta = 0.12
nu = 0.09
kappa=0.5
rho=0.25
call_px=rep(0,24)
implied_vol=rep(0,24)
for (m in 1:24){
params = c(sig,nu,theta,r,m/12,q,S0,kappa,rho)
call_px[[m]]=vg_fft(alpha,n,B,K,params)
implied_vol[[m]]=compute.implied.volatility(r=r,te=m/12,s0= S0,k=K,y=0,call.price= call_px[[m]],lower= 0.001,upper= 0.999)
}
print(call_px)
print(implied_vol)
plot(1:24, implied_vol,main = "original")

#increase theta
S0 = 150.0
r = 0.025
expT = 3/12
K = 150.0
q = 0.0/100.0
sig = 0.4
theta = 0.12+0.12
nu = 0.09
kappa=0.5
rho=0.25
call_px=rep(0,24)
implied_vol=rep(0,24)
for (m in 1:24){
  params = c(sig,nu,theta,r,m/12,q,S0,kappa,rho)
  call_px[[m]]=vg_fft(alpha,n,B,K,params)
  implied_vol[[m]]=compute.implied.volatility(r=r,te=m/12,s0= S0,k=K,y=0,call.price= call_px[[m]],lower= 0.001,upper= 0.999)
}
print(call_px)
print(implied_vol)
plot(1:24, implied_vol,main = "theta from 0.12 to 0.24")
#increase nu
S0 = 150.0
r = 0.025
expT = 3/12
K = 150.0
q = 0.0/100.0
sig = 0.4
theta = 0.12
nu = 0.09+0.09
kappa=0.5
rho=0.25
call_px=rep(0,24)
implied_vol=rep(0,24)
for (m in 1:24){
  params = c(sig,nu,theta,r,m/12,q,S0,kappa,rho)
  call_px[[m]]=vg_fft(alpha,n,B,K,params)
  implied_vol[[m]]=compute.implied.volatility(r=r,te=m/12,s0= S0,k=K,y=0,call.price= call_px[[m]],lower= 0.001,upper= 0.999)
}
print(call_px)
print(implied_vol)
plot(1:24, implied_vol,main = "nu from 0.09 to 0.18")
#increase kappa
S0 = 150.0
r = 0.025
expT = 3/12
K = 150.0
q = 0.0/100.0
sig = 0.4
theta = 0.12
nu = 0.09
kappa=0.5+0.3
rho=0.25
call_px=rep(0,24)
implied_vol=rep(0,24)
for (m in 1:24){
  params = c(sig,nu,theta,r,m/12,q,S0,kappa,rho)
  call_px[[m]]=vg_fft(alpha,n,B,K,params)
  implied_vol[[m]]=compute.implied.volatility(r=r,te=m/12,s0= S0,k=K,y=0,call.price= call_px[[m]],lower= 0.001,upper= 0.999)
}
print(call_px)
print(implied_vol)
plot(1:24, implied_vol,main = "kappa from 0.5 to 0.8")
#increase rho
S0 = 150.0
r = 0.025
expT = 3/12
K = 150.0
q = 0.0/100.0
sig = 0.4
theta = 0.12
nu = 0.09
kappa=0.5
rho=0.25+0.25
call_px=rep(0,24)
implied_vol=rep(0,24)
for (m in 1:24){
  params = c(sig,nu,theta,r,m/12,q,S0,kappa,rho)
  call_px[[m]]=vg_fft(alpha,n,B,K,params)
  implied_vol[[m]]=compute.implied.volatility(r=r,te=m/12,s0= S0,k=K,y=0,call.price= call_px[[m]],lower= 0.001,upper= 0.999)
}
print(call_px)
print(implied_vol)
plot(1:24, implied_vol,main = "rho from 0.25 to 0.5")
#increase theta
S0 = 150.0
r = 0.025
expT = 3/12
K = 150.0
q = 0.0/100.0
sig = 0.4
theta = 0.12+0.24
nu = 0.09
kappa=0.5
rho=0.25
call_px=rep(0,24)
implied_vol=rep(0,24)
for (m in 1:24){
  params = c(sig,nu,theta,r,m/12,q,S0,kappa,rho)
  call_px[[m]]=vg_fft(alpha,n,B,K,params)
  implied_vol[[m]]=compute.implied.volatility(r=r,te=m/12,s0= S0,k=K,y=0,call.price= call_px[[m]],lower= 0.001,upper= 0.999)
}
print(call_px)
print(implied_vol)
plot(1:24, implied_vol, main = "theta from 0.12 to 0.24")

