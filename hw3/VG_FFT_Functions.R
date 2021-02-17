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
   XX=(Eta/3)*exp(-ii*Beta*vj)*Psi_vj*(3+(-1)^J - dirac(J-1))
   ZZ = fft(XX)
   
   #calculate option prices
   Multiplier = exp(-alpha*km)/pi
   ZZ2 <- Multiplier*ZZ
   Km<-exp(km)
   
   #discard strikes that are 0 or infinity to avoid errors in interpolation
   inds <- (!is.infinite(Km)& !is.nan(Km) & (Km>1e-16) & (Km<1e16))
   px_interp <- approx(Km[inds], Re(ZZ2[inds]),method = "linear",xout=K)
   
   fft_price = Re(exp(-r*T)*px_interp$y)
   return(fft_price)
}

#define a dirac delta function
dirac<-function(n){
  y <- rep(0,length(n))
  y[n==0] = 1
  return (y)
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
  
  ii <- complex(real=0,imaginary = 1)
  
  if(nu==0){
    mu=(r-q-theta-0.5*sigma^2)
    y=exp(ii*u*(log(S0)+mu*t))*exp((ii*theta*u-0.5*sigma^2)*t)
  }
  else{
    w = log(1.0-theta*nu-0.5*sigma^2*nu)/nu
    mu = r-q+w
    y=exp(ii*u*(log(S0)+mu*t))*((1.0-ii*nu*theta*u+0.5*nu*sigma^2*u^2))^(-t/nu)
  }
  return (y)
}


plot_change_alpha<-function(alpha, n, B,K, params){
   v=seq(0.01,10,0.01)
   ii<-complex(real=0,imaginary=1)
   #calculate values of characteristic function
   Psi_vj = rep(0,length(v))
   for (zz in 1:length(v)){
      u<-(v[[zz]])
      #numer<-Re(vg_cf(u,params)*exp(-1*K*ii*u))
      numer<-vg_cf(u,params)*exp(-1*(alpha)*u)
      #print(numer)
      
      Psi_vj[[zz]] <- (numer)
   }
   
   
   return(Psi_vj)
}






alpha = 1.5
n = 14
B = 250.0
  
S0 = 100.0
r = 1.25/100
expT = 1.0
K = 100.0
q = 0.0/100.0
sig = 30.0/100.0
theta = -0.25
nu = 1.0
  
params = c(sig,nu,theta,r,expT,q,S0)

K = seq(80,120,by=2.5)
call_px=vg_fft(alpha,n,B,K,params)
print(call_px)

  
price.bsm.option(r = r, te = expT, s0 = S0, k = 80, 
                 sigma =0.1 , y = y)$call



K = 250
a=plot_change_alpha(0.01,n,B,K,params)
b=plot_change_alpha(0.25,n,B,K,params)
c=plot_change_alpha(0.5,n,B,K,params)
d=plot_change_alpha(1,n,B,K,params)
e=plot_change_alpha(1.5,n,B,K,params)
f=plot_change_alpha(5,n,B,K,params)
g=plot_change_alpha(10,n,B,K,params)
h=plot_change_alpha(20,n,B,K,params)

attach(mtcars)
v=seq(0.01,10,0.01)
par(mfrow=c(4,2))
plot(v,a, main="¦Á = 0.01")
plot(v,b, main="¦Á = 0.25")
plot(v,c, main="¦Á = 0.5")
plot(v,d, main="¦Á = 1")
plot(v,e, main="¦Á = 1.5",xlim=c(0,4))
plot(v,f, main="¦Á = 5",xlim=c(0,4))
plot(v,g, main="¦Á = 10",xlim=c(0,1))
plot(v,h, main="¦Á = 20",xlim=c(0,1))














