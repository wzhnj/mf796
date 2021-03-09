a=function(a,b,c,d){
  return (a*b+c*d^2)
}



object = function(x){
  b=x[1]
  d=x[2]
  a=3
  
  return ((a(a,b,4,d)-42)^2+(a(a,b,5,d)-51)^2)
}

optimx(c(1,0),
       object,
       #control=list( save.failures=TRUE, trace=0),
       upper = c(5,5),lower = c(0,0),
       method = "L-BFGS-B"
       )
