library(xtable)

generatedata <- function(n,homo=TRUE,PI=TRUE){
  Z = rbinom(n,1,0.5)
  D0 = rbinom(n,1,0.3)
  D1 = D0 + rbinom(n,1,0.4)
  D1[D1>1] = 1
  S0 = rbinom(n,1,0.3+D0*0.2+(D1-D0)*0.2)
  S1 = rbinom(n,1,0.3+D0*0.3+(D1-D0)*0.3)
  if (PI){
    mu0 = 1
    mu1 = 2
  } else {
    mu0 = 1 + S1
    mu1 = 2 + S0/2
  }
  Y0 = S0*rnorm(n,mu0,0.8)
  Y1 = S1*rnorm(n,mu1,1)
  if (!homo){
    Y0[D1==D0] = (Y0 + rnorm(n,-0.5,0.2))[D1==D0] #nt
    Y1[D1==D0] = (Y1 + rnorm(n,0.8,0.3))[D1==D0] #at
  }
  D = Z*D1 + (1-Z)*D0
  S = D*S1 + (1-D)*S0
  Y = D*Y1 + (1-D)*Y0
  return(list(Z=Z,D=D,S=S,Y=Y,tau=mean((Y1-Y0)[D1>D0&S1*S0==1])))
}


bs = NULL
for (n in c(500,2000,8000)){
bs.n = NULL
for (PI in c(TRUE,FALSE)){
for (homo in c(TRUE,FALSE)){
tau.pace = tau.tsls = se.pace = se.tsls = tau.true = NULL
for (b in 1:1000){
  set.seed(b)
  datb = generatedata(n,homo,PI)
  Z = datb$Z
  D = datb$D
  S = datb$S
  Y = datb$Y
  true = datb$tau
  fit.pace = est(Z,D,S,Y)
  fit.tsls = est.tsls(Z,D,S,Y)
  tau.true = append(tau.true, true)
  tau.pace = append(tau.pace, fit.pace$tau-true)
  tau.tsls = append(tau.tsls, fit.tsls$tau-true)
  se.pace = append(se.pace, fit.pace$se)
  se.tsls = append(se.tsls, fit.tsls$se)
}
cv.pace = mean((tau.pace+1.96*se.pace>=0)*
                 (tau.pace-1.96*se.pace<=0))
cv.tsls = mean((tau.tsls+1.96*se.tsls>=0)*
                 (tau.tsls-1.96*se.tsls<=0))
bs.ab = rbind(c(mean(tau.pace), mean(tau.tsls)), 
                c(sd(tau.pace), sd(tau.tsls)),
                c(mean(se.pace), mean(se.tsls)),
                c(cv.pace, cv.tsls))
bs.n = cbind(bs.n, bs.ab)
}
}
bs = rbind(bs, bs.n)
}
rownames(bs) = rep(c('Bias','SD','SE','Coverage'),3)
round(bs,4)
xtable(bs, digits=3)