est <- function(Z,D,S,Y){
  n = length(Z)
  r = mean(Z==1)
  n1ss = sum(Z==1)
  n11s = sum(Z==1&D==1)
  n0ss = sum(Z==0)
  n00s = sum(Z==0&D==0)
  n10s = sum(Z==1&D==0)
  n01s = sum(Z==0&D==1)
  n111 = sum(Z==1&D==1&S==1)
  n101 = sum(Z==1&D==0&S==1)
  n011 = sum(Z==0&D==1&S==1)
  n001 = sum(Z==0&D==0&S==1)
  b111 = mean(Y[Z*D*S==1])
  b101 = mean(Y[Z*(1-D)*S==1])
  b011 = mean(Y[(1-Z)*D*S==1])
  b001 = mean(Y[(1-Z)*(1-D)*S==1])
  v111 = var(Y[Z*D*S==1])/n111
  v101 = var(Y[Z*(1-D)*S==1])/n101
  v011 = var(Y[(1-Z)*D*S==1])/n011
  v001 = var(Y[(1-Z)*(1-D)*S==1])/n001
  th1 = n11s/n1ss
  th0 = n01s/n0ss
  th11 = n111/n11s
  th10 = n101/n10s
  th01 = n011/n01s
  th00 = n001/n00s
  theta = c(th1, th0, th11, th10, th01, th00, 
            b111, b101, b011, b001)
  theta[is.nan(theta)] = 0
  V = c(r/n/th1/(1-th1), (1-r)/n/th0/(1-th0),
        r*th1/n/th11/(1-th11), r*(1-th1)/n/th10/(1-th10),
        (1-r)*th0/n/th01/(1-th01), (1-r)*(1-th0)/n/th00/(1-th00),
        v111, v101, v011, v001)
  V[is.nan(V)] = 0
  num1 = th1*th11*b111-th0*th01*b011
  den1 = th1*th11-th0*th01
  num0 = (1-th1)*th10*b101-(1-th0)*th00*b001
  den0 = (1-th1)*th10-(1-th0)*th00
  mu1 = num1/den1
  mu0 = num0/den0
  tau = mu1-mu0
  d1 = th11*(b111*den1-num1)/den1^2
  d0 = th01*(-b011*den1+num1)/den1^2
  d11 = th1*(b111*den1-num1)/den1^2
  d10 = 0 
  d01 = th0*(-b011*den1+num1)/den1^2
  d00 = 0
  d111 = th1*th11/den1
  d101 = 0
  d011 = -th0*th01/den1
  d001 = 0
  G1 = c(d1,d0,d11,d10,d01,d00,d111,d101,d011,d001)
  se1 = sqrt(sum(G1*V*G1))
  d1 = th10*(-b101*den0+num0)/den0^2
  d0 = th00*(b001*den0-num0)/den0^2
  d11 = 0
  d10 = (1-th1)*(b101*den0-num0)/den0^2 
  d01 = 0
  d00 = (1-th0)*(-b001*den0+num0)/den1^2
  d111 = 0
  d101 = (1-th1)*th10/den0
  d011 = 0
  d001 = -(1-th0)*th00/den0
  G0 = c(d1,d0,d11,d10,d01,d00,d111,d101,d011,d001)
  se0 = sqrt(sum(G0*V*G0))
  se = sqrt(sum((G1-G0)*V*(G1-G0)))
  return(list(mu1=mu1,mu0=mu0,tau=tau,se1=se1,se0=se0,se=se))
}

ests <- function(Z,D,S){
  n = length(Z)
  r = mean(Z==1)
  n1ss = sum(Z==1)
  n11s = sum(Z==1&D==1)
  n0ss = sum(Z==0)
  n00s = sum(Z==0&D==0)
  n10s = sum(Z==1&D==0)
  n01s = sum(Z==0&D==1)
  n111 = sum(Z==1&D==1&S==1)
  n101 = sum(Z==1&D==0&S==1)
  n011 = sum(Z==0&D==1&S==1)
  n001 = sum(Z==0&D==0&S==1)
  th1 = n11s/n1ss
  th0 = n01s/n0ss
  th11 = n111/n11s
  th10 = n101/n10s
  th01 = n011/n01s
  th00 = n001/n00s
  theta = c(th1, th0, th11, th10, th01, th00)
  V = c(r/n/th1/(1-th1), (1-r)/n/th0/(1-th0),
        r*th1/n/th11/(1-th11), r*(1-th1)/n/th10/(1-th10),
        (1-r)*th0/n/th01/(1-th01), (1-r)*(1-th0)/n/th00/(1-th00))
  theta[is.nan(theta)] = 0
  V[is.nan(V)] = 0
  num = th1*th11+(1-th1)*th10-th0*th01-(1-th0)*th00
  den = th1-th0
  tau = num/den
  d1 = ((th11-th10)*den-num)/den^2
  d0 = ((-th01+th00)*den+num)/den^2
  d11 = th1/den
  d10 = (1-th1)/den
  d01 = -th0/den
  d00 = -(1-th0)/den
  G = c(d1,d0,d11,d10,d01,d00)
  se = sqrt(sum(G*V*G))
  return(list(tau=tau,se=se))
}

est.itt <- function(Z,D,S,Y){
  mu1 = mean(Y[Z==1&S==1])
  mu0 = mean(Y[Z==0&S==1])
  tau = mu1-mu0
  se1 = sd(Y[Z==1&S==1])/sqrt(sum(Z==1&S==1))
  se0 = sd(Y[Z==0&S==1])/sqrt(sum(Z==0&S==1))
  se = sqrt(se1^2+se0^2)
  return(list(mu1=mu1,mu0=mu0,tau=tau,se1=se1,se0=se0,se=se))
}

est.at <- function(Z,D,S,Y){
  mu1 = mean(Y[D==1&S==1])
  mu0 = mean(Y[D==0&S==1])
  tau = mu1-mu0
  se1 = sd(Y[D==1&S==1])/sqrt(sum(D==1&S==1))
  se0 = sd(Y[D==0&S==1])/sqrt(sum(D==0&S==1))
  se = sqrt(se1^2+se0^2)
  return(list(mu1=mu1,mu0=mu0,tau=tau,se1=se1,se0=se0,se=se))
}

est.pp <- function(Z,D,S,Y){
  mu1 = mean(Y[Z==1&Z==D&S==1])
  mu0 = mean(Y[Z==0&Z==D&S==1])
  tau = mu1-mu0
  se1 = sd(Y[Z==1&S==1&Z==D])/sqrt(sum(Z==1&S==1&Z==D))
  se0 = sd(Y[Z==0&S==1&Z==D])/sqrt(sum(Z==0&S==1&Z==D))
  se = sqrt(se1^2+se0^2)
  return(list(mu1=mu1,mu0=mu0,tau=tau,se1=se1,se0=se0,se=se))
}

est.tsls <- function(Z,D,S,Y){
  n = sum(S==1)
  Dhat = predict(lm(D~Z))
  fit3 = lm(Y~Dhat, subset=(S==1))
  tau = fit3$coefficients['Dhat']
  se = summary(fit3)$coefficients['Dhat',2]
  return(list(tau=tau,se=se))
}