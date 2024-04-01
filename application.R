library(causalweight)

dat = causalweight::JC
dim(dat)
Z = dat$assignment
D1 = dat$trainy1
D2 = dat$trainy2
D = ((D1+D2)+abs(D1-D2))/2
S1 = as.numeric(dat$earnq4>0)
S2 = as.numeric(dat$earny2>0)
S3 = as.numeric(dat$earny3>0)
S4 = as.numeric(dat$earny4>0)
Y1 = log(dat$earnq4)
Y2 = log(dat$earny2)
Y3 = log(dat$earny3)
Y4 = log(dat$earny4)
Y1[is.infinite(Y1)] = 0
Y2[is.infinite(Y2)] = 0
Y3[is.infinite(Y3)] = 0
Y4[is.infinite(Y4)] = 0

pa = mean(D[Z==0])
pn = 1-mean(D[Z==1])
pc = 1-pa-pn
round(c(pa,pc,pn), 3)

dat1 = round(c(mean(S1[Z==1&D==1]), mean(S1[Z==1&D==0]), 
               mean(S1[Z==0&D==1]), mean(S1[Z==0&D==0])), 4)
dat2 = round(c(mean(S2[Z==1&D==1]), mean(S2[Z==1&D==0]), 
               mean(S2[Z==0&D==1]), mean(S2[Z==0&D==0])), 4)
dat3 = round(c(mean(S3[Z==1&D==1]), mean(S3[Z==1&D==0]), 
               mean(S3[Z==0&D==1]), mean(S3[Z==0&D==0])), 4)
dat4 = round(c(mean(S4[Z==1&D==1]), mean(S4[Z==1&D==0]), 
               mean(S4[Z==0&D==1]), mean(S4[Z==0&D==0])), 4)
s.tab = rbind(dat1,dat2,dat3,dat4)
dat1 = round(c(mean(Y1[Z==1&D==1&S1==1]), mean(Y1[Z==1&D==0&S1==1]), 
               mean(Y1[Z==0&D==1&S1==1]), mean(Y1[Z==0&D==0&S1==1])), 4)
dat2 = round(c(mean(Y2[Z==1&D==1&S2==1]), mean(Y2[Z==1&D==0&S2==1]), 
               mean(Y2[Z==0&D==1&S2==1]), mean(Y2[Z==0&D==0&S2==1])), 4)
dat3 = round(c(mean(Y3[Z==1&D==1&S3==1]), mean(Y3[Z==1&D==0&S3==1]), 
               mean(Y3[Z==0&D==1&S3==1]), mean(Y3[Z==0&D==0&S3==1])), 4)
dat4 = round(c(mean(Y4[Z==1&D==1&S4==1]), mean(Y4[Z==1&D==0&S4==1]), 
               mean(Y4[Z==0&D==1&S4==1]), mean(Y4[Z==0&D==0&S4==1])), 4)
y.tab = rbind(dat1,dat2,dat3,dat4)
sy.tab = cbind(s.tab, y.tab)
xtable(sy.tab, digits=3)

## effect on S
fit1 = ests(Z,D1,S1)
fit2 = ests(Z,D,S2)
fit3 = ests(Z,D,S3)
fit4 = ests(Z,D,S4)
tau = c(fit1$tau,fit2$tau,fit3$tau,fit4$tau)
se = c(fit1$se,fit2$se,fit3$se,fit4$se)
p = (1-pnorm(abs(tau/se)))*2
ci_l = tau-1.96*se
ci_u = tau+1.96*se
xtable(rbind(tau,se,p), digits=3)
plot(1:4, tau, type='o', lwd=2, ylim=c(-.4,.4), xlab='Year', ylab='Effect',
     main='Effect on employment with year')
points(1:4, ci_l, type='o', lty=5)
points(1:4, ci_u, type='o', lty=5)
abline(h=0, lty=3)

## effect on Y
fit1 = est(Z,D1,S1,Y1)
fit2 = est(Z,D,S2,Y2)
fit3 = est(Z,D,S3,Y3)
fit4 = est(Z,D,S4,Y4)
mu1 = c(fit1$mu1,fit2$mu1,fit3$mu1,fit4$mu1)
mu0 = c(fit1$mu0,fit2$mu0,fit3$mu0,fit4$mu0)
tau = c(fit1$tau,fit2$tau,fit3$tau,fit4$tau)
se1 = c(fit1$se1,fit2$se1,fit3$se1,fit4$se1)
se0 = c(fit1$se0,fit2$se0,fit3$se0,fit4$se0)
se = c(fit1$se,fit2$se,fit3$se,fit4$se)
p = (1-pnorm(abs(tau/se)))*2
ci_l = tau-1.96*se
ci_u = tau+1.96*se
xtable(cbind(mu1,se1,mu0,se0,tau,se,p), digits=3)
plot(1:4, tau, type='o', lwd=2, ylim=c(-.4,.4), xlab='Year', ylab='Effect',
     main='Effect on earnings with year')
points(1:4, ci_l, type='o', lty=5)
points(1:4, ci_u, type='o', lty=5)
abline(h=0, lty=3)

# 2sls, intention-to-treat, as-treated, per-protocol
fit1 = est.tsls(Z,D1,S1,Y1)
fit2 = est.tsls(Z,D,S2,Y2)
fit3 = est.tsls(Z,D,S3,Y3)
fit4 = est.tsls(Z,D,S4,Y4)
mu1 = c(fit1$mu1,fit2$mu1,fit3$mu1,fit4$mu1)
mu0 = c(fit1$mu0,fit2$mu0,fit3$mu0,fit4$mu0)
tau = c(fit1$tau,fit2$tau,fit3$tau,fit4$tau)
se1 = c(fit1$se1,fit2$se1,fit3$se1,fit4$se1)
se0 = c(fit1$se0,fit2$se0,fit3$se0,fit4$se0)
se = c(fit1$se,fit2$se,fit3$se,fit4$se)
p = (1-pnorm(abs(tau/se)))*2
ci_l = tau-1.96*se
ci_u = tau+1.96*se
xtable(cbind(mu1,se1,mu0,se0,tau,se,p), digits=3)
plot(1:4, tau, type='o', lwd=2, ylim=c(-.4,.4), xlab='Year', ylab='Effect',
     main='Effect on earnings with year')
points(1:4, ci_l, type='o', lty=5)
points(1:4, ci_u, type='o', lty=5)
abline(h=0, lty=3)

