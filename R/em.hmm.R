em.hmm <-
function(x, L=2, maxiter=1000, est.null=FALSE)
{
#####################################################################################

## USAGE
 # em.hmm(x, L, maxiter)

## ARGUMENTS
 # x: the observed Z values
 # L: num. of components in the non-null mixture
 # maxiter: the maximum number of iterations

## DETAILS
 # em.hmm calculates the MLE for a HMM model with hidden states being 0/1.
 # the distribution of state 0 is assumed to be normal
 # the distribution of state 1 is assumed to be a normal mixture with L components

## VALUES
 # fuction 'em.hmm' gives the MLE of model parameters and Lfdr estimates for a HMM 
 # pii: the initial state distribution
 # A=(a00 a01\\ a10 a11): transition matrix
 # f0: the null distribution
 # pc, i from 1 to L: 
 # --probability weights of each component in the non-null mixture 
 # f1: an L by 2 matrix
 # --specifying the dist. of each component in the non-null mixture 
 # LIS: the lis statistics
 # BIC: BIC score for the estimated model
 # ni: number of iterations

####################################################################################

NUM<-length(x)
# precision tolerance level
ptol<-1e-4
niter<-0

# Assuming it will converge 
converged=TRUE

### initializing model parameters

########
# L=1
########

if (L==1)
{

pii.new<-c(0.5, 0.5)
A.new<-matrix(c(0.8, 0.2, 0.4, 0.6), 2, 2, byrow=TRUE)
f0.new<-c(2, 1)
f1.new<-c(4, 1)

diff<-1

### The E-M Algorithm

while(diff>ptol && niter<maxiter)
{

niter<-niter+1

pii.old<-pii.new
A.old<-A.new
f0.old<-f0.new
f1.old<-f1.new

## updating the weights and probabilities of hidden states

bwfw.res<-bwfw1.hmm(x, pii.old, A.old, f0.old, f1.old)

# the hidden states probabilities
gamma<-bwfw.res$pr
# the transition variables
dgamma<-bwfw.res$ts
# the rescaling variables
c0<-bwfw.res$rescale

## updating the parameter estimates

# a. initial state distribution

for (i in 0:1)
{
  pii.new[i+1]<-gamma[1, i+1]
}

# b. transition matrix

for (i in 0:1)
{
  for (j in 0:1)
  { 
    q1<-sum(dgamma[i+1, j+1, ])
    q2<-sum(gamma[1:(NUM-1), i+1])
    A.new[i+1, j+1]<-q1/q2  
  }
}

# c. null distribution

# (i). mean
  q5<-sum(gamma[, 1]*x)
  mu0<-q5/sum(gamma[, 1])

# (ii). sd
  q6<-sum(gamma[, 1]*(x-mu0)*(x-mu0))
  sd0<-sqrt(q6/sum(gamma[, 1]))

f0.new<-c(mu0, sd0)

if (!est.null)
{
 f0.new<-c(0, 1)
}

# c. non-null distribution 

q1<-sum(gamma[, 2])
q2<-sum(gamma[, 2]*x)
mu1<-q2/q1
q3<-sum(gamma[, 2]*(x-mu1)*(x-mu1))
sd1<-sqrt(q3/q1)
f1.new<-c(mu1, sd1)

df1<-abs(A.old-A.new)
df2<-abs(f1.old-f1.new)
diff<-max(df1, df2)

if (is.na(diff)) {
   converged=FALSE;
   break;
}

}

# f. the final local fdr statistic
 lfdr<-gamma[, 1]
# g. the loglikelihood 
if (converged) {
 logL<--sum(log(c0))
# let f(L) be the number of parameters needed to be estimated in an HMM.
# if you use N(0,1) as the null distribution (default), f(L)=#(f0)+#(p)+#(f1)+#(A, pi)=0+(L-1)+2*L+3=3L+2
#                                            Otherwise, f(L)=#(f0)+#(p)+#(f1)+#(A, pi)=2+(L-1)+2*L+3=3L+4
# let m be the number of observations
# BIC<-logL-f(L)*log(m)/2 
if (est.null) {
	BIC<-logL-(3*L+4)*log(NUM)/2 
} else {
	BIC<-logL-(3*L+2)*log(NUM)/2 
}

# h. return the results of the E-M algorithm
 em.var<-list(pii=pii.new, A=A.new, f0=f0.new, f1=f1.new, LIS=lfdr, logL=logL, BIC=BIC, ni=niter, converged=converged) 
} else {
 BIC<- logL<- (-Inf)
 em.var<-list(pii=pii.old, A=A.old, f0=f0.old, f1=f1.old, LIS=lfdr, logL=logL, BIC=, ni=niter, converged=converged)
}

}

#######
# L>1
#######

else if (L>1)
{

pii.new<-c(0.5, 0.5)
A.new<-matrix(c(0.95, 0.05, 0.5, 0.5), 2, 2, byrow=TRUE)
pc.new<-rep(1, L)/L
mus<-seq(from=-1, by=1.5, length=L)
sds<-rep(1, L)
f0.new<-c(2, 1)
f1.new<-cbind(mus, sds)

diff<-1

### The E-M Algorithm

while(diff>ptol && niter<maxiter)
{

niter<-niter+1

pii.old<-pii.new
A.old<-A.new
pc.old<-pc.new
f0.old<-f0.new
f1.old<-f1.new

## updating the weights and probabilities of hidden states

bwfw.res<-bwfw.hmm(x, pii.old, A.old, pc.old, f0.old, f1.old)

# the hidden states probabilities
gamma<-bwfw.res$pr
# the transition variables
dgamma<-bwfw.res$ts
# the weight variables
omega<-bwfw.res$wt
# the rescaling variables
c0<-bwfw.res$rescale

## updating the parameter estimates

# a. initial state distribution

for (i in 0:1)
{
  pii.new[i+1]<-gamma[1, i+1]
}

# b. transition matrix

for (i in 0:1)
{
  for (j in 0:1)
  { 
    q1<-sum(dgamma[i+1, j+1, ])
    q2<-sum(gamma[1:(NUM-1), i+1])
    A.new[i+1, j+1]<-q1/q2  
  }
}

# c. null distribution

# (i). mean
  q5<-sum(gamma[, 1]*x)
  mu0<-q5/sum(gamma[, 1])

# (ii). sd
  q6<-sum(gamma[, 1]*(x-mu0)*(x-mu0))
  sd0<-sqrt(q6/sum(gamma[, 1]))

f0.new<-c(mu0, sd0)

if (!est.null)
{
 f0.new<-c(0, 1)
}

# d. non-null mixture distribution

# initializing the vectors of means and variances

mus<-1:L
sds<-1:L

for (c in 1:L)
{

# (i). probability weights
  q1<-sum(omega[, c])
  q2<-sum(gamma[, 2])
  pc.new[c]<-q1/q2
  
# (ii). means
  q3<-sum(omega[, c]*x)
  mus[c]<-q3/q1

# (iii). sds
  q4<-sum(omega[, c]*(x-mus[c])*(x-mus[c]))
  sds[c]<-sqrt(q4/q1)

}
# the non-null mixture distribution
f1.new<-cbind(mus, sds)

df1<-abs(A.old-A.new)
df2<-abs(f1.old-f1.new)
diff<-max(df1, df2)

if (is.na(diff)) {
   converged=FALSE;
   break;
}

}

# f. the final local fdr statistic
 lfdr<-gamma[, 1]
# g. the loglikelihood 
if (converged) {
 logL<--sum(log(c0))
# let f(L) be the number of parameters needed to be estimated in an HMM.
# if you use N(0,1) as the null distribution (default), f(L)=#(f0)+#(p)+#(f1)+#(A, pi)=0+(L-1)+2*L+3=3L+2
#                                            Otherwise, f(L)=#(f0)+#(p)+#(f1)+#(A, pi)=2+(L-1)+2*L+3=3L+4
# let m be the number of observations
# BIC<-logL-f(L)*log(m)/2 
if (est.null) {
	BIC<-logL-(3*L+4)*log(NUM)/2 
} else {
	BIC<-logL-(3*L+2)*log(NUM)/2 
}
# h. return the results of the E-M algorithm
 em.var<-list(pii=pii.new, A=A.new, pc=pc.new, f0=f0.new, f1=f1.new, LIS=lfdr, logL=logL, BIC=BIC, ni=niter, converged=converged) 
} else {
 logL<- (-Inf)
 BIC<- logL<- (-Inf)
 em.var<-list(pii=pii.old, A=A.old, pc=pc.old, f0=f0.old, f1=f1.old, LIS=lfdr, logL=logL, BIC=BIC, ni=niter, converged=converged)
}

}

return (em.var)
}

