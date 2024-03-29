\name{em.hmm}
\Rdversion{1.1}
\alias{em.hmm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
   EM algorithm for HMM to estimate LIS statistic
}
\description{
  em.hmm calculates the MLE for a HMM model with hidden states being 0/1.
  the distribution of observed Z values given state 0 is assumed to be normal and 
  gvien state 1, is assumed to be a normal mixture with L components
}
\usage{
em.hmm(x, L=2, maxiter = 1000, est.null = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
  the observed Z values
}
  \item{L}{
 the number of components in the non-null mixture, default value=2
}
  \item{maxiter}{
 the maximum number of iterations, default value=1000
}
  \item{est.null}{
  logical. If FALSE (the default) set the null distribution as N(0,1), otherwise will estimate the null distribution.
}
}
\details{
  None.
}
\value{
  \item{pii}{the initial state distribution, pii=(prob. of being 0, prob. of being 1)}
  \item{A}{transition matrix, A=(a00 a01\\ a10 a11)}
  \item{f0}{the null distribution}
  \item{pc}{probability weights of each component in the non-null mixture} 
  \item{f1}{an L by 2 matrix, specifying the dist. of each component in the non-null mixture} 
  \item{LIS}{the LIS statistics} 
  \item{ni}{the number of iterations excecuted} 
  \item{logL}{log likelihood} 
  \item{BIC}{BIC score for the estimated model} 
  \item{converged}{Logic, Convergence indicator of the EM procedure}
}
\references{
 Multiple Testing in Genome-Wide Association Studies via Hidden Markov Models, Bioinformatics, 2009
}
\author{
 Wei Z, Sun W, Wang K and Hakonarson H
}
\seealso{
plis
}
\examples{
##(1) Example for analyzing simulated data
grp1.nonNull.loci=c(21:30, 51:60); grp2.nonNull.loci=c(41:60)
grp1.theta<-grp2.theta<-rep(0,200)
grp1.theta[grp1.nonNull.loci]=2; grp2.theta[grp2.nonNull.loci]=2

grp1.zval=rnorm(n=length(grp1.theta),mean=grp1.theta)
grp2.zval=rnorm(n=length(grp2.theta),mean=grp2.theta)
##Group 1
#Use default L=2
grp1.L2rlts=em.hmm(grp1.zval)
#Use true value L=1
grp1.L1rlts=em.hmm(grp1.zval,L=1)
#Choose L by BIC criteria
grp1.Allrlts=sapply(1:3, function(k) em.hmm(grp1.zval,L=k))
BICs=c()
for(i in 1:3) {
  BICs=c(BICs,grp1.Allrlts[[i]]$BIC)
}
grp1.BICrlts=grp1.Allrlts[[which(BICs==max(BICs))]]

rank(grp1.BICrlts$LIS)[grp1.nonNull.loci]
rank(-abs(grp1.zval))[grp1.nonNull.loci]

##Group 2
grp2.Allrlts=sapply(1:3, function(k) em.hmm(grp2.zval,L=k))
BICs=c()
for(i in 1:3) {
  BICs=c(BICs,grp2.Allrlts[[i]]$BIC)
}
grp2.BICrlts=grp2.Allrlts[[which(BICs==max(BICs))]]

rank(grp2.BICrlts$LIS)[grp2.nonNull.loci]
rank(-abs(grp2.zval))[grp2.nonNull.loci]

##PLIS: control global FDR
states=plis(c(grp1.BICrlts$LIS,grp2.BICrlts$LIS),fdr=0.1,adjust=FALSE)$States
#0 accept; 1 reject under fdr level 0.1

##(2) Example for analyzing Genome-Wide Association Studies (GWAS) data
#Information in GWAS.SampleData can be obtained by using PLINK
#http://pngu.mgh.harvard.edu/~purcell/plink/

#not running
#please uncomment to run
#
#data(GWAS.SampleData)
#
#chr1.data=GWAS.SampleData[which(GWAS.SampleData[,"CHR"]==1),]
#chr6.data=GWAS.SampleData[which(GWAS.SampleData[,"CHR"]==6),]
#
##Make sure SNPs in the linear physical order
#chr1.data<-chr1.data[order(chr1.data[,"BP"]),]
#chr6.data<-chr6.data[order(chr6.data[,"BP"]),]
#
##convert p values by chi_sq test to z values; odds ratio (OR) is needed.
#chr1.zval<-rep(0, nrow(chr1.data))
#chr1.ors=(chr1.data[,"OR"]>1)
#chr1.zval[chr1.ors]<-qnorm(chr1.data[chr1.ors, "P"]/2, 0, 1, lower.tail=FALSE)
#chr1.zval[!chr1.ors]<-qnorm(chr1.data[!chr1.ors, "P"]/2, 0, 1, lower.tail=TRUE)
#chr1.L2rlts=em.hmm(chr1.zval)
#
#chr6.zval<-rep(0, nrow(chr6.data))
#chr6.ors=(chr6.data[,"OR"]>1)
#chr6.zval[chr6.ors]<-qnorm(chr6.data[chr6.ors, "P"]/2, 0, 1, lower.tail=FALSE)
#chr6.zval[!chr6.ors]<-qnorm(chr6.data[!chr6.ors, "P"]/2, 0, 1, lower.tail=TRUE)
#chr6.L2rlts=em.hmm(chr6.zval)
#
##Note that for analyzing a chromosome in real GWAS dataset, em.hmm can take as long as 10+ hrs
##L=2 or 3 is recommended for GWAS based on our experience
##em.hmm can be run in parallel for different chromosomes before applying the PLIS procedure
#plis.rlts=plis(c(chr1.L2rlts$LIS,chr6.L2rlts$LIS),fdr=0.01)
#all.Rlts=cbind(rbind(chr1.data,chr6.data), LIS=c(chr1.L2rlts$LIS,chr6.L2rlts$LIS), 
#gFDR=plis.rlts$aLIS, fdr001state=plis.rlts$States)
#all.Rlts[order(all.Rlts[,"LIS"])[1:10],]

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
