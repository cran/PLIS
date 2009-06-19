\name{plis}
\Rdversion{1.1}
\alias{plis}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A multiple testing procedure based on pooled LIS statistics
}
\description{
It controls the global FDR for the pooled hypotheses from different groups
}
\usage{
plis(lis, fdr = 0.001, adjust = TRUE)
}
\arguments{
  \item{lis}{
   pooled LIS statistics estimated from different groups
}
  \item{fdr}{
   nominal fdr level you want to control 
 }
  \item{adjust}{
   logical. If TRUE (the default), will calculate and return "adjusted" LIS value-- the corresponding global FDR if using the LIS statistic as the significance cutoff. It may take hours if you have hundreds of thousands LISs to adjust.
}
}
\value{
 \item{States}{state sequence indicating if the hypotheses should be rejected or not: 0 accepted , 1 rejected}
 \item{aLIS}{the corresponding global FDR if using the LIS statistic as the significance cutoff }	
}
\references{
Multiple Testing in Genome-Wide Association Studies via Hidden Markov Models, Bioinformatics, 2009
}
\author{
Wei Z, Sun W, Wang K and Hakonarson H
}

\seealso{
see em.hmm for examples
}
\keyword{ htest }
