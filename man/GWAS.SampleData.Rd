\name{GWAS.SampleData}
\Rdversion{1.1}
\alias{GWAS.SampleData}
\docType{data}
\title{
Sample GWAS Dataset
}
\description{
Sample GWAS Dataset with 400 SNPs from Chromosome 1 and 6 (200 SNPs each).
}
\usage{data(GWAS.SampleData)}
\format{
  A data frame with 400 observations on the following 6 variables.
  \describe{
    \item{\code{CHR}}{Chromosome ID}
    \item{\code{SNP}}{rs Id}
    \item{\code{BP}}{Phisical Position}
    \item{\code{OR}}{Odds Ratio}
    \item{\code{CHISQ}}{1 d.f. Chi Square test Statistic}
    \item{\code{P}}{P value of 1 d.f. Chi Square test Statistic}
  }
}
\details{
The required values (Odds ratio and P value) can be calculated by using PLINK
}
\references{
Supplementary Material of Multiple Testing in Genome-Wide Association Studies via Hidden Markov Models, Bioinformatics, 2009
}
\examples{
data(GWAS.SampleData)
}
\keyword{datasets}
