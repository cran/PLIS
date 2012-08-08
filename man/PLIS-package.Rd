\name{PLIS-package}
\Rdversion{1.1}
\alias{PLIS-package}
\alias{PLIS}
\docType{package}
\title{
PLIS
}
\description{
PLIS is a multiple testing procedure for testing several groups of hypotheses. Linear dependency is expected from the hypotheses within the same group and is modeled by hidden Markov Models.
It is noted that, for PLIS, a smaller p value does not necessarily imply more significance because of dependency among the hypotheses. A typical applicaiton of PLIS is to analyze genome wide association studies datasets, where SNPs from the same chromosome are treated as a group and exhibit strong linear genomic dependency. 
}
\details{
\tabular{ll}{
Package: \tab PLIS\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2012-08-08\cr
License: \tab GPL-3\cr
LazyLoad: \tab yes\cr
}
main functions: em.hmm & plis
}
\author{
Wei Z, Sun W, Wang K and Hakonarson H\cr
Maintainer: Zhi Wei <zhiwei04@gmail.com>
}
\references{
 Wei Z, Sun W, Wang K and Hakonarson H, Multiple Testing in Genome-Wide Association Studies via Hidden Markov Models, Bioinformatics, 2009
}
\keyword{ package }
\seealso{p.adjust(), in which the traditional procedures are implemented. The adjustment made by p.adjust will not change the original ranking based on the given p values. However, taking into account dependency, PLIS may generate a ranking different from that by p value.
}
