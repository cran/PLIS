\name{bwfw1.hmm}
\Rdversion{1.1}
\alias{bwfw1.hmm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
backward and forward inferences
}
\description{
When L=1, calculate values for backward, forward variables, probabilities of hidden states. A supporting function called by em.hmm.
}
\usage{
bwfw1.hmm(x, pii, A, f0, f1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
the observed Z values
}
  \item{pii}{
(prob. of being 0, prob. of being 1), the initial state distribution 
}
  \item{A}{
A=(a00 a01\\\\a10 a11), transition matrix
}
  \item{f0}{
(mu, sigma), the parameters for null distribution
}
  \item{f1}{(mu[1], sigma[1]\\\\...\\\\mu[L], sigma[L])--an L by 2 matrix, the parameter set for the non-null distribution
}
}
\details{
 calculates values for backward, forward variables, probabilities of hidden states,\cr
  --the lfdr variables and etc. \cr
  --using the forward-backward procedure (Rabiner 89) \cr
  --based on a sequence of observations for a given hidden markov model M=(pii, A, f) \cr
  --see Sun and Cai (2009) for a detailed instruction on the coding of this algorithm
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
 \item{alpha}{rescaled backward variables}
 \item{beta}{rescaled forward variables}
 \item{lfdr}{lfdr variables}
 \item{gamma}{probabilities of hidden states}
 \item{dgamma}{rescaled transition variables}
 \item{omega}{rescaled weight variables}
}
\references{
Multiple Testing in Genome-Wide Association Studies via Hidden Markov Models, Bioinformatics, 2009\cr
Large-scale multiple testing under dependence, Sun W and Cai T (2009), JRSSB, 71, 393-424\cr
A Tutorial on Hidden Markov Models and Selected Applications in Speech Recognition, Rabiner L (1989), Procedings of the IEEE, 77, 257-286.
}
\author{
 Wei Z, Sun W, Wang K and Hakonarson H
}
\keyword{ models }
