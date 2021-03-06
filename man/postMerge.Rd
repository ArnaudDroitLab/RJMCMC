% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rjmcmcMethodsIntern.R
\name{postMerge}
\alias{postMerge}
\title{A internal post treatment function to merge closely positioned
nucleosomes, from the same chromosome,
identified by the \code{\link{rjmcmc}} function}
\usage{
postMerge(startPosForwardReads, startPosReverseReads, resultRJMCMC,
  extendingSize, chrLength, minReads = 5)
}
\arguments{
\item{startPosReverseReads}{a \code{vector} of \code{numeric}, the
start position of all the reverse reads. Beware that the start position of
a reverse read is always higher that the end positition.}

\item{resultRJMCMC}{an object of class 'rjmcmcNucleosomes' or
'rjmcmcNucleosomesMerge' containing informations about nucleosomes.}

\item{extendingSize}{a positive \code{numeric} or a positive \code{integer}
indicating the size of the consensus region used to group closeley
positioned nucleosomes.The minimum size of the consensus region is equal to
twice the value of the \code{extendingSize} parameter. The numeric will
be treated as an integer.}

\item{chrLength}{a positive \code{numeric} or a positive \code{integer}
indicating the lenght of the current chromosome. The length of the
chromosome is used to ensure that the consensus positions are all
located inside the chromosome.}

\item{minReads}{a positive \code{integer} or \code{numeric}, the minimum
number of reads in a potential canditate region. Non-integer values
of \code{minReads} will be casted to \code{integer} and truncated towards
zero. Default: 5.}

\item{startPosFrowardReads}{a \code{vector} of \code{numeric}, the
start position of all the forward reads.}
}
\value{
a \code{array} of \code{numeric}, the updated values of the
nucleosome positions.
}
\description{
A internal helper function which merges closely positioned
nucleosomes to rectify the over splitting and provide a more conservative
approach. Beware that each chromosome must be treated separatly.

The function uses the Bioconductor \code{package} \code{consensusSeeker} to
group closely positioned nucleosomes.
}
\examples{

## Loading dataset
data(RJMCMC_result)
data(reads_demo)

## Results before post-treatment
RJMCMC_result$mu

## Post-treatment function which merged closely positioned nucleosomes
postResult <- RJMCMC:::postMerge(startPosForwardReads =
reads_demo$readsForward, startPosReverseReads = reads_demo$readsReverse,
resultRJMCMC = RJMCMC_result, extendingSize = 80, chrLength = 73500)

## Results after post-treatment
postResult

}
\author{
Pascal Belleau, Astrid Deschenes
}
\keyword{internal}
