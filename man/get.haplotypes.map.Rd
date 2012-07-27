\name{get.haplotypes.map}
\alias{get.haplotypes.map}
\title{Access to the marker map of the haplotypes}
\usage{
    get.haplotypes.map(chrom.name=NULL)
}
\note{
    This function can only be run from inside a clusthaplo context as returned by clusthaplo(...).
}
\arguments{
    \item{chrom.name}{
        A valid chromosome name, or \code{NULL}.
    }
}
\description{
    Retrieves the marker map for the haplotypes from the given clusthaplo context.
}
\value{
    If \code{chrom.name} is \code{NULL}, returns the map as a list for all the chromosomes.
    Otherwise, returns a data.frame containing the marker map for the given chromosome.
}
\seealso{
    \code{\link{get.mcqtl.consensus.map}}
    \code{\link{get.scan.map}}
    \code{\link{clusthaplo}}
}
\examples{
    data(parents.map, parents.gen)
    clu <- clusthaplo(parents.map, NULL, parents.gen)
    length(clu$get.haplotypes.map())
    names(clu$get.haplotypes.map())
}
