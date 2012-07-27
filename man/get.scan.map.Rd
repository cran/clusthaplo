\name{get.scan.map}
\alias{get.scan.map}
\title{Access to the scan map used to compute the pairwise similarities}
\usage{
    get.scan.map(chrom.name=NULL)
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
    Retrieves the scan map used for the similarity scoring from the given clusthaplo context as returned by clusthaplo(...).
}
\value{
    If \code{chrom.name} is \code{NULL}, returns the map as a list for all the chromosomes.
    Otherwise, returns a data.frame containing the scan map data for the given chromosome.
}
\seealso{
    \code{\link{get.mcqtl.consensus.map}}
    \code{\link{get.haplotypes.map}}
    \code{\link{clusthaplo}}
}
\examples{
    data(parents.map, parents.gen)
    clu <- clusthaplo(parents.map, NULL, parents.gen)
    length(clu$get.scan.map())
    names(clu$get.scan.map())
    clu$get.scan.map('chr1')
}
