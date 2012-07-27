\name{get.mcqtl.consensus.map}
\alias{get.mcqtl.consensus.map}
\title{Access to the MCQTL consensus marker map}
\usage{
    get.mcqtl.consensus.map(chrom.name=NULL)
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
    Retrieves the MCQTL consensus marker map from the given clusthaplo context as returned by clusthaplo(...).
}
\value{
    If \code{chrom.name} is \code{NULL}, returns the map as a list for all the chromosomes.
    Otherwise, returns a data.frame containing the marker map data for the given chromosome.
}
\seealso{
    \code{\link{get.haplotypes.map}}
    \code{\link{get.scan.map}}
    \code{\link{clusthaplo}}
}
\examples{
    data(parents.map, parents.gen)
    clu <- clusthaplo(parents.map, NULL, parents.gen)
    # in this case, the MCQTL consensus map is also the haplotypes map.
    length(clu$get.mcqtl.consensus.map())
    names(clu$get.mcqtl.consensus.map())
}
