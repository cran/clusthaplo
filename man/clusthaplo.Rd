\name{clusthaplo}
\alias{clusthaplo}
\title{clusthaplo}
\description{
    Creates a clusthaplo context, from which you can invoke the various
    computations this package provides.
}
\usage{
clusthaplo(haplotypes.map, mcqtl.consensus.map, marker.data,
           kinship.matrix=NULL, na.strings=c('.', 'NA', '-'),
           discard.unknown.markers=T,
           ...)
}
\arguments{
    \item{haplotypes.map}{
        filename or named list of matrices or data.frames (one for each chromosome) for the genetic maps of the haplotypes
    }
    \item{mcqtl.consensus.map}{
        filename or named list of matrices or data.frames (one for each chromosome) for the MCQTL consensus maps
    }
    \item{marker.data}{
        filename or matrix or data.frame for the genotype data (one column per parent)
    }
    \item{kinship.matrix}{
        filename or matrix or data.frame giving the pairwise kinship coefficients between all the haplotypes
    }
    \item{na.strings}{
        defines which NA strings we may encounter while reading the marker data
    }
    \item{discard.unknown.markers}{
        if \code{TRUE}, all the markers that don't have genotype data are removed from the map that is used to compute the similarities.
        Otherwise, those markers are added to the genotype data with \code{NA} observations for all individuals.
    }
    \item{...}{
        optional arguments to pass to config() before returning the clusthaplo context.
    }
}
\value{
    Returns a list containing all the relevant methods provided by this package.
}
\examples{
    data(parents.map, parents.gen)
    is.list(parents.map)
    clu <- clusthaplo(parents.map, NULL, parents.gen)
    names(clu)
    clu$get.chromosome.names()
}
\keyword{haplotype}
\keyword{clustering}
