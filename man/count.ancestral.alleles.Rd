\name{count.ancestral.alleles}
\alias{count.ancestral.alleles}
\title{Count ancestral alleles along a chromosome}
\usage{
count.ancestral.alleles(clusters)
}
\note{
    This function can only be run from inside a clusthaplo context as returned by clusthaplo(...).
}
\arguments{
    \item{clusters}{
        pairwise similarity matrix as returned by \code{\link{pairwise.similarities}}
    }
}
\description{
    Counts the total number of ancestral alleles per scanned locus.
}
\value{
    A data.frame of class \code{clusthaplo.alleles} with columns
    \code{locus} and \code{count}.
    It also has the following attributes:
    \itemize{
        \item \code{chromosome} the name of the selected chromosome
    }
    A special plot function is defined for this datatype (see example).
}
\seealso{
    \code{\link{transitive.closure}}
    \code{\link{select.chromosome}}
    \code{\link{plot.clusthaplo.alleles}}
    \code{\link{config}}
}
\examples{
    data(parents.map, parents.gen)
    clu <- clusthaplo(parents.map, NULL, parents.gen)
    clu$select.chromosome('chr1')
    clu$train()
    tc <- clu$transitive.closure(clu$pairwise.similarities())
    clu$count.ancestral.alleles(tc)
}
