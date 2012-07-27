\name{plot.clusthaplo.alleles}
\alias{plot.clusthaplo.alleles}
\title{Graphic representation of the count of ancestral alleles along a chromosome}
\usage{
\method{plot}{clusthaplo.alleles}(x, ...)
}
\arguments{
    \item{x}{
        A count of ancestral alleles per locus as returned by
        \code{\link{count.ancestral.alleles}}.
    }
    \item{...}{
        Optional arguments include:
        \itemize{
            \item{tick.interval}{
                The interval in centi-Morgans for the horizontal axis.
            }
        }
        The rest of the arguments are passed to \code{\link{plot}}.
    }
}
\description{
    Plots a histogram-like view of the count of ancestral alleles for each scanned locus.
}
\seealso{
    \code{\link{count.ancestral.alleles}}
}
\examples{
    data(parents.map, parents.gen)
    clu <- clusthaplo(parents.map, NULL, parents.gen)
    clu$select.chromosome('chr1')
    clu$train()
    tc <- clu$transitive.closure(clu$pairwise.similarities())
    plot(clu$count.ancestral.alleles(tc))
}
