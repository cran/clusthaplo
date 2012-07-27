\name{plot.clusthaplo.haplotypes}
\alias{plot.clusthaplo.haplotypes}
\title{Graphic representation of the clusterings along a chromosome}
\usage{
\method{plot}{clusthaplo.haplotypes}(x, ...)
}
\arguments{
    \item{x}{
        The result of the transitive closure of a complete set of pairwise clustering as returned by \code{transitive.closure}.
    }
    \item{...}{
        Optional arguments include:
        \itemize{
            \item{ord}{
                The order or subset of the haplotypes to use in drawing the chromosome.
            }
            \item{tick.interval}{
                The interval in centi-Morgans for the horizontal axis.
            }
            \item{expand.cliques}{
                Whether to try and merge colors as much as possible. Makes a much clearer
                representation without modifying the information displayed.
            }
        }
        The rest of the arguments are passed to \code{\link{plot}}.
    }
}
\description{
    Draws a representation of the chromosome.
}
\seealso{
    \code{\link{transitive.closure}}
}
\examples{
    data(parents.map, parents.gen)
    clu <- clusthaplo(parents.map, NULL, parents.gen)
    clu$select.chromosome('chr1')
    clu$train()
    tc <- clu$transitive.closure(clu$pairwise.similarities())
    plot(tc)
}
