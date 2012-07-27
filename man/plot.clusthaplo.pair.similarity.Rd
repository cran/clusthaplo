\name{plot.clusthaplo.pair.similarity}
\alias{plot.clusthaplo.pair.similarity}
\title{Graphic representation of a single similarity signal along a chromosome}
\usage{
\method{plot}{clusthaplo.pair.similarity}(x, ...)
}
\arguments{
    \item{x}{
        A pairwise similarity signal with clustering information
        as returned by \code{\link{cluster.pair}}.
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
    Plots the similarity signal, the threshold value that was used if defined, and the
    clustered segments.
}
\seealso{
    \code{\link{cluster.pair}}
}
\examples{
    data(parents.map, parents.gen)
    clu <- clusthaplo(parents.map, NULL, parents.gen)
    clu$select.chromosome('chr1')
    clu$train()
    clustering12 <- clu$cluster.pair(1, 2)
    plot(clustering12)
}
