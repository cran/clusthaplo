\name{plot.clusthaplo.density}
\alias{plot.clusthaplo.density}
\title{Graphic representation of the window density along a chromosome}
\usage{
\method{plot}{clusthaplo.density}(x, ...)
}
\arguments{
    \item{x}{
        The window density as returned by
        \code{\link{get.window.density}}.
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
    Plots a histogram-like view of the window density along a chromosome.
}
\seealso{
    \code{\link{get.window.density}}
}
\examples{
    data(parents.map, parents.gen)
    clu <- clusthaplo(parents.map, NULL, parents.gen)
    clu$select.chromosome('chr1')
    plot(clu$get.window.density())
}
