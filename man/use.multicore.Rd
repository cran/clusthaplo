\name{use.multicore}
\alias{use.multicore}
\title{Enable parallelization}
\usage{
    use.multicore()
}
\note{
    This function can only be run from inside a clusthaplo context as returned by clusthaplo(...).
}
\description{
    If the package \code{multicore} is loaded, applies parallelization on ALL loops in clusthaplo.
}
\examples{
    data(parents.map, parents.gen)
    clu <- clusthaplo(parents.map, NULL, parents.gen)
    \dontrun{
    clu$use.multicore()
    clu$train()
    }
}
