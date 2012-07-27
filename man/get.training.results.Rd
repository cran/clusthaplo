\name{get.training.results}
\alias{get.training.results}
\title{Access to the last training results}
\usage{
    get.training.results()
}
\note{
    This function can only be run from inside a clusthaplo context as returned by clusthaplo(...).
}
\description{
    Retrieves the results of the last call to train().
}
\value{
    If \code{clustering.method} is set to \code{"hmm"}, the result is a list
    of the 2-state, 3-state and 4-state HMMs that were fitted to the simulation data.
    If \code{clustering.method} is set to \code{"threshold"}, the result is
    the vector of the selected quantiles.
}
\seealso{
    \code{\link{train}}
    \code{\link{config}}
}
\examples{
    data(parents.map, parents.gen)
    clu <- clusthaplo(parents.map, NULL, parents.gen, simulation.Nrep=1, simulation.Np=5)
    clu$select.chromosome('chr1')
    clu$train()
    clu$get.training.results()
}
