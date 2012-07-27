\name{get.chromosome.names}
\alias{get.chromosome.names}
\title{Access to the list of chromosome names}
\usage{
    get.chromosome.names()
}
\note{
    This function can only be run from inside a clusthaplo context as returned by clusthaplo(...).
}
\description{
    Retrieves the chromosome names from a clusthaplo context as returned by clusthaplo(...).
}
\value{
    A character vector containing all the chromosome names defined in the
    given clusthaplo context as returned by clusthaplo(...).
    For use in \code{\link{select.chromosome}}.
}
\seealso{
    \code{\link{select.chromosome}}
    \code{\link{get.selected.chromosome}}
    \code{\link{clusthaplo}}
}
\examples{
    data(parents.map, parents.gen)
    clu <- clusthaplo(parents.map, NULL, parents.gen)
    clu$get.chromosome.names()
}
