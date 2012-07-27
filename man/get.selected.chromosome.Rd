\name{get.selected.chromosome}
\alias{get.selected.chromosome}
\title{Access to the name of the currently selected chromosome}
\usage{
    get.selected.chromosome()
}
\note{
    This function can only be run from inside a clusthaplo context as returned by clusthaplo(...).
}
\description{
    Retrieves the selected chromosome from a clusthaplo context as returned by clusthaplo(...).
}
\value{
    The name of the currently selected chromosome, or \code{NA} if there is none.
}
\seealso{
    \code{\link{clusthaplo}}
    \code{\link{select.chromosome}}
    \code{\link{get.chromosome.names}}
}
\examples{
    data(parents.map, parents.gen)
    clu <- clusthaplo(parents.map, NULL, parents.gen)
    clu$get.selected.chromosome()
    clu$select.chromosome(clu$get.chromosome.names()[1])
    clu$get.selected.chromosome()
    clu$select.chromosome(NA)
    clu$get.selected.chromosome()
}
