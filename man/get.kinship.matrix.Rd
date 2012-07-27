\name{get.kinship.matrix}
\alias{get.kinship.matrix}
\title{Access to the kinship matrix}
\usage{
    get.kinship.matrix()
}
\note{
    This function can only be run from inside a clusthaplo context as returned by clusthaplo(...).
}
\description{
    Retrieves the kinship matrix used by a clusthaplo context.
}
\value{
    A matrix giving the pairwise kinship coefficients between individuals.
}
\seealso{
    \code{\link{clusthaplo}}
}
\examples{
    data(parents.map, parents.gen)
    clu <- clusthaplo(parents.map, NULL, parents.gen)
    clu$get.kinship.matrix() # by default, diag(ncol(parents.gen))
    kmat <- matrix(c(1, .5, 0, .5, 1, .125, 0, .125, 1), ncol=3)
    print(kmat)
    clu <- clusthaplo(parents.map, NULL, parents.gen, kinship=kmat)
    clu$get.kinship.matrix()
}
