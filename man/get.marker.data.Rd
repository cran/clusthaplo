\name{get.marker.data}
\alias{get.marker.data}
\title{Access to the marker data}
\usage{
    get.marker.data()
}
\note{
    This function can only be run from inside a clusthaplo context as returned by clusthaplo(...).
}
\description{
    Retrieves the marker data from a clusthaplo context as returned by clusthaplo(...).
}
\value{
    The marker data matrix with markers as rows and haplotypes as columns.
}
\seealso{
    \code{\link{clusthaplo}}
}
\examples{
    data(parents.map, parents.gen)
    clu <- clusthaplo(parents.map, NULL, parents.gen)
    is.matrix(clu$get.marker.data())
    dim(clu$get.marker.data())
    colnames(clu$get.marker.data())
    rownames(clu$get.marker.data())
}
