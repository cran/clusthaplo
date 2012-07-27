\name{get.window.density}
\alias{get.window.density}
\title{Access to the window density}
\usage{
    get.window.density()
}
\note{
    This function can only be run from inside a clusthaplo context as returned by clusthaplo(...).
    Also, it requires that a chromosome is currently selected.
}
\description{
    Retrieves the window density (number of markers in the window for each
    locus) for the currently selected chromosome.
}
\value{
    A list of class \code{clusthaplo.density} containing the vector
    \code{locus} of the scanned loci and the vector \code{density} of the
    window densities.
    It also has the following attributes:
    \itemize{
        \item \code{window.length} the length of the window in centi-Morgans
        \item \code{chromosome} the name of the selected chromosome
    }
    A special plot function is defined for this datatype (see example).
}
\seealso{
    \code{\link{clusthaplo}}
    \code{\link{select.chromosome}}
    \code{\link{config}}
    \code{\link{plot.clusthaplo.density}}
}
\examples{
    data(parents.map, parents.gen)
    clu <- clusthaplo(parents.map, NULL, parents.gen)
    clu$select.chromosome(clu$get.chromosome.names()[1])
    plot(clu$get.window.density())
}
