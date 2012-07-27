\name{get.config}
\alias{get.config}
\title{Query the configuration of a clusthaplo context.}
\usage{
    get.config(...)
}
\note{
    This function can only be run from inside a clusthaplo context as returned by clusthaplo(...).
}
\arguments{
    \item{...}{
        Any number of character strings. Each string must correspond to a formal argument of \code{config}.
    }
}
\description{
    Queries the current configuration.
}
\value{
    If no argument is provided, returns a named list containing all the current configuration values.
    If exactly one argument is provided, returns the corresponding configuration value.
    If at least two arguments are provided, returns a named list of the corresponding configuration values.
}
\seealso{
    \code{\link{clusthaplo}}
    \code{\link{config}}
}
\examples{
    data(parents.map, parents.gen)
    clu <- clusthaplo(parents.map, NULL, parents.gen)
    print(clu$get.config('na.replace'))
    print(clu$get.config('w1', 'w2'))
    print(clu$get.config())
}
