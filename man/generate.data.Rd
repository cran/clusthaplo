\name{generate.data}
\alias{generate.data}
\title{Simulated dataset generation}
\usage{
generate.data(Nc, Nancestors, Nhaplotypes,
              Nm=80+sample.int(40, Nc, replace=T),
              Nmissing=0, mosaic.Ng=50,
              filename.prefix="generated")
}
\arguments{
    \item{Nc}{
        The number of chromosomes to simulate.
    }
    \item{Nancestors}{
        The number of ancestors to start from.
    }
    \item{Nhaplotypes}{
        The number of haplotypes to simulate.
    }
    \item{Nm}{
        The number of markers per chromosome. Must be an int vector of
        length \code{Nc}.
    }
    \item{Nmissing}{
        The number of missing data in the simulated genotypes matrix.
    }
    \item{mosaic.Ng}{
        The number of intermediary generations to simulate (the \code{Ng}
        parameter of the \code{mosaic} simulation scheme).
    }
    \item{filename.prefix}{
        The path and prefix to use in naming the generated files, which will be:
        \itemize{
            \item \code{prefix}.map
            \item \code{prefix}.geno.txt
        }
    }
}
\description{
    Simulates a genetic map and a set of genotypes.
}
\seealso{
    \code{\link{clusthaplo}}
}
