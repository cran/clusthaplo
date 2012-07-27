\name{run.raw.simulation}
\alias{run.raw.simulation}
\title{Run a simulation without any post-processing}
\usage{
    run.raw.simulation(n.replicates=3, geno.sim)
}
\note{
    This function can only be run from inside a clusthaplo context as returned by clusthaplo(...).
}
\arguments{
    \item{n.replicates}{
        The number of times to repeat the simulation.
    }
    \item{geno.sim}{
        A function taking arguments (genetic.map, number.of.haplotypes) and returning a simulated data set in the
        form of a matrix with haplotypes in columns and markers in rows.
    }
}
\description{
    Computes similarity signals over \code{n.replicates} simulated haplotype datasets.
    If a chromosome is currently selected, the simulation will only run on it. Otherwise,
    it runs on all chromosomes in the scan map.
}
\value{
    \itemize{
        \item If a chromosome is selected: a list of length \code{n.replicates} containing matrices of similarity signals in columns, with as many columns as \code{\link{combn}(number.of.haplotypes, 2)} returns (that is, a list of matrices).
        \item Otherwise: a list of length \code{n.replicates} containing lists of one similarity signal matrix per chromosome in the scan map (that is, a list of lists of matrices).
    }
}
\seealso{
    \code{\link{config}}
    \code{\link{select.chromosome}}
}
\examples{
    data(parents.map, parents.gen)
    clu <- clusthaplo(parents.map, NULL, parents.gen)
    clu$select.chromosome('chr1')
    my.raw.data <- clu$run.raw.simulation(n.replicates=1, geno.sim=function(map, gen) gen[, 1:6])
    print(length(my.raw.data))
    print(dim(my.raw.data[[1]]))
    print(ncol(combn(6, 2)))
}
