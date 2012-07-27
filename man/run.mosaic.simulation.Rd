\name{run.mosaic.simulation}
\alias{run.mosaic.simulation}
\title{Run a simulation using the 'Mosaic' haplotype simulator}
\usage{
    run.mosaic.simulation(n.haplotypes=10,
                          n.generations=50,
                          n.replicates=3,
                          quantiles=seq(.8, 1, by=.01))
}
\note{
    This function can only be run from inside a clusthaplo context as returned by clusthaplo(...).
}
\arguments{
    \item{n.haplotypes}{
        The number of haplotypes to simulate.
    }
    \item{n.generations}{
        The number of intermediary generations to simulate.
    }
    \item{n.replicates}{
        The number of times to repeat the simulation.
    }
    \item{quantiles}{
        A vector of quantiles to give to the function \code{\link{quantile}} in the post-processing stage.
    }
}
\description{
    Computes similarity signals over \code{n.replicates} simulated haplotype datasets.
    If a chromosome is currently selected, the simulation will only run on it. Otherwise,
    it runs on all chromosomes in the scan map.
    Once done, the relevant quantiles are returned.
}
\value{
    A vector of thresholds corresponding to each quantile specified in \code{quantiles}.
}
\seealso{
    \code{\link{config}}
    \code{\link{select.chromosome}}
}
\examples{
    data(parents.map, parents.gen)
    clu <- clusthaplo(parents.map, NULL, parents.gen)
    clu$select.chromosome('chr1')
    my.quantiles <- clu$run.mosaic.simulation(quantiles=seq(.7, 1, by=.01))
    print(my.quantiles)
}
