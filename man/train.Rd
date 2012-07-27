\name{train}
\alias{train}
\title{Perform automated training of clusthaplo}
\usage{
    train(quantiles=seq(.8, 1, by=.01), ...)
}
\note{
    This function can only be run from inside a clusthaplo context as returned by clusthaplo(...).
}
\arguments{
    \item{quantiles}{
        A numeric vector of quantiles to pass to the postprocessing stage.
    }
    \item{...}{
        Supplemental arguments passed to HMMFit from the package RHmm.
    }
}
\description{
    According to the selected \code{clustering.method}, either fits an HMM or estimates a threshold value
    based on the selected \code{threshold.quantile} after computing similarity signals on simulated genotype
    data over \code{simulation.Np} haplotypes, repeated \code{simulation.Nrep} times.

    If a chromosome is currently selected, simulation will run only on this chromosome. Otherwise, simulation will
    run on ALL chromosomes in the loaded dataset.
}
\seealso{
    \code{\link{config}}
    \code{\link{select.chromosome}}
    \code{\link[RHmm]{HMMFit}}
}
\examples{
    data(parents.map, parents.gen)
    clu <- clusthaplo(parents.map, NULL, parents.gen)
    clu$select.chromosome('chr1')
    clu$get.config('clustering.threshold')
    clu$train()
    clu$get.config('clustering.threshold')
}
