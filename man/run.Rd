\name{run}
\alias{run}
\title{Runs the whole clusthaplo pipeline}
\usage{
run(chromosomes=get.chromosome.names(), train.on.all.chromosomes=T,
    output.XML=T, plot.device=png,
    plot.alleles=T, plot.chromosomes=T, plot.pair.similarities=T, plot.densities=T,
    ...)
}
\note{
    This function can only be run from inside a clusthaplo context as returned by clusthaplo(...).
}
\arguments{
    \item{chromosomes}{
        Character vector of chromosome names to run on.
    }
    \item{train.on.all.chromosomes}{
        Logical. If \code{TRUE}, training is done once on all chromosomes, otherwise, clusthaplo will re-train
        on each chromosome.
    }
    \item{output.XML}{
        Logical. \code{TRUE} makes clusthaplo output the MCQTL-LD XML files.
    }
    \item{plot.device}{
        This function will be invoked during the process to open new graphic devices to output the plots.
    }
    \item{plot.alleles}{
        Logical. Controls plotting the count of ancestral alleles for each locus.
    }
    \item{plot.chromosomes}{
        Logical. Controls plotting the haplotypes.
    }
    \item{plot.pair.similarities}{
        Logical. Controls plotting ALL INDIVIDUAL pairwise similarity signals.
    }
    \item{plot.densities}{
        Logical. Controls plotting the window density for each locus.
    }
    \item{...}{
        Supplemental arguments passed to \code{\link{config}}.
    }
}
\description{
    Configures, trains, and runs clusthaplo over all or the given subset of chromosomes.
}
\seealso{
    \code{\link{config}}
    \code{\link{train}}
    \code{\link{select.chromosome}}
    \code{\link{pairwise.similarities}}
    \code{\link{transitive.closure}}
    \code{\link{plot.clusthaplo.haplotypes}}
    \code{\link{plot.clusthaplo.alleles}}
    \code{\link{plot.clusthaplo.density}}
    \code{\link{write.xml}}
}
\examples{
    \dontrun{
        data(parents.map, parents.gen)
        clu <- clusthaplo(parents.map, NULL, parents.gen)
        clu$run('chr1', output.XML=FALSE)
    }
}
