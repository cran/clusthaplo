\name{cluster.pair}
\alias{cluster.pair}
\title{Computes pairwise similarity and clustering for a given pair of haplotypes for the currently selected chromosome.}
\usage{
    cluster.pair(i, j, simi=pair.similarity(i, j))
}
\note{
    This function can only be run from inside a clusthaplo context as returned by clusthaplo(...).
    Also, it requires that a chromosome is currently selected.
    Also, a threshold must be set or an HMM fit prior to invoking this function.
}
\arguments{
    \item{i}{
        The index or name of the first parent.
    }
    \item{j}{
        The index or name of the second parent.
    }
    \item{simi}{
        If the similarity scores over the scanned loci have already been computed, the vector of similarity scores (mainly for internal use).
    }
}
\description{
    Computes the similarity along the selected chromosome for the given pair of haplotypes and estimate clustering on this signal.
}
\value{
    An object of class \code{clusthaplo.pair.similarity} which behaves as a list containing:
    \itemize{
        \item \code{locus} the vector of scanned loci
        \item \code{similarity} the similarity score for each scanned locus
        \item \code{clustering} the clustering value (TRUE/FALSE) for each scanned locus
    }
    It also has the following attributes:
    \itemize{
        \item \code{chromosome} the name of the selected chromosome
        \item \code{i} the index of the first parent
        \item \code{j} the index of the second parent
        \item \code{scoring.method} the scoring method used (see \link{config})
        \item \code{threshold} the threshold value if clustering was based on threshold, NULL otherwise
        \item \code{simulation.type} the simulation type used to train clusthaplo
    }
    A special plot function is defined for this datatype (see example).
}
\seealso{
    \code{\link{plot.clusthaplo.pair.similarity}}
    \code{\link{config}}
    \code{\link{train}}
    \code{\link{select.chromosome}}
}
\examples{
    data(parents.map, parents.gen)
    clu <- clusthaplo(parents.map, NULL, parents.gen)
    clu$select.chromosome('chr1')
    clu$train()
    clustering12 <- clu$cluster.pair(1, 2)
    plot(clustering12)
}
