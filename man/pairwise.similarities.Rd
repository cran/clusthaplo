\name{pairwise.similarities}
\alias{pairwise.similarities}
\title{Computes pairwise similarity on all pairs of haplotypes for the currently selected chromosome.}
\usage{
    pairwise.similarities(hook.pair.simi = function(i, j, ps) ps)
}
\note{
    This function can only be run from inside a clusthaplo context as returned by clusthaplo(...).
}
\arguments{
    \item{hook.pair.simi}{
        Optional callback that will be called after each pairwise similarity computation. Arguments are:
        \itemize{
            \item \code{i} index of the first parent
            \item \code{j} index of the second parent
            \item \code{ps} the similarity signal as returned by \code{\link{pair.similarity}}
        }
    }
}
\description{
    Computes pairwise similarity on all pairs of haplotypes for the currently selected chromosome.
}
\value{
    A matrix of logicals denoting similarity between all pairs of haplotypes (in columns) for each scanned locus.
}
\seealso{
    \code{\link{config}}
    \code{\link{train}}
    \code{\link{select.chromosome}}
}
\examples{
    data(parents.map, parents.gen)
    clu <- clusthaplo(parents.map, NULL, parents.gen)
    clu$select.chromosome('chr1')
    clu$train()
    pws <- clu$pairwise.similarities()
    dim(pws)
}
