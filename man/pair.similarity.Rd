\name{pair.similarity}
\alias{pair.similarity}
\title{Computation of a pairwise similarity signal}
\usage{
pair.similarity(i, j)
}
\note{
    This function can only be run from inside a clusthaplo context as returned by clusthaplo(...).
    Also, it requires that a chromosome is currently selected.
}
\arguments{
    \item{i}{
        The index or name of the first parent.
    }
    \item{j}{
        The index or name of the second parent.
    }
}
\description{
    Computes the similarity between the two given haplotypes along the currently
    selected chromosome.
}
\value{
    A list containing the vectors:
    \itemize{
        \item \code{loci} of the scanned loci,
        \item \code{similarity} of the similarity score for each locus,
        \item \code{density} of the window density for each locus.
    }
}
\seealso{
    \code{\link{cluster.pair}}
    \code{\link{select.chromosome}}
    \code{\link{config}}
}
\examples{
    data(parents.map, parents.gen)
    clu <- clusthaplo(parents.map, NULL, parents.gen)
    clu$select.chromosome('chr1')
    as.data.frame(clu$pair.similarity(1, 2))
}
