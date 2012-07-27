\name{write.clusters}
\alias{write.clusters}
\title{Output a clustering to file}
\usage{
    write.clusters(tc, file)
}
\arguments{
    \item{tc}{
        The result of a transitive closure of pair similarities as returned by \code{\link{transitive.closure}}.
    }
    \item{file}{
        A filename.
    }
}
\description{
    Writes the clustering data to a file.
}
\seealso{
    \code{\link{transitive.closure}}
    \code{\link{read.clusters}}
}
\examples{
    data(parents.map, parents.gen)
    clu <- clusthaplo(parents.map, NULL, parents.gen)
    clu$select.chromosome('chr1')
    clu$train()
    tc <- clu$transitive.closure(clu$pairwise.similarities())
    \dontrun{
    write.clusters(tc, "myclustering.txt")
    }
}
