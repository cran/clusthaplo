\name{read.clusters}
\alias{read.clusters}
\title{Read clustering data from a file}
\usage{
    read.clusters(file)
}
\arguments{
    \item{file}{
        A filename.
    }
}
\description{
    Reads clustering data from a file.
}
\value{
    The result of a transitive closure of pair similarities as returned by \code{\link{transitive.closure}}.
}
\seealso{
    \code{\link{transitive.closure}}
    \code{\link{write.clusters}}
}
\examples{
    data(parents.map, parents.gen)
    clu <- clusthaplo(parents.map, NULL, parents.gen)
    clu$select.chromosome('chr1')
    clu$train()
    tc <- clu$transitive.closure(clu$pairwise.similarities())
    \dontrun{
    write.clusters(tc, "myclustering.txt")
    print(all.equal(read.clusters("myclustering.txt"), tc))
    }
}
