\name{transitive.closure}
\alias{transitive.closure}
\title{Performs transitive closure.}
\usage{
    transitive.closure(pair.simi, Npar=ncol(.gen))
}
\note{
    This function can only be run from inside a clusthaplo context as returned by clusthaplo(...).
}
\arguments{
    \item{pair.simi}{The result of a call to \code{\link{pairwise.similarities}}.}
    \item{Npar}{The number of haplotypes for which the pairwise similarities were computed.}
}
\description{
    Performs transitive closure on the whole set of pairwise similarities for one chromosome.
    Transitive closure means that, at a given locus, if haplotypes A and B match, and A and C match, then
    B and C are assumed to match also.
}
\value{
    A matrix with one column per parent and one row per scanned locus, giving, at each locus and
    for each parent, the smallest parent index belonging in their clique. For instance, if haplotypes 1 3 4 look
    alike amidst 5 haplotypes, and 2 and 5 are singletons, the row will be \code{c(1, 2, 1, 1, 5)}.
}
\seealso{
    \code{\link{pairwise.similarities}}
}
\examples{
    data(parents.map, parents.gen)
    clu <- clusthaplo(parents.map, NULL, parents.gen)
    clu$select.chromosome('chr1')
    clu$train()
    tc <- clu$transitive.closure(clu$pairwise.similarities())
    print(head(tc))
}
