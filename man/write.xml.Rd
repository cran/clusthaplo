\name{write.xml}
\alias{write.xml}
\title{Output for MCQTL-LD}
\usage{
    write.xml(closure)
}
\arguments{
    \item{closure}{
        The result of a transitive closure of pair similarities as returned by \code{\link{transitive.closure}}.
    }
}
\description{
    Writes output suitable for use with MCQTL\_LD. One XML file is written per parent with path \code{xml.output.dir/chromosome.name/LD_file_<parent.name>.XML}.
}
\seealso{
    \code{\link{config}}
    \code{\link{transitive.closure}}
}
\examples{
    data(parents.map, parents.gen)
    clu <- clusthaplo(parents.map, NULL, parents.gen)
    clu$select.chromosome('chr1')
    clu$train()
    tc <- clu$transitive.closure(clu$pairwise.similarities())
    \dontrun{
    write.xml(tc)
    }
}
