\name{select.chromosome}
\alias{select.chromosome}
\title{Select a chromosome to operate on}
\usage{
    select.chromosome(chrom.name)
}
\note{
    This function can only be run from inside a clusthaplo context as returned by clusthaplo(...).
}
\arguments{
    \item{chrom.name}{
        A valid chromosome name, or \code{NA}.
    }
}
\description{
    Selects a chromosome and precomputes some relevant data. \code{NA} stands for no selection.
}
