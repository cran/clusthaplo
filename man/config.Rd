\name{config}
\alias{config}
\title{Configure a clusthaplo context.}
\usage{
    config(w1=NULL, w2=NULL, step.size=NULL, window.length=NULL, na.replace=NULL,
           scoring.method=NULL,
           kinship.threshold=NULL,
           clustering.method=NULL,
           clustering.threshold=NULL,
           simulation.type=NULL, simulation.Np=NULL, simulation.Ng=NULL,
           simulation.Nrep=NULL,
           threshold.quantile=NULL,
           xml.output.dir=NULL,
           plot.output.dir=NULL, plot.height=NULL, plot.width=NULL,
           config=NULL)
}
\note{
    This function can only be run from inside a clusthaplo context as returned by clusthaplo(...).
}
\arguments{
    \item{w1}{
        The w1 weight function in the Li&Jyang similarity score. Possible
        values are \code{"kernel.const"}, \code{"kernel.exp"}, \code{"kernel.gauss"},
        \code{"kernel.unif"}, \code{"kernel.laplace"} or \code{"kernel.null"}.
    }
    \item{w2}{
        The w2 weight function in the Li&Jyang similarity score. Possible
        values are \code{"kernel.const"}, \code{"kernel.exp"}, \code{"kernel.gauss"},
        \code{"kernel.unif"}, \code{"kernel.laplace"} or \code{"kernel.null"}.
    }
    \item{step.size}{
        The step size in centi-Morgans between scanned loci.
    }
    \item{window.length}{
        The total window length in centi-Morgans.
    }
    \item{na.replace}{
        The policy regarding NA's when computing the Identity-by-State function for a pair of haplotypes.
        Can be \code{TRUE}, \code{FALSE}, or \code{NA}.
        \itemize{
            \item \code{TRUE}: \code{NA} is the same allele.
            \item \code{FALSE}: \code{NA} is always a new (different) allele.
            \item \code{NA}: \code{NA} is treated as actually unknown.
        }
    }
    \item{scoring.method}{
        One of \code{"raw"}, \code{"normalized"}, or \code{"kinship"}.
    }
    \item{clustering.method}{
        One of \code{"threshold"}, \code{"hmm"}.
    }
    \item{kinship.threshold}{
        A positive integer representing the number of markers in a window below which the kinship data will be used.
    }
    \item{clustering.threshold}{
        The value above which a pair of haplotypes is considered similar at a given locus,
        if the \code{clustering.method} is \code{"threshold"}.
    }
    \item{simulation.type}{
        The type of simulation used for the training. One of \code{"equi"} or \code{"mosaic"}.
    }
    \item{simulation.Np}{
        The number of parents to simulate for the training.
    }
    \item{simulation.Ng}{
        The number of intermediary generations to simulate for the training (only relevant for \code{"mosaic"}).
    }
    \item{simulation.Nrep}{
        The number of replicates to simulate for the training.
    }
    \item{threshold.quantile}{
        The quantile to use to select the threshold automatically. It must be either a character string in the form "xx\%" or a plain integer xx with 80 <= xx <= 100.
    }
    \item{xml.output.dir}{
        Where to output the XML files for MCQTL.
    }
    \item{plot.output.dir}{
        Where to output the generated graphics.
    }
    \item{plot.width}{
        Default width for the generated graphics.
    }
    \item{plot.height}{
        Default height for the generated graphics.
    }
    \item{config}{
        Use this parameter to recall a complete configuration as returned by a call to \code{config()}.
        It must be a data.frame or a filename pointing to a file containing a data.frame.
    }
}
\description{
    Configures a clusthaplo context.
}
\value{
    If no argument is provided, returns a data.frame containing all the current configuration values.
}
\note{
    \code{"hmm"} method is only available if the package \code{RHmm} is installed.
}
\seealso{
    \code{\link{clusthaplo}}
    \code{\link{get.config}}
}
