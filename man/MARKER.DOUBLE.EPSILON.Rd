\name{MARKER.DOUBLE.EPSILON}
\alias{MARKER.DOUBLE.EPSILON}
\title{package constant: MARKER.DOUBLE.EPSILON}
\usage{
MARKER.DOUBLE.EPSILON
}
\description{
    This variable defines how close two markers should be to be considered the same
    locus by clusthaplo. Often double markers will be exactly at distance zero in
    an input genetic map, but after merging and realigning maps, there might be some
    slight numeric differences.
}
\examples{
    print(MARKER.DOUBLE.EPSILON)
}
