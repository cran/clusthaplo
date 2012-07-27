.kernel.const <- function(window.length, k=1, name='kernel.const') {
    # This function returns a constant kernel function
    ret <- function(xlist) rep(k, length(xlist))
    attributes(ret) <- list(window.length=window.length, k=k, name=name)
    ret
}


.kernel.exp <- function(window.length, lambda=(2 * qexp(.95 / 2)) / window.length) {
    # This function returns an exponentially decreasing kernel function

    # FIXME: only to compare with old ClusthaploV6 code
    #lambda <- 4 / window.length

    ret <- function(xlist) dexp(abs(xlist), lambda)
    attributes(ret) <- list(window.length=window.length, lambda=lambda, name="kernel.exp")
    ret
}


.kernel.laplace <- function(window.length, lambda=qexp(.975) / window.length) {
    # This function returns a Laplace kernel function
    ret <- function(xlist) dexp(abs(xlist), lambda) / 2
    attributes(ret) <- list(window.length=window.length, lambda=lambda, name="kernel.laplace")
    ret
}


.kernel.gauss <- function(window.length, std.dev=window.length / (2 * qnorm(.975))) {
    # This function returns a gaussian kernel function
    ret <- function(xlist) dnorm(xlist, 0, std.dev)
    attributes(ret) <- list(window.length=window.length, std.dev=std.dev, name="kernel.gauss")
    ret
}


.kernel.unif <- function(window.length) {
    # This function returns a uniform kernel function
    .kernel.const(window.length, .95 / window.length, "kernel.unif")
}


.kernel.null <- function(window.length, k=1) {
    # This function returns a null kernel function
    .kernel.const(window.length, 0, "kernel.null")
}


.all.kernels = list(kernel.const=.kernel.const,
                    kernel.exp=.kernel.exp,
                    kernel.gauss=.kernel.gauss,
                    kernel.unif=.kernel.unif,
                    kernel.laplace=.kernel.laplace,
                    kernel.null=.kernel.null)
