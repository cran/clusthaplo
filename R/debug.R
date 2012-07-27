.. <- new.env(hash=TRUE)

.debug <- F
.verbose <- F

..$.debug <- F
..$.verbose <- F


set.debug <- function(dbg) {
    .debug <<- !!dbg
}


set.verbose <- function(vbse) {
    .verbose <<- !!vbse
}


.dumper <- function(flag, ...) {
    if (flag) {
        argv <- list(...)
        if (length(argv) == 1 && is.list(argv[[1]])) {
            print(argv[[1]])
        } else {
            cat(...)
        }
    }
}


.debug.out <- function(...) .dumper(.debug, ...)
.info.out <- function(...) .dumper(.verbose, ...)


environment(set.debug) <- ..
environment(set.verbose) <- ..
environment(.debug.out) <- ..
environment(.info.out) <- ..
