.cliq.graph <- function(cl, R=cl$repr) {
    graph <- matrix(F, ncol=length(cl$index), nrow=length(cl$index))
    #R <- cl$repr

    cliqz <- unique(R[1, ])

    graph[cliqz, cliqz] <- T

    for(pos in 2:nrow(R)) {
        new.cliqz <- unique(R[pos, ])
        graph[new.cliqz, new.cliqz] <- T
        graph[new.cliqz, cliqz] <- T
        graph[cliqz, new.cliqz] <- T
        cliqz <- new.cliqz
    }
    graph[col(graph) == row(graph)] <- F
    graph
}


.singletons <- function(cl) {
    unlist(lapply(1:length(cl$index),
                  function(cli.num) {
                      if (length(cl$index[[cli.num]]) == 1) {
                          cli.num
                      }
                  }))
}



.DSAT <- function(cl, R=cl$repr) {
    graph <- .cliq.graph(cl, R)
    color.vec <- rep(NA, ncol(graph))
    singlz <- .singletons(cl)
    color.vec[singlz] <- 0  # forcefully set this color.

    degrees <- colSums(graph)

    dsat.value <- function(v) {
        color.nei <- unique(na.omit(color.vec[graph[, v]]))
        color.nei <- color.nei[color.nei != 0]
        if (length(color.nei) == 0) {
            return(degrees[v])
        } else {
            return(sum(color.nei != 0))  # modified saturation
        }
    }

    pick.vertex <- function() {
        candidates <- is.na(color.vec)
        if(!any(candidates)) {
            return(NULL)
        }
        candidates <- which(candidates)
        dsat.vec <- sapply(candidates, dsat.value)
        strongest <- candidates[dsat.vec == max(dsat.vec)]
        #cat("candidate vertices", candidates, ".DSAT", max(dsat.vec), "\n")
        if (length(strongest) > 1) {
            ret <- strongest[which.max(degrees[strongest])]
        } else {
            ret <- strongest
        }
        #cat("selecting vertex", ret, "\n")
        ret
    }

    pick.color <- function(v) {
        color.nei <- unique(na.omit(color.vec[graph[, v]]))
        if(max(c(0, color.nei)) == 0) {
            return(1)
        }
        candidates <- 1:(1+length(color.nei))
        ret <- which(!(candidates %in% color.nei))[1]
        #cat("selecting color", ret, "\n")
        ret
    }

    while (!is.null(v <- pick.vertex())) {
        color.vec[v] <- pick.color(v)
        #cat("color[", v, "] = ", color.vec[v], "\n", sep="")
        #cat(color.vec, "\n")
    }
    color.vec
}

