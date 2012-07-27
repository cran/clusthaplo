#source('dsat.R')

.cliques <- function(effet) {
    one_clique <- function(l, i) which(effet[l, ] == i)
    clique_list <- function(l) {
        x <- list()
        xi <- 1
        #for(i in 1:dim(effet)[2]) {
        for(i in sort(unique(effet[l, ]))) {
            o <- one_clique(l, i)
            if(length(o) > 0) {
                #cat(l, i, ":", o, "\n")
                x[[xi]] <- o
                xi <- xi + 1
            }
        }
        x
    }  
    lapply(1:dim(effet)[1], clique_list)
}


.clique.index <- function(cliqlist) {
    unique(unlist(cliqlist, recursive=F))
}


.index.of.clique <- function(cliq, cliqindex) {
    for(i in 1:length(cliqindex)) {
        ci <- cliqindex[[i]]
        if(length(ci) == length(cliq) && all(cliq == ci)) {
            return(i)
        }
    }
}


.cliq.table <- function(idx, cli) {
    n.ind <- max(sapply(idx, max))
    ret <- matrix(ncol=n.ind, nrow=length(cli))
    for(l in 1:length(cli)) {
        for(c in cli[[l]]) {
            ret[l, idx[[c]]] <- c
        }
    }
    ret
}


.build.cliq.struc <- function(cliqlist) {
    ci <- .clique.index(cliqlist)
    cliqlist <- lapply(cliqlist, function(c) { sapply(c, .index.of.clique, ci) })
    list(index=ci,
         cliques=cliqlist,
         n.ind=max(sapply(ci, max)),
         repr=.cliq.table(ci, cliqlist))
}


.useable.index <- function(cli) {
    vec <- 1:cli$n.ind
    lapply(cli$index, function(k) vec %in% k)
}

.select <- function(ui, ord) {
    unlist(lapply(ui, function(v) ord[v[ord]][1]))
}


.cliq.include.matrix <- function(ui) {
    ret <- matrix(F, nrow=length(ui) + 1, ncol=length(ui) + 1)
    for(i in 1:(length(ui) - 1)) {
        csup <- ui[[i]]
        if(sum(csup) == 1) {
            next
        }
        for(j in (i + 1):length(ui)) {
            cinf <- ui[[j]]
            if(sum(cinf) == 1) {
                next
            }
            ret[i, j] <- all((csup & cinf) == cinf) || all((cinf & csup) == csup)
            ret[j, i] <- ret[i, j]
            #ret[i, j] <- all((csup & cinf) == cinf) || all((csup & cinf) == csup)
            #ret[j, i] <- ret[i, j]
        }
    }
    ret
}


.inc.test <- function(cim, r1, r2) {
    u1 <- unique(r1)
    u2 <- unique(r2)
    m1 <- matrix(0, nrow=length(u1), ncol=length(r1))
    m2 <- matrix(0, ncol=length(u2), nrow=length(r2))
    for(i in 1:nrow(m1)) {
        if(sum(r1 == u1[i]) == 1) {
            next
        }
        m1[i, r1 == u1[i]] <- 1
    }
    for(i in 1:ncol(m2)) {
        if(sum(r2 == u2[i]) == 1) {
            next
        }
        m2[r2 == u2[i], i] <- 1
    }
    #print(m1)
    #print(m2)
    ret <- m1 %*% m2
    #print(ret)
    ret <- ret * cim[u1, u2]
    rownames(ret) <- u1
    colnames(ret) <- u2
    ret
}


.select.inc <- function(inc) {
    #cat("--------------------------------------------\n")
    #print(dim(inc))
    if(ncol(inc) == 1 || nrow(inc) == 1) {
        x <- inc == max(inc)
        if(sum(x) > 1) {
            ret <- rep(F, length(x))
        } else {
            ret <- x
        }
    } else {
        ret.col <- t(apply(inc, 1,
                           function(r) {
                               x <- r == max(r)
                               if(sum(x) > 1) {
                                   rep(F, length(r))
                               } else {
                                   x
                               }
                           }))
        ret.row <- apply(inc, 2,
                         function(r) {
                             x <- r == max(r)
                             if(sum(x) > 1) {
                                 rep(F, length(r))
                             } else {
                                 x
                             }
                         })
        #print(ret.row)
        #print(ret.col)
        ret <- ret.row & ret.col
    }
    ret <- matrix(which(ret, arr.ind=T), ncol=2)
    ret[, 1] <- as.numeric(rownames(inc)[ret[, 1]])
    ret[, 2] <- as.numeric(colnames(inc)[ret[, 2]])
    rownames(ret) <- NULL
    ret
}


.maximize.repr <- function(cl) {
    cim <- .cliq.include.matrix(.useable.index(cl))
    sel.inc <- lapply(1:(nrow(cl$repr)-1),
                      function(i) {
                          x <- .inc.test(cim, cl$repr[i,], cl$repr[i+1,])
                          if(sum(x)==0) NULL else .select.inc(x)
                      })

    repr <- cl$repr
    output <- repr
    shift <- length(cl$index)
    for(.prev in 1:(nrow(repr) - 1)) {
        .next <- .prev + 1
        inc <- sel.inc[[.prev]]
        w <- repr[.next, ] == repr[.prev, ]
        output[.next, w] <- output[.prev, w]
        if(is.null(inc) || length(inc) == 0) {
            next
        }
        for(i in 1:nrow(inc)) {
            w <- repr[.next, ] == inc[i, 2]
            prev.cliq <- output[.prev, repr[.prev, ] == inc[i, 1]][1]
            if(prev.cliq %in% output[.next, ]) {
                next
            }
            output[.next, w] <- prev.cliq
        }
    }
    output
}



.draw.chromosome.OBSOLETE <- function(cl,
                                      names=as.character(1:ncol(cl$repr)),
                                      ord=1:ncol(cl$repr),
                                      loci=1:nrow(cl$repr),
                                      from=1,
                                      to=nrow(cl$repr),
                                      tick.interval=25,
                                      bg='gray40',
                                      expand.cliques=F,
                                      ...) {
    if(length(ord) != ncol(cl$repr)) {
        # need to rebuild cliques
        reord <- 1:ncol(cl$repr)
        reord[ord] <- 1:length(ord)
        sel <- reord[.select(.useable.index(cl), unique(c(ord, 1:ncol(cl$repr))))]
        repr <- matrix(sel[cl$repr], nrow=nrow(cl$repr), ncol=ncol(cl$repr))[, ord]
        cl <- .build.cliq.struc(.cliques(repr))
    }
    if(expand.cliques == 'compare') {
        tmp <- .maximize.repr(cl)[, ord]
        R <- matrix(0, ncol=2 * ncol(tmp), nrow=nrow(tmp))
        R[, (1:ncol(tmp)) * 2] <- tmp
        R[, (1:ncol(tmp)) * 2 - 1] <- cl$repr
        tmp <- names
        names[(1:length(tmp)) * 2] <- paste(tmp, ".m", sep="")
        names[(1:length(tmp)) * 2 - 1] <- tmp
        tmp <- ord
        ord[(1:length(tmp)) * 2] <- tmp * 2
        ord[(1:length(tmp)) * 2 - 1] <- tmp * 2 - 1
    } else if(expand.cliques) {
        R <- .maximize.repr(cl)[, ord]
    } else {
        R <- cl$repr[, ord]
    }
    plot.new()
    locus.max <- loci[to]
    sup.args <- list(...)
    pw.args <- list(xlim=c(loci[from], locus.max), ylim=c(length(ord) - .5, 0), mar=c(0, 0, 0, 0), oma=c(0, 0, 0, 0))
    for(n in names(sup.args)) {
        pw.args[[n]] <- sup.args[[n]]
    }
    do.call(plot.window, pw.args)

    if (!is.null(bg)) {
        rect(xleft=loci[from], xright=locus.max, ytop=length(ord) + .1, ybottom=-.6, col=bg)
    }

    .colors <- .DSAT(cl)
    repr <- matrix(.colors[R], ncol=ncol(R), nrow=nrow(R))

    #if(length(palette()) < max(repr)) {
        palette(rainbow(max(repr) + 1))
        #palette(topo.colors(max(repr)))
    #}

    rect.limits <- lapply(1:ncol(repr),
                          function(i) {
                              r <- repr[from:to, i]
                              wd <- which(diff(r) != 0)
                              rbind(c(1, wd), c(wd, length(r)), r[c(1, 1 + wd)])
                          })
    .info.out("Using order", ord, "\n")

    .ticks <- unique(c(seq(0, loci[length(loci)], by=tick.interval),
                       loci[length(loci)]))
    .labels <- as.character(round(.ticks, 1))
    axis(3,  # top
         labels=.labels,
         line=.1,
         at=.ticks)

    axis(1,  # bottom
         labels=.labels,
         line=.1,
         at=.ticks)

    .ticks <- (1:length(ord)) - 0.77
    .labels <- names[ord]
    axis(2,  # left
         las=1,
         line=-1,
         lwd=0,
         lwd.ticks=0,
         labels=.labels,
         hadj=1,
         at=.ticks)

    breaks <- loci[apply(diff(R)!=0, 1, any)]
    abline(v=breaks, col='gray80')

    #col.ord <- sort.int(ord, index.return=T)$ix
    all.rects <- lapply(1:length(ord),
                        function(ind) {
                            r <- rect.limits[[ord[ind]]]
                            color <- r[3, ]
                            ybottom <- ind - 1 + .22 * (color == 0)
                            ytop <- ind - .5 - .22 * (color == 0)
                            color[color == 0] <- 'black'
                            list(xleft=loci[from - 1 + r[1,]],
                                 xright=loci[from - 1 + r[2,]],
                                 color=color,
                                 ybottom=ybottom,
                                 ytop=ytop)
                        })

    color <- unlist(lapply(1:length(ord), function(i) all.rects[[i]]$color))
    ytop <- unlist(lapply(1:length(ord), function(i) all.rects[[i]]$ytop))
    ybottom <- unlist(lapply(1:length(ord), function(i) all.rects[[i]]$ybottom))
    xleft <- unlist(lapply(1:length(ord), function(i) all.rects[[i]]$xleft))
    xright <- unlist(lapply(1:length(ord), function(i) all.rects[[i]]$xright))
    rect(xleft=xleft, ybottom=ybottom, xright=xright, ytop=ytop, col=color, border=color)

    invisible(all.rects)
}


.compare.repr <- function(r1, r2) {
    cl1 <- .build.cliq.struc(.cliques(r1))
    cl2 <- .build.cliq.struc(.cliques(r2))
    if(any(dim(r1) != dim(r2))) {
        stop("Different sizes")
    }
    ret <- T
    for(i in 1:nrow(r1)) {
        for(j in 1:ncol(r1)) {
            cat(i, j, "       \r")
            flush(stdout())
            c1 <- cl1$index[[cl1$repr[i, j]]]
            c2 <- cl2$index[[cl2$repr[i, j]]]
            if(length(c1) != length(c2) || any(c1 != c2)) {
                warning("Cliques differ at locus (", i, ", ", j, ") left=(", paste(c1, collapse=", "), ") right=(", paste(c2, collapse=", "), ")")
                ret <- F
            }
            if(! j %in% c1) {
                warning("Wrong clique assignment in first at locus (", i, ", ", j, "); assigned clique contains (", paste(c1, collapse=", "), ")")
                ret <- F
            }
            if(! j %in% c2) {
                warning("Wrong clique assignment in second at locus (", i, ", ", j, "); assigned clique contains (", paste(c2, collapse=", "), ")")
                ret <- F
            }
        }
    }
    cat("          \n")
    ret
}

.ascii.repr <- function(repr) {
    L <- c(LETTERS)
    L <- c(L, unlist(lapply(L, function(l) lapply(L, function(m) paste(l, m, sep="")))))
    matrix(L[repr], ncol=ncol(repr), nrow=nrow(repr))
}
