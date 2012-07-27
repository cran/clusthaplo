library(ggplot2)
library(grid)

.gp <- function(background='gray80') {
    (ggplot()
     + theme_bw()
     + theme(panel.background=element_rect(fill=background, size=0),
             panel.grid=element_blank(),
             panel.grid.major=element_blank(),
             panel.grid.minor=element_blank(),
             panel.margin=element_blank(),
             axis.ticks.length=unit(.1, "cm")
       )
    )
}

.rect.mat <- function(locus.vec, col.mat, step=1) {
    x <- NULL  # make Rcheck happy
    y <- NULL  # make Rcheck happy
    height <- NULL  # make Rcheck happy
    width <- c(diff(locus.vec), step)
    make.height <- function(col.vec) {
        height <- rep(.8, length(width))
        height[col.vec==0] <- .01
        height
    }
    xvec <- rep(locus.vec + width * .5, ncol(col.mat))
    wvec <- rep(width, ncol(col.mat))
    hvec <- as.vector(apply(col.mat, 2, make.height))
    yvec <- rep(1:ncol(col.mat), each=length(locus.vec))
    fillvec <- 1 + as.integer(as.vector(col.mat))
    d <- data.frame(x=xvec, y=yvec, width=wvec, height=hvec)
    geom_tile(aes(x=x, y=y, width=width, height=height), fill=fillvec, data=d)
}

.draw.chrom <- function(tc, ord=1:ncol(tc), expand.cliques=F, step.size=1, from=min(attr(tc, 'loci')), to=max(attr(tc, 'loci')), background='gray40') {
    x <- NULL  # make Rcheck happy
    y <- NULL  # make Rcheck happy
    cl <- .build.cliq.struc(.cliques(tc))
    if(length(ord) != ncol(cl$repr)) {
        # need to rebuild cliques
        reord <- 1:ncol(cl$repr)
        reord[ord] <- 1:length(ord)
        sel <- reord[.select(.useable.index(cl), unique(c(ord, 1:ncol(cl$repr))))]
        repr <- matrix(sel[cl$repr], nrow=nrow(cl$repr), ncol=ncol(cl$repr))[, ord]
        cl <- .build.cliq.struc(.cliques(repr))
    }
    names <- attr(tc, "ind.names")
    if(expand.cliques == 'compare') {
        tmp <- .maximize.repr(cl)[, ord]
        R <- matrix(0, ncol=2 * ncol(tmp), nrow=nrow(tmp))
        R[, (1:ncol(tmp)) * 2] <- tmp
        R[, (1:ncol(tmp)) * 2 - 1] <- cl$repr
        tmp <- names
        #print(tmp)
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

    .colors <- matrix(.DSAT(cl, R)[R], ncol=ncol(R))
    #print(dim(.colors))

    loci <- attr(tc, 'loci')

    to.keep <- loci <= to & loci >= from

    .colors <- .colors[to.keep, ]
    loci <- loci[to.keep]

    ydf <- data.frame(y=1:ncol(.colors))
    #geom.chromlines <- geom_hline(yintercept=1:ncol(.colors))

    cldf <- data.frame(x=rep(c(from, to), ncol(R)),
                       y=rep(1:ncol(R), each=2))
    geom.chromlines <- geom_line(aes(x=x, y=y, group=y),
                                 data=cldf)
    breaks <- 1 + which(apply(diff(R)!=0, 1, any))

    xl <- xlim(from, to)
    yl <- #ylim(0, 1 + ncol(.colors))
          ylim(names)

    ncolors <- max(max(.colors) + 1, 2)

    .info.out("Using", ncolors, "colors.\n")

    ret <- (.gp(background=background)
            + xl + yl
            + xlab("Map (cM)")
            #+ ylab("Haplotypes")
            + ylab("")
            + geom.chromlines
           )
    if(length(breaks) > 0) {
        geom.loci <- geom_vline(xintercept=loci[-breaks],
                                  colour='#404040C0', size=.1)
        geom.chrom <- .rect.mat(loci, .colors, step.size)
        geom.breaks <- geom_vline(xintercept=loci[breaks],
                                  colour='#000000C0', size=.2)
        ret <- ret + geom.loci + geom.breaks + geom.chrom
    } else {
        geom.loci <- geom_vline(xintercept=loci,
                                  colour='#404040C0', size=.1)
        ret <- ret + geom.loci
    }
    ret + scale_colour_manual(values=rainbow(ncolors))
}
