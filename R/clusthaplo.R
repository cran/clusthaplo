# The following declarations only aim at appeasing R CMD CHECK. Stupid command
# thinks it's smart.
.RHmm.available <- NULL
.clustering.method <- NULL
.clustering.threshold <- NULL
.clusterize <- function(...) stop("Don't use this.")
.desc.map <- NULL
.gen <- NULL
.kinship <- NULL
.na.replace <- NULL
.par.map <- NULL
.plot.height <- NULL
.plot.output.dir <- NULL
.plot.width <- NULL
.Pt <- NULL
.scoring.method <- NULL
.selected.chromosome <- NULL
.simi.data <- NULL
.simulation.Ng <- NULL
.simulation.Np <- NULL
.simulation.Nrep <- NULL
.simulation.type <- NULL
.Smap <- NULL
.Smap.max <- NULL
.step.size <- NULL
.test.map <- NULL
.threshold.quantile <- NULL
.training.results <- NULL
.w1 <- NULL
.w1.fun <- NULL
.w2 <- NULL
.w2.fun <- NULL
.window.length <- NULL
.xml.output.dir <- NULL
.is.trained <- NULL
.kinship.threshold <- NULL

.wrong.context <- "Please use this function from inside a clusthaplo context, as returned by clusthaplo(...)."


count.ancestral.alleles <- function(clusters) {
    if(! 'clusthaplo.haplotypes' %in% class(clusters)) {
        warning("Please provide the result of a transitive closure in order to count ancestral alleles.")
        return(NA)
    }
    ret <- data.frame(locus=attr(clusters, 'loci'),
                      count=apply(clusters, 1, function(r) length(unique(r))))
    attr(ret, 'chromosome') <- attr(clusters, 'chromosome')
    class(ret) <- c('clusthaplo.alleles', class(ret))
    ret
}


.as.index <- function(i) {
    if (!exists('.clusthaplo.tag')) {
        stop(.wrong.context)
    }
    if(length(i) != 1) {
        stop('Only single values are allowed as haplotype index/name')
    }
    if(is.character(i)) {
        which(colnames(.gen) == i)
    } else {
        i
    }
}

pair.similarity <- function(i, j) {
    if (!exists('.clusthaplo.tag')) {
        stop(.wrong.context)
    }
    if(is.na(.selected.chromosome)) {
        stop("Please select a chromosome first.")
    }
    i <- .as.index(i)
    j <- .as.index(j)
    raw.simi <- .compute.similarity(.test.map[[.selected.chromosome]], .simi.data, .gen,
                                    i, j, .w1.fun, .w2.fun,
                                    .na.replace)
    if(.scoring.method == 'raw') {
        return(raw.simi)
    } else {
        if(is.null(.Smap)) {
            .Smap <<-  .compute.similarity(.test.map[[.selected.chromosome]], .simi.data, .gen,
                                           1, 1, .w1.fun, .w2.fun,
                                           .na.replace)
            small.windows <- .simi.data$density <= .kinship.threshold
            #print(.kinship.threshold)
            #print(data.frame(dens=.simi.data$density, small=small.windows))
            if(any(small.windows)) {
                .Smap.max <<- max(.Smap$similarity[small.windows], na.rm=T)
                .Pt <<- .Smap$similarity / .Smap.max
                .Pt[!small.windows] <<- 1
            } else {
                .Smap.max <<- 1
                .Pt <<- rep(1, length(.Smap$similarity))
            }
        }
        # normalized
        raw.simi$similarity <- raw.simi$similarity / .Smap$similarity
        raw.simi$similarity[.simi.data$density == 0] <- 0   # prevent spurious NaN's
        if(.scoring.method == 'kinship') {
            raw.simi$similarity <- .Pt * raw.simi$similarity + (1 - .Pt) * .kinship[i, j]
        }
    }
    if(any(is.nan(raw.simi$similarity) | is.na(raw.simi$similarity))) {
        print(data.frame(simi=raw.simi$similarity, dens=raw.simi$density))
    }
    raw.simi
}


select.chromosome <- function(chrom.name=NULL) {
    if (!exists('.clusthaplo.tag')) {
        stop(.wrong.context)
    }
    #.clu <- parent.env(environment())
    if (is.na(chrom.name) || is.null(chrom.name)) {
        .selected.chromosome <<- NA
        .simi.data <<- NULL
    } else if (chrom.name %in% names(.test.map)) {
        .selected.chromosome <<- chrom.name
        .simi.data <<- .prepare.similarity.data(.desc.map[[chrom.name]],
                                                .test.map[[chrom.name]],
                                                .step.size, .window.length)
        .Smap <<- NULL
        .Smap.max <<- NULL
        .Pt <<- NULL
    } else {
        stop('Invalid chromosome name "', chrom.name, '". Expected one of ',
             paste(names(.test.map), collapse=", "))
    }
}


get.window.density <- function() {
    if (!exists('.clusthaplo.tag')) {
        stop(.wrong.context)
    }
    if(is.na(.selected.chromosome)) {
        stop("Please select a chromosome first.")
    }
    ret <- list(locus=.simi.data$loci, density=.simi.data$density)
    attr(ret, 'window.length') <- .window.length
    attr(ret, 'chromosome') <- .selected.chromosome
    class(ret) <- c('clusthaplo.density', class(ret))
    ret
}


.generic.get.map <- function(map, chrom.name) {
    if (is.null(chrom.name)) {
        lapply(map, function(m) m[, c('locus', 'markers')])
    } else if (chrom.name %in% names(map)) {
        map[[chrom.name]][, c('locus', 'markers')]
    } else {
        stop('Invalid chromosome name "', chrom.name, '". Expected one of ',
             paste(names(map), collapse=", "))
    }
}


get.haplotypes.map <- function(chrom.name=NULL) {
    if (!exists('.clusthaplo.tag')) {
        stop(.wrong.context)
    }
    .generic.get.map(.par.map, chrom.name)
}


get.mcqtl.consensus.map <- function(chrom.name=NULL) {
    if (!exists('.clusthaplo.tag')) {
        stop(.wrong.context)
    }
    .generic.get.map(.desc.map, chrom.name)
}


get.scan.map <- function(chrom.name=NULL) {
    if (!exists('.clusthaplo.tag')) {
        stop(.wrong.context)
    }

    mk.scan <- function(name) {
        p <- data.frame(locus=.test.loci(.desc.map[[name]], .step.size))
        merge(p, .generic.get.map(.desc.map, name), all.x=T)
    }

    if(is.null(chrom.name)) {
        lapply(get.chromosome.names(), mk.scan)
    } else if(chrom.name %in% names(.desc.map)) {
        mk.scan(chrom.name)
    } else {
        stop('Invalid chromosome name "', chrom.name, '". Expected one of ',
             paste(names(.desc.map), collapse=", "))
    }
}


get.config <- function(...) {
    if (!exists('.clusthaplo.tag')) {
        stop(.wrong.context)
    }

    vars <- as.character(unlist(list(...)))
    if(length(vars) == 0) {
        vars <- .legal.names
    }
    illeg <- !vars %in% .legal.names
    if(any(illeg)) {
        stop('Invalid names supplied : ', paste(vars[illeg], sep=', '))
    }
    .get <- function(n) get(paste('.', n, sep=''))
    if(length(vars) == 1) {
        ret <- .get(vars)
    } else {
        ret <- lapply(vars, .get)
        names(ret) <- vars
    }
    ret
}


get.selected.chromosome <- function() {
    if (!exists('.clusthaplo.tag')) {
        stop(.wrong.context)
    }
    .selected.chromosome
}

get.chromosome.names <- function() {
    if (!exists('.clusthaplo.tag')) {
        stop(.wrong.context)
    }
    names(.test.map)
}

get.marker.data <- function() {
    if (!exists('.clusthaplo.tag')) {
        stop(.wrong.context)
    }
    .gen
}

get.kinship.matrix <- function() {
    if (!exists('.clusthaplo.tag')) {
        stop(.wrong.context)
    }
    .kinship
}



.generic.run.simulation <- function(n.replicates, geno.simulator,
                                    collect.replicate=function(r) r,
                                    collect.all=function(a) a) {
    if (!exists('.clusthaplo.tag')) {
        stop(.wrong.context)
    }
    do.a.chrom <- function(chrom.name, sim.clu) {
        .info.out("Simulating on chromosome", chrom.name, "\n")
        sim.clu$select.chromosome(chrom.name)
        apply(combn(ncol(sim.clu$.$.gen), 2), 2,
              function(r) sim.clu$pair.similarity(r[1], r[2])$similarity)
    }
    if(is.null(.selected.chromosome) || is.na(.selected.chromosome)) {
        collect.all(lapply(1:n.replicates,
                           function(i) {
                               sim.clu <- clusthaplo(.par.map, .desc.map,
                                                     do.call(rbind, lapply(.test.map, geno.simulator, .gen)))
                               sim.clu$config(config=config())
                               collect.replicate(lapply(sim.clu$get.chromosome.names(),
                                                        do.a.chrom, sim.clu))
                           }))
    } else {
        chrom.p <- list(.par.map[[.selected.chromosome]])
        chrom.d <- list(.desc.map[[.selected.chromosome]])
        names(chrom.p) <- .selected.chromosome
        names(chrom.d) <- .selected.chromosome
        collect.all(lapply(1:n.replicates,
                           function(i) {
                               sim.clu <- clusthaplo(chrom.p, chrom.d,
                                                     geno.simulator(chrom.p[[1]], .gen))
                               sim.clu$config(config=config())
                               collect.replicate(do.a.chrom(.selected.chromosome, sim.clu))
                           }))
    }
}





run.raw.simulation <- function(n.replicates=3,
                               geno.sim) {
    if (!exists('.clusthaplo.tag')) {
        stop(.wrong.context)
    }
    .generic.run.simulation(n.replicates, geno.sim)
}


run.equi.simulation <- function(n.haplotypes=10, n.replicates=3,
                                quantiles=seq(.8, 1, by=.01)) {
    if (!exists('.clusthaplo.tag')) {
        stop(.wrong.context)
    }
    .generic.run.simulation(n.replicates, .gen.geno.equi(n.haplotypes),
                            collect.all=function(a) quantile(unlist(a), quantiles, na.rm=T))
}


run.mosaic.simulation <- function(n.haplotypes=10, n.generations=50,
                                  n.replicates=3,
                                  quantiles=seq(.8, 1, by=.01)) {
    if (!exists('.clusthaplo.tag')) {
        stop(.wrong.context)
    }
    .generic.run.simulation(n.replicates, .gen.geno.mosaic(n.haplotypes, n.generations),
                            collect.all=function(a) quantile(unlist(a), quantiles, na.rm=T))
}


train <- function(quantiles=seq(.8, 1, by=.01), ...) {
    if (!exists('.clusthaplo.tag')) {
        stop(.wrong.context)
    }
    if(.simulation.type == 'equi') {
        geno.sim <- .gen.geno.equi(.simulation.Np)
    } else if(.simulation.type == 'mosaic') {
        geno.sim <- .gen.geno.mosaic(.simulation.Np, .simulation.Ng)
    }
    if(.clustering.method == 'hmm') {
        # don't simulate, run with the actual marker data.
        train.set <- .generic.run.simulation(1, function(...) .gen,
                                             collect.all=function(a) as.vector(na.omit(unlist(a))))
        .debug.out("length of train set : ", length(train.set), "\n", sep="")
        hmms <- lapply(2:4, function(ns) RHmm::HMMFit(train.set, nStates=ns, dis="NORMAL", ...))
        .debug.out(hmms)
        BICs <- sapply(hmms, function(h) h$BIC)
        hmm4 <- RHmm::HMMFit(train.set, nStates=4, dis="NORMAL", ...)
        .debug.out("2-BIC =", BICs[1], "3-BIC =", BICs[2], "4-BIC =", BICs[3], "\n")
        if(all(is.nan(BICs))) {
            stop("Couldn't fit any HMM. Please select different settings or the 'threshold' clustering method.")
        }
        best.h <- which.min(BICs)
        hmm <- hmms[[best.h]]
        .info.out("Using ", best.h + 1, "-state HMM.\n", sep="")
        #if(hmm3$BIC <= hmm4$BIC) {
        #    hmm <- hmm3
        #    .info.out("Using 3-state HMM.\n")
        #} else {
        #    hmm <- hmm4
        #    .info.out("Using 4-state HMM.\n")
        #}
        means <- hmm$HMM$distribution$mean
        sort.idx <- sort.int(means, index.return=T, decreasing=T)$ix
        state.to.cluster <- sapply(1:length(sort.idx),
                                   function(i) which(sort.idx == i,
                                                     arr.ind=T))
        .training.results <<- list(state.2=hmms[[1]], state.3=hmms[[2]], state.4=hmms[[3]], best=hmm)
        .clusterize <<- function(similarity) {
            v <- RHmm::viterbi(hmm, similarity)
            state.to.cluster[v$states] == 1
        }
    } else {
        .training.results <<- .generic.run.simulation(.simulation.Nrep, geno.sim,
                                                      collect.all=function(a) quantile(unlist(a), quantiles, na.rm=T))
        .clustering.threshold <<- .training.results[.threshold.quantile]
        .clusterize <<- function(similarity) similarity >= .clustering.threshold
    }

}


get.training.results <- function() {
    if (!exists('.clusthaplo.tag')) {
        stop(.wrong.context)
    }
    .training.results
}


pairwise.similarities <- function(hook.pair.simi=function(i, j, ps) ps) {
    if (!exists('.clusthaplo.tag')) {
        stop(.wrong.context)
    }

    if(is.null(.clusterize)) {
        stop("Please train clusthaplo in order to setup the clustering scheme.")
    }

    apply(combn(ncol(.gen), 2, function(r) hook.pair.simi(r[1], r[2], pair.similarity(r[1], r[2]))$similarity),
          2, .clusterize)
}


cluster.pair <- function(i, j, simi=pair.similarity(i, j)) {
    if (!exists('.clusthaplo.tag')) {
        stop(.wrong.context)
    }
    i <- .as.index(i)
    j <- .as.index(j)

    if(is.null(.clusterize)) {
        stop("Please train clusthaplo in order to setup the clustering scheme.")
    }

    simi$clustering <- .clusterize(simi$similarity)
    attr(simi, 'chromosome') <- .selected.chromosome
    attr(simi, 'i') <- i
    attr(simi, 'j') <- j
    attr(simi, 'scoring.method') <- .scoring.method
    attr(simi, 'threshold') <- .clustering.threshold
    attr(simi, 'simulation.type') <- .simulation.type
    class(simi) <- c('clusthaplo.pair.similarity', class(simi))
    simi
}


transitive.closure <- function(pair.simi, Npar=ncol(.gen)) {
    if (!exists('.clusthaplo.tag')) {
        stop(.wrong.context)
    }
    Npos <- nrow(pair.simi)
    mat.cliq <- sapply(1:Npar, function(i) rep(i, Npos))
    clusters <- which(pair.simi, arr.ind=T)
    # clusters contains sequence of (locus, pair.num)
    pair.list <- combn(Npar, 2)
    if(length(clusters) > 0) {
        for(n in 1:nrow(clusters)) {
            pos <- clusters[n, 1]
            pair <- pair.list[, clusters[n, 2]]
            mat.cliq[pos, pair] <- min(mat.cliq[pos, pair])
        }
    }
    class(mat.cliq) <- c(class(mat.cliq), 'clusthaplo.haplotypes')
    attr(mat.cliq, 'chromosome') <- .selected.chromosome
    attr(mat.cliq, 'loci') <- .simi.data$loci
    attr(mat.cliq, 'ind.names') <- colnames(.gen)
    mat.cliq
}


write.clusters <- function(tc, file) {
    tab <- tc
    colnames(tab) <- attr(tc, "ind.names")
    rownames(tab) <- attr(tc, "loci")
    cat("Chromosome: ", attr(tc, 'chromosome'), "\n", sep="", file=file)
    write.table(tab, file, append=T, row.names=T, col.names=T)
}


read.clusters <- function(file) {
    header <- scan(file, what=character(), nlines=1)
    if(length(header)!=2 || header[1] != "Chromosome:") {
        stop("Couldn't read transitive closure from file ", file)
    }
    mat.cliq <- as.matrix(read.table(file, skip=1, header=T))
    class(mat.cliq) <- c(class(mat.cliq), 'clusthaplo.haplotypes')
    attr(mat.cliq, 'chromosome') <- header[2]
    attr(mat.cliq, 'loci') <- as.numeric(rownames(mat.cliq))
    attr(mat.cliq, 'ind.names') <- colnames(mat.cliq)
    attr(mat.cliq, "dimnames") <- NULL
    mat.cliq
}


.run.chromosome <- function(chrom, train.on.each.chromosome, plot.device, plot.chromosomes,
                            plot.pair.similarities, plot.alleles, plot.densities,
                            mk.filename) {
    .info.out("Running on chromosome ", chrom, "...\n", sep="")
    select.chromosome(chrom)
    if(train.on.each.chromosome) { train() }

    hook <- function(i, j, ps) {
        if(plot.pair.similarities) {
            p1 <- colnames(.gen)[i]
            p2 <- colnames(.gen)[j]
            plot.device(filename=mk.filename(chrom,
                                             "chromosome_", chrom,
                                             "_similarity_", p1, "_", p2),
                        width=.plot.width, height=.plot.height)
            plot(cluster.pair(i, j, ps))
            dev.off()
        }
        ps
    }

    tc <- transitive.closure(pairwise.similarities(hook.pair.simi=hook))

    write.clusters(tc, mk.filename(chrom, "clustering.txt", add.ext=F))

    if(plot.chromosomes) {
        plot.device(filename=mk.filename(chrom, "chromosome_", chrom),
                        width=.plot.width, height=.plot.height)
        plot(tc)
        dev.off()
    }

    if(plot.alleles) {
        plot.device(filename=mk.filename(chrom, "chromosome_", chrom, "_ancestral_alleles"),
                        width=.plot.width, height=.plot.height)
        plot(count.ancestral.alleles(tc))
        dev.off()
    }

    if(plot.densities) {
        plot.device(filename=mk.filename(chrom, "chromosome_", chrom, "_marker_density"),
                        width=.plot.width, height=.plot.height)
        plot(get.window.density())
        dev.off()
    }

    tc
}


run <- function(chromosomes=get.chromosome.names(), train.on.all.chromosomes=T,
                output.XML=T, plot.device=png,
                plot.alleles=T, plot.chromosomes=T, plot.pair.similarities=T, plot.densities=T,
                ...) {
    if (!exists('.clusthaplo.tag')) {
        stop(.wrong.context)
    }

    prefix <- file.path(getwd(), .plot.output.dir)
    if(!file.exists(prefix)) {
        dir.create(prefix)
    }

    hack.filename <- rev(as.character(formals(plot.device)$filename))[1]
    hack.ext <- substr(hack.filename, nchar(hack.filename) - 3, nchar(hack.filename))
    mk.filename <- function(chrom.name, ..., add.ext=T) {
        prefix <- file.path(prefix, chrom.name)
        if(!file.exists(prefix)) {
            dir.create(prefix)
        }
        if(add.ext) {
            ret <- file.path(prefix, paste(..., hack.ext, sep=""))
        } else {
            ret <- file.path(prefix, paste(..., sep=""))
        }
        cat("*", ret, "\n")
        ret
    }

    config(...)

    do.train <- is.null(.training.results) || is.na(.training.results) || .training.results != "manual setting"

    if(do.train && train.on.all.chromosomes) {
        select.chromosome(NA)
        train()
    }

    all.tc <- lapply(chromosomes,
                     function(chrom) .run.chromosome(chrom, do.train && !train.on.all.chromosomes,
                                                     plot.device, plot.chromosomes,
                                                     plot.pair.similarities,
                                                     plot.alleles, plot.densities,
                                                     mk.filename))


    if(output.XML) {
        write.xml(all.tc)
    }

    all.tc
}


use.multicore <- function() {
    if (!exists('.clusthaplo.tag')) {
        stop(.wrong.context)
    }
    if(!exists('mclapply')) {
        stop("Please load the 'multicore' package in order to benefit from parallelization")
    }
    lapply <<- get('mclapply')
}


.check.kinship <- function(kinship.matrix) {
    # check symmetry
    if(!all(kinship.matrix == t(kinship.matrix))) {
        warning("Kinship matrix is not symmetric. Using upper triangle.")
        kinship.matrix[lower.tri(kinship.matrix)] <- kinship.matrix[upper.tri(kinship.matrix)]
    }
    # force eigen values to be positive
    x <- eigen(kinship.matrix)
    corrected.values <- sapply(x$values, function(k) max(k, 0))
    if(!all(x$values == corrected.values)) {
        warning("Forcing Kinship matrix to be positive-semidefinite")
        kinship.matrix <- x$vectors %*% diag(corrected.values) %*% solve(x$vectors)
    }
    kinship.matrix
}


.compute.ais.kinship <- function(map.list, marker.data) {
    warning("Computing AIS kinship matrix.")
    kinship <- function(ij) {
        mean(unlist(lapply(map.list,
                    function(chrom) {
                        .compute.I(marker.data, ij[1], ij[2], chrom, na.replace=F)$I
                    })))
    }

    ret <- diag(ncol(marker.data))
    tri <- combn(ncol(marker.data), 2, kinship)
    ret[upper.tri(ret)] <- tri
    ret[lower.tri(ret)] <- t(ret)[lower.tri(ret)]
    ret
}


clusthaplo <- function(haplotypes.map,
                       mcqtl.consensus.map,
                       marker.data,
                       kinship.matrix=NULL,
                       na.strings=c(".", "NA", "-"),
                       discard.unknown.markers=T,
                       ...) {
    .clu <- new.env(hash=TRUE)
    .clu$.clusthaplo.tag <- 'clusthaplo'

    .clu$.RHmm.available <- F

    .clu$.selected.chromosome <- NA
    .clu$.scoring.method <- 'kinship'

    .clu$.w1 <- 'kernel.exp'
    .clu$.w2 <- 'kernel.unif'
    .clu$.step.size <- 1
    .clu$.window.length <- 20
    .clu$.na.replace <- T

    .clu$.w1.fun <- .kernel.exp(.clu$.window.length)
    .clu$.w2.fun <- .kernel.unif(.clu$.window.length)

    .clu$.kinship.threshold <- 10

    .clu$.is.trained <- F

    .clu$.clustering.method <- 'threshold'
    .clu$.simulation.type <- 'equi'
    .clu$.simulation.Nrep <- 3
    #.clu$.simulation.Np <- ncol(.clu$.gen)    # see below, after reading the genotypes...
    .clu$.simulation.Ng <- 50
    .clu$.threshold.quantile <- '95%'
    .clu$.clustering.threshold <- NULL

    .clu$.xml.output.dir <- 'XML'
    .clu$.plot.output.dir <- 'img'

    .clu$.plot.width <- 1024
    .clu$.plot.height <- 768

    # internals (precomputed data relevant to chromosome selection)

    .clu$.simi.data <- NULL
    .clu$.Smap <- NULL
    .clu$.Smap.max <- NULL
    .clu$.Pt <- NULL
    .clu$.training.results <- NULL
    .clu$.clusterize <- NULL

    # read data
    if(class(haplotypes.map) == 'list') {
        .clu$.par.map <- haplotypes.map
    } else if(class(haplotypes.map) == 'character' && length(haplotypes.map) == 1) {
        .clu$.par.map <- .read.map(haplotypes.map)
    } else {
        stop('Parents map MUST be a filename OR a named list of chromosome maps, not ', typeof(haplotypes.map), '.')
    }
    if(is.null(mcqtl.consensus.map)) {
        .clu$.desc.map <- .clu$.par.map
    } else if(class(mcqtl.consensus.map) == 'list') {
        .clu$.desc.map <- mcqtl.consensus.map
    } else if(class(mcqtl.consensus.map) == 'character' && length(mcqtl.consensus.map) == 1) {
        .clu$.desc.map <- .read.map(mcqtl.consensus.map)
    } else {
        stop('Descendants map MUST be a filename OR a named list of chromosome maps OR NULL, not ', typeof(mcqtl.consensus.map), '.')
    }
    if(is.matrix(marker.data) || is.data.frame(marker.data)) {
        .gen <- marker.data
    } else if(class(marker.data) == 'character' && length(marker.data) == 1) {
        .gen <- .read.gen(marker.data)
    } else {
        stop('Marker data MUST be a filename OR a matrix OR a data.frame, not ', typeof(marker.data), '.')
    }

    .clu$.simulation.Np <- ncol(.gen)

    if(is.matrix(kinship.matrix)) {
        .kinship <- kinship.matrix
    } else if(is.character(kinship.matrix)) {
        .kinship <- as.matrix(read.table(kinship.matrix))
    } else if(is.null(kinship.matrix)) {
        #.kinship <- diag(ncol(.gen))
        .kinship <- .compute.ais.kinship(.clu$.par.map, .gen)
    } else {
        stop('Kinship matrix MUST be a filename OR a matrix OR NULL, not ', typeof(kinship.matrix), '.')
    }

    .clu$.kinship <- .check.kinship(.kinship)

    # compute merged map
    .test.map <- .create.test.maps(.clu$.par.map, .clu$.desc.map)

    # Create or discard markers which have NOT been genotyped.
    all.markers <- unique(unlist(lapply(.test.map, function(m) m$markers)))
    full.na <- rownames(.gen)[apply(.gen, 1, function(r) all(is.na(r)))]
    .debug.out("Have", length(all.markers), "markers in map. Have", nrow(.gen) - length(full.na), "genotyped markers.\n")
    markers.without.geno <- as.character(all.markers[!all.markers %in% rownames(.gen)])
    if(!(length(markers.without.geno) == 0 && length(full.na) == 0)) {
        warning(length(markers.without.geno) + length(full.na), " markers don't have genotype data.")
        if(length(markers.without.geno) > 0 && !discard.unknown.markers) {
            nmark <- matrix(NA, ncol=ncol(.gen), nrow=length(markers.without.geno))
            rownames(nmark) <- markers.without.geno
            .gen <- rbind(.gen, nmark)
            warning(length(markers.without.geno), " markers were added to the genotype data (all observations are NA).")
        } else {
            markers.without.geno <- c(markers.without.geno, full.na)
            .clu$.test.map.full <- .test.map
            .test.map <- lapply(.test.map, .remove.markers, markers.without.geno)
            .clu$.markers.without.geno <- markers.without.geno
            .info.out(length(markers.without.geno), "markers were removed from the test map.")
        }
    } else {
        .info.out("All markers have genotype data.\n")
    }
    .clu$.gen <- .gen
    .clu$.test.map <- .test.map

    # Create "bound" methods
    environment(config) <- .clu
    environment(pair.similarity) <- .clu
    environment(get.marker.data) <- .clu
    environment(get.chromosome.names) <- .clu
    environment(get.kinship.matrix) <- .clu
    environment(get.scan.map) <- .clu
    environment(get.mcqtl.consensus.map) <- .clu
    environment(get.haplotypes.map) <- .clu
    environment(select.chromosome) <- .clu
    environment(get.selected.chromosome) <- .clu
    environment(run.equi.simulation) <- .clu
    environment(run.mosaic.simulation) <- .clu
    environment(run.raw.simulation) <- .clu
    environment(.generic.run.simulation) <- .clu
    environment(pairwise.similarities) <- .clu
    environment(transitive.closure) <- .clu
    environment(cluster.pair) <- .clu
    environment(get.window.density) <- .clu
    environment(write.xml) <- .clu
    environment(train) <- .clu
    environment(get.training.results) <- .clu
    environment(run) <- .clu
    environment(.run.chromosome) <- .clu
    environment(get.config) <- .clu
    environment(use.multicore) <- .clu
    environment(.as.index) <- .clu
    .clu$.as.index <- .as.index
    .clu$.run.chromosome <- .run.chromosome
    .clu$.generic.run.simulation <- .generic.run.simulation
    .clu$pair.similarity <- pair.similarity
    .clu$get.window.density <- get.window.density
    .clu$transitive.closure <- transitive.closure
    .clu$pairwise.similarities <- pairwise.similarities
    .clu$write.xml <- write.xml
    .clu$config <- config
    .clu$select.chromosome <- select.chromosome
    .clu$train <- train
    .clu$cluster.pair <- cluster.pair
    environment(lapply) <- .clu
    environment(sapply) <- .clu
    environment(apply) <- .clu
    .clu$lapply <- lapply
    .clu$sapply <- sapply
    .clu$apply <- apply


    # Create object-like list
    ret <- list(get.marker.data=get.marker.data,
                get.kinship.matrix=get.kinship.matrix,
                pair.similarity=pair.similarity,
                get.chromosome.names=get.chromosome.names,
                get.scan.map=get.scan.map,
                get.config=get.config,
                get.haplotypes.map=get.haplotypes.map,
                get.mcqtl.consensus.map=get.mcqtl.consensus.map,
                select.chromosome=select.chromosome,
                get.selected.chromosome=get.selected.chromosome,
                config=config,
                run.equi.simulation=run.equi.simulation,
                run.mosaic.simulation=run.mosaic.simulation,
                run.raw.simulation=run.raw.simulation,
                pairwise.similarities=pairwise.similarities,
                transitive.closure=transitive.closure,
                cluster.pair=cluster.pair,
                write.xml=write.xml,
                train=train,
                get.training.results=get.training.results,
                run=run,
                get.window.density=get.window.density,
                count.ancestral.alleles=count.ancestral.alleles,  # this one doesn't need the environment hack.
                use.multicore=use.multicore,
                .=.clu)
    ret$config(...)
    ret
}


