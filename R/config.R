.check <- function(needle, haystack, pretty.name) {
    if(!needle %in% haystack) {
        stop(pretty.name, " must be one of ", paste(haystack, collapse=", "))
    }
}


config <- function(w1=NULL, w2=NULL, step.size=NULL, window.length=NULL,
                   na.replace=NULL,
                   scoring.method=NULL,
                   kinship.threshold=NULL,
                   clustering.method=NULL,
                   clustering.threshold=NULL,
                   simulation.type=NULL,
                   simulation.Np=NULL,
                   simulation.Ng=NULL,
                   simulation.Nrep=NULL,
                   threshold.quantile=NULL,
                   xml.output.dir=NULL,
                   plot.output.dir=NULL,
                   plot.height=NULL,
                   plot.width=NULL,
                   config=NULL) {

    if (!exists('.clusthaplo.tag')) {
        stop("Please use this function from inside a Clusthaplo context.")
    }

    if(!is.null(config)) {
        if(is.character(config)) {
            config <- read.table(config)
        } else if(!is.data.frame(config)) {
            stop("The 'config' parameter MUST be a filename OR a data.frame.")
        }
        v <- as.character(config$Current.Value)
        # 1             2               3           4             5             6             7             8
        # ""            "kernel.exp"  "kernel.unif" "1"           "20"          "TRUE"        "kinship"     "10"
        # 9             10            11
        # ""            "hmm"         NA
        # 12            13            14            15            16            17
        # ""            "equi"        "3"           "12"          "50"          "95%"
        # 18            19            20            21            22
        # ""            "XML"         "img"         "1024"        "768"        

        if(is.null(w1)) w1 <- v[2]
        if(is.null(w2)) w2 <- v[3]
        if(is.null(step.size)) step.size <- as.numeric(v[4])
        if(is.null(window.length)) window.length <- as.numeric(v[5])
        if(is.null(na.replace)) na.replace <- list(`TRUE`=TRUE, `FALSE`=FALSE, `NA`=NA)[[v[6]]]
        if(is.null(scoring.method)) scoring.method <- v[7]
        if(is.null(kinship.threshold)) kinship.threshold <- v[8]

        if(is.null(clustering.method)) clustering.method <- v[10]
        if(is.null(clustering.threshold)) clustering.threshold <- as.numeric(v[11])

        if(is.null(simulation.type)) simulation.type <- v[13]
        if(is.null(simulation.Nrep)) simulation.Nrep <- as.integer(v[14])
        if(is.null(simulation.Np)) simulation.Np <- as.integer(v[15])
        if(is.null(simulation.Ng)) simulation.Ng <- as.numeric(v[16])
        if(is.null(threshold.quantile)) threshold.quantile <- v[17]

        if(is.null(xml.output.dir)) xml.output.dir <- v[19]
        if(is.null(plot.output.dir)) plot.output.dir <- v[20]
        if(is.null(plot.width)) plot.width <- as.integer(v[21])
        if(is.null(plot.height)) plot.height <- as.integer(v[22])

        return(config(w1, w2, step.size, window.length,
                      na.replace,
                      scoring.method,
                      kinship.threshold,
                      clustering.method,
                      clustering.threshold,
                      simulation.type,
                      simulation.Np,
                      simulation.Ng,
                      simulation.Nrep,
                      threshold.quantile,
                      xml.output.dir,
                      plot.output.dir,
                      plot.height,
                      plot.width))
    }

    kernel.funcs <- names(.all.kernels)
    scoring.methods <- c('raw', 'normalized', 'kinship')
    clustering.methods <- c('threshold', 'hmm')

    if(!is.null(clustering.method)) {
        .check(clustering.method, clustering.methods, 'Clustering method')
        if(clustering.method == 'hmm') {
            # Update .RHmm.available just in case the user did comply and install RHmm only after a first failed attempt.
            .RHmm.available <<- as.logical(require('RHmm', quietly=T))

            #tryCatch({require(RHmm); .RHmm.available <<- T}, error=function(e) .RHmm.available <<- F)
            if(!.RHmm.available) {
                stop("Sorry, 'hmm' method requires the package RHmm to be installed. Try install.packages(\"RHmm\")")
            }
        }
    }
    if(!is.null(w1)) .check(w1, kernel.funcs, 'Weight kernel w1')
    if(!is.null(w2)) .check(w2, kernel.funcs, 'Weight kernel w2')
    if(!is.null(scoring.method)) .check(scoring.method, scoring.methods, 'Scoring method')

    if(!is.null(threshold.quantile)) {
        stop.tq <-"threshold.quantile MUST be a character string in the form \"X%\", X being an integer between 80 and 100, or a plain integer between 80 and 100."
        if(length(threshold.quantile) != 1) {
            threshold.quantile <- NA
        }
        if(is.character(threshold.quantile)) {
            if(substring(threshold.quantile, nchar(threshold.quantile)) != "%") {
                stop(stop.tq)
            }
            x <- as.integer(substring(threshold.quantile, 1, nchar(threshold.quantile) - 1))
            if(is.na(x) || x < 80 || x > 100) {
                stop(stop.tq)
            }
        } else if(is.numeric(threshold.quantile)) {
            threshold.quantile <- as.integer(threshold.quantile)
            if(threshold.quantile < 80 || threshold.quantile > 100) {
                stop(stop.tq)
            }
            threshold.quantile <- paste(threshold.quantile, "%", sep="")
        } else {
            threshold.quantile <- NA
        }
        if(is.na(threshold.quantile)) {
            stop(stop.tq)
        }
    }

    if(is.null(unlist(lapply(names(formals()), function(i) get(i))))) {
        return(data.frame(Current.Value=c(..Similarity..="",
                                          w1=.w1, w2=.w2, step.size=.step.size,
                                          window.length=.window.length, na.replace=.na.replace,
                                          scoring.method=.scoring.method,
                                          kinship.threshold=.kinship.threshold,
                                          ..Clustering..="",
                                          clustering.method=.clustering.method,
                                          clustering.threshold=ifelse(is.null(.clustering.threshold),
                                                                      NA, .clustering.threshold),
                                          ..Training..="",
                                          simulation.type=.simulation.type,
                                          simulation.Nrep=.simulation.Nrep,
                                          simulation.Np=.simulation.Np,
                                          simulation.Ng=.simulation.Ng,
                                          threshold.quantile=.threshold.quantile,
                                          ..Output..="",
                                          xml.output.dir=.xml.output.dir,
                                          plot.output.dir=.plot.output.dir,
                                          plot.width=.plot.width,
                                          plot.height=.plot.height)))
    }

    reconf.chrom <- F
    reconf.thres <- F

    upd.rc <- function(a, b) {
        eq <- is.na(a) && is.na(b) || !is.na(a) && !is.na(b) && a==b
        return(reconf.chrom || !eq)
    }

    upd.rt <- function(a, b) {
        eq <- is.na(a) && is.na(b) || !is.na(a) && !is.na(b) && a==b
        return(reconf.thres || !eq)
    }

    if(!is.null(w1))                   { reconf.chrom <- upd.rc(.w1, w1);                       .w1 <<- w1 }
    if(!is.null(w2))                   { reconf.chrom <- upd.rc(.w2, w2);                       .w2 <<- w2 }
    if(!is.null(window.length))        { reconf.chrom <- upd.rc(.window.length, window.length); .window.length <<- window.length }
    if(!is.null(step.size))            { reconf.chrom <- upd.rc(.step.size, step.size);         .step.size <<- step.size }

    if(!is.null(kinship.threshold))    { reconf.chrom <- upd.rc(.kinship.threshold, kinship.threshold)
                                         .kinship.threshold <<- kinship.threshold
                                       }

    .w1.fun <<- get(paste('.', .w1, sep=''))(.window.length)
    .w2.fun <<- get(paste('.', .w2, sep=''))(.window.length)

    if(!is.null(na.replace))           { reconf.thres <- upd.rt(.na.replace, na.replace);           .na.replace <<- na.replace }
    if(!is.null(simulation.type))      { reconf.thres <- upd.rt(.simulation.type, simulation.type); .simulation.type <<- simulation.type }
    if(!is.null(scoring.method))       { reconf.thres <- upd.rt(.scoring.method, scoring.method);   .scoring.method <<- scoring.method }
    if(!is.null(simulation.Np))        { reconf.thres <- upd.rt(.simulation.Np, simulation.Np);     .simulation.Np <<- simulation.Np }
    if(!is.null(simulation.Nrep))      { reconf.thres <- upd.rt(.simulation.Nrep, simulation.Nrep); .simulation.Nrep <<- simulation.Nrep }
    if(!is.null(threshold.quantile))   {
        if(.clustering.method == 'threshold') {
            reconf.thres <- upd.rt(.threshold.quantile, threshold.quantile)
        }
        .threshold.quantile <<- threshold.quantile
    }
    if(!is.null(simulation.Ng))        {
        if(.simulation.type == 'mosaic') {
            reconf.thres <- upd.rt(.simulation.Ng, simulation.Ng)
        }
        .simulation.Ng <<- simulation.Ng
    }

    if(!is.null(xml.output.dir))       { .xml.output.dir <<- xml.output.dir }
    if(!is.null(plot.output.dir))       { .plot.output.dir <<- plot.output.dir }
    if(!is.null(plot.width))            { .plot.width <<- plot.width }
    if(!is.null(plot.height))           { .plot.height <<- plot.height }
    if(!is.null(clustering.threshold)) {
        .clustering.threshold <<- clustering.threshold
        if(.clustering.method == 'threshold') {
            reconf.thres <- F
            .is.trained <<- T
            .training.results <<- "manual setting"
            .clusterize <<- function(similarity) similarity >= .clustering.threshold
        }
    }


    if(!is.null(clustering.method)) {
        if(.clustering.method != clustering.method) {
            reconf.thres <- T
            .clustering.method <<- clustering.method
        }
    }

    if(!.is.trained || .training.results != "manual setting") {
        if(reconf.chrom || reconf.thres) {
            if(!(is.null(.clusterize))) {
                warning("Configuration modified. Please train Clusthaplo again.")
            }
            .clustering.threshold <<- NULL
            .clusterize <<- NULL
        }
    }

    if(reconf.chrom) {
        select.chromosome(.selected.chromosome)
    }
}


.legal.names <- names(formals(config))
.legal.names <- .legal.names[.legal.names != 'config']
