.output.matrix <- function(prob) {
    paste('                <MATRIX_P number_rows="1" number_columns="1">\n',
          '                    <LINE number="1">\n',
          '                        <PROBABILITY value="', as.numeric(prob), '" />\n',
          '                    </LINE>\n',
          '                </MATRIX_P>\n', sep='')
}

write.xml <- function(closure) {
    if (!exists('.clusthaplo.tag')) {
        stop("Please use this function from inside a Clusthaplo context.")
    }

    if('clusthaplo.haplotypes' %in% class(closure)) {
        closure <- list(closure)
    }

    #loci <- .simi.data$loci
    #chrom.name <- .selected.chromosome
    parent.names <- colnames(.gen)
    prefix <- file.path(getwd(), .xml.output.dir)
    if(!file.exists(prefix)) {
        dir.create(prefix)
    }
    #prefix <- file.path(prefix, chrom.name)
    #if(!file.exists(prefix)) {
    #    dir.create(prefix)
    #}
    make.pos <- function(pa, i, pos, clos) {
        loci <- attr(clos, "loci")
        chrom.name <- attr(clos, "chromosome")
        scan.map <- .desc.map[[chrom.name]]
        mrk <- which(scan.map$locus == loci[pos])
        #.debug.out(scan.map$locus, "\n")
        #mrk <- which(abs(scan.map$locus - loci[pos]) < 1.e-5)
        #.debug.out("Markers at ", loci[pos], ":", paste(mrk, collapse=","), "\n")

        prob <- clos[pos, pa] == i

        mprob <- .output.matrix(prob)
        if (length(mrk) == 0) {
            list('            <POSITION value="', loci[pos], '" marker="no" name=" " weight="-1">\n',
                 mprob, '            </POSITION>\n')
        } else {
            lapply(mrk, function(m) {
                   list('            <POSITION value="', loci[pos],
                     #'" marker="yes" name="', as.character(.desc.map[[chrom.name]]$map$markers[m]),
                     '" marker="yes" name="', as.character(scan.map$markers[m]),
                     '" weight="1">\n', mprob, '            </POSITION>\n')
                 })
        }
    }

    for(pa in 1:ncol(.gen)) {
        #cat(file=paste(pa,'.',i,'.txt', sep=''))
        parent.name <- parent.names[pa]
        #filename <- paste(prefix, 'LD_file_', parent.name, '.xml', sep='')
        filename <- file.path(prefix, paste('LD_file_', parent.name, '.xml', sep=''))
        cat("Writing ", filename, "...\n", sep="")
        .cat <- function(...) cat(..., file=filename, sep="", append=T)
        # create empty file.
        cat(file=filename, append=F)
        .cat('<?xml version="1.0" encoding="ISO-8859-1" ?>\n',
             '<POPULATION name="', parent.name, '" type="parents" nbfondateur="', ncol(.gen), '">\n', collapse="")
        for(i in 1:ncol(.gen)) {
            .cat(unlist(list('    <ANCESTRAL_ALLELE name="F', i, '" number="', i, '">\n',
                             lapply(1:length(closure),
                                    function(ic) {
                                        clos <- closure[[ic]]
                                        list('        <CHROMOSOME name="', attr(clos, "chromosome"), '" number="', ic, '">\n',
                                             lapply(1:length(attr(clos, "loci")), function(pos) make.pos(pa, i, pos, clos)),
                                             '        </CHROMOSOME>\n')
                                    }),
                             '    </ANCESTRAL_ALLELE>\n')), collapse="")
            gc()
        }
        .cat('</POPULATION>\n')
#        .cat(unlist(list('<?xml version="1.0" encoding="ISO-8859-1" ?>\n',
#                         '<POPULATION name="', parent.name, '" type="parents" nbfondateur="', ncol(.gen), '">\n',
#                         lapply(1:ncol(.gen),
#                                function(i) {
#                                    #cat('\n', file=paste(pa,'.',i,'.txt', sep=''), append=T)
#                                    list('    <ANCESTRAL_ALLELE name="F', i, '" number="', i, '">\n',
#                                         lapply(1:length(closure),
#                                                function(ic) {
#                                                    clos <- closure[[ic]]
#                                                    list('        <CHROMOSOME name="', attr(clos, "chromosome"), '" number="', ic, '">\n',
#                                                         lapply(1:length(attr(clos, "loci")), function(pos) make.pos(pa, i, pos, clos)),
#                                                         '        </CHROMOSOME>\n')
#                                                }),
#                                         '    </ANCESTRAL_ALLELE>\n')
#                                }),
#                         '</POPULATION>\n')), collapse="")
    }
}
