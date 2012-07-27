source('../R/clusthaplo.R', chdir=T)
source('../R/debug.R', chdir=T)
source('../R/map.R', chdir=T)
source('../R/similarity.R', chdir=T)
source('../R/config.R', chdir=T)
source('../R/kernel.R', chdir=T)
source('../R/writexml.R', chdir=T)

done <- 0
good <- 0

do.test <- function(simi, file) {
    write.table(simi, file=file)
    test <- read.table(paste('REF.', file, sep=""))
    if(!all(simi == test, na.rm=T)) {
        warning("Test failed: ", file)
    } else {
        good <<- good + 1
    }
    done <<- done + 1
}



clu <- clusthaplo('map1.txt', 'map1.txt', 'test1.txt')
clu$config(w1='kernel.const', w2='kernel.const', na.replace=F, step.size=2, scoring.method='raw')
clu$select.chromosome('Chr1')
do.test(as.data.frame(clu$pair.similarity(1, 2)), 'signal.step2.NA-FALSE.Chr1.txt')
clu$select.chromosome('Chr2')
do.test(as.data.frame(clu$pair.similarity(1, 2)), 'signal.step2.NA-FALSE.Chr2.txt')
clu$config(na.replace=NA)
clu$select.chromosome('Chr1')
do.test(as.data.frame(clu$pair.similarity(1, 2)), 'signal.step2.NA-IGNORED.Chr1.txt')
clu$select.chromosome('Chr2')
do.test(as.data.frame(clu$pair.similarity(1, 2)), 'signal.step2.NA-IGNORED.Chr2.txt')
#sf <- clu$make.similarity.func(w1=kernel.const(1), w2=kernel.const(1), step.size=2, na.replace=T)

#do.test(sf('Chr1', 1, 2), file="signal.step2.NA-FALSE.Chr1.txt")
#do.test(sf('Chr2', 1, 2), file="signal.step2.NA-FALSE.Chr2.txt")

#sf <- clu$make.similarity.func(step.size=2, na.replace=T)
#do.test(sf('Chr1', 1, 2), file="signal.step2.NA-IGNORED.Chr1.txt")
#do.test(sf('Chr2', 1, 2), file="signal.step2.NA-IGNORED.Chr2.txt")

#do.test(clu$get.I('Chr1', 1, 2, na.replace=NA), file="I.1.2.Chr1.txt")
#do.test(clu$get.I('Chr2', 1, 2, na.replace=NA), file="I.1.2.Chr2.txt")

cat(good, "/", done, " tests passed\n", sep="")

quit(status=good-done)
