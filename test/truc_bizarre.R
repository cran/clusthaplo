library(clusthaplo)
library(multicore)
eur <- clusthaplo('misc/geniteurs_euralis.map',
                  'misc/descendant_euralis.map',
                  'misc/geniteurs_euralis_D2007.pedgen',
                  window.length=23)

eur$select.chromosome('chr01')

ploum <- function(w, sco="kinship") {
    eur$config(w1=w, w2=w, scoring.method=sco)
    combn(ncol(eur$get.marker.data()), 2, function(r) eur$pair.similarity(r[1], r[2])$similarity)
}

simi.unif <- ploum('kernel.unif')
simi.const <- ploum('kernel.const')

compare <- function(a, b) {
    which(a != b, arr.ind=T)
}

do <- function(w) {
    eur <- clusthaplo('misc/geniteurs_euralis.map',
                      'misc/descendant_euralis.map',
                      'misc/geniteurs_euralis_D2007.pedgen')
    eur$conf(w1=w, w2=w, clustering.method="hmm")
    eur$select.chromosome("chr01")
    eur$train()
    conf=eur$config()
    tc <- eur$transitive.closure(eur$pairwise.similarities())
    nombre.de.breaks <- sum(apply(diff(tc) != 0, 1, any))  # compte le nombre de lignes où au moins un individu change de clique
    nombre.allele= eur$count.ancestral.alleles(tc)
    list(config=conf,nbr.breaks=nombre.de.breaks,nbr.allele.moyen=mean(nombre.allele$count), tc=tc)
}

print(system.time(k <- mclapply(c('kernel.unif', 'kernel.const'), do)))

t1 <- k[[1]]$tc
t2 <- k[[2]]$tc
diff.per.row <- sapply(1:nrow(t1), function(i) sum(t1[i, ] != t2[i, ]))


#comp.thres <- function(thres) {
#    sum(simi.unif >= thres & simi.const < thres | simi.const >= thres & simi.unif < thres)
#}
#
#tvec <- (1:1000)/1000
#tdif <- tdif <- sapply(tvec, comp.thres)
#hist((simi.unif + simi.const) / 2, breaks=c(0, tvec-.0005, 1), ylim=c(0, max(tdif)))
#lines(tvec, tdif, col="red")
