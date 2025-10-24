dyn.unload("src/hiC_expt.so")
dyn.load("src/hiC_expt.so")

pos <- read.table("pos.tsv")
pos.m <- as.matrix(pos) ## and it's integer!
pos.m <- pos.m[ order(pos.m[,1]), ]
pos.1 <- pos.m[1:10000, ]

assess.links <- function(pos, dims, min.d=1000, toLog10=TRUE){
    pos <- pos[ order(pos[,1]), ]
    tmp <- .Call("assess_region_links", pos, as.integer(dims))
    names(tmp) <- c("obs", "exp", "h", "v")
    if(min.d > 0){
        b <- tmp$v[-1] > log(min.d)
        tmp$obs <- tmp$obs[b,]
        tmp$exp <- tmp$exp[b,]
        i <- which(b)
        tmp$v <- tmp$v[ c(i[1], i+1) ]
    }
    if(toLog10)
        tmp$v <- tmp$v / log(10)
    tmp
}

tmp.1 <- assess.links( pos.m, c(201, 101, 20))

min.e <- with(tmp.1, min(exp[ exp > 0 ]))

## I need the offsets of the sequences here
chr.off <- read.table("chr_off.tsv", sep="\t")

x <- seq(1, max(pos[,2]), 100)
y <- do.call(cbind, apply(chr.off, 1, function(o){
    cbind( o[1] - x, o[2] - x )
}, simplify=FALSE))

par(mfrow=c(3,1))
with(tmp.1, image(h, v, t(obs)))
invisible(apply(y, 2, function(yy){ lines(x, log10(yy), col='blue') }))
with(tmp.1, image(h, v, t(exp)))
invisible(apply(y, 2, function(yy){ lines(x, log10(yy), col='blue') }))
with(tmp.1, image(h, v, t(scale(obs/(min.e+exp)))) )
invisible(apply(y, 2, function(yy){ lines(x, log10(yy), col='blue') }))


## with threading
## pars are : number of columns, number of rows, thinner, nthreads
assess.links.mt <- function(pos, pars, min.d=1000, toLog10=TRUE){
    pos <- pos[ order(pos[,1]), ]
    tmp <- .Call("assess_region_links_mt", pos, as.integer(pars))
    names(tmp) <- c("obs", "exp", "h", "v")
    if(min.d > 0){
        b <- tmp$v[-1] > log(min.d)
        tmp$obs <- tmp$obs[b,]
        tmp$exp <- tmp$exp[b,]
        i <- which(b)
        tmp$v <- tmp$v[ c(i[1], i+1) ]
    }
    if(toLog10)
        tmp$v <- tmp$v / log(10)
    tmp
}


dyn.unload("src/hiC_expt.so")
dyn.load("src/hiC_expt.so")

tmp.1 <- assess.links.mt( pos.m, c(201, 101, 20, 20))

x <- seq(1, max(pos[,2]), 100)
y <- do.call(cbind, apply(chr.off, 1, function(o){
    cbind( o[1] - x, o[2] - x )
}, simplify=FALSE))

min.e <- with(tmp.1, min(exp[ exp > 0 ]))

par(mfrow=c(3,1))
with(tmp.1, image(h, v, t(obs)))
invisible(apply(y, 2, function(yy){ lines(x, log10(yy), col='blue') }))
with(tmp.1, image(h, v, t(exp)))
invisible(apply(y, 2, function(yy){ lines(x, log10(yy), col='blue') }))
with(tmp.1, image(h, v, t(scale(obs/(min.e+exp)))) )
invisible(apply(y, 2, function(yy){ lines(x, log10(yy), col='blue') }))

## lets check the thread timing..
run.times <- vector(mode='list', length=32)
for(i in 1:length(run.times)){
    run.times[[i]] <- system.time( tmp.2 <- assess.links.mt( pos.m, c(201, 101, 20, i) ))
}

run.times <- do.call(rbind, run.times)

elapsed.range <- range(run.times[,'elapsed'])
par(mfrow=c(1,2))
plot(1:nrow(run.times), run.times[,'elapsed'],
     main=sprintf("1 to %d threads (%.2f -> %.2f seconds)",
                  nrow(run.times), elapsed.range[2], elapsed.range[1]))
## We have almost linear increase up to 25x;
plot(1:nrow(run.times), run.times[1,'elapsed'] / run.times[,'elapsed'],
     main=sprintf("max speedup: %.2f",
                  max(run.times[1,'elapsed'] / run.times[,'elapsed'])))

tmp.3 <- assess.links.mt( pos.m, c(201, 101, 20, 47))

## I now allow up to 96 threads; lets see how far it will scale;
## 128 threads should be slower
thread.n <- c(2^(0:5), 48, 49, 96)
run.times.2 <- lapply(thread.n, function(t){
    system.time( tmp <- assess.links.mt( pos.m, c(201, 101, 20, t) ))
})
run.times.2 <- do.call(rbind, run.times.2)

plot(thread.n, run.times.2[,'elapsed'],
     main=sprintf("%d to %d threads (%.2f -> %.2f seconds)",
                  min(thread.n), max(thread.n),
                  max(run.times.2[,'elapsed']), min(run.times.2[,'elapsed'])))
plot(thread.n, run.times.2[1,'elapsed'] / run.times.2[,'elapsed'],
     main=sprintf("max speedup: %.2f",
                  max(run.times.2[1,'elapsed'] / run.times.2[,'elapsed'])))

## we might expect that this should take about 200 seconds.
## lets see: (it took 180 second the first time)
system.time(tmp.4 <-  assess.links.mt( pos.m, c(201, 101, 1, 96) ) )

par(mfrow=c(3,1))
with(tmp.4, image(h, v, t(obs)))
invisible(apply(y, 2, function(yy){ lines(x, log10(yy), col='blue') }))
with(tmp.4, image(h, v, t(exp)))
invisible(apply(y, 2, function(yy){ lines(x, log10(yy), col='blue') }))
with(tmp.4, image(h, v, t(log(obs/(min.e+exp)))) )
invisible(apply(y, 2, function(yy){ lines(x, log10(yy), col='blue') }))

## look at the columns for position up to 6 million in more detail:
i <- which(tmp.4$h < 6e6)
r <- with(tmp.4, v[-1])
l <- with(tmp.4, v[2:length(v) - 1])

i <- i[ colSums(tmp.4$exp[,i]) > 400 ]

plot.new()
with(tmp.4, plot.window(ylim=c(0, max(colSums(obs[,i]))), xlim=range(v)))
with(tmp.4, plot.window(ylim=c(0, 1), xlim=range(v)))
axis(1)
axis(2)
for(j in i){
    y <- c(0, cumsum(tmp.4$obs[,j]))
    y <- y / max(y)
    link.pos <- tmp.4$h[j] + 10^((r + l)/2)
    link.chr <- cut(link.pos, breaks=unique(as.numeric(t(chr.off))))
##    with(tmp.4, lines( (l + r) / 2, cumsum(obs[,j]) ))
    segments(l, y[-length(y)], r, y[-1], col=link.chr, lwd=2)
    input <- readline("next: ")
}

plot.new()
with(tmp.4, plot.window(ylim=c(0, max(colSums(obs[,i]))), xlim=range(v)))
with(tmp.4, plot.window(ylim=c(0, 1), xlim=range(v)))
axis(1)
axis(2)
for(j in i){
    y <- c(0, cumsum(tmp.4$obs[,j] / (min.e + tmp.4$exp[,j])))
    y <- y / max(y)
    link.pos <- tmp.4$h[j] + 10^((r + l)/2)
    link.chr <- cut(link.pos, breaks=unique(as.numeric(t(chr.off))))
##    with(tmp.4, lines( (l + r) / 2, cumsum(obs[,j]) ))
    segments(l, y[-length(y)], r, y[-1], col=link.chr, lwd=2)
    input <- readline("next: ")
}
