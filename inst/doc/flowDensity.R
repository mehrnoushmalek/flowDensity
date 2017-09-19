

##Example1
library(flowCore)
library(flowDensity)
library(GEOmap)
data_dir <- system.file("extdata", package = "flowDensity")
load(list.files(pattern = 'sampleFCS_1', data_dir, full = TRUE))
f
sngl <- flowDensity(f,channels = c("FSC-A","FSC-H"),position = c(F,F),
                    percentile =c(.99999,.99999),use.percentile = c(T,T),

plotDens(f,c(1,2))
lines(sngl@filter,type="l")
###Example2
<<>>=
bcell <- flowDensity(sngl, channels=c(9, 3),
                     position=c(FALSE, NA))



plot(getflowFrame(sngl), bcell)

###Example3
CD19pCD20n <- flowDensity(obj=bcell, channels=c(8, 6),
                        position=c(T,F))
plasmablasts <- flowDensity(obj=CD19pCD20n, channels=c(5, 12),
                            position=c(T, T))

plotDens(getflowFrame(CD19pCD20n), plasmablasts@channels, pch=19)
points(plasmablasts@filter, type='l', col=2, lwd=2)

##Removing margin
f <- nmRemove(f, c("FSC-A", "SSC-A"))

##Multiple step gating
load(list.files(pattern = 'sampleFCS_2', data_dir, full = TRUE))
f2
channels <- c("V500-A", "SSC-A")
# First call to flowDensity
tmp.cp1 <- flowDensity(obj=f2, channels=channels,
                      position=c(TRUE, FALSE), percentile=c(0.25, NA))
# Second call to flowDensity
tmp.cp2 <- flowDensity(obj=tmp.cp1, channels=channels,
                       position=c(TRUE, FALSE), gates=c(FALSE, NA), 
                       percentile=c(NA, 0.85))
# Final call to flowDensity
lymph <- flowDensity(obj=f2, channels=channels,
                     position=c(TRUE, FALSE), gates=c(tmp.cp1@gates[1], 
                     tmp.cp2@gates[2]), ellip.gate=TRUE, scale=.99)

plot(f2, tmp.cp1)

plot(f2, tmp.cp2)

plotDens(f2, channels=channels)
points(lymph@filter, type="l", col=2, lwd=2)

##Using FMOs
load(list.files(pattern = 'sampleFCS_3.Rdata', data_dir, full = TRUE))
f3
load(list.files(pattern = 'sampleFCS_3_FMO', data_dir, full = TRUE))
f3.fmo
f3.gated <- flowDensity(obj=f3, channels=c('BV421-A', 'FSC-A'),
                        position = c(TRUE, NA),use.control = c(TRUE, F)
                        , control = c(f3.fmo, NA))
f3.fmo.gated <- flowDensity(obj=f3.fmo, channels=c('BV421-A', 'FSC-A'),
                            position=c(TRUE, NA),
                            gates=c(f3.gated@gates[1], NA))

plot(f3.fmo, f3.fmo.gated)

plot(f3, f3.gated)
##Using percentile

f3.gated.98p <- flowDensity(obj=f3, channels=c('BV421-A', 'FSC-A'),
                            position = c(TRUE, NA),use.percentile = c(TRUE, NA),
                            percentile = 0.98, use.control = c(TRUE, FALSE),
                            control = c(f3.fmo, NA))
f3.fmo.gated.98p <- flowDensity(obj=f3.fmo, channels=c('BV421-A', 'FSC-A'),
                                position = c(TRUE, NA),
                                gates=c(f3.gated.98p@gates[1], NA))

plot(f3.fmo, f3.fmo.gated.98p)

plot(f3, f3.gated.98p)
#Using other arguments

load(list.files(pattern = 'sampleFCS_2', data_dir, full = TRUE))
thresholds <- deGate(obj = f2,channel = 9)
#Percentile default is .95, which can be changed
thresholds.prcnt <- deGate(f2,channel = 9,use.percentile=T,percentile=.3) 
thresholds.lo <- deGate(f2,channel = 9,use.upper=T,upper=F,alpha = .9)
thresholds.hi <- deGate(f2,channel = 9,use.upper=T,upper=T,alpha = .9)

plotDens(f2,c(9,12))
abline(v=c(thresholds,thresholds.prcnt,thresholds.lo,thresholds.hi),col=c(1,2,3,4))

