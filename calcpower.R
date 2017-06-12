devtools::install_github("chr1swallace/TibsPower")
library(TibsPower)
                                        #devtools::load_all("~/RP/TibsPower")
v <- readRDS("v.rds")

## do calculation - this bit is *long*
result <- tibs.power(v=v,
                     fold.diff=c(1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2),
                     sample.sizes=seq(20,200,by=10),
                     ngenes.diff=c(50,100,200), # shouldn't make much difference
                     p.thr=10^(seq(-10,-2,by=0.1)),
                     nsim=1000)
saveRDS(result,file="result.rds")

result <- readRDS("result.rds")

library(data.table)
x <- as.data.table(result)
x <- x[order(x$fdr.med),]

f <- function(x,y,x.target) {
    hi <- which(x>x.target)[1]
    lo <- hi-1
    d <- (x.target - x[lo])/(x[hi]-x[lo])
    (1-d) * y[lo] + d * y[hi]
}

## FDR=5%
y5 <- x[,.(fnr=f(fdr.med,fnr.med,0.05)),by=c("ngroup1","ngenes.diff","fold.diff")]
library(ggplot2)
library(cowplot)
library(randomFunctions)
library(magrittr)
y5$n <- factor(y5$ngroup1,levels=sort(unique(y5$ngroup1),decreasing=TRUE))
p <- ggplot(y5[ngenes.diff==100 & ngroup1 %% 10==0,],aes(x=fold.diff,y=1-fnr,col=n))+ geom_path() + ylab("Power") + xlab("Expression fold difference, R vs NR") + scale_colour_discrete("Samples\nper group") + ggtitle("RNA-seq differential expression") +
geom_vline(xintercept = 1.5,linetype="dashed",col="grey") + geom_hline(yintercept=0.8,linetype="dashed",col="grey")
randomFunctions:::nicecow(p)
save_plot("rnaseq-power.png", randomFunctions:::nicecow(p),base_height=4)

