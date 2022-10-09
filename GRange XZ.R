library(IRanges)
library(GenomicRanges)
library('rtracklayer')
library('AnnotationHub')
library('dbplyr')
library('BiocFileCache')

## 1.Using rtracklayer, BiomaRt, or annotationHub, import into R the human transcripts for GRCh38.
# library("TxDb.Hsapiens.UCSC.hg38.knownGene")
# txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
# txdb
# transcripts(txdb)
# exons(txdb)
# # Keep only the standard chromosomes.
# tx1 <-transcripts(txdb)
# tx1 <- keepStandardChromosomes(tx1, pruning.mode="coarse")


ahub <-AnnotationHub()
ahub
table(ahub$dataprovider)
query(ahub,c("TxDb", "Homo sapiens", "GRCh38"))
ahub['AH75191']
ahub
gtf <-ahub[['AH75191']]
gtf
columns(gtf)
transcript <- transcripts(gtf)
transcript

## 2.Create a Granges object for the promoters of all protein-coding transcripts, defined for this problem set as 1500 bp upstream of the TSS and 500 bp downstream of the TSS.
pr <- promoters(gtf,upstream = 1500, downstream = 500)
pr

## 3.Using the same methods as in 1, create a Granges object of all CpG islands for the human genome.
# The is the supposed way, but due to the version problem, we cannot find CpG islands for hg38, but can find for hg19 and other versions)
cpg <- query(ahub, "CpG","Homo sapiens","hg38")
cpg1 <- ahub[["AH5086"]]
class(cpg1)

# Therefore, I chose to work with this way to get CpG islands data for hg38.
library(rtracklayer)
session <- browserSession()
genome(session) <- "hg38"
cpg <- session[["CpG Islands"]]
cpg
cpg_record = makeGRangesFromDataFrame(cpg)

## 4.Calculate the fraction of CpG island annotations that overlap a promoter.
cpg_pr_data <-findOverlaps(cpg_record,pr)
fraction <- length(unique(queryHits(cpg_pr_data)))/length(pr)
fraction

cpg_pr_data <-findOverlaps(cpg_record,pr)
fraction <- length(queryHits(cpg_pr_data))/length(pr)
fraction

## 5.Plot the length distributions for CpG islands that do and do not overlap a promoter.
plotRanges <- function(x, xlim=x, main=deparse(substitute(x)),
                       col="black", sep=0.5, ...)
{
  height <- 1
  if (is(xlim, "IntegerRanges"))
    xlim <- c(min(start(xlim)), max(end(xlim))) 
  bins <- disjointBins(IRanges(start(x), end(x) + 1)) 
  plot.new() 
  plot.window(xlim, c(0, max(bins)*(height + sep))) 
  ybottom <- bins * (sep + height) - height 
  rect(start(x)-0.5, ybottom, end(x)+0.5, ybottom + height, col=col, ...) 
  title(main) 
  axis(1)
}

intersect <-intersect(ranges(cpg_record),ranges(pr))
par(mfrow=c(2,1))
overlap <-plotRanges(intersect)
not_overlap <-plotRanges(setdiff(ranges(cpg_record),intersect))

# histogram of the distances to nearest promoter
n.ind=nearest(cpg_record,pr)
dists=distanceToNearest(cpg_record, pr,select="arbitrary")
dists
dist2plot=mcols(dists)[,1]
hist(log10(dist2plot),xlab="log10(dist to nearest promoter)",
     main="Distances")


