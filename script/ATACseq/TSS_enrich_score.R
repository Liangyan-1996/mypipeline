library(ATACseqQC)
library(GenomicRanges)

args = commandArgs(T)
bamfile = args[1]
UCSC2Txdb = args[2]
bam_obj <- readBamFile(bamFile=bamfile, tag=character(0))

library(UCSC2Txdb, character.only = TRUE)

x = paste("txdb <-",UCSC2Txdb,sep=" ")

eval(parse(text = x))

txs <- transcripts(txdb)
tsse <- TSSEscore(bam_obj, txs, seqlev = intersect(seqlevels(bam_obj), seqlevels(txs)), upstream = 2000, downstream = 2000, endSize = 100, width = 100, step = 100)
print(tsse)
