# Biostrings library
library(Biostrings)
# one string
dna1 <- DNAString("ACGT-G")
dna1
# a set of string
dna2 <- DNAStringSet(c("ACG", "ACGT", "ACGTT"))
dna2
IUPAC_CODE_MAP
dna1[2:4]
# add column names to set
names(dna2) <- paste0("seq", 1:3)
dna2
# manuplication to set
width(dna2)
sort(dna2)
rev(dna2)
rev(dna1)
reverseComplement(dna2)
translate(dna2)
# frequency of letters
alphabetFrequency(dna2)
letterFrequency(dna2, letters = "GC")
dinucleotideFrequency(dna2)
consensusMatrix(dna2)
# library of genome
library(BSgenome)
available.genomes()
library("BSgenome.Scerevisiae.UCSC.sacCer2")
Scerevisiae
seqnames(Scerevisiae)
seqlengths(Scerevisiae)
Scerevisiae$chrI
letterFrequency(Scerevisiae$chrI,"GC")
letterFrequency(Scerevisiae$chrI,"GC", as.prob = TRUE)
param = new("BSParams", X = Scerevisiae, FUN = letterFrequency)
bsapply(param, "GC")
unlist(bsapply(param, "GC"))
# total GC content
sum(unlist(bsapply(param, "GC")))/sum(seqlengths(Scerevisiae))
# GC content of each chrosome
unlist(bsapply(param, "GC", as.prob = TRUE))

# Biostrings-Matching
dnaseq <- DNAString("ACGTACGT")
dnaseq
# match a pattern to a genome
matchPattern(dnaseq, Scerevisiae$chrI)
# count match
countPattern(dnaseq, Scerevisiae$chrI)
# match to the whole genome
vmatchPattern(dnaseq, Scerevisiae)
# find transcript binding
matchPWM()
# against a short sequence
pairwiseAlignment()
# allowing intel and mismatches
trimLRPatterns()

#BSgenome- views
vi <- matchPattern(dnaseq, Scerevisiae$chrI)
vi
ranges(vi)
Scerevisiae$chrI[57932:57939]
alphabetFrequency(vi)
shift(vi, 10)

gr = vmatchPattern(dnaseq, Scerevisiae)
vi2 = Views(Scerevisiae, gr)
vi2

library(AnnotationHub)
ah <- AnnotationHub()
qh <- query(ah, c("sacCer2", "genes"))
qh
genes <- ah[["AH7048"]]
# promoters
prom <- promoters(genes)
prom
# trim promoter to remove negative regions
prom <- trim(prom)
prom
# view promoter
promView <- Views(Scerevisiae, prom)
promView
# GC content of promoter
gcProm = letterFrequency(promView, "GC", as.prob = TRUE)
gcProm
# plot GC content
plot(density(gcProm))
abline(v = 0.38)
#load genomic rangs
library(GenomicRanges)
rl <- Rle(c(1,1,1,1,1,1,2,2,2,2,2,4,4,2))
rl
# length
runLength(rl)
# value
runValue(rl)
# retriview
as.numeric(rl)

ir <- IRanges(start = c(2,8), width = 4)
ir
# mean of two ranges in ir
aggregate(rl, ir, FUN = mean)
ir <- IRanges(start = 1:5, width = 3)
ir
coverage(ir)
slice(rl,2)

vi <- Views(rl, IRanges(2,8))
vi
vi <- Views(rl, IRanges(c(2,8), width = 2))
vi
mean(vi)

# for genomic ranges
gr <- GRanges(seqnames = "chr1", ranges = IRanges(start = 1:10, width = 3))
rl <-coverage(gr)
gr
rl
vi = Views(rl, GRanges("chr1", ranges = IRanges(3,7)))
vi
vi = Views(rl, as(GRanges("chr1", ranges = IRanges(3,7)), "RangesList"))
vi$chr1
#GenomicRanges - Lists
gr1 = GRanges(seqnames = "chr1", ranges = IRanges(start = 1:4, width = 3))
gr2 = GRanges(seqnames = "chr2", ranges = IRanges(start = 1:4, width = 3))
gL = GRangesList(gr1 = gr1, gr2 = gr2)
gL
start(gL)
elementLengths(gL)
shift(gL, 10)
findOverlaps(gL, gr2)

#Genomic Features
library(GenomicFeatures)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb = TxDb.Hsapiens.UCSC.hg19.knownGene
txdb
# find gene in a given range
gr = GRanges(seqnames = "chr1", strand = "+", ranges = IRanges(start = 11874, end = 14409))
g = genes(txdb)
subsetByOverlaps(g, gr, ignore.strand = TRUE)
# find transcripts in a given range
t = transcripts(txdb)
subsetByOverlaps(t, gr, ignore.strand = TRUE)
# find exons in a given range
e = exons(txdb)
subsetByOverlaps(e, gr, ignore.strand = TRUE)
#find exons and transcripts
subsetByOverlaps(exonsBy(txdb, by = "tx"), gr)
# find cds
subsetByOverlaps(cds(txdb), gr)
subsetByOverlaps(cdsBy(txdb, by = "tx"), gr)
# transcripts
subset(transcriptLengths(txdb, with.cds_len = TRUE), gene_id == "100287102")

#rtracklayer - Data Import
## ----dependencies, warning=FALSE, message=FALSE--------------------------
library(rtracklayer)
library(AnnotationHub)
library(Rsamtools)

## ----biocLite, eval=FALSE------------------------------------------------
## source("http://www.bioconductor.org/biocLite.R")
## biocLite(c("rtracklayer", "AnnotationHub", "Rsamtools"))

## ----help, eval=FALSE----------------------------------------------------
## ?import
## ?BigWigFile

## ----ahub----------------------------------------------------------------
library(AnnotationHub)
ahub <- AnnotationHub()
table(ahub$rdataclass)

## ----granges-------------------------------------------------------------
ahub.gr <- subset(ahub, rdataclass == "GRanges" & species == "Homo sapiens")
gr <- ahub.gr[[1]]
gr
seqinfo(gr)

## ----BigWig--------------------------------------------------------------
ahub.bw <- subset(ahub, rdataclass == "BigWigFile" & species == "Homo sapiens")
ahub.bw
bw <- ahub.bw[[1]]
bw

## ----importBigWig--------------------------------------------------------
gr1 <- gr[1:3]
out.gr <- import(bw, which = gr1)
out.gr

## ----importBigWig2-------------------------------------------------------
out.rle <- import(bw, which = gr1, as = "Rle")
out.rle

## ----importBigWig3-------------------------------------------------------
gr.chr22 <- GRanges(seqnames = "chr22",
                    ranges = IRanges(start = 1, end = seqlengths(gr)["chr22"]))
out.chr22 <- import(bw, which = gr.chr22, as = "Rle")
out.chr22[["chr22"]]

## ----liftOver------------------------------------------------------------
ahub.chain <- subset(ahub, rdataclass == "ChainFile" & species == "Homo sapiens")
query(ahub.chain, c("hg18", "hg19"))
chain <- ahub.chain[ahub.chain$title == "hg19ToHg18.over.chain.gz"]
chain <- chain[[1]]
gr.hg18 <- liftOver(gr, chain)
gr.hg18

## ----liftOver2-----------------------------------------------------------
table(elementLengths(gr.hg18))

## ----tabixIndex----------------------------------------------------------
library(Rsamtools)
from <- system.file("extdata", "ex1.sam", package="Rsamtools",
                    mustWork=TRUE)
from
to <- tempfile()
zipped <- bgzip(from, to)
idx <- indexTabix(zipped, "sam")

## ----sessionInfo, echo=FALSE---------------------------------------------
sessionInfo()
