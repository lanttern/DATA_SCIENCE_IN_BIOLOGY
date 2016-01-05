#load library----------------------------
library(devtools)
library(Biobase)
library(goseq)
library(DESeq2)
library(MatrixEQTL)
library(broom)
library(limma)
library(org.Mm.eg.db)

#Q1:which of the Mouse genome builds they aligned the reads to?-------------------------------------
supportedGenomes() 
# compare publication date and mouse genome release date
# mm10  Mouse  Dec. 2011    Genome Reference Consortium GRCm38
# mm9   Mouse  Jul. 2007    NCBI Build 37
# mm8   Mouse  Feb. 2006    NCBI Build 36
#mm7    Mouse  Aug. 2005    NCBI Build 35

#Q2:How many genes are differentially expressed at the 5% FDR level using Benjamini-Hochberg correction?--------- 
#What is the gene identifier of the first gene differentially expressed at this level-----------------------------
#load data
con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bottomly_eset.RData")
load(file=con)
close(con)
bot = bottomly.eset
pdata_bot=pData(bot)
fdata_bot = featureData(bot)
edata = exprs(bot)
#remove low expression gene
fdata_bot = fdata_bot[rowMeans(edata) > 5]
edata = edata[rowMeans(edata) > 5, ]
#transform data
edata = log(edata + 1)
#fit limma model
mod = model.matrix(~ pdata_bot$strain)
fit_limma = lmFit(edata, mod)
ebayes_limma = eBayes(fit_limma)
#adjust p value
limma_output = topTable(ebayes_limma, number = dim(edata)[1], adjust.method="BH", sort="none")
names(limma_output)
limma_pvals_adj = limma_output$adj.P.Val
limma_pvals_adj[1:10]
hist(limma_pvals_adj, col = 2)
#count p value less than 0.05
sum(limma_pvals_adj < 0.05)
#find the first gene showed significant differential expression
rownames(edata)[34]#the 34th

#Q3-Q4:What is the top category that comes up as over represented?What is the name of the category?--------------
#get differentiatial expressed gene list
genes = as.integer(limma_pvals_adj < 0.05)
names(genes) = rownames(edata)
not_na = !is.na(genes)
genes = genes[not_na]
head(genes)
sum(genes)
#GO cluster analysis
pwf=nullp(genes,"mm9","ensGene")
head(pwf)
GO.wall=goseq(pwf,"mm9","ensGene")
head(GO.wall)
GO.top10 = GO.wall[1:10,1]

#Q5:How many of the top 10 overrepresented categories are the same for the adjusted and unadjusted analysis?-------------
#fit limma model with lane correction
mod_adj = model.matrix(~ pdata_bot$strain + pdata_bot$lane.number)
fit_limma_adj = lmFit(edata, mod_adj)
ebayes_limma_adj = eBayes(fit_limma_adj)
limma_output_adj = topTable(ebayes_limma_adj, number = dim(edata)[1], adjust.method="BH",sort="none")
limma_pvals_adj_adj = limma_output_adj$adj.P.Val
hist(limma_pvals_adj_adj, col = 2)
sum(limma_pvals_adj_adj < 0.05)
#get differential expressed gene list
genes_adj = as.integer(limma_pvals_adj_adj < 0.05)
not_na = !is.na(genes_adj)
names(genes_adj) = rownames(edata)
genes_adj = genes_adj[not_na]
head(genes_adj)
#find common GO category between non-correction and corrected with lanes
pwf_adj=nullp(genes_adj,"mm9","ensGene")
head(pwf)
GO.wall_adj=goseq(pwf_adj,"mm9","ensGene")
GO.top10_adj = GO.wall_adj[1:10,1]
length(intersect(GO.top10, GO.top10_adj))

