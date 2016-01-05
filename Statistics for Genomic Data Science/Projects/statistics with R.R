# Q1:What are the coefficients for the SNP variable?------
#load library
library(snpStats)
library(broom)
#load data
data(for.exercise)
use <- seq(1, ncol(snps.10), 10)
sub.10 <- snps.10[,use]
snpdata = sub.10@.Data
status = subject.support$cc
#the 3rd snp data
snp3 = as.numeric(snpdata[,3])
#replace 0 with NA
snp3[snp3 == 0] = NA
#fit linear model
lm1 = lm(status ~ snp3)
#logistic regression model
glm1 = glm(status ~ snp3, family = "binomial")
#coefficients of liner model and logistic regression model
tidy(lm1)
tidy(glm1)

#Q3:Fit a logistic regression model on a recessive (need 2 copies of minor allele to confer risk)
#and additive scale for the 10th SNP. Make a table of the fitted values versus the
#case/control status. Does one model fit better than the other?-----
#load library
library(snpStats)
library(broom)
#load data
data(for.exercise)
use <- seq(1, ncol(snps.10), 10)
sub.10 <- snps.10[,use]
snpdata = sub.10@.Data
status = subject.support$cc
#get data of the 10th snp
snp10 = as.numeric(snpdata[,10])
#replace 0 with NA
snp10[snp10 == 0] = NA
#fit logistic regression model on a recessive scale
snp10_rec = (snp10==2)
glm_rec = glm(status ~ snp10_rec, family = "binomial")
#check fitted values
glm_rec$fitted.values
tidy(glm_rec)

#fit logistic regression model on a additive scale
glm2 = glm(status ~ snp10, family = "binomial")
#check fitted values
glm2$fitted.values
tidy(glm2)

#Q4:Fit an additive logistic regression model to each SNP. What is the average effect size? What is the max? What is the minimum?-------
#load library
library(snpStats)
library(broom)
#load data
data(for.exercise)
use <- seq(1, ncol(snps.10), 10)
sub.10 <- snps.10[,use]
snpdata = sub.10@.Data
status = subject.support$cc

results = rep(NA, dim(snpdata)[2]) # set up a dummy vector to save all values
for (i in 1:ncol(snpdata)){         # this line opens the loop and defines how often it will be executed
  snpdata_i = as.numeric(snpdata[,i])
  snpdata_i[snpdata_i == 0] = NA
  glm_i = glm(status ~ snpdata_i, family = "binomial")
  results[i] = tidy(glm_i)$statistic[2]####                 # In the final loop code, save the results into the dummy vector
  }                                   # close the loop
#average effect size
mean(results)
#min effect size
min(results)
#max effect size
max(results)

#Q5:What is the correlation with the results from using snp.rhs.tests and chi.squared ? Why does this make sense?----
library(snpStats)
library(broom)
data(for.exercise)
use <- seq(1, ncol(snps.10), 10)
sub.10 <- snps.10[,use]
snpdata = sub.10@.Data
status = subject.support$cc
#fit an additive logistic regression model to each SNP and square the coefficients
results_coeff = rep(NA, dim(snpdata)[2]) # set up a dummy vector to save all values
for (i in 1:ncol(snpdata)){         # this line opens the loop and defines how often it will be executed
  snpdata_i = as.numeric(snpdata[,i])
  snpdata_i[snpdata_i == 0] = NA
  glm_i = glm(status ~ snpdata_i, family = "binomial")
  results_coeff[i] = glm_i$coefficients # In the final loop code, save the results into the dummy vector
}                                   # close the loop
results_coeff_squre =  results_coeff^2

glm_all = snp.rhs.tests(status ~ 1, snp.data = sub.10)
slotNames(glm_all)
cor(results_coeff_squre, chi.squared(glm_all)) #the result is not match to the correct option

#Q6:Do the log2(data + 1) transform and fit calculate F-statistics for the difference between
#studies/populations using genefilter:rowFtests and using genefilter:rowttests. Do you get the same statistic? 
#Do you get the same p-value?-----
#load data
con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/montpick_eset.RData")
load(file=con)
close(con)
mp = montpick.eset
pdata=pData(mp)
edata=as.data.frame(exprs(mp))
fdata = fData(mp)
#transform data
library(genefilter)
edata = log2(as.matrix(edata) + 1)
#row ttest
tstats_obj = rowttests(edata, pdata$population)
hist(tstats_obj$p.value)
tidy(tstats_obj)
#row Ftest
fstats_obj = rowFtests(edata, as.factor(pdata$population))
hist(fstats_obj$p.value)
tidy(fstats_obj)

#Q7:What is the correlation in the statistics between the two analyses? 
#Are there more differences for the large statistics or the small statistics (hint: Make an MA-plot).----
library(Biobase)
con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/montpick_eset.RData")
load(file=con)
close(con)
mp = montpick.eset
pdata=pData(mp)
edata=as.data.frame(exprs(mp))
edata = edata[rowMeans(edata) > 100,]
fdata = fData(mp)
names(pdata)
library(DESeq2)
library(limma)
#fit with DESeq
de = DESeqDataSetFromMatrix(edata, pdata, ~study)
glm_all_nb = DESeq(de)
result_nb = results(glm_all_nb)
hist(result_nb$stat)

#fit with limma
edata = log2(as.matrix(edata) + 1)
mod = model.matrix(~ pdata$study)
fit_limma = lmFit(edata, mod)
ebayes_limma = eBayes(fit_limma)
limma_all = topTable(ebayes_limma, number = dim(edata)[1])
hist(limma_all$t, col = 4)
names(limma_all)
# correlation
length(result_nb$stat)
length(limma_all$t)
plotMA(result_nb, limma_all)

#Q8:How many results are statistically significant at an FDR of 0.05 in each analysis?----
#number with significant false discovery rate at 0.05 for DEseq
names(result_nb)
fp_bh = p.adjust(result_nb$pvalue, method = "BH")
hist(fp_bh, col = 3)
quantile(fp_bh)
sum(fp_bh < 0.05)

#number with significant false discovery rate at 0.05 for limma
fp_bh_limma = p.adjust(limma_all$P.Value, method = "BH")
hist(fp_bh_limma, col = 4)
quantile(fp_bh_limma)
sum(fp_bh_limma < 0.05)
