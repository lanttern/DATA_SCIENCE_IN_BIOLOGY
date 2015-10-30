# Identify DNA polymorphic sites with Galaxy

[Published Page]( https://usegalaxy.org/u/coursera/p/genomic-data-science-with-galaxyidentify-polymorphic-sites)

This page includes description about analyses of DNA polymorphic sites of father-mother-child sequencing samples:

step 1: load data - the data are loaded from local files, set "fastqsanger" format and "hg19" database on the starting page

step 2: check quality of all sequencing files - use FastQC tool (version: 0.63) to check quality of the sequencing

step 3: mapping - use BWA-MEM tool (version: 0.1) to map sequence to reference genome (choose hg19 as reference), paired end

step 4: add or replace read groups - label each group (the mapping file) using AddOrReplaceReadGroup (version: 1.126.0)

step 5: merge 3 individual mapping files - use MergeSamFiles (version: 1.126.0)

step 6: filter - using filter tools: Filter (version: 1.126.0, remove low quality mapping), MarkDuplicates (version: 1.126.0, filter out duplicated mapping), CleanSam (version: 1.126.0)

step 7: identify polymorphic sites - using FreeBayes tool (version: 0.4) to identify polymorphic sites base on hg19 genome

step 8: filter out false positive sites - using VCFfilter (version: 0.0.3) to select sites where the chance of a false positive call is 1 in 10,000 or better.

step 9: extract workflow and download final vcf file for further analyses.

Stage 2 - analyze data of polymorphic sites based on vcf file

step 10: load data - set format as "vcf", genomic database as hg19

step 11: identify number of snp, mnp, del, ins or complex - using VCFfilter tool (version:0.0.3 ) to select different types of polymorphism (for example: -f "TYPE = snp", select snp only), then using Filter tool (version: 1.1.0) to find duplicated polymorphisms

step 12: identify genes with polymorphic sites - using ANNOVAR Annotate VCF tool (version: 0.1) to annotate the  vcf file in step 10

step 13: count polymorphic sites for each gene - using Group tool (version: 2.1.0, by gene name) to count number of polymorphic sites for each gene

step 14: sort results in step 13 using Sort tool (version: 1.0.3, by descending).

[WorkFlow](https://usegalaxy.org/workflow/display_by_username_and_slug?username=coursera&slug=workflow-constructed-from-history-genomic-data-science-with-galaxy-project---completed)

[History for Analyses](https://usegalaxy.org/u/coursera/h/workflow-constructed-from-history-genomic-data-science-with-galaxy-project---completed)
