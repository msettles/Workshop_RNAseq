
# 1. Load the edgeR package and use the utility function
library("edgeR")

# 2. Set up the experiment
samples <- read.table("samples.txt",header=T,as.is=T)
# Set up how you'll print the sample_id and condition variables, based on the samples.txt sheet
sample_id = apply(samples[,c("SAMPLE_ID")],1,function(x) paste(na.exclude(x),collapse=".")) 
condition = apply(samples[,c("treatment")],1,function(x) paste(na.exclude(x),collapse=".")) 

folder = "04-HTseqCounts"
results = "05-edgeR"
dir.create(results)
annof <- ""

# 3. Identify the count files and read them into R using readDGE
countf <- sapply(file.path(folder,samples$SAMPLE_ID),dir,pattern=".counts$",full.names=T)
counts = readDGE(countf)$counts

#### HTSEQ_COUNT RESULTS
# 4. Filter weakly expressed and noninformative (e.g., non-aligned) features:
noint = rownames(counts) %in% c("__no_feature","__ambiguous","__too_low_aQual",
                                "__not_aligned","__alignment_not_unique")

mean(colSums(counts[!noint,])/colSums(counts))
## MEAN % of reads map to features

cpms = cpm(counts)  ## counts per million
# In edgeR, it is recommended to remove features without 
# at least 1 read per million in n of the samples, 
# where n is the size of the smallest group of replicates, 
keep = rowSums(cpms >1) >=3 & !noint
dim(counts) ## the number of features you started with
counts = counts[keep,]
dim(counts) ## count counts of features you have left over after initial filter
colnames(counts) = samples$SAMPLE_ID

# 5. Create a DGEList object (edgeR's container for RNA-seq count data):
d = DGEList(counts=counts, group=condition)
# 6. Estimate normalization factors using:
d = calcNormFactors(d)
# 7. Inspect the relationships between samples using a multidimensional scaling (MDS) plot, as shown in Figure 4:
pdf(file.path(results,"MDS-edgeR.pdf"))
plotMDS(d, labels=sample_id,
        col = rainbow(length(levels(factor(condition))))[factor(condition)],cex=0.6, main="MDS")
dev.off()

# 8. Estimate tagwise dispersion (simple design) using:
d = estimateCommonDisp(d)
d = estimateTagwiseDisp(d)
# 9. plot the mean-variance relationship:
pdf(file.path(results,"mean.variance-edgeR.pdf"))
plotMeanVar(d, show.tagwise.vars=TRUE, NBline=TRUE,main="MeanVar")
plotBCV(d,main="BCV")
dev.off()

# edgeR —c omplex design
# 10. Create a design matrix to specify the factors that are expected to affect
# expression levels:
design = model.matrix( ~ treatment, samples) ## samples is your sample sheet, treatment is a column in the sample sheet
design
### Here it is pH6 - pH2.5 (so positive fold change indicates higher expression in Normal, negative is higher in Affected)

# 11. Estimate dispersion values, relative to the design matrix, using the Cox-Reid (CR)-adjusted likelihood
d2 = estimateGLMTrendedDisp(d, design)
d2 = estimateGLMTagwiseDisp(d2, design)

# 12. Given the design matrix and dispersion estimates, fit a GLM to each feature:
f = glmFit(d2, design)

# 13. Perform a likelihood ratio test, specifying the difference of interest
de = glmLRT(f, coef=2) ## Treatment coefficient

# 14. Use the topTags function to present a tabular summary of the differential expression statistics
tt = topTags(de, n=nrow(d)) ## all tags, sorted
head(tt$table) ## Check result
table(tt$table$FDR< 0.05) ## the number of "Statistically Differentially Expressed Genes" at an FDR of 0.05

# 15. Inspect the depth-adjusted reads per million for some of the top differentially expressed genes:
nc = cpm(d, normalized.lib.sizes=TRUE)
rn = rownames(tt$table)
head(nc[rn,order(samples$treatment)],5)

# 16. Plot the M (log-fold change) versus A (log-average expression)
deg = rn[tt$table$FDR < .05]
pdf(file.path(results,"smear-edgeR.pdf"))
plotSmear(d, de.tags=deg,main="Smear")
dev.off()


# 17. Save the result table as a CSV file:

write.table(tt$table,file=file.path(results,"toptags_edgeR.annotated.txt"),sep="\t",row.names=T,col.names=T,quote=F)

### if you have annotation
# anno <- read.table(annof,sep="\t",header=T,comment.char="",quote="",as.is=T)
# write.table(data.frame(tt$table,anno[match(rownames(tt$table),anno$Ensembl.Gene.ID),]),file=file.path(results,"toptags_tibia_edgeR.annotated.txt"),sep="\t",row.names=T,col.names=T,quote=F)

### END OF DIFFERENTIAL EXPRESSION ANALYSIS, ADDITIONAL PROCESSING ON THE TABLE CAN BE PERFORMED

