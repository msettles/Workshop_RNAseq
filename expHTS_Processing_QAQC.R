
outputfile <- "QAQC_preprocessing_mds.pdf"
inputfile <- "02-Cleaned/Preprocessing_Summary.log"

samples <- read.table("samples.txt",sep="\t",header=T)

## by default all points are 'black', however we color by the treatment to look for a treatment pattern
col = rep('black',nrow(samples))
if ("treatment" %in% colnames(samples)) col = as.numeric(as.factor(samples$treatment))+1
names(col) <- samples$SAMPLE_ID

# open an output to pdf file
pdf(outputfile)
# read in the summary table
tb <- read.table(inputfile,sep="\t",header=T,row.names=1,as.is=T, quote="", comment.char="",fill=T)
# scale the data
scaled_tb <- apply(tb,2,function(x) scale(as.numeric(sub("%","",x))))
# perform multi-dimentional scaling
mds <- cmdscale(dist(scaled_tb))
# plot the result
plot(mds[,1], mds[,2], type = "n", xlab = "MDS1", ylab = "MDS2", asp = 1, axes = TRUE,
     main = "cmdscale (scaled results)")
text(mds[,1], mds[,2], rownames(tb), cex = 0.6, col = col[match(rownames(tb),names(col))])

# close the connection to the pdf
dev.off()


