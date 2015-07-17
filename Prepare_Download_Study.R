### download data from the article
### http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4059230/
### and SRA project
### http://www.ncbi.nlm.nih.gov/sra/?term=SRP006291
### subselect a random 5M reads (average size 27M too large), and edit read id to be CASAVA 1.8 and output as gzipped

library(parallel)
library(ShortRead)
library(tools)

file2get <- "http://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?save=efetch&rettype=runinfo&db=sra&WebEnv=NCID_1_16454837_130.14.18.48_5555_1435842911_675372944_0MetA0_S_HStore&query_key=1"
download.file(file2get,"SraRunInfo.csv")

samples = c("SRR358569", "SRR358570", "SRR358573", "SRR358574", "SRR358575", "SRR358576")
sample_names = c("NIH8656_pH2.5_3", "NIH8656_pH6_3", "NIH8656_pH6_1", "NIH8656_pH2.5_1", "NIH8656_pH6_2", "NIH8656_pH2.5_2")
treatment = c("pH2.5", "pH6", "pH6", "pH2.5", "pH6", "pH2.5")

cores <- min(detectCores(),length(samples))

write.table(data.frame("SEQUENCE_ID"= samples, "SAMPLE_ID" = sample_names, "treatment" = treatment),"samples.txt",sep="\t",row.names=F,col.names=T,quote=F)

"getSRR" = function(x) paste("ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra",substring(x,1,3),substring(x,1,6),x,paste0(x,".sra"),sep="/")
dir.create("00-RawData")

output <- mclapply(samples,
                   function(sample) {
                        dir.create(file.path("00-RawData",sample))
                        file2get = getSRR(sample)
                        download.file(file2get,file.path("00-RawData",sample,paste0(sample,".sra")))
                        cmd <- paste("fastq-dump --split-3 -F --gzip -O", file.path("00-RawData",sample), file.path("00-RawData", sample, paste0(sample,".sra")),">", file.path("00-RawData", sample, paste0(sample,".sra_extract.txt")))
                        system(cmd)
                        ## read in and subsample to 5M per sample
                        r1 <- grep ('_1.fastq',dir(file.path("00-RawData",sample),full.names = T),value=T)
                        r2 <- grep ('_2.fastq',dir(file.path("00-RawData",sample),full.names = T),value=T)
                        fq1 <- readFastq(r1)
                        fq2 <- readFastq(r2)
                        samp <- sort(sample(1:length(fq1),5e6))
                        fq1_r <- fq1[samp]
                        fq2_r <- fq2[samp]
                        fq1_r@id <- BStringSet(sapply(strsplit(as.character(id(fq1_r)),split=":"),function(x) paste(paste(x[1],"136","FC706VJ",paste(x[2:5],collapse=":"),sep=":"),"1:N:0:",sep=" "))) ## made up run and flowcell id
                        fq2_r@id <- BStringSet(sapply(strsplit(as.character(id(fq2_r)),split=":"),function(x) paste(paste(x[1],"136","FC706VJ",paste(x[2:5],collapse=":"),sep=":"),"2:N:0:",sep=" "))) ## made up run and flowcell id
                        out_r1 <- sub("_1.fastq","_5M_R1.fastq",r1) ## convert to R1,R2 format for compatability with expHTS
                        out_r2 <- sub("_2.fastq","_5M_R2.fastq",r2)
                        if (file_ext(out_r1) != "gz") out_r1 <- paste0(out_r1,".gz")
                        if (file_ext(out_r2) != "gz") out_r2 <- paste0(out_r2,".gz")
                        if (file.exists(out_r1)) unlink(out_r1)
                        if (file.exists(out_r2)) unlink(out_r2)
                        writeFastq(fq1_r,out_r1)
                        writeFastq(fq2_r,out_r2)
                    }, mc.cores = cores)


dir.create("Reference")
## Get Reference Genome:
file2get = "ftp://ftp.ensemblgenomes.org/pub/fungi/release-27/fasta/exophiala_dermatitidis_nih_ut8656/dna/Exophiala_dermatitidis_nih_ut8656.GCA_000230625.1.27.dna.toplevel.fa.gz"
download.file(file2get,file.path("Reference",basename(file2get)))
system(paste("gunzip",file.path("Reference",basename(file2get))))
## Get Reference GTF transcript file:
file2get = "ftp://ftp.ensemblgenomes.org/pub/fungi/release-27/gtf/exophiala_dermatitidis_nih_ut8656/Exophiala_dermatitidis_nih_ut8656.GCA_000230625.1.27.gtf.gz"
download.file(file2get,file.path("Reference",basename(file2get)))
system(paste("gunzip",file.path("Reference",basename(file2get))))
## Get Pfam Annotation:
file2get = "http://www.broadinstitute.org/annotation/genome/Black_Yeasts/download/?sp=EAProteinFamilytoGenes&sp=SExop_derm_V1&sp=S.zip"
download.file(file2get,file.path("Reference",basename(file2get)))
unzip(file.path("Reference",basename(file2get)),exdir = "Reference")
file.remove(file.path("Reference",basename(file2get)))

