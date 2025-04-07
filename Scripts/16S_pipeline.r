library("Rcpp")
library("dada2")
library("fastqcr")
library("ShortRead")
library("Biostrings")
library("phyloseq")
library("MicrobiotaProcess")
library("ggplot2")
library("ranacapa")
library("stringr")


library(devtools)

devtools::install_github("gauravsk/ranacapa")

# set global variables
nthreads=10
nsamples=14
# Set your personal folder as working directory
my_project_directory <- "/home/baccus/Desktop/16S_Zebra"
setwd(my_project_directory) 
# 
# Specify raw reads directory
path <- "/home/baccus/Desktop/16S_Zebra/reads"
list.files(path)
#
#### Check libraries
fastqc(fq.dir = path, # FASTQ files directory - same as before (path)
       qc.dir = "/home/baccus/Desktop/16S_Zebra/fastqc/before_cutadapt", # Results directory in your personal folder
       threads = nthreads# Number of threads
       )
##
### Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq --- check the format of your raw reads
fnFs <- sort(list.files(path, pattern="_R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq.gz", full.names = TRUE))
##
FWD <- "^GTGYCAGCMGCCGCGGTAA...ATTAGAWACCCBNGTAGTCC" # CHANGE - if needed - according to your primers: for 16S 515F-806R = ^GTGYCAGCMGCCGCGGTAA...ATTAGAWACCCBNGTAGTCC --- for ITS ITS3_KYO2-ITS4r = ^GATGAAGAACGYAGYRAA...GCATATCAATAAGCGGAGGA
REV <- "^GGACTACNVGGGTWTCTAAT...TTACCGCGGCKGCTGRCAC" # CHANGE - if needed - according to your primers: for 16S 515F-806R = ^GGACTACNVGGGTWTCTAAT...TTACCGCGGCKGCTGRCAC --- for ITS ITS3_KYO2-ITS4r = ^TCCTCCGCTTATTGATATGC...TTYRCTRCGTTCTTCATC
##
cutadapt <- "/home/baccus/miniconda3/envs/cutadapt/bin/cutadapt" # Cutadapt environment path --- do not change it
path.cut <- file.path(dirname(path), "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))
##
### Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-a", FWD) 
### Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-A", REV) 
### Run Cutadapt
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags, R2.flags,"-j",nthreads, "-n", 2, # -n 2 required to remove FWD and REV from reads
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                             fnFs[i], fnRs[i] # input files
                            )) # If after cutadapt you lose almost ALL of your reads, remove the command "--discard-untrimmed"
}
##
############################################################
#### remove the last n bases from R2 in case of overlap, skip this part if you do not want to trim.
#### it takes around ten seconds per sample because of decompression and compression
###pathtopy <- "/home/bioinfo/bin/trim_reads.py" # path to python script that performs the trimming (on loki it is there), change in case.
###ntocut=9 # change 9 with the number of bases you want to trim from R2 3’
###system(sprintf("python3 %s %s %s",pathtopy,path.cut,ntocut)) # Results are stored in “unoverlapped” -  folder, syster of “cutadapt”.
###path.cut<-file.path(dirname(path), "unoverlapped")
###################################################################
##
# Another quality check - this time after cutadapt
fastqc(fq.dir = path.cut, # FASTQ files directory
       qc.dir = "fastqc/after_cutadapt", # Results directory
       threads = nthreads                    # Number of threads
)
##
### Now we will work with trimmed reads - after cutadapt. Again, check the format of your reads
cutFs <- sort(list.files(path.cut, pattern = "R1_001.fastq.gz", full.names = TRUE))
cutRs <- sort(list.files(path.cut, pattern = "R2_001.fastq.gz", full.names = TRUE))
##
## Extract sample names (the part of the string that matches the pattern of a single uppercase letter followed by exactly two digits).:
###get.sample.name <- function(fname) {
###  result <- str_extract(basename(fname), "[A-Z]\\d{2}")
###  return(result)
###}
get.sample.name <- function(fname) {
  result <- str_split_i(str_split_i(fname, '/', -1),'_',1)
  return(result)
}
sample.names <- unname(sapply(cutFs, get.sample.name))
##### Old version, removes name part starting from "_".  The underscore is referred to the one before R1_001/2.fastq.gz, that is already removed. If your reads names are too long, change the basename indicating what you want to remove from the name.
###get.sample.name <- function(fname) strsplit(basename(fname), "_")[[1]][1] 
###sample.names <- unname(sapply(cutFs, get.sample.name))
##
plotQualityProfile(cutFs[1:2])
plotQualityProfile(cutRs[1:2])
##
filtFs <- file.path(path.cut, "filtered", basename(cutFs))
filtRs <- file.path(path.cut, "filtered", basename(cutRs))
names(filtFs) <- sample.names
names(filtRs) <- sample.names
length_to_keep <- c(132,130)
out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs,truncLen=length_to_keep,
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=nthreads) 
# --discard_untrimmed ?
### Error Rates
errF <- learnErrors(filtFs, multithread=nthreads) # Set multithread=TRUE to use all available cores - if nobody else is using the server!
errR <- learnErrors(filtRs, multithread=nthreads) # Set multithread=TRUE to use all available cores - if nobody else is using the server!
##
plotErrors(errF, nominalQ=TRUE)

#
## ASVs are generated one sample by one. pool=FALSE is recommended. If pool=TRUE, ASVs are generated pulling together all the samples - in this case it will take more time and you will have a much higher number of ASVs.
dadaFs <- dada(filtFs, err=errF, pool=FALSE, multithread=nthreads) # Set multithread=TRUE to use all available cores - if nobody else is using the server!
dadaRs <- dada(filtRs, err=errR, pool=FALSE, multithread=nthreads) # Set multithread=TRUE to use all available cores - if nobody else is using the server!

# Merge paired reads (justConcatenate function: concatenate the paired instead of merging, by putting 10 Ns between them). Set FALSE if your reads can overlap.
#mergers <- modified_mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE, justConcatenate=TRUE)
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE, justConcatenate=FALSE,minOverlap=8)
# Construct sequence table
seqtab <- makeSequenceTable(mergers)
### Remove Chimeras ###
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=nthreads, verbose=TRUE)
# Track reads
getN <- function(x) sum(getUniques(x))
# For one sample:
# track <- cbind(length_to_keep[1],length_to_keep[2], out, getN(dadaFs), getN(mergers), rowSums(seqtab), rowSums(seqtab.nochim))
# For more than one sample:
save.image("~workspace_aftermerging.RData")

track <- cbind(length_to_keep[1],length_to_keep[2], sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("length_FWD","length_REV", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
write.csv(as.data.frame(track), 
          file="quality_info.csv", quote=F, row.names = TRUE)

dir.create("dada2_results")

# Choose path to the appropriate database file:
# SILVA --> "/home/baccus/Desktop/Prestazioni/Pipeline_16S/db/silva_nr99_v138.1_train_set.fa.gz"
# UNITE --> "/home/bioinfo/db/amplicon_db/UNITE_release_29.11.2022.fasta"
taxa <- assignTaxonomy(seqtab.nochim, "/home/baccus/Desktop/Prestazioni/Pipeline_16S/db/silva_nr99_v138.1_train_set.fa.gz", multithread=nthreads) # Set multithread=TRUE to use all available cores - if nobody else is using the server!

# Creating file with Taxonomy
head(taxa)

# You can save the workspace up to now:
save.image("~workspace_afterTAXA.RData")


write.csv(as.data.frame(taxa),
          file="dada2_results/taxa.csv")

# Creating file with Reads
seqtab.transpose <- t(seqtab.nochim)
write.csv(as.data.frame(seqtab.transpose),
          file="dada2_results/count_reads.csv")

# Creating file taxa final
ASV <- paste(paste("ASV_",(1:nrow(taxa))), sep = "")
taxa <- cbind(ASV,taxa)
rownames(taxa) <- NULL
##View(taxa)

# Creating file ASV final
seqtab.transpose <- t(seqtab.nochim)
ASV <- paste(paste("ASV_",(1:nrow(seqtab.transpose))), sep = "")
asv <- cbind(ASV,seqtab.transpose)
rownames(asv) <- NULL
##View(asv)

# merged_2 contains ASVs and Taxonomy
merged_2 <- merge(taxa, asv, by = "ASV")


e <- cbind(asv, rownames(seqtab.transpose))
# merged_3 contains ASVs, Taxonomy, Sequences and reads counts per sample
merged_3 <- merge(taxa, e, by = "ASV")
##View(merged_3)

# If you want to save final tables that include singletons, un-# and run the following commands:
write.csv(as.data.frame(merged_3),
          file="dada2_results/final_table_seq_repres_WITH_SING.csv", quote=F, row.names = FALSE)
write.csv(as.data.frame(merged_2),
          file="dada2_results/final_table_WITH_SING.csv", quote=F, row.names = FALSE)

table_for_phyloseq <- subset(merged_2, apply(merged_2, 1, function(x) any(as.numeric (x) > 5))) # table_for_phyloseq does not contain sequence strings
write.csv(as.data.frame(merged_2),
          file="dada2_results/table_for_phyloseq.csv", quote=F, row.names = FALSE) # Does not have sequence string - needed for being imported to phyloseq

# Filtering keepimg only ASVs with >5 counts in at least one sample - removal of singletons
final_table_NO_SING <- subset(merged_3, apply(merged_3, 1, function(x) any(as.numeric (x) > 5)))

# Creating final file with filtered reads - no singletons
write.csv(as.data.frame(final_table_NO_SING),
          file="dada2_results/final_table.csv", quote=F, row.names = FALSE)


#Phyloseq diversity analysis
df <- read.csv("dada2_results/table_for_phyloseq.csv", sep=",")
head(df)

# Check first with: ncol(df) --- WARNING: for ITS analysis, UNITE is also assigning the species taxonomy, so the file will have 1 more column, so you have to increase the numbers of columns by 1.
nsamples=33 ########################
MySampleNames <- colnames(df[8:(7+nsamples)]) # and MODIFY "---" according to the number of samples: 8:7+YourSampleNumber for 16S, or 9:8+YourSampleNumber for ITS)
MyClusterNames <- c(df[,1]) # This refers to the ASVs name column
Mytaxanames <- colnames(df[2:7]) # This refers to the taxonomy columns --- 2:7 for 16S, 2:8 for ITS
TAXtab <- df[2:7] # 2:7  for 16S, 2:8 for ITS
colnames(TAXtab) <- Mytaxanames
row.names(TAXtab) <- MyClusterNames
# Check first with: ncol(df) --- Same WARNING as before, for ITS analysis
OTUtab <- df[8:(7+nsamples)] # MODIFY "---" according to the number of samples - as before --- 8:7+YourSampleNumber for 16S, 9:8+YourSampleNumber for ITS
colnames(OTUtab) <- MySampleNames
row.names(OTUtab) <- MyClusterNames
TAXtab <- as.matrix(TAXtab)
OTUtab <- as.matrix(OTUtab)

Name <- sample.names
Experiment <- sample.names
metadata <- data.frame(Name, Experiment)
write.csv(metadata, "metadata.csv")

Samples <- as.character(metadata$Name)
Name <- as.character(metadata$Experiment) 

SAMPLEtab <- as.data.frame(matrix(c(Samples,Name), nrow=nsamples, ncol=2)) # MODIFY "---" according to the number of samples in the metadata file - as before
colnames(SAMPLEtab) <- colnames(metadata)
row.names(SAMPLEtab) <- Samples

obj_all_raw <- phyloseq(otu_table(OTUtab, taxa_are_rows=TRUE), phyloseq::tax_table(TAXtab))
sample_data(obj_all_raw) <- SAMPLEtab

head(phyloseq::tax_table(obj_all_raw))

save(obj_all_raw, file = "phyloseq_object.RData")

### Alpha Diversity ###

dir.create("phyloseq")

diversity <- estimate_richness(obj_all_raw, split = TRUE)

write.csv(as.data.frame(diversity),
          file="phyloseq/alpha_diversity.csv", quote=F, row.names = TRUE)

obj <- rarefy_even_depth(obj_all_raw)
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
## BiocManager::install("BiocUpgrade") ## you may need this
BiocManager::install("MicrobiotaProcess")
set.seed(1024)
rareres <- get_rarecurve(obj=obj, chunks=50) # You can increase the number of chunks, if you want more reads subsampling

rare <- rareres[["data"]]
rare_observe <- rare %>% filter(Alpha %in% "Observe")
rare_shannon <- rare %>% filter(Alpha %in% "Shannon")

write.csv(as.data.frame(rare_observe),
          file="phyloseq/rarefaction_observe.csv", quote=F, row.names = FALSE)

write.csv(as.data.frame(rare_shannon),
          file="phyloseq/rarefaction_shannon.csv", quote=F, row.names = FALSE)

# Plot the rarefaction curve
p <- ggrare(obj, step = 100, color = "Name", se = FALSE)
p <- p + geom_line(linewidth = 1.5)
p <- p + theme(axis.title.x = element_text(color="black", size=14, face="bold"))
p <- p + theme(axis.title.y = element_text(color="black", size=14, face="bold"))
p <- p + theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size=14,color="black")) # x axis values
p <- p + theme(axis.text.y = element_text(angle = 0, hjust = 0.5, size=14,color="black")) # y axis values
pf <- p + theme(legend.text = element_text(size = 10, colour = "black"))
pf

ggsave("Rarefaction_curve.png", plot=pf)

### Beta Diversity ###

bray <- phyloseq::distance(obj_all_raw, method = "bray")
ord_bray = ordinate(obj_all_raw, method = "PCoA", distance = bray)

bray_vectors <- ord_bray[["vectors"]]
bray_values <- ord_bray[["values"]]

write.csv(as.data.frame(bray_vectors),
          file="phyloseq/beta_diversity_vectors.csv", quote=F, row.names = TRUE)

write.csv(as.data.frame(bray_values),
          file="phyloseq/beta_diversity_values.csv", quote=F, row.names = TRUE)

# length to keep 215,200 | 132,130

