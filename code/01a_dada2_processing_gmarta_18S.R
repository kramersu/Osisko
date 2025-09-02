library(dada2)
packageVersion("dada2")
library(ShortRead)
packageVersion("ShortRead")
library(Biostrings)
packageVersion("Biostrings")
library(reticulate)
library(DECIPHER)
library(dplyr)

###18S analysis
setwd("D:\\marta/18S/raw")
path<-c("D:\\marta/18S")

myfiles <- list.files(pattern = "fastq.gz$")

fnFs <- sort(list.files(".", pattern = "_R1.fastq.gz$", full.names = TRUE))
fnRs <- sort(list.files(".", pattern = "_R2.fastq.gz$", full.names = TRUE))

FWD <- "ACACTCTTTCCCTACACGACGCTCTTCCGATCTGGCTTAATTTGACTCAACRCG" # 960F
REV <- "GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTGGGCATCACAGACCTGTTAT" # 1438R

allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)


primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}

rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs[[10]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs[[10]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs[[10]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs[[10]]))#lists how many primers in which orientations are found

cutadapt <- 'C:/Users/kramers/AppData/Local/Programs/Python/Python313/Scripts/cutadapt.exe'
system2(cutadapt, args = "--version")
path.cut <- file.path(path, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))

FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)


R1.flags <- paste("-g", FWD, "-a", REV.RC) 
R2.flags <- paste("-G", REV, "-A", FWD.RC) 

# Run Cutadapt
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, "-m", 10, # -n 2 required to remove FWD and REV from reads
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                             fnFs[i], fnRs[i])) # input files, remove empty reads at this stage (min length 10)
}

rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[10]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[10]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[10]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[10]]))

cutFs <- sort(list.files(path.cut, pattern = "_R1.fastq.gz", full.names = TRUE))
cutRs <- sort(list.files(path.cut, pattern = "_R2.fastq.gz", full.names = TRUE))

get.sample.name <- function(fname) strsplit(basename(fname), "_960_")[[1]][1]
sample.names <- unname(sapply(cutFs, get.sample.name))
head(sample.names)

plotQualityProfile(cutFs[2:7]) # Quality drop-off around 225bp
plotQualityProfile(cutRs[2:7]) # Quality drop-off around 225bp

filtFs <- file.path(path, "filtered", basename(cutFs))
filtRs <- file.path(path, "filtered", basename(cutRs))

out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, maxN = 0, maxEE = c(2, 2),
                     truncQ = 2, minLen = 50, rm.phix = TRUE, compress = TRUE, multithread = FALSE)  # on windows, set multithread = FALSE
head(out)

#Remove reads with no or very limited surviving reads:
exists<-file.exists(filtFs)#files that have reads after filtering

errF <- learnErrors(filtFs[exists], multithread = TRUE)
errR <- learnErrors(filtRs[exists], multithread = TRUE)

plotErrors(errF, nominalQ = TRUE)# error estimates a bit shady
plotErrors(errR, nominalQ = TRUE)# error estimates a bit shady


derepFs <- derepFastq(filtFs[exists], verbose = TRUE)
derepRs <- derepFastq(filtRs[exists], verbose = TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names[exists]
names(derepRs) <- sample.names[exists]

dadaFs <- dada(derepFs, err = errF, multithread = TRUE,pool="pseudo")
dadaRs <- dada(derepRs, err = errR, multithread = TRUE,pool="pseudo")

mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
names(mergers)<-sample.names[exists]
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
table(nchar(getSequences(seqtab)))#lots of short fragments in the 220 range, then another peak at 360 and 404 bp, target is supposedly 260 bp long...
seqtab2<-seqtab[,nchar(colnames(seqtab)) %in% 218:405]

seqtab.nochim <- removeBimeraDenovo(seqtab2, method="consensus", multithread=TRUE, verbose=TRUE)

sum(seqtab.nochim)/sum(seqtab2) # 95% of reads remain
table(nchar(getSequences(seqtab.nochim)))

getN <- function(x) sum(getUniques(x))
out_ex<-row.names(out)[exists]
out2<-subset(out,row.names(out) %in% out_ex)
track <- cbind(out2, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers,
                                                                        getN), rowSums(seqtab.nochim))


colnames(track)<-c("input","filtered","denoisedF",
                   "denoisedR","merged","nochim")                                                                                                                                   getN), rowSums(seqtab.nochim))

save.image(file = "dada2_pro_18S_workspace1.RData")

#####
dna <- DNAStringSet(getSequences(seqtab.nochim)) # Create a DNAStringSet from the ASVs
load(paste0(path,"/DB/SILVA_SSU_r138_2_2024.RData")) # Training set from decipher
ids <- IdTaxa(dna, trainingSet, strand="top", processors=NULL, verbose=FALSE)
ranks <- c("domain","kingdom",
           "phylum","class","order","family","genus") 
# Convert the output object of class "Taxa" to a matrix analogous to the output from assignTaxonomy
taxid <- t(sapply(ids, function(x) {
  m <- match(ranks, x$rank)
  taxa <- x$taxon[m]
  taxa[startsWith(taxa, "unclassified_")] <- NA
  taxa
}))
colnames(taxid) <- ranks; rownames(taxid) <- getSequences(seqtab.nochim)
taxid[is.na(taxid)] <- "unassigned" #replace NA with unassigned for later filtering in phyloseq

#handover to phyloseq
library(phyloseq)
ps_base<-phyloseq(otu_table(seqtab.nochim,taxa_are_rows=FALSE),tax_table(taxid))#1745 taxa
#store sequence in refseq(ps_base)
dna<-Biostrings::DNAStringSet(taxa_names(ps_base))
names(dna)<-taxa_names(ps_base)
ps_base<-merge_phyloseq(ps_base,dna)
taxa_names(ps_base)<-paste0("ASV",seq(ntaxa(ps_base)))
ps_base %>%refseq() %>%Biostrings::writeXStringSet("ASV_18S.fna", append=FALSE,compress=FALSE, compression_level=NA, format="fasta")

save.image(file = "dada2_pro_18S_workspace1.RData")


#remove  taxa without phylum assignment
ps_clean<-subset_taxa(ps_base,phylum!="unassigned") #611 taxa remain


save(ps_clean, file ='ps_clean_18S.RData')


