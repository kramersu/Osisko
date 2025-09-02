library(dada2)
packageVersion("dada2")
library(ShortRead)
packageVersion("ShortRead")
library(Biostrings)
packageVersion("Biostrings")
library(reticulate)
library(DECIPHER)
library(dplyr)

#coi analysis
setwd("D:\\marta/coi/raw")
path<-c("D:\\marta/coi")

myfiles <- list.files(pattern = "fastq.gz$")

fnFs <- sort(list.files(".", pattern = "_R1.fastq.gz$", full.names = TRUE))
fnRs <- sort(list.files(".", pattern = "_R2.fastq.gz$", full.names = TRUE))

FWD <- "GGWACWGGWTGAACWGTWTAYCCYCC" #miCOIF
REV <- "TADACYTCDGGRTGDCCRAARAAYCA" # jgHCO2198R


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
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs[[10]]))

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

cutFs <- sort(list.files(path.cut, pattern = "_R1.fastq.gz", full.names = TRUE)) # Pattern name may change!
cutRs <- sort(list.files(path.cut, pattern = "_R2.fastq.gz", full.names = TRUE))

get.sample.name <- function(fname) strsplit(basename(fname), "_coi_")[[1]][1]
sample.names <- unname(sapply(cutFs, get.sample.name))
head(sample.names)

plotQualityProfile(cutFs[2:5]) # Quality drop-off around 280bp
plotQualityProfile(cutRs[2:5]) # Quality drop-off around 280bp

filtFs <- file.path(path, "filtered", basename(cutFs))
filtRs <- file.path(path, "filtered", basename(cutRs))

out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, maxN = 0, maxEE = c(2, 2),
                     truncQ = 2, minLen = 50, rm.phix = TRUE, compress = TRUE, multithread = FALSE)  # on windows, set multithread = FALSE
head(out)

#Remove sample with no or very limited surviving reads:
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
table(nchar(getSequences(seqtab)))#by far most fragments are 313 bp long
seqtab2<-seqtab[,nchar(colnames(seqtab)) %in% 300:336]

seqtab.nochim <- removeBimeraDenovo(seqtab2, method="consensus", multithread=TRUE, verbose=TRUE)
sum(seqtab.nochim)/sum(seqtab2) # 97.4% of reads remain
table(nchar(getSequences(seqtab.nochim)))

getN <- function(x) sum(getUniques(x))
out_ex<-row.names(out)[exists]
out2<-subset(out,row.names(out) %in% out_ex)
track <- cbind(out2, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers,
                                                                        getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged",
                     "nonchim")     
save.image(file = "dada2_pro_coi_workspace1.RData")

dna <- DNAStringSet(getSequences(seqtab.nochim)) # Create a DNAStringSet from the ASVs

library(phyloseq)
ps_base<-phyloseq(otu_table(seqtab.nochim,taxa_are_rows=FALSE))
dna<-Biostrings::DNAStringSet(taxa_names(ps_base))
names(dna)<-taxa_names(ps_base)
ps_base<-merge_phyloseq(ps_base,dna)
taxa_names(ps_base)<-paste0("ASV",seq(ntaxa(ps_base)))
ps_base %>%refseq() %>%Biostrings::writeXStringSet("ASV_coi.fna", append=FALSE,compress=FALSE, compression_level=NA, format="fasta")

taxa_names(ps_base)
tax<-read.csv("coi/results/coi_lca_tax_no_pseudo.csv")
names<-tax$qseqid
tax$qseqid<-NULL
tax<-tax_table(tax)
taxa_names(tax)<-names
colnames(tax)<-c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
tax_table(ps_base)<-tax

save.image(file = 'dada2_pro_coi_workspace1.RData')

####updated tax assignment
taxa<-assignTaxonomy(seqtab.nochim,"../DB/MIDORI2_UNIQ_NUC_GB264_CO1_DADA2.fasta.gz")
#does not work well either: about half of ASVs are assigned to arthropoda, insecta, but non below that level. Everything else is NA

                                                     
save.image(file = '../results/dada2_pro_coi_workspace1.RData')
                                                     
                                                     
                                                     
                                                     
                                                     
