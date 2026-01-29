library(ape)
library(phangorn)
libraRY(DECIPHER)

dna<-readDNAStringSet("trnL_tree.fa")
aligned_DNA<-AlignSeqs(dna)
D<-DistanceMatrix(aligned_DNA)
my_dist_mat<-as.dist(D)
my_nj <- ape::nj(my_dist_mat)
plot(my_nj, "unrooted",cex=0.8)
plot(my_nj,cex=0.5)

my_upgma <- phangorn::upgma(my_dist_mat)
                            