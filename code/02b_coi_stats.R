library(tidyr)
library(dplyr)
library(tidyverse)
library(vegan)
library(pairwiseAdonis)
library(pheatmap)
library(phyloseq)
library(phytools)
library(ggtree)
library(pcaMethods)
library(lavaan)
library(lme4)
library(ggord)
library(microViz)
library(TITAN2)
library(cowplot)
library(ggpubr)

path<-c("D:\\Marta/coi/results")
setwd("D:/Marta/coi/scripts")
load(paste0(path,"/ps_base.RData"))
tax<-read.csv(paste0(path,"/coi_tax_combo.csv"))
names<-tax$qseqid
tax$qseqid<-NULL
tax<-tax_table(tax)
taxa_names(tax)<-names
colnames(tax)<-c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
tax_table(ps_base)<-tax

plot_bar(ps_base,fill="Class",title="All ASVs")#lots of hsa contamination in lake Osisko


ps_clean1<-subset_taxa(ps_base,!(Order=="Primates_9443"))#remove hsa contamination, 120 ASVs remain
ps_cont<-subset_taxa(ps_clean1,Order=="Primates_9443")

plot_bar(ps_clean1,fill="Class",title="cleaned data")

ps_clean2<-prune_samples(sample_sums(ps_clean1)>0,ps_clean1)#24 samples remain

meta<-read.csv("../DB/eDNA_metadata.csv",h=T)
row.names(meta)<-meta$Name

sample_data(ps_clean2)<-meta
ps_clean3<-subset_samples(ps_clean2,Name !="BlankpcrCES3")
ps_clean4<-subset_samples(ps_clean3,Name !="OSIL_2023_D_13b")

#alpha diversity analysis
set.seed(67387)
psrare<-rarefy_even_depth(ps_clean4,sample.size=3830)#94 taxa in 22 samples left
ps_alpha_div <- estimate_richness(psrare, split = TRUE, measure = "Shannon")
meta<-merge(meta,ps_alpha_div,by=0)

boxplot(Shannon~Lake+Site,meta)
fit<-aov(Shannon~Lake*Site,meta)
anova(fit)#lake significant: Dufay significantly higher

fit<-aov(Shannon~Combo,meta)
anova(fit)#marginal

fit<-lm(Shannon~Min_Depth,meta)
summary(fit)#ns

ggplot(meta,aes(x=Min_Depth,y=Shannon,color=Lake))+geom_point()+geom_smooth(method="lm")

ggplot(meta,aes(x=Min_Depth,y=Shannon,color=Lake,shape=Site))+geom_point()+geom_smooth(method="lm")

fit<-lm(Shannon~Treatment,meta)
summary(fit)#significant, exposure has higher diversity
boxplot(Shannon~Treatment,meta)

DD<-subset(meta,Combo=="Dufay_Deep")
DD<-DD[order(DD$Year),]
DL<-subset(meta,Combo=="Dufay_Littoral")
DL<-DL[order(DL$Year),]
OD<-subset(meta,Combo=="Osisko_Deep")
OL<-subset(meta,Combo=="Osisko_Littoral")


fit<-lm(Shannon~Min_Depth,OD)
summary(fit)#ns

fit<-lm(Shannon~Min_Depth,OL)
summary(fit)#ns

fit<-lm(Shannon~Year,DD)
summary(fit)#ns
plot(Shannon~Year,DD,type='b')

fit<-lm(Shannon~Year,DL)
summary(fit)#ns
plot(Shannon~Year,DL,type='b')

fit<-lm(Shannon~Treatment,DL)
summary(fit)#ns

#beta diversity analysis
ps_clean5<-tax_filter(ps_clean4,min_prevalence=3,prev_detection_threshold=5)#26 taxa remain
exclude<-c("DUFD_2023_B_13b","OSID_2023_A_11b","OSID_2023_A_19b","OSID_2023_A_9b")#samples with 0 or very low counts to remove
ps_clean6<-subset_samples(ps_clean5,!(Name %in% exclude))
beta<-transform_sample_counts(ps_clean6,function(x) x/sum(x))#18 samples remain

tabred<-otu_table(beta)
hellinger<-decostand(tabred,method="hellinger")
meta2<-subset(meta,meta$Name %in% row.names(tabred))
adonis_model <- adonis2(hellinger~Combo,data=meta2, permutations = 999)
adonis_model#significant

adonis_model <- adonis(hellinger~Lake*Site,data=meta2, permutations = 999)
adonis_model$aov.tab#all significant

adonis_model <- adonis(tabred~Lake*Min_Depth,data=meta2, permutations = 999)
adonis_model$aov.tab#lake significant

adonis_model <- adonis(hellinger~Treatment,data=meta2, permutations = 999)
adonis_model$aov.tab#significant

dis<-vegdist(tabred,method="bray")
pcoa<-cmdscale(dis,eig=TRUE)
ordiplot(pcoa)
m<-merge(meta2,pcoa$points,by.x="Row.names",by.y=0)
m<-m[order(m$Combo,m$Min_Depth),]
ggplot(m,aes(V1,V2,label=Name2))+geom_point(aes(color=Combo,size=Min_Depth,shape=Treatment))+theme_classic()+geom_text(hjust=1,vjust=1)
ggplot(m,aes(V1,V2,label=Name2))+geom_point(aes(color=Combo,size=Min_Depth))+theme_classic()+geom_path(aes(group=Combo))

#Individual cores
OD<-c("OSID_2023_A_15b","OSID_2023_A_1b","OSID_2023_A_23b","OSID_2023_A_5b")
OD<-subset_samples(beta,Name %in% OD)
dis<-vegdist(otu_table(OD),method="bray")
pcoa<-cmdscale(dis,eig=TRUE)
ordiplot(pcoa)
m<-merge(meta,pcoa$points,by.x="Row.names",by.y=0)
ggplot(m,aes(V1,V2,label=Name2))+geom_point(aes(size=Min_Depth),col="blue")+theme_classic()+geom_text(hjust=1,vjust=1)
sample_names(OD)<-sample_data(OD)$Name2
p1<-plot_bar(OD,fill="Family",title="Osisko Deep")

OL<-c("OSIL_2023_D_13b","OSIL_2023_D_17b","OSIL_2023_D_1b","OSIL_2023_D_5b","OSIL_2023_D_9b")
OL<-subset_samples(beta,Name %in% OL)
dis<-vegdist(otu_table(OL),method="bray")
pcoa<-cmdscale(dis,eig=TRUE)
ordiplot(pcoa)
m<-merge(meta,pcoa$points,by.x="Row.names",by.y=0)
ggplot(m,aes(V1,V2,label=Name2))+geom_point(aes(size=Min_Depth),col="lightblue")+theme_classic()+geom_text(hjust=1,vjust=1)
sample_names(OL)<-sample_data(OL)$Name2
p2<-plot_bar(OL,fill="Family",title="Osisko Littoral")
  
DD<-c("DUFD_2023_B_17b","DUFD_2023_B_1b","DUFD_2023_B_21b","DUFD_2023_B_5b","DUFD_2023_B_9b")
DD<-subset_samples(beta,Name %in% DD)
dis<-vegdist(otu_table(DD),method="bray")
pcoa<-cmdscale(dis,eig=TRUE)
ordiplot(pcoa)
m<-merge(meta,pcoa$points,by.x="Row.names",by.y=0)
ggplot(m,aes(V1,V2,label=Name2))+geom_point(aes(size=Min_Depth),col="green")+theme_classic()+geom_text(hjust=1,vjust=1)
sample_names(DD)<-sample_data(DD)$Name2
p3<-plot_bar(DD,fill="Family",title="Dufay Deep")

DL<-c("DUFL_2023_C_17b","DUFL_2023_C_1b","DUFL_2023_C_29b","DUFL_2023_C_5b","DUFL_2023_C_9b")
DL<-subset_samples(beta,Name %in% DL)
dis<-vegdist(otu_table(DL),method="bray")
pcoa<-cmdscale(dis,eig=TRUE)
ordiplot(pcoa)
m<-merge(meta,pcoa$points,by.x="Row.names",by.y=0)
ggplot(m,aes(V1,V2,label=Name2))+geom_point(aes(size=Min_Depth,shape=Treatment),col="darkgreen")+theme_classic()+geom_text(hjust=1,vjust=1)
sample_names(DL)<-sample_data(DL)$Name2
p4<-plot_bar(DL,fill="Family",title="Dufay Littoral")

ggarrange(p1,p2,p3,p4,labels=c("A","B","C","D"),ncol=2,nrow=2)

#Titan analysis
D<-subset_samples(ps_clean6,sample_data(ps_clean6)$Lake=="Dufay")
Drare<-rarefy_even_depth(D,sample.size=21580)

D_clean<-tax_filter(Drare,min_prevalence=3,prev_detection_threshold=3)

txa_obj<-otu_table(D_clean)
env<-sample_data(D_clean)$Year

to<-titan(env,txa_obj,pur.cut=0.95,rel.cut=0.95,numPerm=1000)#number of observations too low
plot_taxa_ridges(to,z1=FALSE)
plot_sumz(to,filter=FALSE,cumfrq=FALSE,xmin=1920,xmax=2022)

O<-subset_samples(ps_clean6,sample_data(ps_clean6)$Lake=="Osisko")
Orare<-rarefy_even_depth(O,sample.size=20500)
O_clean<-tax_filter(Orare,min_prevalence=3,prev_detection_threshold=3)

txa_obj<-otu_table(O_clean)
env<-sample_data(O_clean)$Min_Depth

to<-titan(env,txa_obj,pur.cut=0.95,rel.cut=0.95,numPerm=1000)#number of observations too small
plot_sumz(to,filter=FALSE,cumfrq=FALSE)
plot_taxa_ridges(to)#doesn't work because no taxa pass filter
