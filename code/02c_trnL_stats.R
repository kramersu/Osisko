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
library(ggpubr)
library(cowplot)

path<-c("D:\\marta/trnL/results")
setwd("D:/marta/trnL/scripts")
load(paste0(path,"/ps_base.RData"))

ps_clean1<-subset_taxa(ps_base,!(Class=="NA"))#15 taxa in 10 samples, 9 samples have reads
ps_clean3<-prune_samples(sample_sums(ps_base)>0,ps_base)

ps_clean2<-prune_samples(sample_sums(ps_clean1)>0,ps_clean1)# only 3 samples remain, so we go ahead with ps_clean3

meta<-read.csv("../DB/eDNA_metadata.csv",h=T)
row.names(meta)<-meta$Name

sample_data(ps_clean3)<-meta
ps_clean4<-subset_samples(ps_clean3,Name!="OSID_2023_A_21")#remove low count sample


#alpha diversity analysis
set.seed(687887)
psrare<-rarefy_even_depth(ps_clean4,sample.size=159039)#48 taxa in 8 samples left
ps_alpha_div <- estimate_richness(psrare, split = TRUE, measure = "Shannon")
meta<-merge(meta,ps_alpha_div,by=0)

boxplot(Shannon~Site,meta)
fit<-aov(Shannon~Site,meta)
anova(fit)#ns

fit<-lm(Shannon~Min_Depth,meta)
summary(fit)#ns

ggplot(meta,aes(x=Min_Depth,y=Shannon,color=Lake,shape=Site))+geom_point()+geom_smooth(method="lm")#increasing diversity with increasing depth

OD<-subset(meta,Combo=="Osisko_Deep")
OL<-subset(meta,Combo=="Osisko_Littoral")


fit<-lm(Shannon~Min_Depth,OD)
summary(fit)#ns

fit<-lm(Shannon~Min_Depth,OL)
summary(fit)#ns

#beta diversity analysis
beta<-transform_sample_counts(ps_clean4,function(x) x/sum(x))

tabred<-otu_table(beta)
hellinger<-decostand(tabred,method="hellinger")
meta2<-subset(meta,meta$Name %in% row.names(tabred))
adonis_model <- adonis2(hellinger~Site,data=meta2, permutations = 999)
adonis_model#ns

adonis_model <- adonis(tabred~Min_Depth,data=meta2, permutations = 999)
adonis_model$aov.tab#ns

dis<-vegdist(tabred,method="bray")
pcoa<-cmdscale(dis,eig=TRUE)
ordiplot(pcoa)
m<-merge(meta2,pcoa$points,by.x="Row.names",by.y=0)
m<-m[order(m$Combo,m$Min_Depth),]
ggplot(m,aes(V1,V2,label=Name2))+geom_point(aes(color=Combo,size=Min_Depth,shape=Treatment))+theme_classic()+geom_text(hjust=1,vjust=1)

#Individual cores
OD<-c("OSID_2023_A_13","OSID_2023_A_1","OSID_2023_A_5","OSID_2023_A_9")
OD<-subset_samples(beta,Name %in% OD)
dis<-vegdist(otu_table(OD),method="bray")
pcoa<-cmdscale(dis,eig=TRUE)
ordiplot(pcoa)
m<-merge(meta,pcoa$points,by.x="Row.names",by.y=0)
ggplot(m,aes(V1,V2,label=Name2))+geom_point(aes(size=Min_Depth),col="blue")+theme_classic()+geom_text(hjust=1,vjust=1)
ggplot(m,aes(V1,V2,label=Name2))+geom_point(aes(size=Min_Depth),col="blue")+theme_classic()+geom_path(aes(group=1))
sample_names(OD)<-sample_data(OD)$Name2
p1<-plot_bar(OD,fill="Family",title="Osisko Deep")

OL<-c("OSIL_2023_D_17","OSIL_2023_D_1","OSIL_2023_D_29","OSIL_2023_D_9")
OL<-subset_samples(beta,Name %in% OL)
dis<-vegdist(otu_table(OL),method="bray")
pcoa<-cmdscale(dis,eig=TRUE)
ordiplot(pcoa)
m<-merge(meta,pcoa$points,by.x="Row.names",by.y=0)
ggplot(m,aes(V1,V2,label=Name2))+geom_point(aes(size=Min_Depth),col="lightblue")+theme_classic()+geom_text(hjust=1,vjust=1)
ggplot(m,aes(V1,V2,label=Name2))+geom_point(aes(size=Min_Depth),col="lightblue")+theme_classic()+geom_path(aes(group=Combo))
sample_names(OL)<-sample_data(OL)$Name2
p2<-plot_bar(OL,fill="Family",title="Osisko Littoral")

ggarrange(p1,p2,labels=c("A","B"),ncol=2)
