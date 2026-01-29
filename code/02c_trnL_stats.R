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
load(paste0(path,"/ps_base.RData"))#10 samples, 8 have reads

meta<-read.csv("../DB/eDNA_metadata.csv",h=T)
row.names(meta)<-meta$Id_sample

sample_data(ps_base)<-meta
ps_clean1<-subset_taxa(ps_base,!(Phylum=="NA"))#50 taxa in 10 samples, 8 samples have reads, only one no mine samples
ps_clean2<-prune_samples(sample_sums(ps_clean1)>0,ps_clean1)

set.seed(8375)
psrare<-rarefy_even_depth(ps_clean2,sample.size=471)#41 taxa in 8 samples left
ps_alpha_div <- estimate_richness(psrare, split = TRUE, measure = "Shannon")
ps_alpha_chao <- estimate_richness(ps_clean2, split = TRUE, measure = "Chao1")
ps_alpha_simp <- estimate_richness(ps_clean2, split = TRUE, measure = "Simpson")
alpha_div<-cbind(ps_alpha_div,ps_alpha_chao,ps_alpha_simp)
meta<-merge(meta,alpha_div,by=0)

sh<-ggplot(data=meta,aes(mine_open,Shannon))+geom_boxplot()+theme_bw()
ch<-ggplot(data=meta,aes(mine_open,Chao1))+geom_boxplot()+theme_bw()
si<-ggplot(data=meta,aes(mine_open,Simpson))+geom_boxplot()+theme_bw()
plot_grid(sh,ch,si,labels=c("A","B","C"))
#no test possible, since there is only one no-mine sample
sh<-ggplot(data=meta,aes(Env,Shannon))+geom_boxplot()+theme_bw()
ch<-ggplot(data=meta,aes(Env,Chao1))+geom_boxplot()+theme_bw()
si<-ggplot(data=meta,aes(Env,Simpson))+geom_boxplot()+theme_bw()
plot_grid(sh,ch,si,labels=c("A","B","C"))

ggplot(meta,aes(x=Year,y=Shannon,color=ID_Site,shape=mine_open))+
  geom_point()+theme_bw()+geom_smooth(method="lm",se=FALSE)+
  facet_wrap(~ID_Site,scales="free")

ggplot(meta,aes(x=Year,y=Shannon,color=ID_Site,shape=mine_open))+
  geom_point()+theme_bw()+geom_smooth(method="lm",se=FALSE)

ggplot(meta,aes(x=Year,y=Chao1,color=ID_Site,shape=mine_open))+
  geom_point()+theme_bw()+geom_smooth(method="lm",se=FALSE)+
  facet_wrap(~ID_Site,scales="free")

ggplot(meta,aes(x=Year,y=Chao1,color=ID_Site,shape=mine_open))+
  geom_point()+theme_bw()+geom_smooth(method="lm",se=FALSE)

ggplot(meta,aes(x=Year,y=Simpson,color=ID_Site,shape=mine_open))+
  geom_point()+theme_bw()+geom_smooth(method="lm",se=FALSE)+
  facet_wrap(~ID_Site,scales="free")

ggplot(meta,aes(x=Year,y=Simpson,color=ID_Site,shape=mine_open))+
  geom_point()+theme_bw()+geom_smooth(method="lm",se=FALSE)

fit<-lm(Shannon~Env,meta)
summary(fit)#marginal

fit<-lm(Chao1~Env,meta)
summary(fit)#sig

fit<-lm(Simpson~Env,meta)
summary(fit)#marginal


OD<-subset(meta,ID_Site=="OSI_DEEP")#4 samples
OD<-OD[order(OD$Year),]
OL<-subset(meta,ID_Site=="OSI_LIT")#4 samples
OL<-OL[order(OL$Year),]

fit<-lm(Shannon~Year,OD)
summary(fit)#ns
fit<-lm(Chao1~Year,OD)
summary(fit)#marginal
fit<-lm(Simpson~Year,OD)
summary(fit)#ns

fit<-lm(Shannon~Year,OL)
summary(fit)#ns
fit<-lm(Chao1~Year,OL)
summary(fit)#ns
fit<-lm(Simpson~Year,OL)
summary(fit)#ns

#beta diversity analysis
beta<-transform_sample_counts(ps_clean2,function(x) x/sum(x))

tabred<-otu_table(beta)
bc<-vegdist(tabred,method="bray")

pcoa<-cmdscale(bc,eig=TRUE)
ordiplot(pcoa)
m<-merge(meta,pcoa$points,by.x="Id_sample",by.y=0)
ggplot(m,aes(V1,V2))+geom_point(aes(color=mine_open,shape=ID_Site))+theme_classic()
ggplot(m,aes(V1,V2))+geom_point(aes(color=pre_post,shape=ID_Site))+theme_classic()

env_sub<-meta[c(16,20,23,24,47,25)]
row.names(env_sub)<-meta$Row.names
ev<-envfit(pcoa,env_sub,permutations=999)#ns

ggplot(m,aes(V1,V2))+
  geom_point(aes(color=mine_open,shape=ID_Site))+
  theme_classic()+
  geom_hline(yintercept=0,linetype="dashed")+
  geom_vline(xintercept=0,linetype="dashed")


#rel abundance plot of different divisions
group_ps = tax_glom(beta, taxrank="Order")# tax glom to order level

barplot<-plot_bar(group_ps,x="Year",fill="Order")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  facet_wrap(~ID_Site)+theme_bw()