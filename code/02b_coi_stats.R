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
library(lme4)
library(microViz)
library(ggplot2)
library(randomForest)
library(cowplot)

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
ps_clean1<-prune_samples(sample_sums(ps_clean1)>600,ps_clean1)#22 samples remain

plot_bar(ps_clean1,fill="Class",title="cleaned data")
plot_bar(ps_clean1,fill="Order",title="cleaned data")

meta<-read.csv("../DB/eDNA_metadata.csv",h=T)
row.names(meta)<-meta$Id_sample
sample_data(ps_clean1)<-meta

#alpha diversity analysis
set.seed(684)
psrare<-rarefy_even_depth(ps_clean1,sample.size=min(sample_sums(ps_clean1)))#87 taxa in 22 samples left
ps_alpha_div <- estimate_richness(psrare, split = TRUE, measure = "Shannon")
ps_alpha_chao <- estimate_richness(ps_clean1, split = TRUE, measure = "Chao1")
ps_alpha_simp <- estimate_richness(ps_clean1, split = TRUE, measure = "Simpson")
alpha_div<-cbind(ps_alpha_div,ps_alpha_chao,ps_alpha_simp)
meta<-merge(meta,alpha_div,by=0)

boxplot(Shannon~mine_open,meta)
boxplot(Shannon~pre_post,meta)
wilcox.test(Shannon~mine_open,meta)#sig lower in yes category
wilcox.test(Shannon~pre_post,meta)#ns
wilcox.test(Shannon~Site,meta)#sig
boxplot(Shannon~Site,meta)#Dufay has higher diversity
wilcox.test(Shannon~Env,meta)#ns
sh<-ggplot(data=meta,aes(mine_open,Shannon))+geom_boxplot()+theme_bw()
ch<-ggplot(data=meta,aes(mine_open,Chao1))+geom_boxplot()+theme_bw()
si<-ggplot(data=meta,aes(mine_open,Simpson))+geom_boxplot()+theme_bw()
plot_grid(sh,ch,si,labels=c("A","B","C"))

wilcox.test(Chao1~mine_open,meta)#sig.
wilcox.test(Chao1~pre_post,meta)#marginal
wilcox.test(Chao1~Site,meta)#sig
wilcox.test(Chao1~Env,meta)#ns
boxplot(Chao1~mine_open,meta)
boxplot(Chao1~Site,meta)

wilcox.test(Simpson~mine_open,meta)#sig
wilcox.test(Simpson~pre_post,meta)#ns
wilcox.test(Simpson~Site,meta)#sig
wilcox.test(Simpson~Env,meta)#ns
boxplot(Simpson~mine_open,meta)

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

fit<-lm(Shannon~Year*ID_Site,meta)
summary(fit)#ns

fit<-lm(Chao1~Year*ID_Site,meta)
summary(fit)#ns

fit<-lm(Simpson~Year*ID_Site,meta)
summary(fit)#ns

DD<-subset(meta,ID_Site=="DUF_DEEP")
DD<-DD[order(DD$Year),]
DL<-subset(meta,ID_Site=="DUF_LIT")
DL<-DL[order(DL$Year),]
OD<-subset(meta,ID_Site=="OSI_DEEP")
OD<-OD[order(OD$Year),]
OL<-subset(meta,ID_Site=="OSI_LIT")
OL<-OL[order(OL$Year),]

fit<-lm(Shannon~Year,OD)
summary(fit)#ns
fit<-lm(Chao1~Year,OD)
summary(fit)#ns
fit<-lm(Simpson~Year,OD)
summary(fit)#ns

fit<-lm(Shannon~Year,OL)
summary(fit)#ns
fit<-lm(Chao1~Year,OL)
summary(fit)#ns
fit<-lm(Simpson~Year,OL)
summary(fit)#ns

fit<-lm(Shannon~Year,DD)
summary(fit)#ns
fit<-lm(Chao1~Year,DD)
summary(fit)#ns
fit<-lm(Simpson~Year,DD)
summary(fit)#ns

fit<-lm(Shannon~Year,DL)
summary(fit)#ns
fit<-lm(Chao1~Year,DL)
summary(fit)#ns
fit<-lm(Simpson~Year,DL)
summary(fit)#ns

#beta diversity analysis
beta<-transform_sample_counts(ps_clean1,function(x) x/sum(x))

tabred<-otu_table(beta)
bc<-vegdist(tabred,method="bray")
meta$Id_sample==row.names(tabred)#assume order of samples is the same
adonis_model <- adonis2(bc~mine_open,data=meta, permutations = 999)
adonis_model#marginal

anos_mod<-anosim(bc,meta$mine_open,permutations=999)
summary(anos_mod)#sig

adonis_model <- adonis2(bc~pre_post,data=meta, permutations = 999)
adonis_model#ns

anos_mod<-anosim(bc,meta$pre_post,permutations=999)
summary(anos_mod)#ns

pcoa<-cmdscale(bc,eig=TRUE)
ordiplot(pcoa)
m<-merge(meta,pcoa$points,by.x="Id_sample",by.y=0)
ggplot(m,aes(V1,V2))+geom_point(aes(color=mine_open,shape=ID_Site))+theme_classic()
ggplot(m,aes(V1,V2))+geom_point(aes(color=pre_post,shape=ID_Site))+theme_classic()

env_sub<-meta[c(16,20,23,24,47,25)]
row.names(env_sub)<-meta$Row.names
ev<-envfit(pcoa,env_sub,permutations=999)#year significant. Cu marginal
arrow.df<-data.frame(ev$vectors$arrows)
arrow.df<-arrow.df[1:2,]
arrowmu<-ordiArrowMul(ev)
arrow_map = aes(xend = Dim1*arrowmu, yend = Dim2*arrowmu, x = 0, y = 0, shape = NULL, color = NULL, 
                label = row.names(arrow.df))
label_map = aes(x = 1.2 * Dim1*arrowmu, y = 1.2 * Dim2*arrowmu, shape = NULL, color = NULL, 
                label = row.names(arrow.df))
arrowhead=arrow(length=unit(0.03,"npc"))

ggplot(m,aes(V1,V2))+
  geom_point(aes(color=mine_open,shape=ID_Site))+
  theme_classic()+
  geom_segment(arrow_map,size=0.5,data=arrow.df,arrow=arrowhead)+
  geom_text(label_map,size=4,data=arrow.df)+
  geom_hline(yintercept=0,linetype="dashed")+
  geom_vline(xintercept=0,linetype="dashed")

#rel abundance plot of different phyla
group_ps = tax_glom(beta, taxrank="Phylum")# tax glom to phylum level, only chordata present across most samples

barplot<-plot_bar(group_ps,x="Sampling_Order",fill="Phylum")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  facet_wrap(~ID_Site)+theme_bw()

chord_beta<-subset_taxa(beta,Phylum=="Chordata_7711")#77 taxa
chord_genus<-tax_glom(chord_beta, taxrank="Genus")# tax glom to genus level
plot_bar(chord_genus,x="Sampling_Order",fill="Genus")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  facet_wrap(~ID_Site)+theme_bw()

#focus on fish
acti<-subset_taxa(ps_clean1,Class=="Actinopteri_186623")#76 taxa
acti<-subset_taxa(acti,Species!="Ptychocheilus_oregonensis_71769")#remove pike minnow not present in Quebec
acti1<-prune_samples(sample_sums(acti)>2,acti)#16 samples remain
set.seed(8234)
actirare<-rarefy_even_depth(acti1,sample.size=min(sample_sums(acti1)))#72 taxa in 16 samples left
acti_alpha_div <- estimate_richness(actirare, split = TRUE, measure = "Shannon")
acti_alpha_chao <- estimate_richness(acti1, split = TRUE, measure = "Chao1")
acti_alpha_simp <- estimate_richness(acti1, split = TRUE, measure = "Simpson")
acti_meta<-sample_data(actirare)
acti_meta<-merge(acti_meta,acti_alpha_div,by=0)
row.names(acti_meta)<-acti_meta$Row.names
acti_meta$Row.names<-NULL
acti_meta<-merge(acti_meta,acti_alpha_chao,by=0)
row.names(acti_meta)<-acti_meta$Row.names
acti_meta$Row.names<-NULL
acti_meta<-merge(acti_meta,acti_alpha_simp,by=0)
row.names(acti_meta)<-acti_meta$Row.names
acti_meta$Row.names<-NULL


sh<-ggplot(data=acti_meta,aes(mine_open,Shannon))+geom_boxplot()+theme_bw()
ch<-ggplot(data=acti_meta,aes(mine_open,Chao1))+geom_boxplot()+theme_bw()
si<-ggplot(data=acti_meta,aes(mine_open,Simpson))+geom_boxplot()+theme_bw()
plot_grid(sh,ch,si,labels=c("A","B","C"))

wilcox.test(Shannon~mine_open,acti_meta)#sig
wilcox.test(Chao1~mine_open,acti_meta)#sig
wilcox.test(Simpson~mine_open,acti_meta)#sig


ggplot(acti_meta,aes(x=Year,y=Shannon,color=ID_Site,shape=mine_open))+
  geom_point()+theme_bw()+geom_smooth(method="lm",se=FALSE)+
  facet_wrap(~ID_Site,scales="free")

fit<-lm(Shannon~Year,acti_meta)
summary(fit)#ns
fit<-lm(Chao1~Year,acti_meta)
summary(fit)#sig negative
fit<-lm(Simpson~Year,acti_meta)
summary(fit)#ns

beta_acti<-transform_sample_counts(acti1,function(x) x/sum(x))
bc_acti<-vegdist(otu_table(beta_acti),method="bray")
meta_acti<-sample_data(beta_acti)
adonis_model <- adonis2(bc_acti~meta_acti$mine_open, permutations = 999)
adonis_model#sig

#fish pcoa
pcoa<-cmdscale(bc_acti,eig=TRUE)
ordiplot(pcoa)
m<-merge(acti_meta,pcoa$points,by.x="Id_sample",by.y=0)
ggplot(m,aes(V1,V2))+geom_point(aes(color=mine_open,shape=ID_Site))+theme_classic()
ggplot(m,aes(V1,V2))+geom_point(aes(color=pre_post,shape=ID_Site))+theme_classic()

env_sub<-meta[c(16,20,23,24,47,25)]
row.names(env_sub)<-meta$Row.names
env_sub<-subset(env_sub,row.names(env_sub) %in% m$Id_sample)
row.names(env_sub)==m$Id_sample

ev<-envfit(pcoa,env_sub,permutations=999)#Cu, Al, Ca significant, year marginal
arrow.df<-data.frame(ev$vectors$arrows)
arrow.df<-arrow.df[1:4,]
arrowmu<-ordiArrowMul(ev)
arrow_map = aes(xend = Dim1*arrowmu, yend = Dim2*arrowmu, x = 0, y = 0, shape = NULL, color = NULL, 
                label = row.names(arrow.df))
label_map = aes(x = 1.2 * Dim1*arrowmu, y = 1.2 * Dim2*arrowmu, shape = NULL, color = NULL, 
                label = row.names(arrow.df))
arrowhead=arrow(length=unit(0.03,"npc"))

ggplot(m,aes(V1,V2))+
  geom_point(aes(color=mine_open,shape=ID_Site))+
  theme_classic()+
  geom_segment(arrow_map,size=0.5,data=arrow.df,arrow=arrowhead)+
  geom_text(label_map,size=4,data=arrow.df)+
  geom_hline(yintercept=0,linetype="dashed")+
  geom_vline(xintercept=0,linetype="dashed")

anos_mod<-anosim(bc_acti,meta_acti$mine_open,permutations=999)
summary(anos_mod)#sig

###analyze fish by trophic level
tax_group<-read.csv("../results/acti1_tax.csv")
row.names(tax_group)<-tax_group$X
tax_group$X<-NULL
tax_tax<-tax_table(tax_group)
row.names(tax_tax)<-row.names(tax_group)
colnames(tax_tax)<-colnames(tax_group)
acti2<-acti1
tax_table(acti2)<-tax_tax
plot_bar(acti_spec2,x="Sampling_Order",fill="Species")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  facet_wrap(~ID_Site)+theme_bw()


plot_bar(acti_spec2,x="Sampling_Order",fill="group")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  facet_wrap(~ID_Site)+theme_bw()

####correlations with qPCR data
acti_spec<-tax_glom(acti1,taxrank="Species")
acti_cor<-merge(data.frame(otu_table(acti_spec)),meta,by.x=0,by.y="Id_sample")
#Sander
cor.test(acti_cor$ASV1,acti_cor$eSAVI2)#sig! 0.73 cor
plot(acti_cor$ASV1,acti_cor$eSAVI2)
#Catostomus
cor.test(acti_cor$ASV3,acti_cor$eCACO4)#sig! 0.77 cor
plot(acti_cor$ASV3,acti_cor$eCACO4)
#Perca
cor.test(acti_cor$ASV4,acti_cor$ePEFL1)#ns
plot(acti_cor$ASV4,acti_cor$ePEFL1)
#Esox
cor.test(acti_cor$ASV8,acti_cor$eESLU1)#sig, 0.998 cor!!!
plot(acti_cor$ASV8,acti_cor$eESLU1)

#random forest for indicator taxa (on species level)
forest_mat<-as.data.frame(otu_table(acti_spec))
samp_class<-data.frame(sample_data(acti_spec))
forest_mat$mine<-samp_class$mine_open
tax_class<-data.frame(tax_table(acti_spec))
acti.rf<-randomForest(as.factor(mine)~.,data=forest_mat,importance=TRUE)#OOB error: 6.25
imp<-importance(acti.rf)
imp<-merge(tax_class,imp,by=0)
imp<-imp[rev(order(imp[,12])),]
#Sander most important, followed by Perca

#export otu table for stamp analysis
beta_acti_otus<-data.frame(t(otu_table(beta_acti)))
beta_tax<-data.frame(tax_table(beta_acti))
beta_m<-merge(beta_tax,beta_acti_otus,by=0)
beta_m$Row.names<-NULL
beta_m$group<-NULL
beta_samp<-data.frame(sample_data(beta_acti))
beta_samp2<-cbind(beta_samp$Id_sample,beta_samp$Env,beta_samp$mine_open,beta_samp$pre_post)
write.table(beta_m,"../results/beta_acti_otus.txt",quote=FALSE,sep="\t",row.names=FALSE)
write.table(beta_samp2,"../results/beta_acti_meta_data.txt",quote=FALSE,sep="\t",row.names=FALSE)
#stamp: no significant differences for any species between exposed and non-exposed
#also make stamp file for groups #no significant stamp differences

save.image(file="coi_stats.RData")

