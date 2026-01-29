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
library(microViz)
library(randomForest)
library(cowplot)

path<-c("D:\\marta/18S/results")
load(paste0(path,"/ps_clean_18S_pr2.RData"))#1176 taxa

new_tax<-read.csv("../results/tax_pr2_group.csv")
row.names(new_tax)<-new_tax$X
new_tax$X<-NULL
taxid<-tax_table(new_tax)
row.names(taxid)<-row.names(new_tax)
colnames(taxid)<-c("supergroup","division","class","order","family","genus","species","group")
tax_table(ps_clean)<-taxid #557 taxa are euk microbes

ps_clean1<-prune_samples(sample_sums(ps_clean)>500,ps_clean)#23 samples remain
ps_clean1<-prune_taxa(taxa_sums(ps_clean1)>0,ps_clean1)#546 taxa remain

meta<-read.csv("../DB/eDNA_metadata.csv",h=T)
row.names(meta)<-meta$Id_sample

sample_data(ps_clean1)<-meta

#alpha diversity analysis all microbial euk.
set.seed(572686837)
psrare<-rarefy_even_depth(ps_clean1,sample.size=540)#342 taxa in 23 samples left
ps_alpha_div <- estimate_richness(psrare, split = TRUE, measure = "Shannon")
ps_alpha_chao <- estimate_richness(ps_clean1, split = TRUE, measure = "Chao1")
ps_alpha_simp <- estimate_richness(ps_clean1, split = TRUE, measure = "Simpson")
alpha_div<-cbind(ps_alpha_div,ps_alpha_chao,ps_alpha_simp)
meta<-merge(meta,alpha_div,by=0)

sh<-ggplot(data=meta,aes(mine_open,Shannon))+geom_boxplot()+theme_bw()
ch<-ggplot(data=meta,aes(mine_open,Chao1))+geom_boxplot()+theme_bw()
si<-ggplot(data=meta,aes(mine_open,Simpson))+geom_boxplot()+theme_bw()
plot_grid(sh,ch,si,labels=c("A","B","C"))

wilcox.test(Shannon~mine_open,meta)#ns
wilcox.test(Chao1~mine_open,meta)#ns
wilcox.test(Simpson~mine_open,meta)#ns

wilcox.test(Shannon~pre_post,meta)#ns
wilcox.test(Chao1~pre_post,meta)#ns
wilcox.test(Simpson~pre_post,meta)#ns


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

fit<-lm(Shannon~Year,meta)
summary(fit)#sig positive

fit<-lm(Chao1~Year,meta)
summary(fit)#ns

fit<-lm(Simpson~Year,meta)
summary(fit)#sig. positive

DD<-subset(meta,ID_Site=="DUF_DEEP")#five samples
DD<-DD[order(DD$Year),]
DL<-subset(meta,ID_Site=="DUF_LIT")#just a single sample
DL<-DL[order(DL$Year),]
OD<-subset(meta,ID_Site=="OSI_DEEP")#13 samples
OD<-OD[order(OD$Year),]
OL<-subset(meta,ID_Site=="OSI_LIT")#four samples
OL<-OL[order(OL$Year),]

fit<-lm(Shannon~Year,OD)
summary(fit)#sig
fit<-lm(Chao1~Year,OD)
summary(fit)#ns
fit<-lm(Simpson~Year,OD)
summary(fit)#sig

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

#beta diversity analysis
beta<-transform_sample_counts(ps_clean1,function(x) x/sum(x))

tabred<-otu_table(beta)
bc<-vegdist(tabred,method="bray")
meta$Id_sample==row.names(tabred)#assume order of samples is the same
adonis_model <- adonis2(bc~mine_open,data=meta, permutations = 999)
adonis_model#sig

anos_mod<-anosim(bc,meta$mine_open,permutations=999)
summary(anos_mod)#sig

adonis_model <- adonis2(bc~pre_post,data=meta, permutations = 999)
adonis_model#marginal

anos_mod<-anosim(bc,meta$pre_post,permutations=999)
summary(anos_mod)#ns

pcoa<-cmdscale(bc,eig=TRUE)
ordiplot(pcoa)
m<-merge(meta,pcoa$points,by.x="Id_sample",by.y=0)
ggplot(m,aes(V1,V2))+geom_point(aes(color=mine_open,shape=ID_Site))+theme_classic()
ggplot(m,aes(V1,V2))+geom_point(aes(color=pre_post,shape=ID_Site))+theme_classic()

env_sub<-meta[c(16,20,23,24,47,25)]
row.names(env_sub)<-meta$Row.names
ev<-envfit(pcoa,env_sub,permutations=999)#year, cu, al, S significant
arrow.df<-data.frame(ev$vectors$arrows)
arrow.df<-arrow.df[c(1,2,3,6),]
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


#rel abundance plot of different divisions
group_ps = tax_glom(beta, taxrank="division")# tax glom to division level

barplot<-plot_bar(group_ps,x="Order",fill="division")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  facet_wrap(~ID_Site)+theme_bw()

#by trophic group
#heterotrophs
het1<-subset_taxa(ps_clean1,group=="heterotrophs")#198 taxa
set.seed(57237)
hetrare<-rarefy_even_depth(het1,sample.size=208)#139 taxa in 11 samples left
het_alpha_div <- estimate_richness(hetrare, split = TRUE, measure = "Shannon")
het_alpha_chao <- estimate_richness(het1, split = TRUE, measure = "Chao1")
het_alpha_simp <- estimate_richness(het1, split = TRUE, measure = "Simpson")

het_meta<-sample_data(het1)
het_meta<-merge(het_meta,het_alpha_div,by=0,all.x=TRUE)
row.names(het_meta)<-het_meta$Row.names
het_meta$Row.names<-NULL
het_meta<-merge(het_meta,het_alpha_chao,by=0)
row.names(het_meta)<-het_meta$Row.names
het_meta$Row.names<-NULL
het_meta<-merge(het_meta,het_alpha_simp,by=0)
row.names(het_meta)<-het_meta$Row.names
het_meta$Row.names<-NULL

wilcox.test(Shannon~mine_open,het_meta)#ns
wilcox.test(Chao1~mine_open,het_meta)#ns
wilcox.test(Simpson~mine_open,het_meta)#ns

sh<-ggplot(data=het_meta,aes(mine_open,Shannon))+geom_boxplot()+theme_bw()
ch<-ggplot(data=het_meta,aes(mine_open,Chao1))+geom_boxplot()+theme_bw()
si<-ggplot(data=het_meta,aes(mine_open,Simpson))+geom_boxplot()+theme_bw()
plot_grid(sh,ch,si,labels=c("A","B","C"))

ggplot(het_meta,aes(x=Year,y=Shannon,color=ID_Site,shape=mine_open))+
  geom_point()+theme_bw()+geom_smooth(method="lm",se=FALSE)+
  facet_wrap(~ID_Site,scales="free")

ggplot(het_meta,aes(x=Year,y=Shannon,color=ID_Site,shape=mine_open))+
  geom_point()+theme_bw()+geom_smooth(method="lm",se=FALSE)

fit<-lm(Shannon~Year,het_meta)
summary(fit)#ns
fit<-lm(Chao1~Year,het_meta)
summary(fit)#ns
fit<-lm(Simpson~Year,het_meta)
summary(fit)#ns

het1<-prune_samples(sample_sums(het1)>0,het1)
beta_het<-transform_sample_counts(het1,function(x) x/sum(x))
bc_het<-vegdist(otu_table(beta_het),method="bray")
meta_het<-sample_data(beta_het)
adonis_model <- adonis2(bc_het~meta_het$mine_open, permutations = 999)
adonis_model#sig

anos_mod<-anosim(bc_het,meta_het$mine_open,permutations=999)
summary(anos_mod)#sig

pcoa<-cmdscale(bc_het,eig=TRUE)
ordiplot(pcoa)
m<-merge(het_meta,pcoa$points,by.x="Id_sample",by.y=0)

env_sub<-het_meta[c(15,19,22,23,46,24)]
ev<-envfit(pcoa,env_sub,permutations=999)#Al, S significant, Cu marginal
arrow.df<-data.frame(ev$vectors$arrows)
arrow.df<-arrow.df[c(2,3,6),]
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

#mixotrophs
mix1<-subset_taxa(ps_clean1,group=="mixotrophs")#11 taxa
set.seed(3769)
mixrare<-rarefy_even_depth(mix1,sample.size=105)#11 taxa in 7 samples left
mix_alpha_div <- estimate_richness(mixrare, split = TRUE, measure = "Shannon")
mix_alpha_chao <- estimate_richness(mix1, split = TRUE, measure = "Chao1")
mix_alpha_simp <- estimate_richness(mix1, split = TRUE, measure = "Simpson")

mix_meta<-sample_data(mix1)
mix_meta<-merge(mix_meta,mix_alpha_div,by=0,all.x=TRUE)
row.names(mix_meta)<-mix_meta$Row.names
mix_meta$Row.names<-NULL
mix_meta<-merge(mix_meta,mix_alpha_chao,by=0)
row.names(mix_meta)<-mix_meta$Row.names
mix_meta$Row.names<-NULL
mix_meta<-merge(mix_meta,mix_alpha_simp,by=0)
row.names(mix_meta)<-mix_meta$Row.names
mix_meta$Row.names<-NULL

wilcox.test(Shannon~mine_open,mix_meta)#ns
wilcox.test(Chao1~mine_open,mix_meta)#marg.
wilcox.test(Simpson~mine_open,mix_meta)#marg.

sh<-ggplot(data=mix_meta,aes(mine_open,Shannon))+geom_boxplot()+theme_bw()
ch<-ggplot(data=mix_meta,aes(mine_open,Chao1))+geom_boxplot()+theme_bw()
si<-ggplot(data=mix_meta,aes(mine_open,Simpson))+geom_boxplot()+theme_bw()
plot_grid(sh,ch,si,labels=c("A","B","C"))

ggplot(mix_meta,aes(x=Year,y=Shannon,color=ID_Site,shape=mine_open))+
  geom_point()+theme_bw()+geom_smooth(method="lm",se=FALSE)+
  facet_wrap(~ID_Site,scales="free")

ggplot(mix_meta,aes(x=Year,y=Shannon,color=ID_Site,shape=mine_open))+
  geom_point()+theme_bw()+geom_smooth(method="lm",se=FALSE)

fit<-lm(Shannon~Year,mix_meta)
summary(fit)#ns
fit<-lm(Chao1~Year,mix_meta)
summary(fit)#ns
fit<-lm(Simpson~Year,mix_meta)
summary(fit)#sig.negative

mix1<-prune_samples(sample_sums(mix1)>0,mix1)
beta_mix<-transform_sample_counts(mix1,function(x) x/sum(x))
bc_mix<-vegdist(otu_table(beta_mix),method="bray")
meta_mix<-sample_data(beta_mix)
adonis_model <- adonis2(bc_mix~meta_mix$mine_open, permutations = 999)
adonis_model#ns

anos_mod<-anosim(bc_mix,meta_mix$mine_open,permutations=999)
summary(anos_mod)#ns

pcoa<-cmdscale(bc_mix,eig=TRUE)
ordiplot(pcoa)
m<-merge(data.frame(meta_mix),pcoa$points,by.x="Id_sample",by.y=0)

env_sub<-m[c(15,19,22,23,46,24)]
ev<-envfit(pcoa,env_sub,permutations=999)#Year marginal
#arrow.df<-data.frame(ev$vectors$arrows)
#arrow.df<-arrow.df[c(2,3,6),]
#arrowmu<-ordiArrowMul(ev)
#arrow_map = aes(xend = Dim1*arrowmu, yend = Dim2*arrowmu, x = 0, y = 0, shape = NULL, color = NULL, 
#                label = row.names(arrow.df))
#label_map = aes(x = 1.2 * Dim1*arrowmu, y = 1.2 * Dim2*arrowmu, shape = NULL, color = NULL, 
#                label = row.names(arrow.df))
#arrowhead=arrow(length=unit(0.03,"npc"))

ggplot(m,aes(V1,V2))+
  geom_point(aes(color=mine_open,shape=ID_Site))+
  theme_classic()+
  geom_hline(yintercept=0,linetype="dashed")+
  geom_vline(xintercept=0,linetype="dashed")

#phototrophs
pho1<-subset_taxa(ps_clean1,group=="phototrophs")#208 taxa
set.seed(795)
phorare<-rarefy_even_depth(pho1,sample.size=147)#119 taxa in 23 samples left
pho_alpha_div <- estimate_richness(phorare, split = TRUE, measure = "Shannon")
pho_alpha_chao <- estimate_richness(pho1, split = TRUE, measure = "Chao1")
pho_alpha_simp <- estimate_richness(pho1, split = TRUE, measure = "Simpson")
pho_meta<-sample_data(pho1)
pho_meta<-merge(pho_meta,pho_alpha_div,by=0,all.x=TRUE)
row.names(pho_meta)<-pho_meta$Row.names
pho_meta$Row.names<-NULL
pho_meta<-merge(pho_meta,pho_alpha_chao,by=0)
row.names(pho_meta)<-pho_meta$Row.names
pho_meta$Row.names<-NULL
pho_meta<-merge(pho_meta,pho_alpha_simp,by=0)
row.names(pho_meta)<-pho_meta$Row.names
pho_meta$Row.names<-NULL

wilcox.test(Shannon~mine_open,pho_meta)#ns
wilcox.test(Chao1~mine_open,pho_meta)#ns
wilcox.test(Simpson~mine_open,pho_meta)#ns

sh<-ggplot(data=pho_meta,aes(mine_open,Shannon))+geom_boxplot()+theme_bw()
ch<-ggplot(data=pho_meta,aes(mine_open,Chao1))+geom_boxplot()+theme_bw()
si<-ggplot(data=pho_meta,aes(mine_open,Simpson))+geom_boxplot()+theme_bw()
plot_grid(sh,ch,si,labels=c("A","B","C"))

fit<-lm(Shannon~Year,pho_meta)
summary(fit)#sig positive
fit<-lm(Chao1~Year,pho_meta)
summary(fit)#ns
fit<-lm(Simpson~Year,pho_meta)
summary(fit)#sig positive

pho1<-prune_samples(sample_sums(pho1)>0,pho1)
beta_pho<-transform_sample_counts(pho1,function(x) x/sum(x))
bc_pho<-vegdist(otu_table(beta_pho),method="bray")
meta_pho<-sample_data(beta_pho)
adonis_model <- adonis2(bc_pho~meta_pho$mine_open, permutations = 999)
adonis_model#sig

anos_mod<-anosim(bc_pho,meta_pho$mine_open,permutations=999)
summary(anos_mod)#sig

pcoa<-cmdscale(bc_pho,eig=TRUE)
ordiplot(pcoa)
m<-merge(data.frame(meta_pho),pcoa$points,by.x="Id_sample",by.y=0)

env_sub<-m[c(15,19,22,23,46,24)]
ev<-envfit(pcoa,env_sub,permutations=999)#Year, Cu, Al, S significant
arrow.df<-data.frame(ev$vectors$arrows)
arrow.df<-arrow.df[c(1,2,3,6),]
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

#parasites
par1<-subset_taxa(ps_clean1,group=="parasites")#55 taxa
set.seed(795)
parrare<-rarefy_even_depth(par1,sample.size=73)#44 taxa in 5 samples left
par_alpha_div <- estimate_richness(parrare, split = TRUE, measure = "Shannon")
par_alpha_chao <- estimate_richness(par1, split = TRUE, measure = "Chao1")
par_alpha_simp <- estimate_richness(par1, split = TRUE, measure = "Simpson")
par_meta<-sample_data(par1)
par_meta<-merge(par_meta,par_alpha_div,by=0,all.x=TRUE)
row.names(par_meta)<-par_meta$Row.names
par_meta$Row.names<-NULL
par_meta<-merge(par_meta,par_alpha_chao,by=0)
row.names(par_meta)<-par_meta$Row.names
par_meta$Row.names<-NULL
par_meta<-merge(par_meta,par_alpha_simp,by=0)
row.names(par_meta)<-par_meta$Row.names
par_meta$Row.names<-NULL

wilcox.test(Shannon~mine_open,par_meta)#ns
wilcox.test(Chao1~mine_open,par_meta)#ns
wilcox.test(Simpson~mine_open,par_meta)#ns

sh<-ggplot(data=par_meta,aes(mine_open,Shannon))+geom_boxplot()+theme_bw()
ch<-ggplot(data=par_meta,aes(mine_open,Chao1))+geom_boxplot()+theme_bw()
si<-ggplot(data=par_meta,aes(mine_open,Simpson))+geom_boxplot()+theme_bw()
plot_grid(sh,ch,si,labels=c("A","B","C"))

fit<-lm(Shannon~Year,par_meta)
summary(fit)#ns
fit<-lm(Chao1~Year,par_meta)
summary(fit)#ns
fit<-lm(Simpson~Year,par_meta)
summary(fit)#sig

par1<-prune_samples(sample_sums(par1)>0,par1)
beta_par<-transform_sample_counts(par1,function(x) x/sum(x))
bc_par<-vegdist(otu_table(beta_par),method="bray")
meta_par<-sample_data(beta_par)
adonis_model <- adonis2(bc_par~meta_par$mine_open, permutations = 999)
adonis_model#ns

anos_mod<-anosim(bc_par,meta_par$mine_open,permutations=999)
summary(anos_mod)#ns

pcoa<-cmdscale(bc_par,eig=TRUE)
ordiplot(pcoa)
m<-merge(data.frame(meta_pho),pcoa$points,by.x="Id_sample",by.y=0)

env_sub<-m[c(15,19,22,23,46,24)]
ev<-envfit(pcoa,env_sub,permutations=999)#ns
#arrow.df<-data.frame(ev$vectors$arrows)
#arrow.df<-arrow.df[c(1,2,3,6),]
#arrowmu<-ordiArrowMul(ev)
#arrow_map = aes(xend = Dim1*arrowmu, yend = Dim2*arrowmu, x = 0, y = 0, shape = NULL, color = NULL, 
#                label = row.names(arrow.df))
#label_map = aes(x = 1.2 * Dim1*arrowmu, y = 1.2 * Dim2*arrowmu, shape = NULL, color = NULL, 
#                label = row.names(arrow.df))
#arrowhead=arrow(length=unit(0.03,"npc"))

ggplot(m,aes(V1,V2))+
  geom_point(aes(color=mine_open,shape=ID_Site))+
  theme_classic()+
  geom_hline(yintercept=0,linetype="dashed")+
  geom_vline(xintercept=0,linetype="dashed")


#rel abundance plot of trophic groups
beta2<-tax_filter(ps_clean1,min_prevalence=2)
beta2<-transform_sample_counts(beta2,function(x) x/sum(x))

beta3<-beta2
tax_rel<-tax_table(beta2)
tax_rel1<-tax_rel[,c(1,8)]

tax_table(beta3)<-tax_rel1

group_ps = tax_glom(beta3, taxrank="group")# tax glom to group level

barplot<-plot_bar(group_ps,x="Order",fill="group")+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    facet_wrap(~ID_Site)+theme_bw()

#random forest for indicator taxa (on division level)
beta_division = tax_glom(beta, taxrank="division")
forest_mat<-as.data.frame(otu_table(beta_division))
samp_class<-data.frame(sample_data(beta_division))
forest_mat$mine<-samp_class$mine_open
tax_class<-data.frame(tax_table(beta_division))
division.rf<-randomForest(as.factor(mine)~.,data=forest_mat,importance=TRUE)
imp<-importance(division.rf)
imp<-merge(tax_class,imp,by=0)
imp<-imp[rev(order(imp[,13])),]#Gyrista and Chlorophyta divisions with the highest Gini index

#export otu table for stamp analysis
beta_otus<-data.frame(t(otu_table(beta_division)))
beta_tax<-data.frame(tax_table(beta_division))
beta_tax<-beta_tax[,c(1,2)]
beta_m<-merge(beta_tax,beta_otus,by=0)
beta_m$Row.names<-NULL
beta_m$group<-NULL
beta_samp<-data.frame(sample_data(beta))
beta_samp2<-cbind(beta_samp$Id_sample,beta_samp$Env,beta_samp$mine_open,beta_samp$pre_post)
write.table(beta_m,"../results/beta_otus.txt",quote=FALSE,sep="\t",row.names=FALSE)
write.table(beta_samp2,"../results/beta_meta_data.txt",quote=FALSE,sep="\t",row.names=FALSE)

#STAMP with groups
beta_otus<-data.frame(t(otu_table(beta)))
beta_tax<-data.frame(tax_table(beta))
beta_m<-merge(beta_tax,beta_otus,by=0)
beta_m$Row.names<-NULL
beta_m$supergroup<-NULL
beta_m$division<-NULL
beta_m$class<-NULL
beta_m$order<-NULL
beta_m$family<-NULL
beta_m$genus<-NULL
beta_m$species<-NULL
write.table(beta_m,"../results/beta_otus_group.txt",quote=FALSE,sep="\t",row.names=FALSE)

save.image(file="18S_stats.RData")
