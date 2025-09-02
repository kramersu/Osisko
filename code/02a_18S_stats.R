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

path<-c("D:\\marta/18S/results")
load(paste0(path,"/ps_clean_18S.RData"))
class_exclude<-c("unassigned","Incertae Sedis")#all eukaryotic identified to the class level

ps_clean1<-subset_taxa(ps_clean,!(class %in% class_exclude))#512 ASVs remain

ps_clean2<-prune_samples(sample_sums(ps_clean1)>0,ps_clean1)

meta<-read.csv("18S/DB/eDNA_metadata.csv",h=T)
row.names(meta)<-meta$Name

sample_data(ps_clean2)<-meta

#alpha diversity analysis
set.seed(572686837)
psrare<-rarefy_even_depth(ps_clean2,sample.size=1217)#287 taxa in 37 samples left
ps_alpha_div <- estimate_richness(psrare, split = TRUE, measure = "Shannon")
meta<-merge(meta,ps_alpha_div,by=0)

boxplot(Shannon~Lake+Site,meta)
fit<-aov(Shannon~Lake*Site,meta)
anova(fit)#ns

fit<-aov(Shannon~Combo,meta)
anova(fit)#ns

fit<-lm(Shannon~Min_Depth,meta)
summary(fit)#sig. positive correlation

ggplot(meta,aes(x=Min_Depth,y=Shannon,color=Lake,shape=Site))+geom_point()+geom_smooth(method="lm")+theme_bw()#increasing diversity with increasing depth

fit<-lm(Shannon~Treatment,meta)
summary(fit)#ns
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
summary(fit)#marginal

fit<-lm(Shannon~Year,DD)
summary(fit)#ns
plot(Shannon~Year,DD,type='b')

fit<-lm(Shannon~Year,DL)
summary(fit)#ns
plot(Shannon~Year,DL,type='b')

fit<-lm(Shannon~Treatment,DL)
summary(fit)#ns

#beta diversity analysis
ps_clean3<-tax_filter(ps_clean2,min_prevalence=3,prev_detection_threshold=5)#66 taxa remain
beta<-transform_sample_counts(ps_clean3,function(x) x/sum(x))

tabred<-otu_table(beta)
hellinger<-decostand(tabred,method="hellinger")
adonis_model <- adonis2(hellinger~Combo,data=meta, permutations = 999)
adonis_model#marginal

adonis_model <- adonis(hellinger~Lake*Site,data=meta, permutations = 999)
adonis_model$aov.tab#Lake marginal

adonis_model <- adonis(tabred~Lake*Min_Depth,data=meta, permutations = 999)
adonis_model$aov.tab#min depth significant, lake marginal

adonis_model <- adonis(hellinger~Treatment,data=meta, permutations = 999)
adonis_model$aov.tab#marginal

dis<-vegdist(tabred,method="bray")
pcoa<-cmdscale(dis,eig=TRUE)
ordiplot(pcoa)
m<-merge(meta,pcoa$points,by.x="Row.names",by.y=0)
m<-m[order(m$Combo,m$Min_Depth),]
ggplot(m,aes(V1,V2,label=Name2))+geom_point(aes(color=Combo,size=Min_Depth,shape=Treatment))+theme_classic()+geom_text(hjust=1,vjust=1)
ggplot(m,aes(V1,V2,label=Name2))+geom_point(aes(color=Combo,size=Min_Depth))+theme_classic()+geom_path(aes(group=Combo))

#Individual cores
OD<-subset_samples(beta,Combo=="Osisko_Deep")
dis<-vegdist(otu_table(OD),method="bray")
pcoa<-cmdscale(dis,eig=TRUE)
ordiplot(pcoa)
m<-merge(meta,pcoa$points,by.x="Row.names",by.y=0)
m<-m[order(m$Name2),]
ggplot(m,aes(V1,V2,label=Name2))+geom_point(aes(size=Min_Depth),col="blue")+theme_classic()+geom_text(hjust=2,vjust=2)
sample_names(OD)<-sample_data(OD)$Name2
p1<-plot_bar(OD,fill="class",title="Osisko Deep")
  
OL<-subset_samples(beta,Combo=="Osisko_Littoral")
dis<-vegdist(otu_table(OL),method="bray")
pcoa<-cmdscale(dis,eig=TRUE)
ordiplot(pcoa)
m<-merge(meta,pcoa$points,by.x="Row.names",by.y=0)
ggplot(m,aes(V1,V2,label=Name2))+geom_point(aes(size=Min_Depth),col="lightblue")+theme_classic()+geom_text(hjust=1,vjust=1)
sample_names(OL)<-sample_data(OL)$Name2
p2<-plot_bar(OL,fill="class",title="Osisko Littoral")
  
DD<-subset_samples(beta,Combo=="Dufay_Deep")
dis<-vegdist(otu_table(DD),method="bray")
pcoa<-cmdscale(dis,eig=TRUE)
ordiplot(pcoa)
m<-merge(meta,pcoa$points,by.x="Row.names",by.y=0)
ggplot(m,aes(V1,V2,label=Name2))+geom_point(aes(size=Min_Depth),col="green")+theme_classic()+geom_text(hjust=1,vjust=1)
sample_names(DD)<-sample_data(DD)$Name2
p3<-plot_bar(DD,fill="class",title="Dufay Deep")

DL<-subset_samples(beta,Combo=="Dufay_Littoral")
dis<-vegdist(otu_table(DL),method="bray")
pcoa<-cmdscale(dis,eig=TRUE)
ordiplot(pcoa)
m<-merge(meta,pcoa$points,by.x="Row.names",by.y=0)
ggplot(m,aes(V1,V2,label=Name2))+geom_point(aes(size=Min_Depth,shape=Treatment),col="darkgreen")+theme_classic()+geom_text(hjust=1,vjust=1)
sample_names(DL)<-sample_data(DL)$Name2
p4<-plot_bar(DL,fill="class",title="Dufay Littoral")

ggarrange(p1,p2,p3,p4,labels=c("A","B","C","D"),ncol=2,nrow=2)

#Titan analysis
D<-subset_samples(ps_clean3,sample_data(ps_clean3)$Lake=="Dufay")
Drare<-rarefy_even_depth(D,sample.size=93487)

D_clean<-tax_filter(Drare,min_prevalence=3,prev_detection_threshold=3)

txa_obj<-otu_table(D_clean)
env<-sample_data(D_clean)$Year

to<-titan(env,txa_obj,pur.cut=0.95,rel.cut=0.95,numPerm=1000)#no taxa pass the filter
plot_taxa_ridges(to,z1=FALSE)
plot_sumz(to,filter=FALSE,cumfrq=FALSE,xmin=1920,xmax=2022)

O<-subset_samples(ps_clean3,sample_data(ps_clean3)$Lake=="Osisko")
#remove samples with low counts
O2<-subset_samples(O,Name !="OSID_2023_A_11a" & Name !="OSID_2023_A_17a" & Name != "OSID_2023_A_21a")
Orare<-rarefy_even_depth(O2,sample.size=13949)

O_clean<-tax_filter(Orare,min_prevalence=3,prev_detection_threshold=3)

txa_obj<-otu_table(O_clean)
env<-sample_data(O_clean)$Min_Depth

to<-titan(env,txa_obj,pur.cut=0.95,rel.cut=0.95,numPerm=1000)#no taxa pass the filter
plot_sumz(to,filter=FALSE,cumfrq=FALSE)
plot_taxa_ridges(to)#doesn't work because no taxa pass filter
