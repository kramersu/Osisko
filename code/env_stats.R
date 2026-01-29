library(tidyr)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(vegan)
library(ffaframework)
library(Kendall)


setwd("D:\\Marta")
data<-read.csv("cores_meta_data.csv",h=T)

#many metals are correlated, sort into clusters
metals<-data[c(21:63)]
metals<-na.omit(metals)
mm<-cor(metals)
hie_ave<-hclust(dist(mm),method="average")
dend<-as.dendrogram(hie_ave)
plot(dend)
#five clusters
#cluster1: As, Ag, Ni, Tl, Sb, Cd, Pb, Fe, Co, Zn, Cu, Se #increases in Osisko: choose Cu as a representative
#cluster2: S, Yb, Tm, Lu, Dy, Gd, Tb, Er, Y, Ho, Nd, Eu, Sm #increases, then decresases in Osisko: choose S as representative
#cluster3: Pr, La, Ce #downward trajectory, decreases, then recovers in Osisko Deep: choose Pr
#cluster4: Sc, Cr, Ba, Mg, V, Sr, Ti, Al, K, Th Th #drecreases, then recovers in Osisko: choose Al
#cluster5: Ca, Mn, Pt, Na, U# noisy:, choose Ca

env<-data[c(5,6,8,10,12,14,21,26,24,48,25)]
row.names(env)<-data$Id_sample
env_cur<-na.omit(env)#60 sampels remain
env_cur.pca<-pca(env_cur,scale=TRUE)
summary(env_cur.pca)
ev_cur<-env_cur.pca$CA$eig
ev_cur[ev_cur>mean(ev_cur)]#first 4 PCs explain more than mean
PCs<-scores(env_cur.pca,display="sites",choices=c(1:4))
PC_loadings<-scores(env_cur.pca, choices = 1:3, display = "species", scaling = 0)
arrow_scl<-scores(env_cur.pca,display='species',scaling=1)
arrow.df<-data.frame(arrow_scl,labels=rownames(arrow_scl))
arrow_map = aes(xend = PC1, yend = PC2, x = 0, y = 0, shape = NULL, color = NULL, 
                label = labels)
label_map = aes(x = 1.2 * PC1, y = 1.2 * PC2, shape = NULL, color = NULL, 
                label = labels)
arrowhead=arrow(length=unit(0.03,"npc"))

PCs<-merge(PCs,data,by.x=0,by.y="Id_sample")
ggplot(data=PCs, aes(PC1,PC2,col=ID_Site,shape=mine_open))+geom_point()+
  geom_segment(arrow_map,size=0.5,data=arrow.df,arrow=arrowhead)+
  geom_text(label_map,size=4,data=arrow.df)+theme_classic()+
  geom_hline(yintercept=0,linetype="dashed")+geom_vline(xintercept=0,linetype="dashed")

#qPCR plots
data2<-subset(data,Year>1899)
ggplot(data=data2,aes(Year,eCACO4))+geom_bar(stat="identity",fill="darkgreen")+facet_grid(ID_Site~.,scales="free")+theme_bw()
ggplot(data=data2,aes(Year,eSAVI2))+geom_bar(stat="identity",fill="lightblue")+facet_grid(ID_Site~.,scales="free")+theme_bw()
ggplot(data=data2,aes(Year,ePEFL1))+geom_bar(stat="identity",fill="darkred")+facet_grid(ID_Site~.,scales="free")+theme_bw()
ggplot(data=data2,aes(Year,eESLU1))+geom_bar(stat="identity",fill="orange")+facet_grid(ID_Site~.,scales="free")+theme_bw()

#Trends and change points: OD
OD<-subset(data,ID_Site=="OSI_DEEP")
MannKendall(OD$eSAVI2)#ns
MannKendall(OD$eCACO4)#ns
MannKendall(OD$ePEFL1)#ns
MannKendall(OD$eESLU1)#ns
MannKendall(OD$Cu)#sig, increase with time
MannKendall(OD$S)#sig, decrease with time
MannKendall(OD$Pr)#sig, decrease with time
MannKendall(OD$Al)#ns
MannKendall(OD$Ca)#sig, decrease with time

#Change points
OD2<-OD[1:36,]#remove NA
eda_pettitt_test(OD2$Cu,OD2$Year)#change point: 1931, p=0.023
eda_mks_test(OD2$Cu,OD2$Year)#ns
eda_pettitt_test(OD2$S,OD2$Year)#change point: 1963, p=0.032
eda_mks_test(OD2$S,OD2$Year)#ns
eda_pettitt_test(OD2$Pr,OD2$Year)#change point: 1960
eda_mks_test(OD2$Pr,OD2$Year)#ns
eda_pettitt_test(OD2$Ca,OD2$Year)#change point: 1966, p=0.014
eda_mks_test(OD2$Ca,OD2$Year)#sig., 1983
eda_pettitt_test(OD2$Al,OD2$Year)#change point: 1935, p=0.003
eda_mks_test(OD2$Al,OD2$Year)#sig., 1960, 1955, 1942
eda_pettitt_test(OD2$Ca,OD2$Year)#change point: 1966, p=0.014
eda_mks_test(OD2$Ca,OD2$Year)#sig., 1983

#Trends and change points: OL
OL<-subset(data,ID_Site=="OSI_LIT")
MannKendall(OL$eSAVI2)#ns
MannKendall(OL$eCACO4)#marg. positive
MannKendall(OL$ePEFL1)#marg.positive
MannKendall(OL$eESLU1)#ns
MannKendall(OL$Cu)#sig increase with time
MannKendall(OL$S)#ns
MannKendall(OL$Pr)#sig, decrease with time
MannKendall(OL$Al)#marg. negative
MannKendall(OL$Ca)#ns

#Change points
eda_pettitt_test(OL$Cu,OL$Year)#ns
eda_mks_test(OL$Cu,OL$Year)#ns
eda_pettitt_test(OL$S,OL$Year)#ns
eda_mks_test(OL$S,OL$Year)#ns
eda_pettitt_test(OL$Pr,OL$Year)#ns
eda_mks_test(OL$Pr,OL$Year)#ns
eda_pettitt_test(OL$Al,OL$Year)#ns
eda_mks_test(OL$Al,OL$Year)#ns
eda_pettitt_test(OL$Ca,OL$Year)#ns
eda_mks_test(OL$Ca,OL$Year)#ns

#Trends and change points: DD
DD<-subset(data,ID_Site=="DUF_DEEP")
MannKendall(DD$eSAVI2)#ns
MannKendall(DD$eCACO4)#ns
MannKendall(DD$ePEFL1)#ns
MannKendall(DD$eESLU1)#ns
MannKendall(DD$Cu)#ns
MannKendall(DD$S)#ns
MannKendall(DD$Pr)#marg decrease over time
MannKendall(DD$Al)#ns
MannKendall(DD$Ca)#ns

eda_pettitt_test(DD$Cu,DD$Year)#ns
eda_mks_test(DD$Cu,DD$Year)#ns
eda_pettitt_test(DD$S,DD$Year)#ns
eda_mks_test(DD$S,DD$Year)#ns
eda_pettitt_test(DD$Al,DD$Year)#ns
eda_mks_test(DD$Al,DD$Year)#ns
eda_pettitt_test(DD$Pr,DD$Year)#ns
eda_mks_test(DD$Pr,DD$Year)#ns
eda_pettitt_test(DD$Ca,DD$Year)#ns
eda_mks_test(DD$Ca,DD$Year)#ns

#Trends and change points: DL
DL<-subset(data,ID_Site=="DUF_LIT")
MannKendall(DL$eSAVI2)#ns
MannKendall(DL$eCACO4)#ns
MannKendall(DL$ePEFL1)#ns
MannKendall(DL$eESLU1)#ns
MannKendall(DL$Cu)#sig increase with time
MannKendall(DL$S)#ns
MannKendall(DL$Pr)#marg decrease with time
MannKendall(DL$Al)#sig decrease over time
MannKendall(DL$Ca)#ns

#Change points
DL2<-DL[1:7,]
eda_pettitt_test(DL2$Cu,DL2$Year)#ns
eda_mks_test(DL2$Cu,DL2$Year)#ns
eda_pettitt_test(DL2$S,DL2$Year)#ns
eda_mks_test(DL2$S,DL2$Year)#ns
eda_pettitt_test(DL2$Al,DL2$Year)#ns
eda_mks_test(DL2$Al,DL2$Year)#ns
eda_pettitt_test(DL2$Pr,DL2$Year)#ns
eda_mks_test(DL2$Pr,DL2$Year)#ns
eda_pettitt_test(DL2$Ca,DL2$Year)#ns
eda_mks_test(DL2$Ca,DL2$Year)#ns

#metal graphs
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
n = 4
cols = gg_color_hue(n)

ggplot(data=data,aes(Year,Cu,group=ID_Site))+
  geom_line(aes(col=ID_Site),size=1)+
  geom_point(aes(col=ID_Site))+
  theme_bw()+xlim(1900,2025)+
  geom_vline(xintercept=1927,size=1.2,linetype="dashed", col="grey")+
  geom_vline(xintercept=1931,size=1.2,col="#00BFC4",linetype="dashed")

ggplot(data=data,aes(Year,S,group=ID_Site))+
  geom_line(aes(col=ID_Site),size=1)+
  geom_point(aes(col=ID_Site))+
  theme_bw()+xlim(1900,2025)+
  geom_vline(xintercept=1927,size=1.2,linetype="dashed", col="grey")+
  geom_vline(xintercept=1963,size=1.2,col="#00BFC4",linetype="dashed")

ggplot(data=data,aes(Year,Al,group=ID_Site))+
  geom_line(aes(col=ID_Site),size=1)+
  geom_point(aes(col=ID_Site))+
  theme_bw()+xlim(1900,2025)+
  geom_vline(xintercept=1927,size=1.2,linetype="dashed", col="grey")+
  geom_vline(xintercept=1935,size=1.2,col="#00BFC4",linetype="dashed")

ggplot(data=data,aes(Year,Pr,group=ID_Site))+
  geom_line(aes(col=ID_Site),size=1)+
  geom_point(aes(col=ID_Site))+
  theme_bw()+xlim(1900,2025)+
  geom_vline(xintercept=1927,size=1.2,linetype="dashed", col="grey")+
  geom_vline(xintercept=1960,size=1.2,col="#00BFC4",linetype="dashed")

ggplot(data=data,aes(Year,Ca,group=ID_Site))+
  geom_line(aes(col=ID_Site),size=1)+
  geom_point(aes(col=ID_Site))+theme_bw()+xlim(1900,2025)+
  geom_vline(xintercept=1927,size=1.2,linetype="dashed", col="grey")+
  geom_vline(xintercept=1966,size=1.2,col="#00BFC4",linetype="dashed")+
  geom_vline(xintercept=1983,size=1.2,col="#00BFC4",linetype="dashed")

#qPCR tests
wilcox.test(data$eSAVI2~data$mine_open)#ns
wilcox.test(data$eCACO4~data$mine_open)#sig less in yes category, p=0.031
boxplot(data$eCACO4~data$mine_open)#less in yes than in no
wilcox.test(data$ePEFL1~data$mine_open)#marginal
boxplot(data$ePEFL1~data$mine_open)#less in yes than in no
wilcox.test(data$eESLU1~data$mine_open)#sig, p=0.023
boxplot(data$eESLU1~data$mine_open)#less in yes than in no

wilcox.test(data$eSAVI2~data$pre_post)#ns
wilcox.test(data$eCACO4~data$pre_post)#marginal
boxplot(data$eCACO4~data$pre_post)#less in yes than in no
wilcox.test(data$ePEFL1~data$pre_post)#ns
wilcox.test(data$eESLU1~data$pre_post)#ns


wilcox.test(data$eSAVI2~data$Site)#ns
wilcox.test(data$eCACO4~data$Site)#marginal
wilcox.test(data$ePEFL1~data$Site)#ns
wilcox.test(data$eESLU1~data$Site)#sig

#qPCR: water compared to top sediment
data<-read.csv("w_s_qPCR.csv")
ggplot(data)+
  geom_bar(aes(x=Core,y=eSAV1),stat="identity",fill="skyblue")+
  geom_errorbar(aes(x=Core,ymin=eSAV1-eSAV1_error,ymax=eSAV1+eSAV1_error),width=0.25,size=1)+
  facet_grid(Class~.,scales="free")+theme_bw()

ggplot(data)+
  geom_bar(aes(x=Core,y=eCACO4),stat="identity",fill="darkgreen")+
  geom_errorbar(aes(x=Core,ymin=eCACO4-eCACO4_error,ymax=eCACO4+eCACO4_error),width=0.25,size=1)+
  facet_grid(Class~.,scales="free")+theme_bw()

ggplot(data)+
  geom_bar(aes(x=Core,y=ePELF1),stat="identity",fill="darkred")+
  geom_errorbar(aes(x=Core,ymin=ePELF1-ePELF1_error,ymax=ePELF1+ePELF1_error),width=0.25,size=1)+
  facet_grid(Class~.,scales="free")+theme_bw()

ggplot(data)+
  geom_bar(aes(x=Core,y=eESLU1),stat="identity",fill="orange")+
  geom_errorbar(aes(x=Core,ymin=eESLU1-eESLU1_error,ymax=eESLU1+eESLU1_error),width=0.25,size=1)+
  facet_grid(Class~.,scales="free")+theme_bw()+
  ylim(0,140000)