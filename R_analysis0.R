#UpSetR, Figure 2
classiTbl2=(apply(classiTbl,2,as.numeric))
classiTbl2=as.data.frame(classiTbl2)
library(UpSetR)
windows()
upset(classiTbl2, sets = c("iC1","iC2", "iC3", "iNC", "SimpleFalk", "TaxJanis", "MockCont", "HumanHit", "MouseHit"), sets.bar.color = "#56B4E9",
      order.by = "freq", empty.intersections = "on")

#############################################################################
#Plot Age and PMD vs 16s copies, Suppl. Figure 1
Human=meta2$Source2=="Human"
conc=meta2$cor_16Spermg
sel=!is.na(conc)
conc2=conc[sel]
windows()
plotAp(meta2$PMD[sel], as.numeric(log10(conc2+1)), fit = "1", xlab = "PMD [min]", ylab = "log 16s copies/mg")
windows()
plotAp(meta2$Age[sel], as.numeric(log10(conc2+1)), fit = "1", xlab = "Age [yrs]", ylab = "log 16s copies/mg")

#############################################################################
#deskriptive statistics  with new meta, Suppl. Table 3
meta2 <- read.delim("~/R/meta2.txt")
Cls=meta2$Source2
human=Cls=="Human"
summary(meta2$cor_16Spermg[human])
mouse=Cls=="mouse"
empty=Cls=="Empty"
summary(meta2$cor_16Spermg[empty])
summary(meta2$cor_16Spermg[mouse])

summary(meta2$qPCR16s_neg[human])
summary(meta2$qPCR16s_neg[empty])
summary(meta2$qPCR16s_neg[mouse])

Cls2=meta2$Source1
PD=Cls2=="Parkinson"
CO=Cls2=="Control"
KitUKB=Cls2=="KitomeUKB"
KitQIB=Cls2=="KitomeQIB"
Buffer=Cls2=="DNA_Buffer"
Water=Cls2=="sterileWater"
GRF=Cls2=="Mouse GRF"
SPF=Cls2=="Mouse"
SterilBrain=Cls2=="Non PD Human"
summary(meta2$cor_16Spermg[PD])
summary(meta2$cor_16Spermg[CO])
summary(meta2$cor_16Spermg[GRF])
summary(meta2$cor_16Spermg[SPF])
summary(meta2$cor_16Spermg[KitUKB])
summary(meta2$cor_16Spermgue[KitQIB])
summary(meta2$cor_16Spermg[Buffer])
summary(meta2$cor_16Spermg[SterilBrain])
SW=meta2$SampleName
SW=SW==c("SW1", "SW2", "SW3")
summary(meta2$cor_16Spermg[SW])
summary(meta2$qPCR16s_neg[SW])
summary(meta2$qPCR16s_neg[PD])
summary(meta2$qPCR16s_neg[CO])
summary(meta2$qPCR16s_neg[GRF])
summary(meta2$qPCR16s_neg[SPF])
summary(meta2$qPCR16s_neg[KitUKB])
summary(meta2$qPCR16s_neg[KitQIB])
summary(meta2$qPCR16s_neg[Buffer])
summary(meta2$qPCR16s_neg[SterilBrain])
Mock=meta2$Source1
Mock=Mock=="MockPK"
summary(meta2$qPCR16sper_mgTissue[Mock])
summary(meta2$qPCR16s_neg[Mock])

#statistics
library(dunn.test)
PCR=as.vector(meta2$cor_16Spermg)
factor=as.factor(meta2$Source1)
kruskal.test(PCR, factor)
dunn.test(PCR, factor, method="bh",wrap=T)
factor2=as.factor(meta2$Source2)
kruskal.test(PCR, factor2)
dunn.test(PCR, factor2, method="bh",wrap=T)

Uncl=Hnew[, 1]=="?"
Uncl=as.data.frame(Uncl)
Mouse=as.data.frame(classiTbl[,9])
Human=as.data.frame(classiTbl[,8])
rownames(Mouse)%in%rownames(Uncl)
rownames(Human)%in%rownames(Uncl)
###########################################################################
#Stacked barplot RDP/Blast, Figure 1 A
#load files
setwd("~/R")
load("C:/Users/f/Documents/R/Contaminants.Rdata")
RDP<- read.table("RDP.txt",sep="\t",row.names=1)
#Off targets
HUMAN=as.vector(classiTbl[,8])
MOUSE=as.vector(classiTbl[,9])

#Phylumlevel
#BLAST Phylum
BlastM=as.matrix(Hnew[MOUSE,2])
BlastH=as.matrix(Hnew[HUMAN,2])
BLAST=merge(table(BlastM), table(BlastH), all=TRUE)
#RDP with filter on phylum level
RDPP=RDP[c(5,7)]
RDPM=RDPP[MOUSE,c(1,2)]
RDPH=RDPP[HUMAN,c(1,2)]

library(tidyverse)

RDPHf=RDPH%>% filter(RDPH$V8>= 0.8)
RDPHuncl=RDPH%>% filter(RDPH$V8< 0.8)
RDPMf=RDPM%>% filter(RDPM$V8>= 0.8)
RDPMuncl=RDPM%>% filter(RDPM$V8< 0.8)
RDPMOUSE=merge(sum(table(RDPMf)),sum(table(RDPMuncl)))
RDPHUMAN=merge(sum(table(RDPHf)), sum(table(RDPHuncl)))

##combine all to one table
RDPClassiPhylum<- read.delim("~/R/RDPClassiPhylum.txt", row.names=1)
RDPClassiPhylum=as.matrix(RDPClassiPhylum)

#stacked barplot
library(viridis)
library(viridisLite)
library(RColorBrewer)
windows()
barplot((RDPClassiPhylum),main="Off targets Classifier",beside=FALSE,col=brewer.pal(8,"Set3"), legend.text =TRUE, ylim=range(pretty(c(0,850))), args.legend = list(bty = "n", x = "topleft", ncol = 2, cex=0.65))

#Domainlevel
#load files
#BLAST Domain
BlastM1=as.matrix(Hnew[MOUSE,1])
BlastH1=as.matrix(Hnew[HUMAN,1])
BlAST1=merge(table(BlastM1), table(BlastH1), all=TRUE)
#RDP with filter on Domain level
RDPD=RDP[c(2,4)]
RDPM1=RDPD[MOUSE,c(1,2)]
RDPH1=RDPD[HUMAN,c(1,2)]

library(tidyverse)
RDPHf1=RDPH1%>% filter(RDPH1$V5>= 0.8)
RDPHuncl1=RDPH1%>% filter(RDPH1$V5< 0.8)
RDPMf1=RDPM1%>% filter(RDPM1$V5>= 0.8)
RDPMuncl1=RDPM1%>% filter(RDPM1$V5< 0.8)
RDPMOUSE1=merge(sum(table(RDPMf1)),sum(table(RDPMuncl1)))
RDPHUMAN1=merge(sum(table(RDPHf1)), sum(table(RDPHuncl1)))


####include more pipelines; Domain level
setwd("~/R")
dada2 <- read.delim("~/R/dada2.txt", row.names=1, comment.char="#")
mothur <- read.delim("~/R/mothur.txt", header=FALSE, row.names=1)
qiime1_sortmerna <- read.delim("~/R/qiime1_sortmerna.txt", row.names=1)
qiime1_uclust_tax <- read.delim("~/R/qiime1_uclust_tax.txt", row.names=1)

MothurH=as.matrix(mothur[HUMAN, 1])
MothurM=as.matrix(mothur[MOUSE, 1])
table(MothurH)
table(MothurM)

Dada2H=as.matrix(dada2[HUMAN,1])
Dada2M=as.matrix(dada2[MOUSE,1])
table(Dada2H)
table(Dada2M)

Qiime1H=as.matrix(qiime1_sortmerna[HUMAN, 1])
Qiime1M=as.matrix(qiime1_sortmerna[MOUSE, 1])
table(Qiime1H)
table(Qiime1M)

Qiime2H=as.matrix(qiime1_uclust_tax[HUMAN, 1])
Qiime2M=as.matrix(qiime1_uclust_tax[MOUSE, 1])
table(Qiime2H)
table(Qiime2M)

RDPDomain<- read.delim("~/R/RDPDomain.txt", row.names=1)
RDPDomain=as.matrix(RDPDomain)
library(RColorBrewer)

n <- 25
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

windows()
barplot((RDPDomain),main="Off targets Classifier Domain level",beside=FALSE,col=brewer.pal(12,"Set3"), legend.text =TRUE, ylim=range(c(0,850)), args.legend = list(bty = "n", x = "topright", ncol = 2, cex=1.0), las=3)

#2. Version
levels(RDPDomain)<- c("?", "Eukaryota", "Archaea", "Bacteria")

barplot((RDPDomain),main="Off targets Classifier Domain level",
        beside=FALSE,
        col=c("cadetblue3", "brown3", "gold1","00000000")
        , legend.text =TRUE, ylim=range(c(0,850)), args.legend = list(bty = "n", x = "topright", ncol = 1, cex=1.0), las=3)

####on phylum level

MothurH1=as.matrix(mothur[HUMAN,c(1,2)])
MothurM1=as.matrix(mothur[MOUSE, c(1,2)])
table(MothurH1)
table(MothurM1)

Dada2H1=as.matrix(dada2[HUMAN,c(1,2)])
Dada2M1=as.matrix(dada2[MOUSE, c(1,2)])
table((Dada2H1[,1]=="?") & (Dada2H1[,2]=="?"))
table((Dada2H1[,1]=="Bacteria") &(Dada2H1[,2]=="?"))
table((Dada2H1[,1]=="Archaea") &(Dada2H1[,2]=="?"))
table(Dada2H1[,2])
table((Dada2M1[,1]=="?") & (Dada2M1[,2]=="?"))
table((Dada2M1[,1]=="Bacteria") &(Dada2M1[,2]=="?"))
table((Dada2M1[,1]=="Archaea") &(Dada2M1[,2]=="?"))
table(Dada2M1[,2])

Qiime1H1=as.matrix(qiime1_sortmerna[HUMAN, 2])
Qiime1M1=as.matrix(qiime1_sortmerna[MOUSE, 2])
table(Qiime1H1)
table(Qiime1M1)

Qiime2H1=as.matrix(qiime1_uclust_tax[HUMAN, 2])
Qiime2M1=as.matrix(qiime1_uclust_tax[MOUSE, 2])
table(Qiime2H1)
table(Qiime2M1)


RDPClassiPhylum<- read.delim("~/R/RDPClassiPhylum.txt", row.names=1)
RDPClassiPhylum=as.matrix(RDPClassiPhylum) 

sort(rowSums(RDPClassiPhylum))
#use Top7, pick with Hand/Excel

library(RColorBrewer)
windows()
barplot((RDPClassiPhylum),main="Off targets Classifier Phylum level"
        ,beside=FALSE,col=brewer.pal(12,"Set3"), legend.text =TRUE, ylim=range(c(0,850)),args.legend = list(bty = "n", x = "topright", ncol = 2, cex=1.0), las=3)

#2. Version
barplot((RDPClassiPhylum),main="Off targets Classifier Phylum level"
        ,beside=FALSE,col=c("cadetblue3","coral2","darkseagreen3","deeppink4","darksalmon", "brown3","gold1","00000000"), legend.text =TRUE, ylim=range(c(0,850)),args.legend = list(bty = "n", x = "topright", ncol = 1, cex=1.0), las=3)


###check further 
BlastM.1=as.matrix(Hnew[MOUSE,])
BlastH.1=as.matrix(Hnew[HUMAN,])
SpeciesH=as.matrix(table(BlastH.1[,7]))
SpeciesM=as.matrix(table(BlastM.1[,7]))
table(Hnew[,2]=="?")


#######################################################################
#additional Zotus Classifier, Suppl. Figure 5
#check nachträgliche zotus
zotu= c("Zotu24", "Zotu593", "Zotu590", "Zotu1665", "Zotu643","Zotu1896","Zotu873","Zotu607","Zotu939","Zotu1931","Zotu34","Zotu187","Zotu231","Zotu20","Zotu56","Zotu88","Zotu95","Zotu126","Zotu13","Zotu121","Zotu113","Zotu256","Zotu259","Zotu187","Zotu231","Zotu256","Zotu259","Zotu383","Zotu817","Zotu389","Zotu525","Zotu790","Zotu849","Zotu1096","Zotu1257","Zotu1368","Zotu495","Zotu577","Zotu965","Zotu1085","Zotu1130","Zotu1206","Zotu1369","Zotu1376","Zotu1402","Zotu1535","Zotu1586","Zotu1823","Zotu1873","Zotu2331","Zotu2494","Zotu2685","Zotu2916","Zotu189","Zotu717","Zotu746","Zotu799","Zotu884","Zotu1101","Zotu1247","Zotu1321","Zotu1574","Zotu1727","Zotu1978","Zotu1982","Zotu2000","Zotu2186","Zotu2213","Zotu2315","Zotu2386","Zotu2412","Zotu2445","Zotu2698","Zotu2724","Zotu2917")

#blast domain
BlastZ=as.matrix(Hnew[zotu,1])
#RDP with filter on Domain level
RDPD=RDP[c(2,4)]
RDPZ=RDPD[zotu,c(1,2)]

library(tidyverse)
RDPZf=RDPZ%>% filter(RDPZ$V5>= 0.8)
RDPZuncl1=RDPZ%>% filter(RDPZ$V5< 0.8)

RDPZ1=merge(sum(table(RDPZf)),sum(table(RDPZuncl1)))

#mothur domain
MothurZ=as.matrix(mothur[zotu, 1])
#dada2 domain
Dada2Z=as.matrix(dada2[zotu,1])
#qiime domain
Qiime1Z=as.matrix(qiime1_sortmerna[zotu, 1])
Qiime2Z=as.matrix(qiime1_uclust_tax[zotu, 1])

##phylum
#BLAST Phylum
BlastZ1=as.matrix(Hnew[zotu,2])

#RDP with filter on phylum level
RDPZ2=RDP[c(5,7)]
RDPZZ=RDPZ2[zotu,c(1,2)]
RDPZZf=RDPZZ%>% filter(RDPZZ$V8>= 0.8)
RDPZZuncl=RDPZZ%>% filter(RDPZZ$V8< 0.8)
RDPzotu=merge(sum(table(RDPZZf)),sum(table(RDPZZuncl)))
#mothur phylum
MothurZ1=as.matrix(mothur[zotu,c(1,2)])
#dada2 phylum
Dada2Z1=as.matrix(dada2[zotu,c(1,2)])
#qiime phylum
Qiime1Z1=as.matrix(qiime1_sortmerna[zotu, 2])

Qiime2Z1=as.matrix(qiime1_uclust_tax[zotu, 2])

#combine to one file

ZotuDomain <- read.delim("~/R/ZotuDomain", row.names=1)
ZotuDomain=as.matrix(ZotuDomain)
library(RColorBrewer)
windows()
barplot((ZotuDomain),main="Additional Zotus Classifier Domain level",
        beside=FALSE,
        col=c("cadetblue3", "brown3", "gold1","00000000")
        , legend.text =TRUE, ylim=range(c(0,100)), args.legend = list(bty = "n", x = "topright", ncol = 1, cex=1.0), las=3)


ZotuPhylum <- read.delim("~/R/ZotuPhylum", row.names=1)
ZotuPhylum=as.matrix(ZotuPhylum)

sort(rowSums(ZotuPhylum))
library(RColorBrewer)
windows()
barplot((ZotuPhylum),main="Additional Zotus Classifier Phylum level"
        ,beside=FALSE,col=c("cadetblue3","coral2","darkseagreen3","deeppink4","darksalmon", "brown3","gold1","steelblue","00000000"), legend.text =TRUE, ylim=range(c(0,120)),args.legend = list(bty = "n", x = "topright", ncol = 2, cex=1.0), las=3)

###################################################################
# check Mock for off target
#set mocks 
Mock=M[, c("NR", "Nr4", "Nr6", "Nr8", "Nr10", "Nr20")]
Genus=as.matrix(Hnew[,6])
                write.table(Genus, "Genus.txt")
                write.table(Mock, "Mock.txt")
                #merge Genus and Mock by hand in Excel as TaxaMock
                TaxaMock <- read.delim("~/R/TaxaMock.txt", row.names=1)                
                TaxaMock=as.matrix(TaxaMock)
                sort(table(TaxaMock[,1]), decreasing = T)
                
                #off targets in Mock
                MockH=as.matrix(TaxaMock[HUMAN,])
                MockM=as.matrix(TaxaMock[MOUSE,])
                write.table(MockH, "MockH.txt")
                write.table(MockM, "MockM.txt")
                
                