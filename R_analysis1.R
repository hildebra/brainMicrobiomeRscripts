#script for brain microbiome analysis in Bedarf et al.



M=read.table("C:/Users/falkh/OneDrive/science/data/LT/Brain4n/OTU.txt",header=TRUE,"\t",row.names=1,comment.char="")
H=read.table("C:/Users/falkh/OneDrive/science/data/LT/Brain4n/hiera_BLAST.txt",header=TRUE,"\t",row.names=1,as.is=TRUE)
Met=read.table("C:/Users/falkh/OneDrive/science/data/LT/Brain4n/meta/meta.txt",header=TRUE,"\t",row.names=1,as.is=TRUE,comment.char="")

M=t(t(M)/colSums(M))
#colSums(M)
#create new unclass 16S counts
Hnew = H[dimnames(M)[[1]],]
Hnew[is.na(Hnew)]="?"
uncl=Hnew[,1]=="?" | Hnew[,1]=="Eukaryota"
names(uncl)=dimnames(M)[[1]]
#added to account for all the extensive tests..

#dimnames(M)[[1]]%in%rmZotus

#based on pre done filtering now..
load(file=paste0("C:/Users/falkh/OneDrive/science/data/LT/Brain4n/","Conti1.Rdata"))
rmZotus=names(which(rowSums(classiTbl)>0))
table(uncl[dimnames(M)[[1]]%in%rmZotus])
uncl[dimnames(M)[[1]]%in%rmZotus] = TRUE

#save this list to remove these zotus
rmZotus2=names(uncl[uncl])
length(rmZotus2)
save(rmZotus2,file=paste0("C:/Users/falkh/OneDrive/science/data/LT/Brain4n/","Conti2.Rdata"))


Met=Met[dimnames(M)[[2]],]
Met=as.data.frame(Met)
newC=Met$`qPCR16sper_mgTissue`-  colSums(M[uncl,]) * Met$`qPCR16sper_mgTissue`
plot(sort(newC))

#MetN=cbind(Met,`cor_16Spermg`=newC)
MetN = Met;
MetN[,"cor_16Spermg"]=newC
write.table(MetN,file="C:/Users/falkh/OneDrive/science/data/LT/Brain4n/meta/meta2.txt",sep="\t",quote=FALSE)

#H=metaInfo$hierachy$mat

Mp=as.matrix(t(t(M)/colSums(M)))

plot(colSums(M[uncl,]),log="")

mockS=c("NR2","Nr4","Nr6","Nr8","Nr10","Nr20")
mockE=c("Nr19","Nr9","Nr7","Nr5","Nr3","Nr1")
#head(Mp[,mockE])

contT=names(sort(rowSums(Mp[,mockE],na.rm=TRUE)/6,TRUE)[1:20])

colSums(Mp[contT[1:6],mockE])

colSums(M[,mockS]>0)




if (0){ #load and reduce OTU set to off-targets only
  load(file=paste0(opt$data$path,"Conti1.Rdata"))
  mousHum = names(which(rowSums(classiTbl[,c("HumanHit","MouseHit")])>0))
  mous = names(which(rowSums(classiTbl[,c("MouseHit"),drop=FALSE])>0))
  human = rowSums(classiTbl[,c("HumanHit"),drop=FALSE])>0
  names(which(mousHum))
  library(ape)
  dna=read.FASTA("C:/Users/hildebra/OneDrive/science/data/LT/Brain4n/otus.fa")
  
  
  write.FASTA(dna[mousHum],"C:/Users/hildebra/OneDrive/science/data/LT/Brain4n/offtarget.fa")
  H=read.table("C:/Users/falkh/OneDrive/science/data/LT/Brain4n/hiera_BLAST.txt",header=TRUE,"\t",row.names=1,as.is=TRUE)
  CdiffZotus= rownames(H[mousHum,])[!is.na(H[mousHum,7]) & H[mousHum,7]=="Clostridioides difficile"]
  
  
  BL=read.table("C:/Users/hildebra/OneDrive/science/data/LT/Brain4n/blast/tax.1.blast")
  subB=BL[BL[,1]%in%CdiffZotus,]
  subB1=subB[subB[,3]>95,]
  mean((subB1[,7]+subB1[,8])/2)
}





if (0){
  inP="C:/Users/hildebra/OneDrive/science/data/LT/BrainEGA16S/"
  M=read.table(paste0(inP,"OTU.txt"),header=TRUE,row.names=1,as.is=TRUE)
  tax=read.table(paste0(inP,"hiera_BLAST.txt"),header=TRUE,row.names=1,as.is=TRUE,sep="\t")
  
  contTax=unique(c("Afipia","Aquabacterium","Asticcacaulis","Aurantimonas","Beijerinckia","Bosea","Bradyrhizobium","Brevundimonas","Caulobacter","Craurococcus","Devosia","Hoeflea","Mesorhizobium","Methylobacterium","Novosphingobium","Ochrobactrum","Paracoccus","Pedomicrobium","Phyllobacterium","Rhizobium","Roseomonas","Sphingobium","Sphingomonas","Sphingopyxis","Acidovorax","Azoarcus","Azospira","Burkholderia","Comamonas","Cupriavidus","Curvibacter","Delftia","Duganella","Herbaspirillum","Janthinobacterium","Kingella","Leptothrix","Limnobacter","Massilia","Methylophilus","Methyloversatilis","Oxalobacter","Pelomonas","Polaromonas","Ralstonia","Schlegelella","Sulfuritalea","Undibacterium","Variovorax","Acinetobacter","Enhydrobacter","Enterobacter","Escherichia","Nevskia","Pseudomonas","Pseudoxanthomonas","Psychrobacter","Stenotrophomonas","Xanthomonas","Aeromicrobium","Arthrobacter","Beutenbergia","Brevibacterium","Corynebacterium","Curtobacterium","Dietzia","Geodermatophilus","Janibacter","Kocuria","Microbacterium","Micrococcus","Microlunatus","Patulibacter","Propionibacterium","Rhodococcus","Tsukamurella","Abiotrophia","Bacillus","Brevibacillus","Brochothrix","Facklamia","Paenibacillus","Streptococcus","Chryseobacterium","Dyadobacter","Flavobacterium","Hydrotalea","Niastella","Olivibacter","Pedobacter","Wautersiella","ThermusDeinococcus",
                   "Actinomyces","Corynebacterium","Arthrobacter","Rothia","Propionibacterium","Atopobium","Sediminibacterium","Porphyromonas","Prevotella","Chryseobacterium","Capnocytophaga","Chryseobacterium","Flavobacterium","Pedobacter","TM7","Bacillus","Geobacillus","Brevibacillus","Paenibacillus","Staphylococcus","Abiotrophia","Granulicatella","Enterococcus","Lactobacillus","Streptococcus","Clostridium","Coprococcus","Anaerococcus","Dialister","Megasphaera","Veillonella","Fusobacterium","Leptotrichia","Brevundimonas","Afipia","Bradyrhizobium","Devosia","Methylobacterium","Mesorhizobium","Phyllobacterium","Rhizobium","Methylobacterium","Phyllobacterium","Roseomonas","Novosphingobium","Sphingobium","Sphingomonas","Achromobacter","Achromobacter","Burkholderia","Acidovorax","Comamonas","Curvibacter","Pelomonas","Cupriavidus","Duganella","Herbaspirillum","Janthinobacterium","Massilia","Oxalobacter","Ralstonia","Leptothrix","kingella","Neisseria","Escherichia","Haemophilus","Acinetobacter","Enhydrobacter","Pseudomonas","Stenotrophomonas","Xanthomonas"))
  
  Hnew=Hold = tax
  #Hnew = matrix(NA,dim(M)[[1]],8)
  #dimnames(Hnew) = list(dimnames(M)[[1]],dimnames(Hold)[[2]])
  
  #Hnew[dimnames(Hold)[[1]],] = Hold
  Hnew[is.na(Hnew)]="?"
  
  contL2 = Hnew[,7] %in%  contTax |  Hnew[,6] %in%  contTax | Hnew[,1]=="Eukaryota"
  sum(M[contL2,])/sum(M)
  sum(M[!contL2 & Hnew[,5] != "?",])/sum(M)
  
  
  
}





#second part analysis
M=read.table("OTU.txt")
cls=read.table("metadata.txt")

Mat=list(Normed=t(t(M)/rowSums(M),Feature=M)


if (0){#test contaminants
	
	if (0){#create contaminant lists..
		library(decontam);
		#stores all the different methods
		classiTbl = matrix(NA,dim(mat$Feature)[[1]],9)
		dimnames(classiTbl) = list(dimnames(mat$Feature)[[1]] , c("iC1","iC2","iC3","iNC","SimpleFalk","TaxJanis","MockCont","HumanHit","MouseHit"))
		
		#conc=fit$`16sper_mgTissue`
		conc=fit$qPCR16sper_mgTissue
		neg1=fit$Source2
		neg1[neg1=="mouse"]="Human"
		neg=neg1=="Empty"
		sel=conc>0 & !is.na(conc)# & (fit$ExtrRun==2 & !is.na(fit$ExtrRun))
		#sum(is.na(t(mat$Normed[,sel])))
		#classical isContimant in 3 variants:
		#iC1=isContaminant(t(mat$Feature[,sel]), conc=conc[sel],neg=neg[sel], method="frequency")#,batch=fit$ExtrRun[sel])
		#iC2=isContaminant(t(mat$Feature[,sel]), conc=conc[sel],neg=neg[sel], method="prevalence")#,batch=fit$ExtrRun[sel])
		iC1=isContaminant(t(mat$Feature), neg=neg, method="prevalence",threshold=0.1)#,batch=fit$ExtrRun[sel])
		iC2=isContaminant(t(mat$Feature[,sel]), conc=conc[sel],neg=neg[sel], method="frequency",threshold=0.1)#,batch=fit$ExtrRun[sel])
		iC3=isContaminant(t(mat$Feature[,sel]), conc=conc[sel],neg=neg[sel], method="either",threshold=0.1)#,batch=fit$ExtrRun[sel])
		
		classiTbl[,"iC1"]= iC1$contaminant
		classiTbl[,"iC2"]= iC2$contaminant
		classiTbl[,"iC3"]= iC3$contaminant
		
		#harder isNotCont algo, 0.5 is default
		iNC1=isNotContaminant(t(mat$Normed),neg = neg,threshold=0.5)
		classiTbl[,"iNC"]= !iNC1
		
		#ad hoc var of contaminant detection from me
		M=mat$Normed
		contList=rowSums(M[,fit$Source2=="Empty"]>0.01)>=1
		classiTbl[,"SimpleFalk"]= contList
		
		
		#and check for Janis taxa
		contTax=unique(c("Afipia","Aquabacterium","Asticcacaulis","Aurantimonas","Beijerinckia","Bosea","Bradyrhizobium","Brevundimonas","Caulobacter","Craurococcus","Devosia","Hoeflea","Mesorhizobium","Methylobacterium","Novosphingobium","Ochrobactrum","Paracoccus","Pedomicrobium","Phyllobacterium","Rhizobium","Roseomonas","Sphingobium","Sphingomonas","Sphingopyxis","Acidovorax","Azoarcus","Azospira","Burkholderia","Comamonas","Cupriavidus","Curvibacter","Delftia","Duganella","Herbaspirillum","Janthinobacterium","Kingella","Leptothrix","Limnobacter","Massilia","Methylophilus","Methyloversatilis","Oxalobacter","Pelomonas","Polaromonas","Ralstonia","Schlegelella","Sulfuritalea","Undibacterium","Variovorax","Acinetobacter","Enhydrobacter","Enterobacter","Escherichia","Nevskia","Pseudomonas","Pseudoxanthomonas","Psychrobacter","Stenotrophomonas","Xanthomonas","Aeromicrobium","Arthrobacter","Beutenbergia","Brevibacterium","Corynebacterium","Curtobacterium","Dietzia","Geodermatophilus","Janibacter","Kocuria","Microbacterium","Micrococcus","Microlunatus","Patulibacter","Propionibacterium","Rhodococcus","Tsukamurella","Abiotrophia","Bacillus","Brevibacillus","Brochothrix","Facklamia","Paenibacillus","Streptococcus","Chryseobacterium","Dyadobacter","Flavobacterium","Hydrotalea","Niastella","Olivibacter","Pedobacter","Wautersiella","ThermusDeinococcus",
						 "Actinomyces","Corynebacterium","Arthrobacter","Rothia","Propionibacterium","Atopobium","Sediminibacterium","Porphyromonas","Prevotella","Chryseobacterium","Capnocytophaga","Chryseobacterium","Flavobacterium","Pedobacter","TM7","Bacillus","Geobacillus","Brevibacillus","Paenibacillus","Staphylococcus","Abiotrophia","Granulicatella","Enterococcus","Lactobacillus","Streptococcus","Clostridium","Coprococcus","Anaerococcus","Dialister","Megasphaera","Veillonella","Fusobacterium","Leptotrichia","Brevundimonas","Afipia","Bradyrhizobium","Devosia","Methylobacterium","Mesorhizobium","Phyllobacterium","Rhizobium","Methylobacterium","Phyllobacterium","Roseomonas","Novosphingobium","Sphingobium","Sphingomonas","Achromobacter","Achromobacter","Burkholderia","Acidovorax","Comamonas","Curvibacter","Pelomonas","Cupriavidus","Duganella","Herbaspirillum","Janthinobacterium","Massilia","Oxalobacter","Ralstonia","Leptothrix","kingella","Neisseria","Escherichia","Haemophilus","Acinetobacter","Enhydrobacter","Pseudomonas","Stenotrophomonas","Xanthomonas"))
		
		Hold = metaInfo$hierachy$mat
		Hnew = matrix(NA,dim(M)[[1]],8)
		dimnames(Hnew) = list(dimnames(M)[[1]],dimnames(Hold)[[2]])
		
		Hnew[dimnames(Hold)[[1]],] = Hold
		Hnew[is.na(Hnew)]="?"

		contL2 = Hnew[,6] %in%  contTax |  Hnew[,5] %in%  contTax | Hnew[,1]=="Eukaryota"
		classiTbl[,"TaxJanis"]= contL2
		
		#Zotus in mock community
		mockZotus=c("Zotu10","Zotu6","Zotu15","Zotu9","Zotu28","Zotu32","Zotu2408")
		mockConts=rowSums(M[,fit$Source2=="Mock"])>0
		mockConts[mockZotus] = FALSE
		classiTbl[,"MockCont"]= mockConts
		
		
		michum=read.table(file=paste0(opt$data$path,"/host.txt"),sep="\t",row.names = NULL,header=FALSE)
		humanzotus=as.character(michum[michum[,2]=="human",1])
		mousezotus=as.character(michum[michum[,2]=="mouse",1])
		plot(sort(michum[,3]))
		head(sort(table(c(mousezotus,humanzotus))))
		
		table(Hnew[humanzotus,5])
		table(Hnew[mousezotus,5])
		
		sum(rownames(metaInfo$hierachy$mat)%in% humanzotus)
		sum(!rownames(metaInfo$hierachy$mat)%in% humanzotus)
		
		
		classiTbl[,"HumanHit"]=FALSE
		classiTbl[humanzotus,"HumanHit"]=TRUE
		classiTbl[,"MouseHit"]=FALSE
		classiTbl[mousezotus,"MouseHit"]=TRUE
		
		save(metaInfo,Hnew,classiTbl,fit,clto,mat,file=paste0(opt$data$path,"Contaminants.Rdata"))
		load(file=paste0(opt$data$path,"Contaminants.Rdata"))
		save(classiTbl,file=paste0(opt$data$path,"Conti1.Rdata"))
		
		#plot frequency based
		plot(1,1,type="n",ylim=c(-10,-1),xlim=c(6,11),xlab="16S concentration",ylab="OTU abundance")
		cnt=1
		for (xx in contL4){
		  ssel=(M[xx,])>0
		  oo=order(conc[ssel])
		  points(log(conc[ssel][oo]),log(M[xx,ssel][oo]),col=cnt,pch=18)
		  cnt=cnt+1
		}
		#plot prevalence based
		
		df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
							contaminant=contamdf.prev$contaminant)
		ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
		  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")
	}
}

if (0){ plot contaminant abundances
  load(file=paste0(opt$data$path,"Conti1.Rdata"))
  M = mat$Normed
	contM=cbind(colSums(M[classiTbl[,1],]))
	for (x in 2:dim(classiTbl)[2]){
	  contM=cbind(contM,colSums(M[classiTbl[,x],]))
	}
	contM[,8]=contM[,8]*100
	clsE=as.character(fit$ExtractionMethod)
	clsE[!clsE%in%c("PK1","PK2")]=NA
	clsE=factor(clsE)
	t.test(contM[,8]~clsE)
	par(mfrow=c(2,1))
	abund=fit$qPCR16s_well
	boxplot(contM[,8]~clsE,main="P=8e-8",
			ylab="Off-target read %",xlab="Extraction Method")
	sub=!is.na(clsE)&contM[,8]>0.1
	plotAp(contM[sub,8],abund[sub],log="",fit="1",col=1,cex.lab=1,
		 ylab="bacterial load (16S qPCR, uncorrected)",xlab="Off-target read %")
	cor.test(fit$qPCR16s_well[sub],contM[sub,8],method="spearman")
	
	dimnames(contM)[[2]] = dimnames(classiTbl)[[2]]
#test		for unclassified
	uncl=metaInfo$hierachy$mat[,1]=="?"
	plot(colSums(M[uncl,]))
	
	cls = fit$Source2
	cls[cls=="Mock"]=NA
	cls[cls=="mouse"]="Human"
	cls[cls=="Empty"&fit$Run=="3rd"] = NA
	
	cls = factor(cls)
	
	#newC=fit$`16s_neg`-  colSums(M[uncl,]) * fit$`16sper_mgTissue`
	#newC[!is.na(cls)& cls=="Empty"]=fit$`16s_neg`[!is.na(cls)&cls=="Empty"]
	newC = fit$cor_16Spermg
	
	par(mfrow=c(1,2)) 
	ylims=c(0,1250)
	boxplot((fit$`qPCR16sper_mgTissue`)~cls,ylim=ylims,col=c("red","blue"),ylab="log10 16S copies/mg Tissue",main="Uncorrected\nP=0.0006")
	boxplot((newC+1)~cls,col=c("red","blue"),ylim=ylims,ylab="",main="Corrected for \na) off-targets\nb) contaminants\nP=0.07")
	kruskal.test(newC~cls)
	kruskal.test(fit$`16s_neg`~cls)

}

if (0){ #plot overlap in categories
	if (0 && ijk==6){ #only works with unfiltered OTUs
	normV2 = as.numeric(fit[["qPCR16sper_mgTissue"]])
	names(normV2) = names(fit[[1]])
	normV = log10(normV2+1)
	load(file=paste0(opt$data$path,"Conti1.Rdata"))
	M = mat$Normed[dimnames(classiTbl)[[1]],]
	#contM=cbind(colSums(M))
	contM=matrix(NA,dim(M)[[2]],dim(classiTbl)[2])
	for (x in 1:dim(classiTbl)[2]){
	# contM=cbind(contM,colSums(M[classiTbl[,x],]))
	contM[,x]=colSums(M[classiTbl[,x],])
	}
	dimnames(contM)[[2]] = dimnames(classiTbl)[[2]]
	dimnames(contM)[[1]] = dimnames(M)[[2]]
	contM=cbind(contM, `remains`=colSums(M[rowSums(classiTbl)==0,]),`total`=colSums(M[,]) )
	contCnt=contM*normV
	#unique entries
	mousHum = rowSums(classiTbl[,c("HumanHit","MouseHit")])>0
	conts = rowSums(classiTbl)>0 & !mousHum
	rem=!conts & !mousHum
	contM2=cbind(`remains`=colSums(M[rem,]),`Contaminants`=colSums(M[conts,]), `off-Target`=colSums(M[mousHum,]))
	contCnt2=contM2*normV2

	library(UpSetR)
	classiTbl2=(apply(classiTbl,2,as.numeric))
	classiTbl2=as.data.frame(classiTbl2)
	retUps = upset(classiTbl2, sets = c("iC1","iC2", "iC3", "iNC", "SimpleFalk", "TaxJanis", "MockCont", "HumanHit", "MouseHit"), 
		sets.bar.color = "#56B4E9",order.by = "freq", empty.intersections = "on")

	exclSC=c(paste("Br",seq(79,86),sep=""),paste("Nr",seq(21,26),sep=""),paste("Nr",seq(1,10),sep=""),"Nr20","NR2")
	exclB = !is.na(fit[[normTag]]) & !names(fit[[1]])%in%exclSC

	S2f = factor(fit$Source2)
	levels(S2f) = levels(S2f)[c(1,3,2,4)]

}


#waffle plot to visualize absolute abundance:
if (0){
	sel=!is.na(clto$disc$empty)
	M2=t(t(mat$Normed)*fit$qPCR16sper_mgTissue)
	empM = matrix(NA,dim(mat$Normed)[1],5)
	colnames(empM) = levels(clto$disc$empty)
	rownames(empM) = rownames(M2)
	lvls=levels(clto$disc$empty)
	for (xx in lvls){
	  sel=!is.na(clto$disc$empty) & clto$disc$empty == xx
	  empM[,xx] = rowMeans(M2[,sel,drop=FALSE],na.rm=TRUE)
	}
	
	colSums(empM)
	
	isConta=rowSums(classiTbl[,1:7])>0
	isOff=rowSums(classiTbl[,8:9])>0
	
	par(mfrow=c(1,1))
	barplot( rbind(colSums(empM[isConta,]),  colSums(empM[isOff,])) )
	
	
	
	sub1=sub2=ori = rowMeans(mat$Normed[,sel,drop=FALSE])

	prepie=function(v,n=20){
	  ord=order(v,decreasing=TRUE)
	  
	  vr=c(v[ord[seq(n)]],`other`=sum(v[ord[(n+1):length(v)]],na.rm=TRUE))
	  names(vr)=shortenHierNames(names(vr))
	  vr[vr<0]=0
	  vr
	}
	for (p in c(1,2,3,5)){
	  waffle(prepie(empM[,p],10)*1,colors=rainbow(11),rows=25,keep=FALSE,title =lvls[p],
			 xlab="1 square == 1 bacteria / 100 ul")
	  pdf(file=paste0("Waffle",p,".pdf"),7,5)
	  waffle(prepie(empM[,p],10)*1,colors=rainbow(11),rows=25,keep=FALSE,title =lvls[p],
		   xlab="1 square == 1 bacteria / 100 ul")
	  dev.off()
	}
	#head(classiTbl)
	sub1=ori[rowSums(classiTbl[,1:7])==0]
	sumOff=sum(ori[rowSums(classiTbl[,8:9,drop=FALSE])>0])
	sumConta=sum(ori[rowSums(classiTbl[,1:7,drop=FALSE])>0])
	sub2=ori[rowSums(classiTbl[,1:9])==0]
	
	
	
	
	library(waffle)
	#library(extrafont)
	waffle(prepie(ori)*1000,colors=rainbow(21),rows=10,size=0.5,title="No contaminant control",flip=FALSE)
	
	
	waffle(prepie(ori)*1000*(1-sumConta),colors=rainbow(21),rows=10,size=0.5,title="No contaminant control",flip=TRUE)
	waffle(prepie(ori)*1000*(1-sumOff),colors=rainbow(21),rows=10,size=0.5,title="No contaminant control",flip=TRUE)
	
	par(mfrow=c(1,3))
	pie(prepie(ori),main="No contaminant control")
	pie(prepie(sub1),main="Contaminant removed")
	pie(prepie(sub2),main="Contaminant + off-target\nremoved")
	
waffle(savings/392, rows=7, size=0.5, 
   colors=c("#c7d4b6", "#a3aabd", "#a0d0de", "#97b5cf"), 
   title="Average Household Savings Each Year", 
   xlab="1 square == $392")

}