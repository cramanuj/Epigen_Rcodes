#!/usr/bin/env Rscript

#################################################################################################################
## Analyzing the final results
##
## LAST UPDATE: March 7, 2016
## AUTHOR: Chaitanya Acharya
##
## NOTES:
## 1) This script aggregates all the text files from all the directories and computes significant relationships
##
#################################################################################################################
lib.list = c("qqman","qvalue","ggplot2","reshape2","GEOquery")
for(i in 1:length(lib.list)){
	if(any(installed.packages()[,1]==lib.list[i])){
		library(lib.list[i],character.only=T)}else{
			source("http://bioconductor.org/biocLite.R");
			biocLite(lib.list[i]);
			require(lib.list[i],character.only=T)};
}

## Gene expression annotation
gene_annot = parseGEO("GPL6104.soft")
gene_annot = Table(gene_annot@dataTable)
genes = gene_annot[gene_annot$Chromosome!="X",]
genes = genes[genes$Chromosome!="Y",]
genes$Chromosome = as.numeric(as.character(genes$Chromosome))
genes = subset(genes,!is.na(genes$Chromosome))
genes$Entrez_Gene_ID = as.numeric(genes$Entrez_Gene_ID)
genes = genes[,c(1,6,9,18,20)]
colnames(genes)[1:2]=c("Probes","Gene")
genes$start = unlist(lapply(strsplit(as.character(genes$Probe_Coordinates), "\\-"), "[", 1))
genes$end = unlist(lapply(strsplit(as.character(genes$Probe_Coordinates), "\\-"), "[", 2))
genes=na.omit(genes)
genes = genes[,-5]
dim(genes)
head(genes)

## Methylation expression annotation
cpg_annot = parseGEO("GPL8490.soft")
cpg_annot = Table(cpg_annot@dataTable)
cpg_annot = cpg_annot[,c(1,9,19,20,22,23,29,30,34,35)]
cpg_annot = cpg_annot[cpg_annot$CPG_ISLAND!="FALSE",]
cpg_annot$Gene_ID = gsub("GeneID:","",cpg_annot$Gene_ID)
cpg_annot$RANGE_START = as.numeric(cpg_annot$RANGE_START)
cpg_annot$RANGE_END = as.numeric(cpg_annot$RANGE_END)
cpg_annot = cpg_annot[cpg_annot$Chr != "X",]
cpg_annot = cpg_annot[cpg_annot$Chr != "Y",]
dim(cpg_annot)
head(cpg_annot)

#
############################################
#

dirs = dir()[file.info(dir())$isdir]
dirs = dirs[grep("dir",dirs)]

## Read all the results
fname1 = "JT_cisAnalysis.txt";
fname2 = "DELTA_cisAnalysis.txt";
fname3 = "TBT_eQTL_cis.txt";
fname4 = "TBT_mQTL_cis.txt";
fname5 = "JAG_results.txt";
fname6 = "ANOVA_cisAnalysis.txt";

jt_eqtl = delta_eqtl = tbt_eqtl = tbt_mqtl = jag_eqtl = anova_eqtl = vector(mode="list",length=length(dirs))

for(i in 1:length(dirs)){
	print(i)
	setwd(paste("dir",i,sep=""));
	jt_eqtl[[i]] = read.delim(fname1,header=T)
	delta_eqtl[[i]] = read.delim(fname2,header=T)
	tbt_eqtl[[i]] = read.delim(fname3,header=T)
	tbt_mqtl[[i]] =read.delim(fname4,header=T)
	jag_eqtl[[i]] = read.delim(fname5,header=T)
	anova_eqtl[[i]] = read.delim(fname6,header=T)
	setwd("../")
}

### Joint test
jt_out = do.call("rbind",jt_eqtl)
jt_out$qval = qvalue(jt_out$Upsi)$qval
jt_out = jt_out[,c(1,2,3,9)]
names(jt_out)[4]=c("Upsi")
jt_sig = jt_out[jt_out$Upsi<=0.05,]
jt_triplets = unique(paste(jt_sig$Genes,jt_sig$CpG,jt_sig$SNP,sep="-"))
length(jt_triplets)

fin = jt_out

## ANOVA
anova_out = do.call("rbind",anova_eqtl)
anova_eqtl = cbind(anova_out[,c(1:3)],apply(anova_out[,c(4:7)],2,function(x) qvalue(x)$qval))
apply(anova_eqtl[,-c(1:3)],2,function(x) sum(x<=0.05))
anova_eqtl$MinP = apply(anova_eqtl[,4:7],1,min)
sum(anova_eqtl$MinP<=(0.05/4))
anova_sig = anova_eqtl[anova_eqtl$MinP<=(0.05/4),]
anova_triplets = unique(paste(anova_sig$Gene,anova_sig$CpG,anova_sig$SNP,sep="-"))
length(anova_triplets)

round(100 * sum(anova_triplets %in% jt_triplets)/length(anova_triplets),0)	## Triplets common between JT and ANOVA

## DELTA
delta_out = do.call("rbind",delta_eqtl)
delta_out$qval = qvalue(delta_out$Udel)$qval

## TBT
tbt_eqtl_out = do.call("rbind",tbt_eqtl)
tbt_eqtl_out$Pairs = paste(tbt_eqtl_out$Gene,tbt_eqtl_out$SNP,sep="-")
tbt_eqtl_out = tbt_eqtl_out[!(duplicated(tbt_eqtl_out$Pairs)),-7]
tbt1 = cbind(tbt_eqtl_out[,c(1:2)],apply(tbt_eqtl_out[,-c(1:2)],2,function(x) qvalue(x)$qval))
apply(tbt1[,-c(1:2)],2,function(x) sum(x<=0.05))
tbt1$MinP = apply(tbt1[,-c(1:2)],1,min)
sum(tbt1$MinP<=(0.05/4))
tbt1_sig = tbt1[tbt1$MinP<=(0.05/4),]

tbt_mqtl_out = do.call("rbind",tbt_mqtl)
tbt_mqtl_out$Pairs = paste(tbt_mqtl_out$CpG,tbt_mqtl_out$SNP,sep="-")
tbt_mqtl_out = tbt_mqtl_out[!(duplicated(tbt_mqtl_out$Pairs)),-7]
tbt2 = cbind(tbt_mqtl_out[,c(1:2)],apply(tbt_mqtl_out[,-c(1:2)],2,function(x) qvalue(x)$qval))
apply(tbt2[,-c(1:2)],2,function(x) sum(x<=0.05))
tbt2_sig = tbt2[tbt2$MinP<=(0.05/4),]
sum(tbt2$MinP<=(0.05/4))

genes_tmp = merge(genes,tbt1,by.x=names(genes)[1],by.y=names(tbt1)[1],sort=F)
genes_tmp = genes_tmp[genes_tmp$MinP<=(0.05/4),]
dim(genes_tmp)
cpg_temp = merge(cpg_annot,tbt2,by.x=names(cpg_annot)[1],by.y=names(tbt2)[1],sort=F)
cpg_temp = cpg_temp[cpg_temp$MinP<=(0.05/4),]
dim(cpg_temp)

t = na.omit(cpg_temp[match(cpg_temp$Gene_ID,unique(genes_tmp$Entrez_Gene_ID)),])
t$Gene_ID = as.numeric(t$Gene_ID)
t = t[order(t$Gene_ID),]
t1 = merge(t,genes_tmp,by.x=names(t)[5],by.y=names(genes_tmp)[3])

temp1 =data.frame("Gene"=as.character(t1[,17]),"CpG"=as.character(t1[,2]),"SNP"=as.character(t1$SNP.x))
dim(temp1); head(temp1)
temp2 =data.frame("Gene"=as.character(t1[,17]),"CpG"=as.character(t1[,2]),"SNP"=as.character(t1$SNP.y))
dim(temp2); head(temp2)
temp3 = rbind(temp1,temp2)
tbt_triplets = unique(paste(temp3$Gene,temp3$CpG,temp3$SNP,sep="-"))

## JAGUAR
jag_eqtl_out = do.call("rbind",jag_eqtl)

jag1_eqtl = jag_eqtl_out[,c(2,1,3)]
jag1_eqtl$Pairs = paste(jag1_eqtl$Gene,jag1_eqtl$SNP,sep="-")
jag1_eqtl = jag1_eqtl[!(duplicated(jag1_eqtl$Pairs)),-4]
jag1_eqtl$qvalue = qvalue(jag1_eqtl$pvalue)$qval
jag1_sig = jag1_eqtl[jag1_eqtl$qvalue<=0.05,]
genes_tmp = merge(genes,jag1_sig,by.x=names(genes)[1],by.y=names(jag1_sig)[1],sort=F)
dim(genes_tmp)

jag_mqtl_out = jag_eqtl_out[,c(4,1,5)]
names(jag_mqtl_out)[3]=c("pvalue")
jag_mqtl_out$Pairs = paste(jag_mqtl_out$CpG,jag_mqtl_out$SNP,sep="-")
jag_mqtl_out = jag_mqtl_out[!(duplicated(jag_mqtl_out$Pairs)),-4]
jag_mqtl_out$qvalue = qvalue(jag_mqtl_out$pvalue)$qval
dim(jag_mqtl_out)
jag_mqtl_sig = jag_mqtl_out[jag_mqtl_out$qvalue<=0.05,]
cpg_temp = merge(cpg_annot,jag_mqtl_sig,by.x=names(cpg_annot)[1],by.y=names(jag_mqtl_sig)[1],sort=F)
dim(cpg_temp); head(cpg_temp)

tt = na.omit(cpg_temp[match(cpg_temp$Gene_ID,unique(genes_tmp$Entrez_Gene_ID)),])
tt$Gene_ID = as.numeric(tt$Gene_ID)
tt = tt[order(tt$Gene_ID),]
t2 = merge(tt,genes_tmp,by.x=names(tt)[5],by.y=names(genes_tmp)[3])

temp4 =data.frame("Gene"=as.character(t2[,14]),"CpG"=as.character(t2[,2]),"SNP"=as.character(t2$SNP.x))
dim(temp4); head(temp4)
temp5 =data.frame("Gene"=as.character(t2[,14]),"CpG"=as.character(t2[,2]),"SNP"=as.character(t2$SNP.y))
dim(temp5); head(temp5)
temp6 = rbind(temp4,temp5)
temp6 = temp6[!(duplicated(temp6)),]

jag_triplets = paste(temp6$Gene,temp6$CpG,temp6$SNP,sep="-")

round(100 * sum(tbt_triplets %in% jt_triplets)/length(tbt_triplets),0)		## Triplets common between TBT and JT
round(100 * sum(anova_triplets %in% jt_triplets)/length(anova_triplets),0)	## Triplets common between ANOVA and JT
round(100 * sum(jag_triplets %in% jt_triplets)/length(jag_triplets),0)		## Triplets common between JAG and JT

round(100 * sum(tbt_triplets %in% anova_triplets)/length(anova_triplets),0)	## Triplets common between ANOVA and TBT
round(100 * sum(tbt_triplets %in% jag_triplets)/length(tbt_triplets),0)		## Triplets common between TBT and ANOVA
round(100 * sum(jag_triplets %in% anova_triplets)/length(anova_triplets),0)	## Triplets common between JAG and ANOVA


## Plots

pdf("jt_plots.pdf")
par(mfrow=c(1,2))
plot(-log10(jt_out$Upsi),-log10(anova_eqtl$MinP),ylab="-log10 pvalue from ANOVA",xlab="-log10 pvalue from Joint Test",pch=".")
plot(-log10(jt_out$Upsi),-log10(delta_out$qval),ylab="-log10 pvalue from Variance Component Test",xlab="-log10 pvalue from Joint Test",pch=".")
par(mfrow=c(1,1))
dev.off()

