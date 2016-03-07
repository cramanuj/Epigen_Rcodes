#!/usr/bin/env Rscript

#################################################################################################################
## Preprocess Gibbs et al gene expression and methylation data
##
## LAST UPDATE: March 7, 2016
## AUTHOR: Chaitanya Acharya
##
## NOTES:
## 1) This script makes use of Series_Matrix_File deposited on the Gene Expression Omnibus server
## 2) This script needs methylation batch information text file (since the one accessed from GEO is corrupted)
## 3) Covariate-adjusted gene expression and methylation data are deposited as text files
## 4) The final files are sliced into fragments of 100 genes/CpG sites and deposited into multiple subdirectories
##	on which other scripts are run
## 5) Please DO NOT DELETE the .SOFT FILES. They will be used else where in the analysis.
#################################################################################################################

lib.list = c("GEOquery","JAGUAR")

for(i in 1:length(lib.list)){
	if(any(installed.packages()[,1]==lib.list[i])){
		library(lib.list[i],character.only=T)}else{
			source("http://bioconductor.org/biocLite.R");
			biocLite(lib.list[i]);
			require(lib.list[i],character.only=T)};
}

rm(list=ls())
cwd = getwd()
if(!interactive()) pdf("plots_gibbs_etal.pdf",paper="USr")

############################################################
### Download series matrix files from GEO
##
############################################################

geo_data = getGEO("GSE15745",GSEMatrix=T,destdir = getwd())
class(geo_data)
objects(geo_data)

###########################################################
## Preprocessing gene expression data
##
###########################################################

# Expression data
GeneExp = geo_data$`GSE15745-GPL6104_series_matrix.txt.gz`
dim(GeneExp)

# Phenotype data
pheno_gene_exp = as(GeneExp@phenoData, "data.frame")
dim(pheno_gene_exp)
pheno_gene_exp = pheno_gene_exp[,c(1,8,10:14)]
colnames(pheno_gene_exp) = c("Patient_ID","Tissue","Gender","Age","PMI","TissueBank","Batch")
pheno_gene_exp$Patient_ID = gsub("-CRBLM-mRNA","",pheno_gene_exp$Patient_ID)
pheno_gene_exp$Patient_ID = gsub("-FCTX-mRNA","",pheno_gene_exp$Patient_ID)
pheno_gene_exp$Patient_ID = gsub("-PONS-mRNA","",pheno_gene_exp$Patient_ID)
pheno_gene_exp$Patient_ID = gsub("-TCTX-mRNA","",pheno_gene_exp$Patient_ID)
pheno_gene_exp$Patient_ID = as.factor(pheno_gene_exp$Patient_ID)
pheno_gene_exp$Tissue = as.factor(gsub("Human Brain Tissue: ","",pheno_gene_exp$Tissue))
pheno_gene_exp$Gender = as.factor(gsub("gender: ","",pheno_gene_exp$Gender))
pheno_gene_exp$Age = as.numeric(gsub("age: ","",pheno_gene_exp$Age))
pheno_gene_exp$PMI = as.numeric(gsub("pmi: ","",pheno_gene_exp$PMI))
pheno_gene_exp$TissueBank = as.factor(gsub("tissuebank: ","",pheno_gene_exp$TissueBank))
pheno_gene_exp$Batch = as.factor(gsub("prep_hyb_batch: ","",pheno_gene_exp$Batch))

# Remove sample from different ethnicity
out = c(grep("UMARY-927",pheno_gene_exp$Patient_ID),
		grep("UMARY-4545",pheno_gene_exp$Patient_ID),
		grep("JHU-1344",pheno_gene_exp$Patient_ID))
pheno_gene_exp = pheno_gene_exp[-(out),]
head(pheno_gene_exp)

# PCA
GeneExp = as.data.frame(exprs(GeneExp))
GeneExp = GeneExp[,-(out)]
GeneExp[GeneExp<0]<- 1
GeneExp = as.matrix(log2(GeneExp))
data_pca = prcomp(t(GeneExp))
pca = (t(GeneExp) %*% data_pca$rotation)[,1:3]
b1 = pheno_gene_exp$Tissue; b2 = pheno_gene_exp$TissueBank; b3 = pheno_gene_exp$Batch

# PCA plots
#pdf("pca_gene_exp.pdf")
par(mfrow=c(2,2))
screeplot(data_pca,10,type="lines",main="Principal Components")
plot(pca,col= factor(as.numeric(factor(b1))),cex=1,pch=16,main="Gibbs et al gene expression data");
legend("bottomright",legend=levels(factor(b1)),fill=levels(factor(as.numeric(factor(b1)))),title="Brain Region",bty="n",cex=0.6);
plot(pca,col= factor(as.numeric(factor(b2))),cex=1,pch=16,main="Gibbs et al gene expression data");
legend("bottomright",legend=levels(factor(b2)),fill=levels(factor(as.numeric(factor(b2)))),title="Tissue Bank",bty="n",cex=0.6);
plot(pca,col= factor(as.numeric(factor(b3))),cex=1,pch=16,main="Gibbs et al gene expression data");
legend("bottomright",legend=levels(factor(b3)),fill=levels(factor(as.numeric(factor(b3)))),title="Hybridization Batch",bty="n",cex=0.6);
par(mfrow=c(1,1))
#dev.off()

# Outlier identification
cutoff = median(sort(pca[,1])) + 1.5*IQR(sort(pca[,1]))
outliers = which(pca[,1]>cutoff,arr.ind=T)
GeneExp = GeneExp[,-outliers]
pheno_gene_exp = pheno_gene_exp[-outliers,]

# Gene annotation
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
gene_coords_bed = genes[,c(4,5,6,1)]
write.table(gene_coords_bed,"gene_coords.bed",sep="\t",col.names=F,row.names=F,quote=F)

gene_probes = rownames(GeneExp)[rownames(GeneExp)%in%genes$Probes]
genes = subset(genes,genes$Probes %in% gene_probes)
GeneExp = GeneExp[match(gene_probes,rownames(GeneExp)),]

# Missing samples
miss_samples = names( which(table(pheno_gene_exp$Patient_ID) < 4) )
out = c("UMARY-927","UMARY-4545","JHU-1344")
miss_samples = miss_samples[!(miss_samples %in% out)]
batch_split = split(pheno_gene_exp,pheno_gene_exp$Tissue)
dataF = list()
brain_regions = c("CRBLM","FCTX","PONS","TCTX")
for(i in 1:length(batch_split)){
	print(i)
	b = batch_split[[i]]
	d = GeneExp[,rownames(b)]
	colnames(d) = as.character(b$Patient_ID)
	data_pca = prcomp(t(d))
	pca = (t(d) %*% data_pca$rotation)[,1:50]
	mod=model.matrix(~.,cbind(b[,-c(1:2)],pca))
	fit = lm(t(d)~ mod)
	d = data.frame(t(apply(t(fit$residuals),1,scale)),check.names=F);
	colnames(d) = as.character(b$Patient_ID)
	new_samples = miss_samples[!(miss_samples %in% colnames(d))]
	NA_col = data.frame(matrix(NA,nrow(d),length(new_samples)))
	colnames(NA_col) = new_samples
	d = cbind(d,NA_col)
	dataF[[i]] = d[,order(colnames(d))]
}
GeneExp = do.call("cbind",dataF)
dim(GeneExp)
EXP_SAMPLES = colnames(GeneExp)[1:145]
write.table(GeneExp,"GeneExp_Gibbs.txt",sep="\t",col.names=NA,quote=F)

###########################################################
## Preprocessing methylation expression data
##
###########################################################

# Methylation data
MethExp = geo_data$'GSE15745-GPL8490_series_matrix.txt.gz'
dim(MethExp)
MethExp = as.data.frame(exprs(MethExp))
MethExp = na.omit(MethExp)
dim(MethExp)

# CpG Island annotation
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
cpg_probes = rownames(MethExp)[rownames(MethExp)%in%cpg_annot$ID]
MethExp = MethExp[match(cpg_probes,rownames(MethExp)),]
dim(MethExp)
cpg = subset(cpg_annot,cpg_annot$ID %in% cpg_probes)
dim(cpg)
cpg_bed_file = cpg[,c(2,9,10,1)]
write.table(cpg_bed_file,"cpg_coords.bed",sep="\t",col.names=F,row.names=F,quote=F)

# Phenotype data
pheno_meth_exp = read.delim("GSE15745_batchInfo_methyData.txt",header=T,check.names=F)
dim(pheno_meth_exp)
pheno_meth_exp$Patient_ID = gsub("-CRBLM-CpG","",pheno_meth_exp$Patient_ID)
pheno_meth_exp$Patient_ID = gsub("-FCTX-CpG","",pheno_meth_exp$Patient_ID)
pheno_meth_exp$Patient_ID = gsub("-PONS-CpG","",pheno_meth_exp$Patient_ID)
pheno_meth_exp$Patient_ID = gsub("-TCTX-CpG","",pheno_meth_exp$Patient_ID)
pheno_meth_exp$Patient_ID = as.factor(pheno_meth_exp$Patient_ID)
out = c(grep("UMARY-927",pheno_meth_exp$Patient_ID),
grep("UMARY-4545",pheno_meth_exp$Patient_ID),
grep("JHU-1344",pheno_meth_exp$Patient_ID))
pheno_meth_exp = pheno_meth_exp[-(out),]
head(pheno_meth_exp)

MethExp = MethExp[,colnames(MethExp) %in% pheno_meth_exp$GEO_ID]
dim(MethExp)

# PCA analysis
data_pca = prcomp(t(as.matrix(MethExp)))
pca = (t(MethExp) %*% data_pca$rotation)[,1:3]
b1 = as.factor(pheno_meth_exp$TissueType); b2 = as.factor(pheno_meth_exp$TissueBank); b3 = as.factor(pheno_meth_exp$Batch)

par(mfrow=c(2,2))
screeplot(data_pca,10,type="lines",main="Principal Components")
plot(pca,col= factor(as.numeric(factor(b1))),cex=1,pch=16,main="Gibbs et al methylation data");
legend("topright",legend=levels(factor(b1)),fill=levels(factor(as.numeric(factor(b1)))),title="Brain Region",bty="n",cex=0.6);
plot(pca,col= factor(as.numeric(factor(b2))),cex=1,pch=16,main="Gibbs et al methylation data");
legend("topright",legend=levels(factor(b2)),fill=levels(factor(as.numeric(factor(b2)))),title="Tissue Bank",bty="n",cex=0.6);
plot(pca,col= factor(as.numeric(factor(b3))),cex=1,pch=16,main="Gibbs et al methylation data");
legend("topright",legend=levels(factor(b3)),fill=levels(factor(as.numeric(factor(b3)))),title="Hybridization Batch",bty="n",cex=0.6);
par(mfrow=c(1,1))

batch_split = split(pheno_meth_exp,pheno_meth_exp$TissueType)

# Missing samples
for(i in 1:length(batch_split)){
	print(i)
	b = batch_split[[i]]
	d = MethExp[,as.character(b$GEO_ID)]
	colnames(d) = as.character(b$Patient_ID)
	mod=model.matrix(~.,cbind(b[,-c(1:3)]))
	fit = lm(t(d)~ mod)
	d = data.frame(t(apply(t(fit$residuals),1,scale)),check.names=F);
	colnames(d) = as.character(b$Patient_ID)
	miss_samples = EXP_SAMPLES[!(EXP_SAMPLES %in% b$Patient_ID)]
	NA_col = data.frame(matrix(NA,nrow(d),length(miss_samples)))
	colnames(NA_col) = miss_samples
	d = cbind(d,NA_col)
	d = d[,colnames(d) %in% EXP_SAMPLES]
	dataF[[i]] = d[,order(names(d))]
}
MethExp = do.call("cbind",dataF)
dim(MethExp)
write.table(MethExp,"MethExp_Gibbs.txt",sep="\t",col.names=NA,quote=F)

# Merge gene expression and methylation data based on Entrez Gene ID

tmp = na.omit( merge(genes,cpg,by.x=names(genes)[3],by.y=names(cpg)[5],sort=F) )
dim(tmp); head(tmp)
tab = aggregate(ID ~ Probes,tmp,"length")
dim(tab)
table(tab$ID)

GeneExp_new = GeneExp[match(tmp$Probes,rownames(GeneExp)),]
dim(GeneExp_new)
MethExp_new = MethExp[match(tmp$ID,rownames(MethExp)),]
dim(MethExp_new)

# Slice the methylation expression data
jaguar_slice(MethExp_new,size=100,path=cwd)
dirs = dir()[file.info(dir())$isdir]
dirs = dirs[grep("dir",dirs)]
for(i in 1:length(dirs)){
	setwd(paste(cwd,"/dir",i,sep=""));
	system("cp GeneExp_matrix.txt MethExp_matrix.txt")
}
setwd(cwd)
# Slice the gene expression data
jaguar_slice(GeneExp_new,size=100,path=cwd)

## End of plots
if(!interactive()) dev.off()

## SessionInfo
sessionInfo()

