#!/usr/bin/env Rscript

#################################################################################################################
## Apply the variance component score test to measure Genotype x Methylation x Tissue effect or in other words
##  tissue-specific methylation effect on genotype
##
## LAST UPDATE: March 7, 2016
## AUTHOR: Chaitanya Acharya
##
## NOTES:
## 1) This script should be run on each partition of gene expression/CpG methylation data
## 2) The results from the analysis are stored in text file labeled "DELTA_cisAnalysis.txt"
## 3) Please change the path to the genotype file  and bed files before running this script
#################################################################################################################

lib.list = c("data.table","Rcpp","lme4","plyr")
for(i in 1:length(lib.list)){
	if(any(installed.packages()[,1]==lib.list[i])){
		library(lib.list[i],character.only=T)}else{
			source("http://bioconductor.org/biocLite.R");
			biocLite(lib.list[i]);
			require(lib.list[i],character.only=T)};
}

## Rcpp script for the analysis
sourceCpp("joint_test.cpp")

## RUN
tNAME = c("CRBLM","FCTX","PONS","TCTX")
GeneExp = read.delim("GeneExp_matrix.txt",row.names=1,header=T,check.names=F)
dim(GeneExp)
MethExp = read.delim("MethExp_matrix.txt",row.names=1,header=T,check.names=F)
dim(MethExp)
samples_all = colnames(GeneExp)
gene_probes = unlist(lapply(strsplit(as.character(rownames(GeneExp)), "\\."), "[", 1))
cpg_probes = unlist(lapply(strsplit(as.character(rownames(MethExp)), "\\."), "[", 1))

GenoMat = read.delim("Genotype_matrix.txt",as.is=T,row.names=1,check.names=F)
samples = unique(names(GeneExp))
GenoMat = GenoMat[,samples]
dim(GenoMat)
mean_geno = rowMeans(GenoMat)
GenoMat = GenoMat - mean_geno
nobs = ncol(GenoMat)
rownames(GenoMat) = gsub("_.*","",rownames(GenoMat))

k = length(tNAME)
ind = as.factor(rep(1:nobs,k))
tissue = as.factor(rep(c(1:k),each=nobs))

SNP = GenoMat[,rep(colnames(GenoMat), k)]
colnames(SNP) = rep(colnames(GenoMat),k)
SNP = SNP[,colnames(GeneExp)]
colnames(SNP) = colnames(GeneExp)

## Identify cis components
cisDist = 100000
geneC = read.table("gene_coords.bed",header=F,stringsAsFactors=F)
geneC = geneC[match(gene_probes,geneC[,4]),]
geneC[grep(":",geneC[,3]),3] = unlist(lapply(strsplit(geneC[grep(":",geneC[,3]),3], "\\:"), "[", 2))
snpC = read.table("snp_coords.bed",header=F,stringsAsFactors=F)
snpC[,1] = as.numeric(gsub("chr","",snpC[,1]))
snpC.split = split(snpC,snpC[,1])
geneC.split = split(geneC,geneC[,1])
newGENO = apply(geneC,1,function(x){
	gene_chr = as.numeric(x[1]); gene_tss = as.numeric(x[2])
	snp_chr = match(as.character(gene_chr),names(snpC.split))
	snp_sites = snpC.split[[snp_chr]][,2]
	maxDist = abs(gene_tss-snp_sites)
	keep = which(maxDist<=cisDist)
	cis_snps = snpC.split[[snp_chr]][keep,4]
	geno_mat = as.matrix(SNP[rownames(SNP) %in% cis_snps,])
})
geno_dim = ldply(newGENO,dim)[,2]

out = which(geno_dim==0)
if(length(out)>0){
	GeneExp = GeneExp[-out,]
	MethExp = MethExp[-out,]
	newGENO = newGENO[-out]
	gene_probes = gene_probes[-out]
	cpg_probes = cpg_probes[-out]
}
gene_out = vector(mode="list",length=nrow(GeneExp))
cat("Number of genes: ",nrow(GeneExp),"\n")
for(i in 1:nrow(GeneExp)){
	cat("Gene: ",i,"\n");
	Geno = newGENO[[i]]
	pD = rep(0,nrow(Geno))
	for(s in 1:nrow(Geno)){
		data = data.frame("IND"=as.factor(ind),"Gene"=as.numeric(GeneExp[i,]),"Tissue"=as.factor(tissue),"Meth"=as.numeric(MethExp[i,]),"Geno"=as.numeric(Geno[s,]))
		data = data[order(data$IND),]
		fit =  suppressWarnings(lmer(Gene~0+Tissue+Geno*Meth+(1|IND)+(0+Geno|Tissue)+(0+Meth|Tissue),data,REML=F));
		est.eps = sigma(fit)^2; est.tau = VarCorr(fit)[[1]][1]; est.gamma = VarCorr(fit)[[2]][1]; est.theta = VarCorr(fit)[[3]][1];
		data = na.omit(data)
		Yhat = data$Gene - (model.matrix(fit) %*% fixef(fit));
		A = as.matrix(getME(fit,"Z"))[,1:nobs];
		B = as.matrix(model.matrix(fit)[,1:k] * data$Geno)
		C = B*data$Meth
		D = as.matrix(model.matrix(fit)[,1:k] * data$Meth)
		Z = as.matrix(data$Geno)
		MG = matrix( (Z*data$Meth)-mean(Z*data$Meth),nrow(data),1);
		V = as.matrix( (est.eps * diag(nrow(data))) + (est.tau * tcrossprod(A)) + (est.gamma * tcrossprod(B)) +  (est.theta * tcrossprod(D)) );
		pD[s] = delta_test(as.matrix(Yhat),C,V);
	}
	gene_out[[i]] = pD
}

fin = data.frame(	"Genes"=rep(gene_probes,as.numeric(geno_dim[geno_dim!=0])),
					"CpG" = rep(cpg_probes,as.numeric(geno_dim[geno_dim!=0])),
					"SNP"=unlist(lapply(newGENO,rownames)),
					"Udel" = unlist(gene_out)	)

write.table(fin,"DELTA_cisAnalysis.txt",row.names=F,sep="\t",quote=F)
