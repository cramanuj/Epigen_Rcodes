#!/usr/bin/env Rscript

######################################################################################################
## Tissue-by-Tissue analysis using the following linear model -
##
## Y = J\alpha + G\beta + M\lambda + MG\phi + \Xi 
##
## LAST UPDATE: March 7, 2016
## AUTHOR: Chaitanya Acharya
##
## NOTES:
## 1) This script should be run on each partition of gene expression/CpG methylation data
## 2) The results from the analysis are stored in a text file labeled  "ANOVA_cisAnalysis.txt"
## 3) Please change the path to the genotype file  and bed files before running this script
########################################################################################################

lib.list = c("plyr")

for(i in 1:length(lib.list)){
	if(any(installed.packages()[,1]==lib.list[i])){
		library(lib.list[i],character.only=T)}else{
			source("http://bioconductor.org/biocLite.R");
			biocLite(lib.list[i]);
			require(lib.list[i],character.only=T)};
}

## Function to extract the omnibus p-value
tbt <- function (modelobject) {
    if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
    f <- summary(modelobject)$fstatistic
    p <- pf(f[1],f[2],f[3],lower.tail=F)
    attributes(p) <- NULL
    return(p)
}
##

tNAME = c("CRBLM","FCTX","PONS","TCTX")
GeneExp = read.delim("GeneExp_matrix.txt",row.names=1,header=T,check.names=F)
dim(GeneExp)
MethExp = read.delim("MethExp_matrix.txt",row.names=1,header=T,check.names=F)
dim(MethExp)
samples_all = colnames(GeneExp); samples = unique(samples_all);
gene_probes = rownames(GeneExp); cpg_probes = rownames(MethExp);
gene_probes = unlist(lapply(strsplit(as.character(rownames(GeneExp)), "\\."), "[", 1))
cpg_probes = unlist(lapply(strsplit(as.character(rownames(MethExp)), "\\."), "[", 1))
GenoMat = read.delim("Genotype_matrix.txt",header=T,row.names=1,check.names=F)
dim(GenoMat)
nobs = ncol(GenoMat)
GenoMat = GenoMat[,match(samples,colnames(GenoMat))]
rownames(GenoMat) = gsub("_.*","",rownames(GenoMat))

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
	geno_mat = as.matrix(GenoMat[rownames(GenoMat) %in% cis_snps,])
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

dataF = dataM = list()
dataF[[1]] = GeneExp[,1:145]
dataF[[2]] = GeneExp[,146:290]
dataF[[3]] = GeneExp[,291:435]
dataF[[4]] = GeneExp[,436:580]
dataM[[1]] = MethExp[,1:145]
dataM[[2]] = MethExp[,146:290]
dataM[[3]] = MethExp[,291:435]
dataM[[4]] = MethExp[,436:580]

tissue_eqtls=list()
for(i in 1:length(dataF)){
	cat("Tissue: ",tNAME[i],"\n")
	exp_data = dataF[[i]];  meth_data = dataM[[i]]
	gene_out = list()
	for(g in 1:nrow(exp_data)){
		cat("Gene: ",g,"\n")
		gene = as.numeric(exp_data[g,])
		cpg = as.numeric(meth_data[g,])
		snpdata = newGENO[[g]]
		snp_out = numeric(0)
		for(s in 1:nrow(snpdata)){
			snp = snpdata[s,]
			snp_out[s] = tbt(lm(gene~cpg*snp))
		}
		gene_out[[g]] = snp_out
	}
	tissue_eqtls[[i]] = unlist(gene_out)
}
fin = do.call("cbind",tissue_eqtls)
colnames(fin) = tNAME
fin=cbind( 	"Gene" = as.character(rep(rownames(GeneExp),as.numeric(geno_dim[geno_dim!=0]))),
			"CpG" = as.character(rep(rownames(MethExp),as.numeric(geno_dim[geno_dim!=0]))),
			"SNP" = unlist(lapply(newGENO,rownames)),
			fin)
fin[,1] = unlist(lapply(strsplit(as.character(fin[,1]), "\\."), "[", 1))
fin[,2] = unlist(lapply(strsplit(as.character(fin[,2]), "\\."), "[", 1))

write.table(fin,"ANOVA_cisAnalysis.txt",sep="\t",row.names=F,quote=F)
