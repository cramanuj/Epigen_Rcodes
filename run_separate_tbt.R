#!/usr/bin/env Rscript

#################################################################################################################
## Separate tissue-by-tissue analysis of gene expression and methylation data using the following linear model -
##
## Y = G\beta + \Xi
##	where Y is either gene expression or methylation data
## 
##
## LAST UPDATE: March 7, 2016
## AUTHOR: Chaitanya Acharya
##
## NOTES:
## 1) This script should be run on each partition of gene expression/CpG methylation data
## 2) The results from the analysis are stored in two text files labeled "TBT_eQTL_cis.txt" and "TBT_mQTL_cis.txt"
## 3) Please change the path to the genotype file  and bed files before running this script
#################################################################################################################

lib.list = c("plyr")

for(i in 1:length(lib.list)){
        if(any(installed.packages()[,1]==lib.list[i])){
                library(lib.list[i],character.only=T)}else{
                        source("http://bioconductor.org/biocLite.R");
                        biocLite(lib.list[i]);
                        require(lib.list[i],character.only=T)};
}

tNAME = c("CRBLM","FCTX","PONS","TCTX")
GeneExp = read.delim("GeneExp_matrix.txt",row.names=1,header=T,check.names=F)
dim(GeneExp)
MethExp = read.delim("MethExp_matrix.txt",row.names=1,header=T,check.names=F)
dim(MethExp)
samples_all = colnames(GeneExp); samples = unique(samples_all)
gene_probes = unlist(lapply(strsplit(as.character(rownames(GeneExp)), "\\."), "[", 1))
cpg_probes = unlist(lapply(strsplit(as.character(rownames(MethExp)), "\\."), "[", 1))

GenoMat = read.delim("/home/igsp1/ca31/GeneRegulation/Brain/Genotype_matrix.txt",header=T,row.names=1,check.names=F)
dim(GenoMat)
nobs = ncol(GenoMat)
GenoMat = GenoMat[,match(samples,colnames(GenoMat))]
rownames(GenoMat) = gsub("_.*","",rownames(GenoMat))

cisDist = 100000
geneC = read.table("/hpchome/igsp1/ca31/GeneRegulation/gibbs_brain/gene_coords.bed",header=F,stringsAsFactors=F)
geneC = geneC[match(gene_probes,geneC[,4]),]
geneC[grep(":",geneC[,3]),3] = unlist(lapply(strsplit(geneC[grep(":",geneC[,3]),3], "\\:"), "[", 2))
snpC = read.table("/hpchome/igsp1/ca31/GeneRegulation/Brain/snp_coords.bed",header=F,stringsAsFactors=F)
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
	newGENO = newGENO[-out]
	MethExp = MethExp[-out,]
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

tissue_eqtls=tissue_mqtls=list()
for(i in 1:length(dataF)){
	cat("Tissue: ",tNAME[i],"\n")
	data = dataF[[i]]; meth = dataM[[i]]
	samples = colnames(dataF[[i]])
	gene_out = meth_out = list()
	for(g in 1:nrow(data)){
		cat("Gene: ",g,"\n")
		gene = as.numeric(data[g,])
		cpg =as.numeric(meth[g,])
		geno = newGENO[[g]]
		if(nrow(geno)==1){
			snpdata = as.matrix(t(geno[,match(samples,colnames(geno))]))
		}else{
			snpdata = geno[,match(samples,colnames(geno))]
		}
		eqtl_out = mqtl_out = NULL
		for(s in 1:nrow(snpdata)){
			snp = snpdata[s,]
			eqtl_out[s] = coef(summary(lm(gene~snp)))[2,4]
			mqtl_out[s] = coef(summary(lm(cpg~snp)))[2,4]
		}
		gene_out[[g]] = eqtl_out
		meth_out[[g]] = mqtl_out
	}
	tissue_eqtls[[i]] = as.numeric(unlist(gene_out))
	tissue_mqtls[[i]] = as.numeric(unlist(meth_out))
}
fin_eqtl = do.call("cbind",tissue_eqtls)
fin_mqtl = do.call("cbind",tissue_mqtls)
colnames(fin_eqtl) = colnames(fin_mqtl) = tNAME
out_eqtl = cbind("Gene"=rep(gene_probes,as.numeric(geno_dim[geno_dim!=0])),"SNP"=unlist(lapply(newGENO,rownames)),fin_eqtl)
out_meqtl = cbind("CpG"= rep(cpg_probes,as.numeric(geno_dim[geno_dim!=0])),"SNP"=unlist(lapply(newGENO,rownames)),fin_mqtl)
out_eqtl[,1] = unlist(lapply(strsplit(as.character(out_eqtl[,1]), "\\."), "[", 1))
out_meqtl[,1] = unlist(lapply(strsplit(as.character(out_meqtl[,1]), "\\."), "[", 1))
dim(out_eqtl); dim(out_meqtl)

write.table(out_eqtl[!(duplicated(out_eqtl)),],"TBT_eQTL_cis.txt",sep="\t",row.names=F,quote=F)
write.table(out_meqtl[!(duplicated(out_meqtl)),],"TBT_mQTL_cis.txt",sep="\t",row.names=F,quote=F)
