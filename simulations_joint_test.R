#!/usr/bin/env Rscript

#############################################################################
##	Run simulations testing our joint model
##	Individual components of our score test will be visualized
##
##	Author: Chaitanya Acharya
##	Updated: March 7, 2016
############################################################################

lib.list = c("lme4","plyr","mvtnorm","Rcpp","RcppArmadillo","compiler")
for(i in 1:length(lib.list)){
	if(any(installed.packages()[,1]==lib.list[i])){
		library(lib.list[i],character.only=T)}else{
			source("http://bioconductor.org/biocLite.R");
			biocLite(lib.list[i]);
			require(lib.list[i],character.only=T)};
}

sourceCpp("joint_test.cpp");

## Function to run simulations
joint_test = function(nobs=100,k=5,maf=0.10,mu=5,lambda=5,tau=1,theta=1,eps=1,
						bta=0,gamma=0,delta=0,phi=0){
	snp = rbinom(nobs,2,maf);
	v = rnorm(k,0,gamma);
	w = rnorm(k,0,delta);
	x = rnorm(k,0,theta);
	alpha = rep(mu,k);
	meth_dat = rmvnorm(nobs,rep(0,k),diag(k));
	Y = do.call("rbind",lapply(1:nobs,function(i){ 
		alpha + (bta*snp[i]) + (lambda*meth_dat[i,]) + (phi*snp[i]*meth_dat[i,]) + rnorm(1,0,tau) + (v*snp[i]) + (w*snp[i]*meth_dat[i,]) + (x*meth_dat[i,]) + rnorm(k,0,eps) }
		));
	data = data.frame("IND"=as.factor(rep(1:nobs,each=k)),"Gene"=as.vector(t(Y)),"Geno" = rep(snp,each=k),"Tissue"=as.factor(rep(1:k,nobs)),"Meth"=as.vector(t(meth_dat)));
	fit = suppressWarnings(lmer(Gene~0+Tissue+Meth+(1|IND)+(0+Meth|Tissue),data,REML=F));
	est_eps = sigma(fit)^2; est_tau = VarCorr(fit)[[1]][1]; est_theta = VarCorr(fit)[[2]][1];
	Yhat = data$Gene - model.matrix(fit) %*% fixef(fit);
	A = as.matrix(getME(fit,"Z"))[,1:nobs];
	Z = matrix(rep(snp-mean(snp),each=k),nobs*k,1);
	MG = matrix( (Z*data$Meth)-mean(Z*data$Meth),nobs*k,1);
	geno.list = meth.list = meth.geno.list = list();
	for(i in 1:nobs){
		geno.list[[i]] = (snp[i]-mean(snp))*diag(k);
		meth.list[[i]] = meth_dat[i,]*diag(k);
		meth.geno.list[[i]] = ( (snp[i]-mean(snp))*meth_dat[i,] )*diag(k);
	}
	B = do.call("rbind",geno.list);
	C = do.call("rbind",meth.geno.list);
	D = do.call("rbind",meth.list);
	
	V = as.matrix( (est_eps * diag(nobs*k)) + (est_tau * tcrossprod(A)) + (est_theta * tcrossprod(D)) );
	return(joint_test1(as.matrix(Yhat),Z,MG,B,C,V));
}

test = cmpfun(joint_test)

nsim=1000
nt = c(5)
maf=c(0.30)
bta = phi = c(0,0.25)
gamma = delta = c(0,0.15,0.25)

outMAF = vector(mode="list",length=length(maf))
outGAM = vector(mode="list",length=length(gamma))
outDELTA = vector(mode="list",length=length(delta))
outBTA = vector(mode="list",length=length(bta))
outPHI = vector(mode="list",length=length(phi))
out = vector(mode="list",length=length(nt))

for(i in 1:length(nt)){
	cat("Number of tissues: ",nt[i],"\n")
	for(j in 1:length(maf)){
		cat("\n MAF: ",maf[j],"\t")
		for(b in 1:length(bta)){
			cat("\n BTA: ",bta[b],"\t")
			for(d in 1:length(delta)){
				cat("\n DELTA: ",delta[d],"\t")
				for(g in 1:length(gamma)){
					cat("\n GAM: ",gamma[g],"\t")
					for(p in 1:length(phi)){
						cat("PHI: ",phi[p],"\t")
						FUN = do.call("rbind",rlply(nsim,test(nobs=100,k=nt[i],bta=bta[b],gamma=gamma[g],maf=maf[j],delta=delta[d],phi=phi[p])))
						power = apply(FUN,2,function(x)sum(x<=0.05)/nsim)
						testOUT = c("PHI"=phi[p],power)
						outPHI[[p]] = c("MAF"=maf[j],"BTA"=bta[b],"DELTA"=delta[d],"GAM"=gamma[g],testOUT)
					}
					outGAM[[g]] = do.call("rbind",outPHI)
				}
				outDELTA[[d]] = do.call("rbind",outGAM)
			}
			outBTA[[b]] = do.call("rbind",outDELTA)
		}
		outMAF[[j]] = do.call("rbind",outBTA)
	}
	out[[i]] = do.call("rbind",outMAF)
	names(out)[[i]] = nt[i]
	colnames(out[[i]])[6:10]=c("Ubta","Uphi","Ugam","Udel","Upsi")
	filename = paste("powSims_nt_",nt[i],"_",format(Sys.Date(),format="%b%d.%Y"),".txt",sep="")
	write.table(out[[i]],filename,sep="\t",row.names=F,quote=F)
	cat("\n")
}
out
