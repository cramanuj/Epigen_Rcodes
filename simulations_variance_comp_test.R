#!/usr/bin/env Rscript

###########################################################################################
##	Run simulations testing our joint model against the variance component score test and
##		restricted likelihood ratio test (RLRT)
##
##
##	Author: Chaitanya Acharya
##	Updated: March 7, 2016
###########################################################################################

## Load the following libraries

lib.list = c("lme4","plyr","mvtnorm","Rcpp","RcppArmadillo","compiler","RLRsim")
for(i in 1:length(lib.list)){
	if(any(installed.packages()[,1]==lib.list[i])){
		library(lib.list[i],character.only=T)}else{
			source("http://bioconductor.org/biocLite.R");
			biocLite(lib.list[i]);
			require(lib.list[i],character.only=T)};
}

sourceCpp("joint_test.cpp");


## Main Function

main = function(nobs=100,k=5,maf=0.3,mu=5,lambda=5,eps=1,tau=1,theta=1,
					bta=0, gamma=0, delta=0, phi=0){

	snp = rbinom(nobs,2,maf);
	v = rnorm(k,0,gamma);
	w = rnorm(k,0,delta);
	x = rnorm(k,0,theta);
	alpha = rep(mu,k);
	meth_dat = rmvnorm(nobs,rep(0,k),diag(k));
	Y = do.call("rbind",lapply(1:nobs,function(i){ 
		alpha + bta*snp[i] + lambda*meth_dat[i,] + phi*snp[i]*meth_dat[i,] + rnorm(1,0,tau) + v*snp[i] + w*snp[i]*meth_dat[i,] + x*meth_dat[i,] + rnorm(k,0,eps) }
		));
	data = data.frame("IND"=as.factor(rep(1:nobs,each=k)),"Gene"=as.vector(t(Y)),"Geno" = rep(snp,each=k),"Tissue"=as.factor(rep(1:k,nobs)),"Meth"=as.vector(t(meth_dat)));
	fit =  suppressWarnings(lmer(Gene~0+Tissue+Geno*Meth+(1|IND)+(0+Geno|Tissue)+(0+Meth|Tissue),data,REML=F));
	est.eps = sigma(fit)^2; est.tau = VarCorr(fit)[[1]][1]; est.gamma = VarCorr(fit)[[2]][1]; est.theta = VarCorr(fit)[[3]][1];
	Yhat = data$Gene - (model.matrix(fit) %*% fixef(fit));
	A = as.matrix(getME(fit,"Z"))[,1:nobs];
	geno.list = meth.list = meth.geno.list = list();
	for(i in 1:nobs){
		geno.list[[i]] = (snp[i]-mean(snp))*diag(k);
		meth.list[[i]] = meth_dat[i,]*diag(k);
		meth.geno.list[[i]] = ( (snp[i]-mean(snp))*meth_dat[i,] )*diag(k);
	}
	MG = matrix( (data$Geno*data$Meth)-mean(data$Geno*data$Meth),nobs*k,1);
	B = do.call("rbind",geno.list);
	C = do.call("rbind",meth.geno.list);
	D = do.call("rbind",meth.list);
	Z = matrix(rep(snp-mean(snp),each=k),nobs*k,1);
	V = as.matrix( (est.eps * diag(nobs*k)) + (est.tau * tcrossprod(A)) + (est.gamma * tcrossprod(B)) +  (est.theta * tcrossprod(D)) );
	
	pD = delta_test(as.matrix(Yhat),C,V);

## Joint test
	fitJ = lmer(Gene~0+Tissue+Meth+(1|IND)+(0+Meth|Tissue),data,REML=F);
	est_eps = sigma(fitJ)^2; est_tau = VarCorr(fitJ)[[1]][1]; est_theta = VarCorr(fitJ)[[2]][1];
	Yhat = data$Gene - model.matrix(fitJ) %*% fixef(fitJ);
	V = as.matrix( (est_eps * diag(nobs*k)) + (est_tau * tcrossprod(A)) + (est_theta * tcrossprod(D)) );
	JT = joint_test2(as.matrix(Yhat),Z,MG,B,C,V);
	
## RLR sim
	fit_full = suppressWarnings(lmer(Gene~0+Tissue+Geno*Meth+(1|IND)+(0+Geno|Tissue)+(0+Geno:Meth|Tissue)+(0+Meth|Tissue),data,REML=F));
	fitM = suppressWarnings(lmer(Gene~0+Tissue+Geno*Meth+(0+Geno:Meth|Tissue),data));
	rlrt = exactRLRT(fitM,fit_full,fit)
	return(c("DELTA_TEST"=pD,"RLRT"=rlrt$p.value,"JOINT_TEST"=JT));
}

test = cmpfun(main)

###################
# Power simulations
###################

nsim=1000

nt = c(5)
maf=c(0.30)
bta = phi = c(0,0.25)
gamma = delta =  c(0,0.15,0.25)

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
					cat("GAM: ",gamma[g],"\t")
					for(p in 1:length(phi)){
						cat("PHI: ",phi[p],"\t")
						FUN = do.call("rbind",rlply(nsim,test(nobs=100,k=nt[i],bta=bta[b],phi=phi[b],gamma=gamma[g],maf=maf[j],delta=delta[d])))
						power = apply(FUN,2,function(x)sum(x<=0.05)/nsim)
						testOUT = c("GAM"=gamma[g],power)
						outPHI[[p]] = c("MAF"=maf[j],"BTA"=bta[b],"PHI"=phi[p],"DELTA"=delta[d],"GAM"=gamma[g],testOUT)
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
	filename = paste("deltaTEST_nt_",nt[i],"_",format(Sys.Date(),format="%b%d.%Y"),".txt",sep="")
	write.table(out[[i]],filename,sep="\t",row.names=F,quote=F)
	cat("\n")
}
out
