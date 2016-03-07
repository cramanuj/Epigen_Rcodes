#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <stdint.h>
#include <vector>
#include <functional>
#include <algorithm>
#include <iostream>
#include <progress.hpp>

using namespace Rcpp;
using namespace std;

/*
 
 jag_fun - main JAGUAR algorithm

 Author: Chaitanya Acharya
 Updated on: Aug 27, 2015
 ***RCPP ARMADILLO version of JAGUAR***
 
 Returns a p-value indicating the significance of association between a gene-SNP pair
 Our joint score test statistic is computed as --
 
 U_\psi = Y^t . V^{-1} . (0.5 a_\gamma XX^t + a_\beta GG^t) . V^{-1} . Y
 
 Arguments: 
			Eps		-> \hat{\Epsilon} from the model
			Tau		-> \hat{\Tau} from the model
 			k		-> vector indicating number of tissues per observation
			Y		-> Residuals from the model
			snp		-> vector of mean centered genotypes 
			R		-> Matrix indicating the presence/absence of samples in the study ('1' indicates presence; '0' indicates absence)
 
 Returns: P-value computed from the Satterthwaite method
 
*/
// [[Rcpp::depends(RcppProgress)]]

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
double jag_fun(const double Eps, const double Tau, const arma::vec k, const arma::mat& Y, const arma::vec snp, const arma::mat& R){
	
	const double nobs = Y.n_cols;
	const double kmax = R.n_cols;
	double Ub=0.0, meanB=0.0, meanG=0.0;
	
	arma::vec Yhat(kmax);
	arma::mat Ugam(kmax,nobs,arma::fill::zeros);
	arma::mat varGMAT(kmax,kmax,arma::fill::zeros);
	arma::rowvec Ysum = sum(Y,0);
	
	for(int i =0; i<nobs; i++){
		double G = snp[i]-mean(snp);
		double V1 = (Eps+(k[i]-1)*Tau) / ( (Eps*Eps)+k[i]*Eps*Tau );
		double V2 = - ( (Tau)/(Eps*Eps+k[i]*Eps*Tau) );
		arma::mat V(kmax,kmax); V.fill(V2);
		arma::vec V1_vec = rep(V1,kmax); V.diag() = V1_vec;
		arma::rowvec R_vec = R.row(i);
		arma::mat RR = trans(R.row(i)) * R.row(i);
		varGMAT+= G*G*(RR % V);
		arma::colvec Yhat = Y.col(i);
		Ugam.col(i) = (trans(R_vec)*G) % ( (V1*Yhat) + (V2*(Ysum[i]-Yhat)) );
		Ub+=G * (V1 + (k[i]-1)*V2)*Ysum[i];
		meanB+= G * G * (V1+(k[i]-1)*V2) * k[i];
		meanG+=k[i] * G * G * V1;
	}
	
	const double U2beta = Ub*Ub;
	const double Ugamma = 0.5 * accu( sum(Ugam,1) % sum(Ugam,1) );
	const double varB = 2 * meanB * meanB;
	meanG = 0.5 * meanG;
	const double varG = 0.5 * accu(varGMAT%varGMAT);
	const double cov = accu( sum(varGMAT,0)%sum(varGMAT,0) );
	const double agam = (varB - cov)/(varB+varG-2*cov);
	const double abeta = (varG - cov)/(varB+varG-2*cov);
	const double Upsi = abeta*U2beta + agam*Ugamma;
	const double meanPSI = (abeta*meanB)+ (agam*meanG);
	const double varPSI = (abeta*varB*abeta) + (agam*agam*varG);
	const double a1 = varPSI/(2*meanPSI);
	const double a2 = (2*meanPSI*meanPSI)/varPSI;
	const double scaled_upsi = Upsi/a1;
	double jag_pval = 1 - R::pchisq(scaled_upsi,a2,1,0);
	if(jag_pval==0) jag_pval = 2e-16;
	return jag_pval;
}

/*
 
 GENEapply - Genome-wide application of JAGUAR algorithm
 
 Author: Chaitanya Acharya
 Updated on: Aug 27, 2015
 ***RCPP ARMADILLO version of JAGUAR***
 
 Returns a matrix of p-values (indicating the significance of association between a gene-SNP pair) with genes on rows and SNPs on columns
 
 Arguments: 
			geno	-> Genotype matrix in allele dosage format
			Y		-> List of Genes as matrices of dimension k by n where k is the # of tissues and n is total individuals
			Eps		-> Vector of MLE \hat{\epsilon} from the model for every gene
			Tau		-> Vector of MLE \hat{\tau} from the model for every gene
			k		-> Vector indicating the number of tissues for every individual (allows for missing tissues)
			R		-> Matrix indicating the presence/absence of samples in the study ('1' indicates presence; '0' indicates absence)
 
 Returns: Matrix of joint test p-values computed from the Satterthwaite method
 
 */

// [[Rcpp::export]]
arma::mat GENEapply(const arma::mat& geno, List Y, const arma::vec Eps, const arma::vec Tau, const arma::vec k, const arma::mat& R, bool display_progress=true){
	
	int ngenes = Y.size();
	int nsnps = geno.n_rows;
	
	arma::mat GENEout_pval(ngenes,nsnps);
	arma::vec SNPout_pval(nsnps);
	arma::rowvec snp(R.n_rows);
  Progress p(ngenes, display_progress);
	for(int i =0; i<ngenes; i++){
		Rcpp::checkUserInterrupt();
    p.increment();
		arma::mat Ymat = Y(i);
		double est_eps = Eps(i);
		double est_tau = Tau(i);
		for(int j=0; j<nsnps; j++){
			Rcpp::checkUserInterrupt();
			SNPout_pval(j) = jag_fun(est_eps,est_tau,k,Ymat,vectorise(geno.row(j)),R);
		}
		GENEout_pval.row(i) = vectorise(SNPout_pval,1);
	}
	return GENEout_pval;
}

/*
 
 rowsumscpp - RCPP function to compute row sums of a matrix
			- Used by jagSIM and vcSIM functions
			- Faster than ARMADILLO version of row sums (perhaps slower than R's rowSums())
			- Original code by hadley Wickham
			- Source: http://adv-r.had.co.nz/Rcpp.html
 
 */

// [[Rcpp::export]]
NumericVector rowsumscpp(NumericMatrix x) {
	
	int nrow = x.nrow(), ncol = x.ncol();
	NumericVector out(nrow);
	
	for (int i = 0; i < nrow; i++) {
		Rcpp::checkUserInterrupt();
		double total = 0;
		for (int j = 0; j < ncol; j++) {
			total += x(i, j);
		}
		out[i] = total;
	}
	return out;
}


/*
 jagSIM - JAGUAR code for simulations ONLY
 
 Author: Chaitanya Acharya
 Updated on: Aug 27, 2015
 ***RCPP version of JAGUAR in case of no missing data (balanced design)***
 
 Returns a p-value indicating the significance of association between a gene-SNP pair
 Our joint score test statistic is computed as --
 
 U_\psi = Y^t . V^{-1} . (0.5 a_\gamma XX^t + a_\beta GG^t) . V^{-1} . Y
 
 Arguments: 
				Eps -> \hat{\Epsilon} from the model
				Tau -> \hat{\Tau} from the model
				k	-> number of tissues per observation
				Y	-> Residuals from the model
				snp -> vector of mean centered genotypes 
 
 Returns: P-value computed from the Satterthwaite method
 
 */


// [[Rcpp::export]]
double jagSIM(double Eps, double Tau, double k, NumericMatrix Y, NumericVector snp){	
	
	if(k<2) stop("At least two groups are expected!"); 
	double nobs = Y.ncol();
	if(nobs != snp.size()) stop("Unequal sample sizes detected between Gene expression and Genotype data");
	NumericVector Yhat(k);
	double V1 = (Eps+(k-1)*Tau) / ( (Eps*Eps)+k*Eps*Tau);
	double V2 = - (Tau)/(Eps*Eps+k*Eps*Tau);
	NumericVector Ug(k); 
	NumericMatrix Ugam(k,nobs);
	double mG = mean(snp);
	double Ubta = 0.0;
	for(int i=0; i<nobs; i++){
		Rcpp::checkUserInterrupt();
		Yhat = Y(_,i);
		double G = snp[i]-mG;
		double Ysum = sum(Yhat);
		Ubta+=(G * Ysum)/(Eps + k * Tau);
		for(int t=0; t<k; t++){
			Ug[t] = ( (V1*Yhat[t])+(V2*(Ysum-Yhat[t])) ) * G;
		}
		Ugam(_,i)=Ug;
	}
	double Ugamma = 0.5 * sum(rowsumscpp(Ugam)*rowsumscpp(Ugam));
	double U2beta = Ubta*Ubta;
	double snp2 = sum((snp-mG)*(snp-mG));
	double meanB = snp2 * (V1+(k-1)*V2) * k;
	double varB = 2 * meanB * meanB;
	double meanG = 0.5 * sum( (snp-mean(snp))*(snp-mean(snp))*V1*k );
	double varG = 0.5 * snp2 * snp2 * (V1*V1+(k-1)*V2*V2) * k;
	double cov = (snp2 * (V1+(k-1)*V2))*(snp2 * (V1+(k-1)*V2))*k;
	double agam = (varB - cov)/(varB+varG-2*cov);
	double abeta = (varG - cov)/(varB+varG-2*cov);
	double Upsi = abeta*U2beta + agam*Ugamma;
	double meanPSI = (abeta * meanB) + (agam * meanG); 
	double varPSI = (agam * agam * varG) + (abeta * abeta * varB);
	double a1 = varPSI/(2*meanPSI);
	double a2 = (2*meanPSI*meanPSI)/varPSI;
	double scaled_upsi = Upsi/a1;
	double pval = 1 - R::pchisq(scaled_upsi,a2,1,0);
	if(pval==0) pval = 2e-16;
	return pval;
}

/*
 vcSIM - Variance component score test code for simulations ONLY
 
 Author: Chaitanya Acharya
 Updated on: Aug 27, 2015
 ***RCPP version of the variance component score test in case of no missing data (balanced design)***
 
 Returns a p-value indicating the significance of association between a gene-SNP pair

 Variance component score test statistic is computed as --
 
 U_\gamma = 0.5 * (Y^t . V^{-1} . XX^t . V^{-1} . Y)
 
 Arguments: 
			Eps -> \hat{\Epsilon} from the model
			Tau -> \hat{\Tau} from the model
			k -> vector indicating number of tissues per observation
			Y -> Residuals from the model
			snp -> vector of mean centered genotypes 
 
 Returns: P-value computed from the Satterthwaite method
 
 */


// [[Rcpp::export]]
double vcSIM(double Eps, double Tau, double k, NumericMatrix Y, NumericVector snp){  
	
	if(k<2) stop("At least two groups are expected!"); 
	double nobs = Y.ncol();
	if(nobs != snp.size()) stop("Unequal sample sizes detected between Gene expression and Genotype data");	
	NumericVector Yhat(k);
	double V1 = (Eps+(k-1)*Tau) / ( (Eps*Eps)+k*Eps*Tau);
	double V2 = - (Tau)/(Eps*Eps+k*Eps*Tau);
	NumericVector Ug(k); NumericMatrix Ugam(k,nobs);
	double mG = mean(snp);
	for(int i=0; i<nobs; i++){
		Yhat = Y(_,i);
		double Ysum = sum(Yhat);
		double G = snp[i]-mG;
		for(int t=0; t<k; t++){
			Ug[t] = ( (V1*Yhat[t])+(V2*(Ysum-Yhat[t])) ) * G;
		}
		Ugam(_,i)=Ug;
	}
	double Ugamma = 0.5 * sum(rowsumscpp(Ugam)*rowsumscpp(Ugam));
	double snp2 = sum((snp-mG)*(snp-mG));
	double meanG = 0.5 * sum( (snp-mean(snp))*(snp-mean(snp))*V1*k );
	double varG = 0.5 * snp2 * snp2 * (V1*V1+(k-1)*V2*V2) * k;
	double a1 = varG/(2*meanG); 
	double a2 = (2*meanG*meanG)/varG;
	double pval = 1 - R::pchisq(Ugamma/a1,a2,1,0);
	return pval;
}

/*
 cis_eqtl - JAGUAR code for analysis of the cis regions ONLY
			with validations through permutation-resampling method. 
 
 Author: Chaitanya Acharya
 Updated on: Aug 27, 2015
 
 Returns a p-value indicating the significance of association between a gene-SNP pair
 Our joint score test statistic is computed as --
 
 U_\psi = Y^t . V^{-1} . (0.5 a_\gamma XX^t + a_\beta GG^t) . V^{-1} . Y
 
 Arguments: 
				Eps		-> \hat{\Epsilon} from the model
				Tau		-> \hat{\Tau} from the model
				k		-> vector indicating number of tissues per observation
				Y		-> Residuals from the model
				snp		-> vector of mean centered genotypes 
				R		-> Matrix indicating the presence/absence of samples in the study ('1' indicates presence; '0' indicates absence)
 
 
 Returns: P-value computed from the Satterthwaite method
 
 */

// [[Rcpp::export]]
SEXP cis_eqtl(List SNP,List Y,IntegerMatrix perm_mat,arma::vec Eps,arma::vec Tau,arma::vec k,const arma::mat& R, bool v,bool display_progress=true){
	
	bool verbose = v;
	int ngenes = Y.size();
	NumericVector fdr_adj_pval(ngenes);
	
	if(verbose){
		int nobs = R.n_rows;
		int B = perm_mat.nrow();
		const double invB = 1.0/B;
		NumericVector per_gene_minP(B);
		NumericVector rest_perm_minP(B-1);
		//For every gene
		Progress p(ngenes, display_progress);
		for(int i = 0; i<ngenes; i++){
			Rcpp::checkUserInterrupt();
			p.increment();
			const arma::mat& Ymat = Y[i];
			const double est_eps = Eps[i];
			const double est_tau = Tau[i];
			NumericMatrix geno = SNP[i];
			int nsnps = geno.nrow();
			arma::vec SNPout_perm_pval(nsnps);
			//***** Begin PERMUTATION-RESAMPLING
			for(int p = 0;p<B; p++){
				NumericVector snp_perm(nobs);
				IntegerVector perm_order = perm_mat(p,_);
				//Compute the score test for every SNP
				for(int j=0; j<nsnps; j++){
					Rcpp::checkUserInterrupt();
					NumericVector s = geno(j,_);
					snp_perm = s[ perm_order ];
					arma::vec snp(snp_perm.begin(),snp_perm.size(),false);
					SNPout_perm_pval(j) = jag_fun(est_eps,est_tau,k,Ymat,snp,R);
				}
				per_gene_minP(p) = min(SNPout_perm_pval);
			}//***** End PERMUTATIONS
			double gene_minP = per_gene_minP[0];
			rest_perm_minP = per_gene_minP[ seq(1,per_gene_minP.size()-1) ];
			fdr_adj_pval(i) = (1.0 + count_if(rest_perm_minP.begin(), rest_perm_minP.end(), bind2nd(std::less<double>(), gene_minP)))  * invB;
		}
		return wrap(fdr_adj_pval);
	}else{
		int ngenes = Y.size();
		List GENEout_pval(ngenes);
		Progress p(ngenes, display_progress);
		for(int i =0; i<ngenes; i++){
			Rcpp::checkUserInterrupt();
			p.increment();
			const arma::mat& Ymat = Y(i);
			const double est_eps = Eps[i];
			const double est_tau = Tau[i];
			const arma::mat& geno = SNP[i];
			int nsnps = geno.n_rows;
			arma::vec SNPout_pval(nsnps);
			for(int j=0; j<nsnps; j++){
				Rcpp::checkUserInterrupt();
				SNPout_pval[j] = jag_fun(est_eps,est_tau,k,Ymat,vectorise(geno.row(j)),R);
			}
			GENEout_pval[i] = vectorise(SNPout_pval,1);
		}
		return wrap(GENEout_pval);
	}
}
