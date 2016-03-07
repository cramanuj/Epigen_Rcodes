/*
 Rcpp script to run our joint model and a variance component score test
 
 Author: Chaitanya Acharya
 Updated: March 7, 2016
 */

#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <stdint.h>
#include <vector>
#include <functional>
#include <algorithm>
#include <iostream>

using namespace Rcpp;
using namespace std;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

SEXP joint_test1(const arma::mat Yhat, const arma::mat Z, const arma::mat MG, const arma::mat B, const arma::mat C, const arma::mat V ){

	const arma::mat iV = inv(V);
	double varB = 2*trace(Z.t()*iV*Z*Z.t()*iV*Z);
	double meanB = trace(Z.t()*iV*Z);
	double meanP = trace( MG.t() * iV * MG );
	double varP =  2 * trace(MG.t()*iV*MG * MG.t()*iV*MG);
	double varG = 0.5*trace(B.t()*iV*B * B.t()*iV*B);
	double meanG = 0.5*trace(B.t()*iV*B);
	double varD = 0.5*trace(C.t()*iV*C * C.t()*iV*C);
	double meanD = 0.5*trace(C.t()*iV*C);
	
	double covBG = trace(iV*Z*Z.t()*iV*B*B.t());
	double covBP = 2 * trace(iV*Z*Z.t()*iV*MG*MG.t());
	double covBD = trace(iV*Z*Z.t()*iV*C*C.t());
	double covPG = trace(iV*MG*MG.t()*iV*B*B.t());
	double covPD = trace(iV*MG*MG.t()*iV*C*C.t());
	double covGD = 0.5*trace(iV*B*B.t()*iV*C*C.t());

	arma::mat V_psi(4,4,arma::fill::zeros);
	V_psi(0,0) = varB; V_psi(0,1)=covBP;V_psi(0,2)=covBG;V_psi(0,3)=covBD;
	V_psi(1,0) = covBP; V_psi(1,1)=varP; V_psi(1,2)=covPG;V_psi(1,3)=covPD;
	V_psi(2,0) = covBG; V_psi(2,1)=covPG;V_psi(2,2)=varG;V_psi(2,3)=covGD;
	V_psi(3,0) = covBD; V_psi(3,1)=covPD; V_psi(3,2)=covGD;V_psi(3,3)=varD;
	const arma::mat m1(4,1,arma::fill::ones);
	const arma::mat m2(1,4,arma::fill::ones);
	const arma::mat a = (inv(V_psi)*m1) / as_scalar(m2*inv(V_psi)*m1);

	const double U_bta = as_scalar(Yhat.t()*iV*Z*Z.t()*iV*Yhat);
	double bta_pval = 1 - R::pchisq(U_bta/(varB/(2*meanB)),2*(meanB*meanB)/varB,1,0);
	if(bta_pval==0) bta_pval = 2e-16;

	const double U_phi = as_scalar(Yhat.t()*iV*MG*MG.t()*iV*Yhat);
	double phi_pval = 1 - R::pchisq(U_phi/(varP/(2*meanP)),2*(meanP*meanP)/varP,1,0);
	if(phi_pval==0) phi_pval = 2e-16;
	
	const double U_gam = 0.5 * as_scalar(Yhat.t()*iV*B*B.t()*iV*Yhat);
	double gam_pval = 1 - R::pchisq(U_gam/(varG/(2*meanG)),2*(meanG*meanG)/varG,1,0);
	if(gam_pval==0) gam_pval = 2e-16;
	
	const double U_del = 0.5 * as_scalar(Yhat.t()*iV*C*C.t()*iV*Yhat);
	double del_pval = 1 - R::pchisq(U_del/(varD/(2*meanD)),2*(meanD*meanD)/varD,1,0);
	if(del_pval==0) del_pval = 2e-16;
	
	const double U_psi = a(0,0)*U_bta + a(1,0)*U_phi + a(2,0)*U_gam + a(3,0)*U_del;
	const double meanPSI = a(0,0)*meanB + a(1,0)*meanP + a(2,0)*meanG + a(3,0)*meanD;
	const double varPSI = a(0,0)*a(0,0)*varB + a(1,0)*a(1,0)*varP + a(2,0)*a(2,0)*varG + a(3,0)*a(3,0)*varD;
	double psi_pval = 1 - R::pchisq(U_psi/(varPSI/(2*meanPSI)),2*(meanPSI*meanPSI)/varPSI,1,0);
	if(psi_pval==0) psi_pval = 2e-16;
		
	arma::rowvec out(5,arma::fill::zeros);
	out[0] = bta_pval; out[1] = phi_pval; out[2] = gam_pval; out[3] = del_pval; out[4] = psi_pval;
	return wrap(out);
}

// [[Rcpp::export]]

SEXP joint_test2(const arma::mat Yhat, const arma::mat Z, const arma::mat MG, const arma::mat B, const arma::mat C, const arma::mat V ){
	
	const arma::mat iV = inv(V);
	double varB = 2*trace(Z.t()*iV*Z*Z.t()*iV*Z);
	double meanB = trace(Z.t()*iV*Z);
	double meanP = trace( MG.t() * iV * MG );
	double varP =  2 * trace(MG.t()*iV*MG * MG.t()*iV*MG);
	double varG = 0.5*trace(B.t()*iV*B * B.t()*iV*B);
	double meanG = 0.5*trace(B.t()*iV*B);
	double varD = 0.5*trace(C.t()*iV*C * C.t()*iV*C);
	double meanD = 0.5*trace(C.t()*iV*C);
	
	double covBG = trace(iV*Z*Z.t()*iV*B*B.t());
	double covBP = 2 * trace(iV*Z*Z.t()*iV*MG*MG.t());
	double covBD = trace(iV*Z*Z.t()*iV*C*C.t());
	double covPG = trace(iV*MG*MG.t()*iV*B*B.t());
	double covPD = trace(iV*MG*MG.t()*iV*C*C.t());
	double covGD = 0.5*trace(iV*B*B.t()*iV*C*C.t());

	arma::mat V_psi(4,4,arma::fill::zeros);
	V_psi(0,0) = varB; V_psi(0,1)=covBP;V_psi(0,2)=covBG;V_psi(0,3)=covBD;
	V_psi(1,0) = covBP; V_psi(1,1)=varP; V_psi(1,2)=covPG;V_psi(1,3)=covPD;
	V_psi(2,0) = covBG; V_psi(2,1)=covPG;V_psi(2,2)=varG;V_psi(2,3)=covGD;
	V_psi(3,0) = covBD; V_psi(3,1)=covPD; V_psi(3,2)=covGD;V_psi(3,3)=varD;

	const arma::mat m1(4,1,arma::fill::ones);
	const arma::mat m2(1,4,arma::fill::ones);
	const arma::mat a = (inv(V_psi)*m1) / as_scalar(m2*inv(V_psi)*m1);

	const double U_bta = as_scalar(Yhat.t()*iV*Z*Z.t()*iV*Yhat);
	const double U_phi = as_scalar(Yhat.t()*iV*MG*MG.t()*iV*Yhat);
	const double U_gam = 0.5 * as_scalar(Yhat.t()*iV*B*B.t()*iV*Yhat);
	const double U_del = 0.5 * as_scalar(Yhat.t()*iV*C*C.t()*iV*Yhat);
	const double U_psi = a(0,0)*U_bta + a(1,0)*U_phi + a(2,0)*U_gam + a(3,0)*U_del;
	const double meanPSI = a(0,0)*meanB + a(1,0)*meanP + a(2,0)*meanG + a(3,0)*meanD;
	const double varPSI = a(0,0)*a(0,0)*varB + a(1,0)*a(1,0)*varP + a(2,0)*a(2,0)*varG + a(3,0)*a(3,0)*varD;
	double out = 	1 - R::pchisq(U_psi/(varPSI/(2*meanPSI)),2*(meanPSI*meanPSI)/varPSI,1,0);
	if(out == 0) out=2e-16;
	return(wrap(out));
}

// [[Rcpp::export]]

SEXP beta_gamma_test(const arma::mat Yhat, const arma::mat Z, const arma::mat B, const arma::mat V ){
	
	const arma::mat iV = inv(V);
	double varB = 2*trace(Z.t()*iV*Z*Z.t()*iV*Z);
	double meanB = trace(Z.t()*iV*Z);
	double varG = 0.5*trace(B.t()*iV*B * B.t()*iV*B);
	double meanG = 0.5*trace(B.t()*iV*B);
	double cov = trace(iV*Z*Z.t()*iV*B*B.t());
	
	double a_gam = (varB-cov)/(varB+varG-2*cov); 
	double a_beta = (varG-cov)/(varB+varG-2*cov);
	
	const double U_bta = as_scalar(Yhat.t()*iV*Z*Z.t()*iV*Yhat);
	const double U_gam = 0.5 * as_scalar(Yhat.t()*iV*B*B.t()*iV*Yhat);
	const double U_psi = a_beta*U_bta + a_gam*U_gam;
	const double meanPSI = a_beta*meanB + a_gam*meanG;
	const double varPSI = a_beta*a_beta*varB + a_gam*a_gam*varG;
	return wrap(1 - R::pchisq(U_psi/(varPSI/(2*meanPSI)),2*(meanPSI*meanPSI)/varPSI,1,0));
}

// [[Rcpp::export]]

SEXP delta_test(const arma::mat Yhat, const arma::mat C, const arma::mat V ){
	
	const arma::mat iV = inv(V);
	double U_delta = as_scalar( (Yhat.t() * iV * C * C.t() * iV * Yhat) );
	double varD = trace(iV * C * C.t() * iV * C * C.t() );
	double meanD = trace(C.t() * iV * C);
	double b = 2*(meanD*meanD)/varD; double a = varD/(2*meanD);
	return wrap(1 - R::pchisq(U_delta/a,b,1,0));
	
}
