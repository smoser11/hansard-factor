#include "NBfactor.h"
#include <iostream>

#include "RNG.h"
#include "PolyaGamma.hpp"
#include <exception>

#ifdef USE_R
#include "R.h"
#include "Rmath.h"
#endif

using namespace Rcpp ;
using namespace RcppArmadillo;
using namespace std;

void rpg_devroye(double *x, int *n, double *z, int *num)
{
  
  RNG r;
  PolyaGamma pg;

  #ifdef USE_R
  GetRNGstate();
  #endif

  for(int i=0; i < *num; ++i){
    x[i] = pg.draw(n[i], z[i], r);
  }

  #ifdef USE_R
  PutRNGstate();
  #endif
  
} // rpg

// [[Rcpp::export]]
SEXP NBfactor(SEXP Ys, int ntopic, int Ndrawn, int Nburnn){
  
	Rcpp::NumericMatrix Yr(Ys); 				//	Create Rcpp matrix from SEXP
	int N = Yr.nrow(), P = Yr.ncol();			//	dim(YY)

	arma::mat YY(Yr.begin(), N, P, false);   	// reuses memory and avoids extra copy

	int K = ntopic; 								//	No of topics

	// Parameter Setting
	arma::mat V = arma::eye(K,K);
	arma::mat V0 = arma::eye(K,K);

	double a = 0.01;
	int h0 = 6, niu_a = 10, niu_alpha = 10;

	arma::colvec niu = 0.98 * arma::ones(K);
	arma::colvec q = 0.9 * arma::ones(K);

	double h=1.0, c_s = 1.0, d_s = 1.0;

	// Initialization
	Rcpp::NumericVector drawsB = Rcpp::rnorm(P*K,0,1);
	arma::mat B(drawsB.begin(),P,K,false);

	Rcpp::NumericVector drawsf = Rcpp::rnorm(N*K,0,1);
	arma::mat f(drawsf.begin(),N,K,false);
	
	Rcpp::NumericVector drawsal = Rcpp::rnorm(P,0,niu_alpha);
	arma::rowvec alpha(drawsal.begin(),P,false);
	Rcpp::NumericVector drawga = Rcpp::rnorm(N,0,1);
	arma::colvec gamma(drawga.begin(),N,false);

	arma::mat Psi = arma::zeros(N,P);
	arma::mat Omega = arma::zeros(N,P);
	arma::mat Z = arma::zeros(N,P);
	arma::mat zhat = arma::zeros(N,P);

	//Gibbs Sampler
	//int Ndrawn = 5; 
	//int Nburnn = 0;

	// Store Parameter
	arma::cube BBn = arma::zeros(P,K,Ndrawn);
	arma::cube FFn = arma::zeros(N,K,Ndrawn);

	arma::mat Alfn = arma::zeros(P,Ndrawn);
	arma::mat Gamn = arma::zeros(N,Ndrawn);

	double test;

	//Iteration
	for (int iter = 0; iter < (Ndrawn+Nburnn); ++iter)
	{
		//Calculate Psi
		for (int i = 0; i < N; ++i)
		{
			arma::rowvec temp(P);
			temp.fill(gamma[i]);
			Psi.row(i) = alpha + temp + arma::trans(B*arma::trans(f.row(i))); 
		}

		//Sample omega
		for (int i = 0; i < N; ++i)
		{
			for (int j = 0; j < P; ++j)
			{
        
          double *xx;
          int *nn;
          double *zz;
          int *numm;
    
          double xxx=0.5;
          int nnn=YY(i,j)+h;
          double zzz=Psi(i,j);
          int numnum=1;
    
          xx=&xxx;
          nn = &nnn;
          zz = &zzz;
          numm = &numnum;

          rpg_devroye(xx,nn,zz,numm);     
             
				  Omega(i,j) = *xx;
			}
		}

		//Calculate Z
		for (int i = 0; i < N; ++i)
		{
			for (int j = 0; j < P; ++j)
			{
				Z(i,j) = (YY(i,j)-h)/(2*Omega(i,j));
			}
		}

		//Sample gamma_i
		for (int i = 0; i < N; ++i)
		{
			double Vgamma = 1/(1+sum(Omega.row(i)));

			double Mgamma = Vgamma*(a+arma::as_scalar((Z.row(i)-alpha-arma::trans(B*arma::trans(f.row(i))))*arma::trans(Omega.row(i))));
			//double Mgamma = Vgamma*(a+)
			//double Mgamma = Vgamma*(a+ sum(arma::trans(B*arma::trans(f.row(i)))));
			gamma(i) = as<double>(Rcpp::rnorm(1,Mgamma,sqrt(Vgamma)));
			 
		}

		//Sample a
		double Va = 1/(N+1/niu_a);
		double Ma = Va*sum(gamma);
		a = as<double>(Rcpp::rnorm(1,Ma,sqrt(Va)));

		//Sample alpha_j
		for (int j = 0; j < P; ++j)
		{
			double Valpha = 1/(1/niu_alpha+sum(Omega.col(j)));
			double Malpha = Valpha*arma::as_scalar((arma::trans((Z.col(j)-gamma-f*arma::trans(B.row(j))))*Omega.col(j)));
			//double Malpha = Valpha*sum((Z.col(j)-gamma-f*arma::trans(B.row(j)))%Omega.col(j));
			alpha(j) = as<double>(Rcpp::rnorm(1,Malpha,sqrt(Valpha)));
		}

		//Sample f[i,]
		for (int i = 0; i < N; ++i)
		{
			arma::rowvec temp(P);
			temp.fill(gamma[i]);
			arma::mat Ome = diagmat(Omega.row(i));
			arma::mat Vf = arma::inv(arma::inv(V)+arma::trans(B)*Ome*B);
			arma::colvec Mf = Vf*arma::trans(B)*Ome*arma::trans(Z.row(i)-alpha-temp);
			Rcpp::NumericVector drawtemp = Rcpp::rnorm(K,0,1);
			arma::colvec ftemp(drawtemp.begin(),K,false);
			f.row(i) = arma::trans(Mf+chol(Vf)*ftemp);

		}

		//Sample V

		//Sample B
		for (int s = 0; s < K; ++s)
		{
			arma::mat Bs = B;
			Bs.shed_col(s);
			arma::mat fs = f;
			fs.shed_col(s);
			double Vhat, bhat, prob1,qhat;
			for (int i = 0; i < N; ++i)
			{
				arma::rowvec temp(P);
				temp.fill(gamma[i]);
				zhat.row(i) = Z.row(i)-alpha-temp-arma::trans(Bs*arma::trans(fs.row(i)));
			}
			for (int j = 0; j < P; ++j)
			{
				if (j<s)
				{
					B(j,s) = 0.0;
				}
				else
				{
					Vhat = 1/(1/niu(s)+sum(Omega.col(j)%f.col(s)%f.col(s)));
					bhat = Vhat*sum(f.col(s)%Omega.col(j)%zhat.col(j));
					test = bhat;
					prob1 = Rcpp::dnorm(NumericVector::create(0.0),test*1.0,sqrt(Vhat))[0];
					if(prob1 < 0.00001){
						qhat = 1.0;
					}
					else{
						double rhat = Rcpp::dnorm(NumericVector::create(0.0),0.0,niu(s)*1.0)[0]/prob1;
						qhat = rhat/((1.0-q(s))/q(s)+rhat);
					}

					//sample B_js
					if(j == s){
            
            RNG rr;
            
            B[j,s] = rr.tnorm(0.0, bhat, sqrt(Vhat));
            
					}
					else{
						double rnum = as<double>(Rcpp::runif(1));
						if(rnum < qhat){
							B(j,s) = as<double>(Rcpp::rnorm(1,bhat,sqrt(Vhat)));
						}
						else{
							B(j,s) = 0.0;
						}
					}
					
				}
			}
		}

		//draw niu
		arma::colvec ms = arma::zeros(K);
		for (int s = 0; s < K; ++s)
		{
			for (int j = (s+1); j < P; ++j)
			{
				if(abs(B(j,s))>1e-100){
					ms(s) = ms(s)+1.0;
				}
			}
			double temp = arma::as_scalar(arma::trans(B.col(s))*B.col(s));
			niu(s) = 1/as<double>(Rcpp::rgamma(1,(2+ms(s))/2,(2.0+temp)/2));

			//draw q
			q(s) = as<double>(Rcpp::rbeta(1,1+ms(s),1.0+P-s-ms(s)));

		}

		//Store parameter
		if(iter > (Nburnn-1)){
			BBn.slice(iter-Nburnn) = B;
			FFn.slice(iter-Nburnn) = f;
			Alfn.col(iter-Nburnn) = arma::trans(alpha);
			Gamn.col(iter-Nburnn) = gamma;

		}

		cout<<"\n iter = "<<iter;

		






	

	}
	









	return Rcpp::List::create(
        //Rcpp::Named("testB") = Vgamma,
        Rcpp::Named("BBn") = BBn,
        Rcpp::Named("FFn") = FFn,
        Rcpp::Named("Alfn") = Alfn,
        Rcpp::Named("Gamn") = Gamn
        //Rcpp::Named("stderr")       = stderrest,
        //Rcpp::Named("testFunction") = f(3)
    ) ;	



	

	





}



