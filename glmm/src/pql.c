#include "myheader.h"

/* fn.inner.trust */
void fnInnerTrust (double *mypar, double *Y, double *X, double *Z, double *A, int nbeta, double *cache, int ntrials, int family.mcml)
{
	int beta [nbeta];
	int i;
	for (i = 0; i < nbeta;  i++){
		beta[i] = i;
	}
	
	int s [nbeta]
	for (i = 0; i < nbeta; i++){
		s[i] = -i;
	}

	n_trials = ntrials;	

	/*eta<-X%*%beta+Z%*%A%*%s*/
	double *x = Calloc(*X, double);
	double *y = Calloc(*Y, double);
	double *z = Calloc(*Z, double);
	double *a = Calloc(*A, double);
	
	matvecmult(x, beta, n, nbeta, xbeta);
	matvecmult(z, a, n, n, ZA);
	matvecmult(ZA, s, n, nbeta, ZAs);
	addvec(xbeta, ZAs, n, eta);

	/*family.mcml<-getFamily(family.mcml)
	if(family.mcml$family.glmm=="bernoulli.glmm") {mu<-family.mcml$cp(eta)}
	if(family.mcml$family.glmm=="poisson.glmm"){mu<-family.mcml$cp(eta)}
	if(family.mcml$family.glmm=="binomial.glmm"){mu<-family.mcml$cp(eta, ntrials)}*/
	

	/* Being with family = Bernouli */
	cp3(eta, n, 1, ntrials, mu); 



	/*if(family.mcml$family.glmm=="bernoulli.glmm"){
		value<-.C(C_elc, as.double(Y), as.double(X), as.integer(nrow(X)), as.integer(ncol(X)), as.double(eta),as.integer(1),as.integer(ntrials), value=double(1), gradient=double(ncol(X)), hessian=double((ncol(X)^2)))$value-.5*s%*%s}
	if(family.mcml$family.glmm=="poisson.glmm"){
		value<-.C(C_elc, as.double(Y), as.double(X), as.integer(nrow(X)), as.integer(ncol(X)), as.double(eta), as.integer(2), as.integer(ntrials), value=double(1), gradient=double(ncol(X)), hessian=double((ncol(X)^2)))$value-.5*s%*%s}
	if(family.mcml$family.glmm=="binomial.glmm"){
		value<-.C(C_elc, as.double(Y), as.double(X), as.integer(nrow(X)), as.integer(ncol(X)), as.double(eta), as.integer(3), as.integer(ntrials), value=double(1), gradient=double(ncol(X)), hessian=double((ncol(X)^2)))$value-.5*s%*%s}*/



	double value = elc(y, x, n, nbeta, eta, 1,  n_trials, 1.0, nbeta, nbeta*nbeta) - 0.5*matvecmult(s, nbeta, s, nbeta, val);

	
	/*gradient calculation */
	
	matvecmult(x, y-mu, n, n, db);
	matvecmult(a, z, n, n, az);
	matvecmult(az, y-mu, n, n, azymu);
	double ds = azymu-s;
	double gradient [2] = (db, ds);

	/* hessian calculation */


	/*Bernoulli */
	cpp3(eta, n, 1, ntrials, cdub);
	int kyoo = sizeof(a);
	matvecmult(-x, cdub, n, n, piece);
	matvecumlt(piece, x, n, n, piece1);
	matvecmult(piece, z, n, n, piece2_1);
	matvecmult(piece2_1, a, n, n, piece2);
	matvecmult(-a, z, n, n, piece3_1);
	matvecmult(piece3_1, cdub, n, n, piece3_2);
	matvecmult(piece3_2, x, n, n, piece3);
	matvecmult(-a, z, n, n, piece4_1);
	matvecmult(piece4_1, cdub, n, n, piece4_2);
	matvecmult(piece4_2, z, n, n, piece4_3);
	matvecmult(piece4_3, a, n, n, piece4);
	piece4 = piece4 - kyoo;
	double hessian [4][n]

	int j;
	for (i = 0, i < 4, i++){
		for (j = 0, j <n, j++){
			if(i=0){ hessian [i][j] = piece1[i]; }
			if(i=1){ hessian [i][j] = piece2[i]; }
			if(i=2){ hessian [i][j] = piece3[i]; }
			if(i=3){ hessian [i][j] = piece4[i]; }
		}
	}
}






