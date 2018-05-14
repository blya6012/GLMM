pql <-
function(mod.mcml,family.mcml,cache){
   	if (! missing(cache))
       	    stopifnot(is.environment(cache))
	eek<-getEk(mod.mcml$z)
	
	#need inits for parameters
	sigma<-rep(1,length(mod.mcml$z))
	beta<-rep(0,ncol(mod.mcml$x))
	
	#need inits for random effects
	nrand<-lapply(mod.mcml$z,ncol)
	nrandom<-unlist(nrand)
	totnrandom<-sum(nrandom)
	s<-rep(1,totnrandom)
	
	#need to give this Y, X, etc
	Y = mod.mcml$y
	X = mod.mcml$x
	Z = do.call(cbind,mod.mcml$z)
	ntrials = mod.mcml$ntrials
	
	outer.optim<-suppressWarnings(optim(par=sigma, fn=fn.outer,  beta=beta, s=s, Y=Y ,X=X ,Z=Z ,eek=eek ,family.mcml=family.mcml, cache=cache, ntrials=ntrials))

	list(sigma=outer.optim$par)

}


fn.outer <-
function(par,beta,s,Y,X,Z,eek,family.mcml,cache, ntrials){
	sigma<-par
	Aks<-Map("*",eek,sigma)
	A<-addVecs(Aks) #at this point still a vector
	A<-diag(A) #takes the vector and makes it a diag matrix

 	if (! missing(cache)) {
             stopifnot(is.environment(cache))
             if (exists("s.twid", envir = cache, inherits = FALSE)) {
                 stopifnot(is.numeric(cache$s.twid))
                 stopifnot(is.finite(cache$s.twid))
                 stopifnot(length(cache$s.twid) == length(s))
                 s <- cache$s.twid
             }
#
	     if (exists("beta.twid", envir = cache, inherits = FALSE)) {
                 stopifnot(is.numeric(cache$beta.twid))
                 stopifnot(is.finite(cache$beta.twid))
                 stopifnot(length(cache$beta.twid) == length(beta))
                 beta <- cache$beta.twid
             }
         }

	nbeta<-length(beta)
	#run trust
	inner.optim<-trust(fn.inner.trust,parinit=c(beta,s), rinit=5, rmax=10000, minimize=F, Y=Y, X=X, Z=Z, A=A, nbeta=nbeta, family.mcml=family.mcml, cache=cache, ntrials = ntrials)
	
	#get beta and s
	beta.twid<-cache$beta.twid
	s.twid<-cache$s.twid
	
	#calculate eta.twid
	eta.twid<-(X%*%beta.twid+Z%*%A%*%s.twid)
	
	#calculate el using eta.twid. this is piece1
	family.mcml<-getFamily(family.mcml)

	if(family.mcml$family.glmm=="bernoulli.glmm"){
		piece1<-.C(C_elc, as.double(Y), as.double(X), as.integer(nrow(X)), as.integer(ncol(X)), as.double(eta.twid), as.integer(1), as.integer(ntrials), value=double(1), gradient=double(ncol(X)), hessian=double((ncol(X)^2)))$value}
	if(family.mcml$family.glmm=="poisson.glmm"){
		piece1<-.C(C_elc, as.double(Y), as.double(X), as.integer(nrow(X)), as.integer(ncol(X)), as.double(eta.twid), as.integer(2), as.integer(ntrials), value=double(1), gradient=double(ncol(X)), hessian=double((ncol(X)^2)))$value}
	if(family.mcml$family.glmm=="binomial.glmm"){
		piece1<-.C(C_elc, as.double(Y), as.double(X), as.integer(nrow(X)), as.integer(ncol(X)), as.double(eta.twid), as.integer(3), as.integer(ntrials), value=double(1), gradient=double(ncol(X)), hessian=double((ncol(X)^2)))$value}
	
	#calculate W = c''(eta.twid)
	if(family.mcml$family.glmm=="bernoulli.glmm") {dubya<-family.mcml$cpp(eta.twid)}
	if(family.mcml$family.glmm=="poisson.glmm"){dubya<-family.mcml$cpp(eta.twid)}
	if(family.mcml$family.glmm=="binomial.glmm"){dubya<-family.mcml$cpp(eta.twid, ntrials)}
	W<-diag(as.vector(dubya))
	
	#calculate thatthing = A t(Z)WZA+I
	thatthing<-A%*%t(Z)%*%W%*%Z%*%A+diag(nrow(A))
	L<-chol(thatthing)
	
	#calculate piece2a = .5 log det(thatthing)
	piece2a<-sum(log(diag(L)))
	
	
	#calculate piece 2b =  .5 s.twid's.twid 
	piece2b <- .5*s%*%s
	
	#then minvalue= piece2a+piece2b-piece1
	minvalue<-piece2a+piece2b-piece1
	
	as.vector(minvalue)
}

fn.inner.trust <-
function(mypar,Y,X,Z,A,family.mcml,nbeta,cache, ntrials )
{
	stuff <- .C(C_fnInnerTrust, as.double(mypar), as.double(Y), as.double(X), as.double(Z), as.double(A), as.integer(nbeta), as.integer(ntrials), as.integer(family.mcml))

	cache$s.twid<-s
	cache$beta.twid<-beta
	
	list(value=value,gradient=gradient,hessian=hessian)
}
