




iter = strtoi(Sys.getenv("ITER"))
i = strtoi(Sys.getenv("SGE_TASK_ID"))

source("helper_library.r")



# Always get it running on your home computer first, from a freshly started R session.
#    (But then edit it out of the code thereafter...)

# Dropbox 
# setwd("C:/Users/danielmc/Dropbox/Work/phD File/Stat Research/Publicly shared stuff/Grid")
# source("helper_library.r")


epsilon = .00001


# Make separate function to handle the creation of X and Y data, because this is repeated so frequently.
XYfunc = function(betas, efunc, n){
	# Inputs: 
	# efunc: 1 - 0     , 2 - N(0,1) , 3 - X_*^2 N(0,1)    [3 cases]
	# X_*^2 is the sum of the squared covariate values, which were in the model.
	# n    : We restrict ourselves to n = {25, 100}
	# betas: this stores the true slope coefficients for the 4 covariates we use, {X, X^2, log(X), exp(X)}.
	# Output: matrix with x_i's in 1st 4 columns (not including column of 1's!), y's thereafter (in the 5th).
	
	# Generate x_i's, eps's and y's.
	X_mat = matrix(rnorm(n=4*n,mean=0,sd=1),nrow=n,ncol=4)
	X_mat[,3] = abs(X_mat[,3])
	betas_repl = many_rows(betas,n)
	
	if(efunc == 1){eps = rep(0,n)}
	if(efunc == 2){eps = rnorm(n)}
	if(efunc == 3){eps = apply(X_mat^2*betas_repl,1,sum)*rnorm(n)}
	
	fX_mat = X_mat
	fX_mat[,2] = fX_mat[,2]^2
	fX_mat[,3] = log(fX_mat[,3])
	fX_mat[,4] = exp(fX_mat[,4])
	y = apply(fX_mat * betas_repl,1,sum) + eps
	
	Xy = cbind.data.frame(X_mat,y)
	colnames(Xy) = c("x1","x2","x3","x4","y")
	return(Xy)
	}
	
# Make quick function to get coefficients beta.hat given X and Y data.
bhat <- function(Xy, epsilon=.00001){      
	# input : Xy matrix, x_i's in the 1st 4 columns, y's thereafter.
	# output: beta coefficients.  Assumes 4 predictors (excluding intercept).
	Xm = as.matrix(cbind(rep(1,n),Xy[,1:4]))
	betahat <- solve(t(Xm)%*%Xm + diag(epsilon,5))%*%t(Xm)%*%Xy[,"y"]     # add just a smidge to the diagonal to prevent matrix instability, ridge regression style.
	return(betahat)
	}

 
resids = function(Xy,betahat){
	# input: Xy matrix *without 1's column in X*, betahat coefficients  (b/c we already have - to speed up run-time)
	# output: nx1 vector of residuals.
	Xm = as.matrix(cbind(rep(1,n),Xy[,1:4]))
	yhat = Xm %*% betahat
	resids = Xy[,"y"] - yhat
	return(resids)
	}


# Get sandwich estimate of standard error for beta hat
sand = function(Xm, resids){
	# Input: X matrix *inclusive of column of 1's*, vector of residuals.
	# Output: 5 element vector with sand est of SE for icept/beta's (no covariances)
	n=length(resids)
	# input the residuals and calculated the RMSE weighted sum on X'X
	tmp1 = 1/n*(t(Xm)%*%diag(as.vector(resids^2))%*%Xm)             
	tmp3 = solve(1/n*(t(Xm)%*%Xm+diag(epsilon,5)))           # add just a smidge to the diagonal to prevent matrix instability.
	sand.cov = 1/n*(tmp3%*%tmp1%*%tmp3)
	return(sqrt(diag(sand.cov)))
	}

	
se.bhat = function(Xm, resids){
	# Input: X matrix *inclusive of column of 1's*, vector of residuals.
	# Output: 2x1 traditional est of SE for beta's (no covariances) - (X'X)^(-1) s, s = RSS/(N-2)
	s.sq = sum(resids^2)/(n-5)
	return(sqrt(diag(solve(t(Xm)%*%Xm+diag(epsilon,5))*s.sq)))    # add just a smidge to the diagonal to prevent matrix instability.
	}

	

	
tfunc  = function(betas, efunc, n, alpha, G = 30, B = 1000, B2 = 1000){
	# Inputs: 
	# betas: 1 through 15, some combination of {X, X^2, log(X), exp(X)} (coefficients are what are specified, vector of length 4) [15 cases]
	# efunc: 1 - 0     , 2 - N(0,1)     , 3 - X_*^2*N(0,1)      [3 cases]
	# n    : We restrict ourselves to n = {25, 100, 400}  [3 cases]
	# alpha: We restrict ourselves to alpha = {0.05} for target coverage of 90%.
	# Total number of cases = 15*3*3 = 135. I will work in chunks.
	# G    : # of grid points to evaluate calibrated percentile method on each side. (=20)
	# B    : # of 1st stage bootstrap samples
	# B2   : # of 2nd stage bootstrap samples
	
	# Output: intermediate outputs in list format  - 
	# Intervals: 90% for following: z-interval, z-sandwich interval (sand in denominator), simple percentile interval, and perc-cal interval.
	# Parameters: {beta_1, ..., beta_4}    (the slope coefficients)
	
	# Initialize all the matrices and vectors to be used below.
	
	print(Sys.time())

	thetas.hat.boot = matrix(NA,nrow=B,ncol=4)
	thetas.hat.dboot = matrix(NA,nrow=B2,ncol=4)    # we are just pre-allocating for one bootstrap run each time (not all runs)

	perc.int = matrix(NA,nrow=4,ncol=2)
	perc.cal.int = matrix(NA,nrow=4,ncol=2)

	theta.qtl.lgrid.lo = array(NA,dim=c(B,G,4))       # these are what I use to track when we get the desired coverage of thetas.hat.
	theta.qtl.lgrid.hi = array(NA,dim=c(B,G,4))
	log.lgrid.lo = array(NA,dim=c(B,G,4))
	log.lgrid.hi = array(NA,dim=c(B,G,4))
	l.grid.lo = seq(.0001, alpha*3, length.out=G)     # these are the quantile grid points I evaluate over.
	l.grid.hi = 1-l.grid.lo
	percs.lo = matrix(NA,nrow=4,ncol=G)
	percs.hi = matrix(NA,nrow=4,ncol=G)
	
	
	nbig=1000000
	xy.big = XYfunc(betas, efunc, nbig)
	thetas.true = bhat(xy.big)[2:5,1]
	rm(xy.big)
	
	# Start main loop
	# t1 = Sys.time()
	Xy = XYfunc(betas, efunc=efunc, n=n)
	Xm = as.matrix(cbind(rep(1,n),Xy[,1:4]))
	bhats = bhat(Xy)
	thetas.hat = bhats[2:5,1]
	for(j in 1:B){
		Xy.boot = Xy[sample(1:n,n,replace=TRUE),]
		Xm.boot = as.matrix(cbind(rep(1,n),Xy.boot[,1:4]))
		bhats.boot = bhat(Xy.boot)
		thetas.hat.boot[j,] = bhats.boot[2:5,1]
		for(k in 1:B2){
			Xy.boot2 = Xy.boot[sample(1:n,n,replace=TRUE),]
			thetas.hat.dboot[k,] = bhat(Xy.boot2)[2:5,1]
			}
		for(coefs in 1:4){
			theta.qtl.lgrid.lo[j,,coefs] = quantile(thetas.hat.dboot[,coefs], l.grid.lo)
			theta.qtl.lgrid.hi[j,,coefs] = quantile(thetas.hat.dboot[,coefs], l.grid.hi)		
			}
		}
	Resid = resids(Xy,bhats)
	theta.se = se.bhat(Xm, Resid)[2:5]	    # this is the traditional SE estimate.
	theta.se.sand = sand(Xm, Resid)[2:5]      # this is the sandwich SE estimate.
	
	z.int = cbind(thetas.hat+qnorm(alpha)*theta.se, thetas.hat+qnorm(1-alpha)*theta.se)    # first the simple intervals.
	z.int.sand = cbind(thetas.hat+qnorm(alpha)*theta.se.sand, thetas.hat+qnorm(1-alpha)*theta.se.sand)
	
	perc.int = cbind(apply(thetas.hat.boot,2,quantile,alpha),apply(thetas.hat.boot,2,quantile,1-alpha))     # next the simple percentile interval.
	
	for(coefs in 1:4){
		log.lgrid.lo[,,coefs] = theta.qtl.lgrid.lo[,,coefs] < thetas.hat[coefs]        # finally the double BS calibrated percentile interval.
		log.lgrid.hi[,,coefs] = theta.qtl.lgrid.hi[,,coefs] > thetas.hat[coefs]        
		percs.lo[coefs,] = 1-colSums(log.lgrid.lo[,,coefs])/B
		percs.hi[coefs,] = 1-colSums(log.lgrid.hi[,,coefs])/B
		lo.log = (percs.lo[coefs,1:(G-1)]<alpha)*(percs.lo[coefs,2:G]>=alpha)
		hi.log = (percs.hi[coefs,1:(G-1)]<alpha)*(percs.hi[coefs,2:G]>=alpha)
		if(sum(lo.log)>0){l.alpha.lo = l.grid.lo[lo.log*(1:(G-1))]}
		if(sum(lo.log)<=0){l.alpha.lo = .0001}
		if(sum(hi.log)>0){l.alpha.hi = l.grid.hi[hi.log*(1:(G-1))]}
		if(sum(hi.log)<=0){l.alpha.hi = 1-.0001}
		perc.cal.int[coefs,] = quantile(thetas.hat.boot[,coefs],c(l.alpha.lo, l.alpha.hi))
		}
				
	# t2 = Sys.time()
	
	print(Sys.time())
	
	return(list(thetas.hat=thetas.hat, theta.se=theta.se, theta.se.sand=theta.se.sand, z.int=z.int, z.int.sand=z.int.sand,
		perc.int=perc.int, perc.cal.int=perc.cal.int, thetas.true=thetas.true))
		
	}


	
	
	





set.seed(1600*(iter-1)+i)      # reproducibility / error checking


# Global settings

alpha = .05
G = 50
B = 100						# originally 1000
B2 = 100 					# originally 1000

# Run over all cases.

cases_mat = indexes.n.p(15,2,3)
colnames(cases_mat) = c("betas","ns","efuncs")
ns = c(25,100)
betas_mat = indexes.n.p(2,2,2,2, num.params=4)-1
betas_mat = betas_mat[-1,]
case = cases_mat[iter,]
betas = betas_mat[case[1],]
n = ns[case[2]]
output = tfunc(betas, efunc=case[3], n=n, alpha=alpha, G = G, B = B, B2 = B2)
save(output, file=paste("output/output", iter, i, "mat", sep="."))



