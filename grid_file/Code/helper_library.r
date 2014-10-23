


append_intercept = function(x){
	# Give me a design matrix (doesn't need to be a vector), and I'll stick a column of 1's on the LHS.
	x = as.matrix(x)
	n = nrow(x)
	with_int = cbind(rep(1, n), x)
	with_int
	}

beta_hat_from_X_and_y = function(X, y, epsilon=.00001){
	# Give me a design matrix X and the response y, and I'll give you back the beta_hat vector.
	# (I add in an epsilon of error to ensure invertability. 
	beta_hat = solve(t(X)%*%X+diag(epsilon,2))%*%t(X)%*%as.matrix(y)
	beta_hat
	}




	
	
	
	
#################################################################################################################
	
RMSE = function(y, y_hat, dfs = 2, boot.SE = FALSE){
	# Give me predictions and actuals and I'll give you back the RMSE.
	# I assume the DF adjustment is 2 (simple linear regression) but allow for others too (ie. for the mean, it'd be 1)
	n = length(y)
	rmse = sqrt(sum((y-y_hat)^2) / (n - dfs))
	if(boot.SE==FALSE){return(rmse)}
	if(boot.SE==TRUE){
		B = 10000
		rmse_boot = rep(NA,B)
		for(i in 1:B){
			boot_ind = sample(1:n, size=n,replace=TRUE)
			rmse_boot[i] = sqrt(sum((y[boot_ind]-y_hat[boot_ind])^2) / (n - dfs))
			}
		se_rmse_boot = sd(rmse_boot)
		return(list(rmse=rmse,se_rmse_boot=se_rmse_boot))
		}
	}
	
	
	
#################################################################################################################

gamma_ratio = function(a,b){
	# inputs: a and b for gamma(a+b)/gamma(a)
	# output: evaluation of ratio - when elements are big enough, switches to approximation.
	# accepts vector valued output - assumes either none, one or both elements of ratio are scalar values
	# THIS DOES NOT HELP WHEN a+b LARGE, BUT a IS NOT LARGE - ASYmPTOTICS DO NOT KICK IN, SO USE OLD METHOD.
	n1 = length(a); n2 = length(b)
	if(n1!=n2){
		if(n1 > n2){
			b = rep(b,n1)
			}
		if(n1 < n2){
			a = rep(a,n2)
			}
		}
	n = length(a)
	logvec = (a>14)*((a+b)>14)*(1:n)
	output = rep(NA,n)
	output[logvec] = sqrt((a[logvec]+b[logvec]-1)/(a[logvec]-1))*
						(1+b[logvec]/(a[logvec]-1))^(a[logvec]-1)*((a[logvec]+b[logvec]-1)/exp(1))^b[logvec]
	output[logvec==0] = gamma(a[logvec==0]+b[logvec==0])/gamma(a[logvec==0])
	return(output)
	}

	
	
#################################################################################################################
	
MAPE = function(y.hat, y){
	# give me your actual and expected, and I'll return mean absolute prediction error.
	is.na.logvec = (1-is.na(y.hat)|is.na(y))*(1:length(y))  # subset out any NA's.
	sqrt( mean( abs(y.hat[is.na.logvec] - y[is.na.logvec])) )
	}

	
	


#################################################################################################################

pop_least_squares_beta_func = function(nsim=10000, Xfunc=1, sample_size){
	xtx_mat_1 = matrix(0,nrow =  2, ncol = nsim)
	xtx_mat_2 = matrix(0,nrow =  2, ncol = nsim)
	xty_mat = matrix(0,nrow =  2, ncol = nsim)
	for(i in 1 : nsim){
    
		if(Xfunc==1){x = rnorm(n = sample_size)}
		if(Xfunc==2){x = abs(rnorm(n = sample_size))}
		if(Xfunc==3){x = exp(rnorm(n = sample_size))}
		y = x^2
		X = append_intercept(x)

		xtx = t(X)%*%X
		xty = t(X)%*%y

		xtx_mat_1[,i] = xtx[,1]
		xtx_mat_2[,i] = xtx[,2]
		xty_mat[,i] = xty

		# if(i %%100 == 0){print(i)}
		}  
	E_xty = apply(xty_mat, 1, mean)
	E_xtx_1 = apply(xtx_mat_1, 1, mean)
	E_xtx_2 = apply(xtx_mat_2, 1, mean)

	E_xtx_inv = solve(cbind(E_xtx_1, E_xtx_2))

	beta =   E_xtx_inv %*%    E_xty
	beta
	}


	
	
#################################################################################################################
	
many_cols = function(vec, nreps){
	# you give me a vector and the number of replications you want, stacked together as columns side by side, 
	# and I'll spit you back that matrix.

	return(t(rep(1,nreps) %*% t.default(vec)))
	}


#################################################################################################################

	
many_rows = function(vec, nreps){
	# you give me a vector and the number of replications you want, stacked together as rows on top of one another,
	# and I'll spit you back that matrix.

	return(rep(1,nreps) %*% t.default(vec))
	}
	



#################################################################################################################

# Error function
erf     <- function(x) {  # 2*pnorm(sqrt(2)*x)-1
    pchisq(2*x^2,1)*sign(x)
}

# Inverse error function
erfinv  <- function(y) {
        y[abs(y) > 1] <- NA
        sqrt(qchisq(abs(y),1)/2) * sign(y)
}

# Complementary error function
erfc    <- function(x) {  # 1 - erf(x)
    2 * pnorm(-sqrt(2) * x)
}

# Inverse complementary error function
erfcinv <- function(y) {
    y[y < 0 | y > 2] <- NA
    -qnorm(y/2)/sqrt(2)
}

# Scaled complementary error function
erfcx   <- function(x) {
    exp(x^2) * erfc(x)
}

# Complex error function
erfz    <- function(z)
{
    if (is.null(z)) return( NULL )
    else if (!is.numeric(z) && !is.complex(z))
        stop("Argument 'z' must be a numeric or complex scalar or vector.")

    a0 <- abs(z)
    c0 <- exp(-z * z)
    z1 <- ifelse (Re(z) < 0, -z, z) 

	i <- a0 <= 5.8
	work.i <- i
	cer <- rep(NA, length = length(z))
    if ( sum(work.i) > 0) {
        cs <- z1
        cr <- cs
        for (k in 1:120) {
            cr[work.i] <- cr[work.i] * z1[work.i] * z1[work.i]/(k + 0.5)
            cs[work.i] <- cs[work.i] + cr[work.i]
            work.i <- work.i & (abs(cr/cs) >= 1e-15)
	    if (sum(work.i) == 0) break
        }
        cer[i] <- 2 * c0[i] * cs[i]/sqrt(pi)
    }
	work.i <- !i
    if( sum(work.i) > 0) {
        cl <- 1/z1
        cr <- cl
        for (k in 1:13) {
            cr[work.i] <- -cr[work.i] * (k - 0.5)/(z1[work.i] * z1[work.i])
            cl[work.i] <-  cl[work.i] + cr[work.i]
            work.i <- work.i & (abs(cr/cl) >= 1e-15)
	    if (sum(work.i) == 0) break
        }
        cer[!i] <- 1 - c0[!i] * cl[!i]/sqrt(pi)
    }
	cer[ Re(z) < 0] <- -cer[ Re(z) < 0]
    return(cer)
}


# Imaginary error function
erfi <- function(z) {
	-1i * erfz(1i * z)
}








#################################################################################################################

match_mat = function(mat1, mat2){
	# let mat2 be the dictionary to match against (r2 by c2), let mat1 be the big matrix with dupes (r1 by c1).  
	# return a vector of length r1 whose elements are the row numbers in mat2 corresponding to each of the rows
	# in mat1 (thus, something between 1 and r2).

	match(data.frame(t(mat1)), data.frame(t(mat2)))
	}


	
	
	
	
#################################################################################################################

match_vec_vs_mat = function(vec, mat){
	# You give me a vector of length p, and a matrix of dimension n by p, and I'll give you back all the rows
	# in the matrix which exactly match that vector.

	row_inds = which(rowSums(mat != rep(vec, each=nrow(mat))) == 0)
	return(row_inds)
	}



	
#################################################################################################################
	
rep_by_element = function(vec, n){
	# You give me a vector and the number of times you want things repeated and I'll give you back a vector
	# where each of the elements are repeated, element by element.  ie: rep_by_element(c(1,3,2),2) = 113322.	

	unlist(lapply(vec,rep,n))
	}



	
#################################################################################################################
	
choose2 = function(n, k){
	# choose2(n,k) is the same as n choose k built into R, except if n and k are negative, this 
	# evaluates to 1.
	
	max(choose(n,k),1)
	
	# ifelse(n<0,ifelse(k<0,1,choose(n,k)),choose(n,k))
	}




#################################################################################################################

smooth.lin <- function(x, y, bw=1/5, nloc=999) {
	# You give me x and y, and I'll give you a local linear kernel bandwith smooth.
	# Input: x and y, the bandwith and the number of grid points to split x up into.
	# Output: x grid and y.hats's evaluated at those x's.
	xmin <- min(x);   xmax <- max(x)   # needed for the bandwidth
	xs <- seq(xmin, xmax, length=nloc) # the grid points
	fx <- rep(NA, nloc)                # allocate space for the fits
	act.bw <- bw * (xmax-xmin)         # the actual bandwidth
	for(i in 1:nloc) {                 # for all grid points ...
		w <- dnorm(x=x-xs[i],sd=act.bw)  # normal or Gaussian kernel weights
		w <- w/sum(w)                    # force sum(w)=1
		xbar <- sum(w*x)                 # see the formulas above...
		ybar <- sum(w*y)
		x.y  <- sum(w*(x-xbar)*(y-ybar))
		x.x  <- sum(w*(x-xbar)*(x-xbar))
		fx[i] <- x.y/x.x*(xs[i]-xbar) + ybar # evaluate the line at the grid point
		}
	xf <- cbind(xs,fx)     # Return the grid and its fitted values
	#return(xf)
	return(list(x=xs,y=fx))
	}

	

	




#################################################################################################################
# (Before I learned about expand.grid() )


# For enumerating all scenarios
indexes.n.p <- function(n1=2,n2=2,n3=2,n4=2,n5=2,n6=2,n7=2,n8=2,n9=2,n10=2,n11=2,num.params=3){
	A <- rep(1:n1,n2)
	B <- 0
	for(i in 1:n2){
		B <- c(B,rep(i,n1))
		}
	B <- B[-1]
	AB <- cbind(A,B)

	if(num.params>2){
		AB.n <- c(0,0)	  
		for(i in 1:n3){
			AB.n <- rbind(AB.n,AB)
			}
		AB.n <- AB.n[-1,]
		C <- 0
		for(i in 1:n3){
			C <- c(C,rep(i,n1*n2))
			}
		C <- C[-1]
		ABC <- cbind(AB.n,C)
		if(num.params>3){
			ABC.n <- c(0,0,0)
			for(i in 1:n4){
				ABC.n <- rbind(ABC.n,ABC)
				}
			ABC.n <- ABC.n[-1,]
			D <- 0
			for(i in 1:n4){
				D <- c(D,rep(i,n1*n2*n3))
				}
			D <- D[-1]
			ABCD <- cbind(ABC.n,D)		
			
			if(num.params>4){
				ABCD.n <- c(0,0,0,0)
				for(i in 1:n5){
				  ABCD.n <- rbind(ABCD.n,ABCD)
				  }
				ABCD.n <- ABCD.n[-1,]

				E <- 0
				for(i in 1:n5){
				  E <- c(E,rep(i,n1*n2*n3*n4))
				  }
				E <- E[-1]
				ABCDE <- cbind(ABCD.n,E)

				if(num.params>5){
					ABCDE.n <- c(0,0,0,0,0)
					for(i in 1:n6){
					  ABCDE.n <- rbind(ABCDE.n,ABCDE)
					  }
					ABCDE.n <- ABCDE.n[-1,]

					F <- 0
					for(i in 1:n6){
					  F <- c(F,rep(i,n1*n2*n3*n4*n5))
					  }
					F <- F[-1]
					ABCDEF <- cbind(ABCDE.n,F)
					
					if(num.params>6){	
						ABCDEF.n <- c(0,0,0,0,0,0)
						for(i in 1:n7){
						  ABCDEF.n <- rbind(ABCDEF.n,ABCDEF)
						  }
						ABCDEF.n <- ABCDEF.n[-1,]

						G <- 0
						for(i in 1:n7){
						  G <- c(G,rep(i,n1*n2*n3*n4*n5*n6))
						  }
						G <- G[-1]
						ABCDEFG <- cbind(ABCDEF.n,G)
						
						if(num.params>7){
							ABCDEFG.n <- c(0,0,0,0,0,0,0)
							for(i in 1:n8){
							  ABCDEFG.n <- rbind(ABCDEFG.n,ABCDEFG)
							  }
							ABCDEFG.n <- ABCDEFG.n[-1,]

							H <- 0
							for(i in 1:n8){
							  H <- c(H,rep(i,n1*n2*n3*n4*n5*n6*n7))
							  }
							H <- H[-1]
							ABCDEFGH <- cbind(ABCDEFG.n,H)
							
							if(num.params>8){
								ABCDEFGH.n <- c(0,0,0,0,0,0,0,0)
								for(i in 1:n9){
								  ABCDEFGH.n <- rbind(ABCDEFGH.n,ABCDEFGH)
								  }
								ABCDEFGH.n <- ABCDEFGH.n[-1,]

								I <- 0
								for(i in 1:n9){
								  I <- c(I,rep(i,n1*n2*n3*n4*n5*n6*n7*n8))
								  }
								I <- I[-1]
								ABCDEFGHI <- cbind(ABCDEFGH.n,I)

								if(num.params>9){
									ABCDEFGH.n <- c(0,0,0,0,0,0,0,0,0)
									for(i in 1:n10){
									  ABCDEFGHI.n <- rbind(ABCDEFGHI.n,ABCDEFGHI)
									  }
									ABCDEFGHI.n <- ABCDEFGHI.n[-1,]

									J <- 0
									for(i in 1:n10){
									  J <- c(J,rep(i,n1*n2*n3*n4*n5*n6*n7*n8*n9))
									  }
									J <- J[-1]
									ABCDEFGHIJ <- cbind(ABCDEFGHI.n,J)

									if(num.params>10){
										ABCDEFGHI.n <- c(0,0,0,0,0,0,0,0,0,0)
										for(i in 1:n11){
										  ABCDEFGHIJ.n <- rbind(ABCDEFGHIJ.n,ABCDEFGHIJ)
										  }
										ABCDEFGHIJ.n <- ABCDEFGHIJ.n[-1,]

										K <- 0
										for(i in 1:n11){
										  K <- c(K,rep(i,n1*n2*n3*n4*n5*n6*n7*n8*n9*n10))
										  }
										K <- K[-1]
										ABCDEFGHIJK <- cbind(ABCDEFGHIJ.n,K)
										}
									
									}
								}
							}
						}
					}				
				}
			}
		}
		
	if(num.params==2){return(AB)}
	if(num.params==3){return(ABC)}
	if(num.params==4){return(ABCD)}
	if(num.params==5){return(ABCDE)}
	if(num.params==6){return(ABCDEF)}
	if(num.params==7){return(ABCDEFG)}
	if(num.params==8){return(ABCDEFGH)}
	if(num.params==9){return(ABCDEFGHI)}
	if(num.params==10){return(ABCDEFGHIJ)}
	if(num.params==11){return(ABCDEFGHIJ)}
	}


	
	
	
	
	
# For enumerating binary sequences - you give me the length, lenth = {1, 2, 3, ...} and I'll spit you back the 
# listing of them all.

enum_binary = function(len){
	output = c(0,1)
	if(len>1){
		for(i in 2:len){
			output0 = cbind(output,rep(0,2^(i-1)))
			output1 = cbind(output,rep(1,2^(i-1)))
			output = rbind(output0,output1)
			}
		}
	return(output)
	}
	
	
		
	