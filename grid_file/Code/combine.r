

# For a particular iter, this recombines all associated i's, then returns the desired matrix.
# Take output/output.iter.i for all i, load them in, then combine them.




# Dropbox 
# setwd("C:/Users/Dan/Dropbox/Work/phD File/Stat Research/Publicly shared stuff/Grid")

source("helper_library.r")

epsilon = .00001


# Now we need to convert this into a format that throws all into one matrix
output.matrix = function(output){
	# unpack the output into what we need.
	thetas.true = output$thetas.true
	z.int = output$z.int
	z.int.sand = output$z.int.sand
	perc.int = output$perc.int
	perc.cal.int = output$perc.cal.int
	thetas.hat = output$thetas.hat
	
	
	# make indicator vectors for coverage.
	logvec.z.int = (many_cols(thetas.true,2) > z.int[,1])*(many_cols(thetas.true,2) < z.int[,2])
	logvec.z.int.sand = (many_cols(thetas.true,2) > z.int.sand[,1])*(many_cols(thetas.true,2) < z.int.sand[,2])
	logvec.perc.int = (many_cols(thetas.true,2) > perc.int[,1])*(many_cols(thetas.true,2) < perc.int[,2])
	logvec.perc.cal.int = (many_cols(thetas.true,2) > perc.cal.int[,1])*(many_cols(thetas.true,2) < perc.cal.int[,2])

	}

	





# Get the specifications for all cases.  

cases_mat = indexes.n.p(15,2,3)
colnames(cases_mat) = c("betas","ns","efuncs")
ns = c(25,100)
betas_mat = indexes.n.p(2,2,2,2, num.params=4)-1
betas_mat = betas_mat[-1,]
colnames(betas_mat) = c("b1","b2","b3","b4")

	
	

ncases = 90
output_mat = matrix(NA,nrow=1,ncol=16)
output_mat_avg = matrix(NA,nrow=1,ncol=4)
method_vec = c("z","z-sand","perc","perc-cal")		
colnames(output_mat) = rep_by_element(c(paste("covg%:",method_vec,sep="")),4)
colnames(output_mat_avg) = c(paste("covg%:",method_vec,sep=""))

# Run over all i's for a particular value of iter.

# for(iter in 1:ncases){
for(iter in 1:2){

	nsim = 2
	
	for(i in 1:nsim){
	
		load(paste("output/output",iter,i,"mat",sep="."))
		
		thetas.true = output$thetas.true
		z.int = output$z.int
		z.int.sand = output$z.int.sand
		perc.int = output$perc.int
		perc.cal.int = output$perc.cal.int
		thetas.hat = output$thetas.hat

		# make indicator vectors for coverage.
		
		new_logvec.z.int = as.numeric((thetas.true > z.int[,1])+(thetas.true < z.int[,2]))==2
		new_logvec.z.int.sand = as.numeric((thetas.true > z.int.sand[,1])+(thetas.true < z.int.sand[,2]))==2
		new_logvec.perc.int = as.numeric((thetas.true > perc.int[,1])+(thetas.true < perc.int[,2]))==2
		new_logvec.perc.cal.int = as.numeric((thetas.true > perc.cal.int[,1])+(thetas.true < perc.cal.int[,2]))==2
		
		new_logvec = c(new_logvec.z.int, new_logvec.z.int.sand, new_logvec.perc.int, new_logvec.perc.cal.int)
		
		if(i==1){
			logvec = new_logvec
			}
		if(i>1){
			logvec = rbind(logvec,new_logvec)
			}	
		}

	colnames(logvec) = rep_by_element(method_vec, 4)       # Reporting results coefficient by coefficient
	new_output_covg = apply(logvec,2,mean)        # this is the average coverage for the 4 slope coefficients
	new_output_SE = apply(logvec,2,sd)/sqrt(nsim)        # this is the average coverage for the 4 slope coefficients
	output_mat = rbind(output_mat, new_output_covg)
	}

	
	
output_mat = output_mat[-1,]
	
results_full = cbind(1:ncases,cases_mat[1:ncases,],output_mat)	
results_avg = cbind(1:ncases,cases_mat[1:ncases,],output_mat_avg)	
rownames(results_full) = rep("",ncases)
rownames(results_avg) = rep("",ncases)
results_full
round(results_avg,2)
apply(output_mat_avg,2,mean)

write.csv(output_mat,file="output_mat.csv")
write.csv(output_mat_avg,file="output_mat_avg.csv")
write.csv(results_avg,file="results_avg.csv")







