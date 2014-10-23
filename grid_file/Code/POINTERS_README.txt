#Kilroy wuz here


# HANDY COMMANDS
# 1. Delete all jobs?					qdel -u danielmc
# 2. See status of current jobs?		qstat
# 3. Don't output all e and o logs		-o /dev/null
	# i.e.,  echo "export ITER=$ITER; echo "running ITER: $ITER TASK: \$SGE_TASK_ID"; R --no-save < loop.r" | qsub -N moms-$ITER -o /dev/null -j y -t 1-1000 -cwd
    # Only do this once you are perfectly confident the code runs without error! No error checking, just output.	


# MY PROCESS FOR RUNNING A JOB
# 1. CREATE loop.r, loop.sh and combine.r FILES
# 		-- loop.r - this is the R code which is looped over, many times, saving output at the bottom
#		-- loop.sh - this is a shell script which says "run loop.r many times with this variable in it"
#		-- combine.r - this is the R code which takes the output from all loops and combines it together again.
# 2. RUN loop.sh 
# 3. CHECK INTERMEDIATE OUTPUT FREQUENTLY
#   	-- Look at e and o files
#		-- Run intermediate output
# 4. PUT ALL OUTPUT INTO A LOCAL FOLDER AND RUN ANALYSIS ON OUTPUT FROM THERE.

	
	
# POINTERS
# -- *Save intermediate output.*  Trust me, tremendous increase in error checking, making sure results as expected.
# -- *Ask Hugh if you want more nodes.*  Default is 50, I've gotten 230 at least a few times.
# -- *Run small jobs first.* Be humble, it will probably bomb somewhere the first few times.
# -- *Check your log files often.* It is of paramount importance to make sure things on track.
# -- *Save output to separate folder.* The main network drive will slow down otherwise.
# -- *Don't use for loops except in tiny amounts.* It will crash the job scheduler.
# -- *Don't crash the grid.* Be super careful when you need to do huge jobs...





# EXECUTING A JOB
# -- Go to PuTTY
# 	-- Host Name:			 	danielmc@unix.wharton.upenn.edu
# 	-- (Port: 					22)
# 	-- (Connection type: 		SSH)
































