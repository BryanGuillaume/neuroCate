This folder contains the shell and R scripts/functions that were used for the Monte Carlo simulations of the "neuroCate paper" currently in submission.

They were intended to work on a cluster with a PBS scheduler using the following commands:

	qsub -J 1-10000 simuCompaCateStrategy_0Z_000alpha.sh 
	qsub -J 1-10000 simuCompaCateStrategy_0Z_010alpha.sh 
	qsub -J 1-10000 simuCompaCateStrategy_0Z_030alpha.sh 
	qsub -J 1-10000 simuCompaCateStrategy_2Z_000alpha.sh 
	qsub -J 1-10000 simuCompaCateStrategy_2Z_010alpha.sh 
	qsub -J 1-10000 simuCompaCateStrategy_2Z_030alpha.sh 

Each shell script was paired with the R script "simuCompaCateStrategy.R". In order to work, this R script required the R package "cate" in version 1.1.1 and some custom R functions that are also in this folder:

	myFaPC.R
	myFaPXEM.R
	myCate.R

Each run of the script "simuCompaCateStrategy.R" produced the results for the 16 noise distributions described in the paper. It was taking 3 inputs from the shell scripts, being the index of the simulation realisation, the number of unknown covariates modelled and the level of confounding of the unknow covariates. 



