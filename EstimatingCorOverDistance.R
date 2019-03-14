
# Estimating the correlation over distance


rm(list=ls())


# libraries & single.functions
	#
	library(zoo)
	source("~/additionalFunctions.R")



### data => c("legend_meta","M13_meta","M13_Prx","M13_MPI6k","M13_T21k","M13_rawMPI6k","M13_rawT21k")
	#
	load("~/M13.RData")
	
	

### parameter
	#
	tsc.in <- 1000
	detrending = FALSE		# TRUE for tsc.in <- 400
	mean.res <- tsc.in/2
	int.method = "linear"
	filter.dt <- 10
	length.timser <- 6000
	beta.powerlaw <- 1
	rep <- 1000
	#
	bin <- seq(from=0, to=20000, by=2000)
	binseq <- seq(from=1000, to=19000, by=2000)

	
	
### pre-calculations
	#
	naList <- which(M13_meta[[4]] > mean.res)
	#
	splitvec <- list()
		splitvec[[1]] <- which(M13_meta[[5]] == "Uk37")
		splitvec[[2]] <- which(M13_meta[[5]] == "Mg/Ca")
		splitvec[[3]] <- which(M13_meta[[5]] == "TEX86")
		splitvec[[4]] <- which(M13_meta[[5]] == "bioindicator")
		splitvec[[5]] <- which(M13_meta[[5]] == "ice")
		splitvec[[6]] <- which(M13_meta[[5]] == "other")
	
	
	
### estimations  => Output: $dist,$CorMat,$qCorMat,$CorFkt,$qCorFkt,$QuantCorFkt,$numCorFkt
	#
	Prx.split <- spatial.cor.analysis.splitvec(
					data = M13_Prx, 
					metadat = M13_meta,
					detrending = detrending,
					tsc.in = tsc.in,
					filt.position = filter.dt,
					int.method = int.method,
					startTimser = filter.dt,
					endTimser = length.timser,
					null.hyp = TRUE,
					quant = c(0.05,0.95),
					powlaw.beta = beta.powerlaw,
					rep = rep,
					split = TRUE,
					splitlist = splitvec,
					binning = TRUE,
					bin = bin,
					binseq = binseq,
					na = TRUE,
					nalist = naList) 
	#
	MPI6k.split <- spatial.cor.analysis.splitvec(
					data = M13_MPI6k, 
					metadat = M13_meta,
					detrending = detrending,
					tsc.in = tsc.in,
					filt.position = filter.dt,
					int.method = int.method,
					startTimser = filter.dt,
					endTimser = length.timser,
					null.hyp = FALSE,
					split = TRUE,
					splitlist = splitvec,
					binning = TRUE,
					bin = bin,
					binseq = binseq,
					na = TRUE,
					nalist = naList) 
	#
	T21k.split <- spatial.cor.analysis.splitvec(
					data = M13_T21k, 
					metadat = M13_meta,
					detrending = detrending,
					tsc.in = tsc.in,
					filt.position = filter.dt,
					int.method = int.method,
					startTimser = filter.dt,
					endTimser = length.timser,
					null.hyp = FALSE,
					split = TRUE,
					splitlist = splitvec,
					binning = TRUE,
					bin = bin,
					binseq = binseq,
					na = TRUE,
					nalist = naList) 


