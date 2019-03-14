
# Estimating time uncertainty and signal-to-noise ratio (UnSNR)


rm(list=ls())


# libraries & single.functions
	#
	library(zoo)
	source("~/additionalFunctions.R")



### data => c("legend_meta","M13_meta","M13_Prx","M13_MPI6k","M13_T21k","M13_rawMPI6k","M13_rawT21k")
	#
	load("~/M13.RData")
		metadat <- M13_meta
		prxdat <- M13_Prx
		moddat <- M13_MPI6k
		modtemp <- M13_rawMPI6k
	
	

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
	bin <- c(0,5000)
	binseq <- c(2500)
	#
	c.timeuncertainty <- seq(from=0,to=400,by=50) 
	c.snr <- c(0.01,0.02,0.05,0.08,0.1,0.2,0.5,0.8,1,2,5,8,10,20,50,80,100,0)
	points <- c(2000,4000,6000)
	unsnr.dist <- 5000
	


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
	#
	outTime <- seq(from=filter.dt, to=length.timser, by=filter.dt)



### estimations
	#
	#CorMatrix for real proxy data
		#
		Prx.split <- spatial.cor.analysis.splitvec(data=prxdat, metadat=metadat, detrending=detrending, tsc.in=tsc.in, filt.position=filter.dt, int.method=int.method, startTimser=filter.dt, endTimser=length.timser, null.hyp=FALSE, quant=c(0.05,0.95), powlaw.beta=beta.powerlaw, rep, split=TRUE, splitlist=splitvec, binning=TRUE, bin=bin, binseq=binseq, na=TRUE, nalist=naList)
	#
	#CorMatrix of modified model data
		#
		modified.modCorMat <- modified.modCorMat.min <- modified.modCorMat.max <- array(NA, c(length(prxdat), length(prxdat), length(c.snr), length(c.timeuncertainty)))
		#
		for (i.timeuncertainty in 1:length(c.timeuncertainty)){
   			print(paste("Working on time uncertainty ",i.timeuncertainty))
   			for (i.snr in 1:length(c.snr)){
   				print(paste("Working on snr ",i.snr))
    			    tmp.CorMat <- array(NA, c(length(prxdat), length(prxdat), rep))
    			    for (i.rep in 1:rep){
        			print(paste("rep= ",i.rep,sep=""))
        			tmp.modified.moddat <- matrix(NA, length(outTime), length(prxdat))
            			for (i.Proxy in 1:length(prxdat)){
            			tmp.subsampled.temp <- SubsampleTimeseriesBlock(ts=modtemp[[i.Proxy]], timepoints=index(prxdat[[i.Proxy]]))
            			tmp.modified.moddat[,i.Proxy] <- modifyIrregData(irreg.temperature=tmp.subsampled.temp, irreg.age=index(prxdat[[i.Proxy]]), time.un=c.timeuncertainty[i.timeuncertainty], snr=c.snr[i.snr], points=points, tsc.in=tsc.in, detrending=detrending, filt.position=filter.dt, int.method=int.method, startTimser=filter.dt, endTimser=length.timser)
            			}
            			tmp.CorMat[,,i.rep] <- corMatrix(tmp.modified.moddat)
            			for(u in 1:length(naList)){
						tmp.CorMat[,naList[u],i.rep] <- NA
						tmp.CorMat[naList[u],,i.rep] <- NA
					}     
        			}
			modified.modCorMat[,,i.snr,i.timeuncertainty] <- apply(tmp.CorMat, c(1,2), mean, na.rm=TRUE) 
			modified.modCorMat.min[,,i.snr,i.timeuncertainty] <- apply(tmp.CorMat, c(1,2), min, na.rm=TRUE)
			modified.modCorMat.max[,,i.snr,i.timeuncertainty] <- apply(tmp.CorMat, c(1,2), max, na.rm=TRUE) 
		    }
		}



### post-calculations
	#
	dist.matrix <- Prx.split$dist
	unsnr.dist.pos <- which(dist.matrix <= unsnr.dist)
	#
	splitvec.modified.modCorMat <- list()
		for(l in 1:length(splitvec)){
			splitvec.modified.modCorMat[[l]] <- array(NA, dim=c(length(metadat[[2]]), length(metadat[[2]]), length(c.snr), length(c.timeuncertainty)))
		}
		for(i.timeuncertainty in 1:length(c.timeuncertainty)){
	 		print(paste("Working on time uncertainty ",i.timeuncertainty))
			for(i.snr in 1:length(c.snr)){
				print(paste("Working on snr ",i.snr))
				tmp <- splitMat(modified.modCorMat[,,i.snr,i.timeuncertainty], splitvec)
				#
				splitvec.modified.modCorMat[[1]][,,i.snr,i.timeuncertainty] <- tmp[[1]]
				splitvec.modified.modCorMat[[2]][,,i.snr,i.timeuncertainty] <- tmp[[2]]
				splitvec.modified.modCorMat[[3]][,,i.snr,i.timeuncertainty] <- tmp[[3]]
				splitvec.modified.modCorMat[[4]][,,i.snr,i.timeuncertainty] <- tmp[[4]]
				splitvec.modified.modCorMat[[5]][,,i.snr,i.timeuncertainty] <- tmp[[5]]
				splitvec.modified.modCorMat[[6]][,,i.snr,i.timeuncertainty] <- tmp[[6]]
			}
		}
	#
	misfit.unsnrMat.all <- misfit.unsnrMat.Uk37 <- misfit.unsnrMat.MgCa <- misfit.unsnrMat.TEX86 <- misfit.unsnrMat.bioind <- misfit.unsnrMat.ice <- misfit.unsnrMat.other <- array(NA, dim=c(length(c.snr),length(c.timeuncertainty)))
#
for(i.timeuncertainty in 1:length(c.timeuncertainty)){
	print(paste("Working on time uncertainty ",i.timeuncertainty))
	for(i.snr in 1:length(c.snr)){
		print(paste("Working on snr ",i.snr))
		misfit.unsnrMat.Uk37[i.snr,i.timeuncertainty] <- abs((mean(splitvec.modified.modCorMat[[1]][,,i.snr,i.timeuncertainty][unsnr.dist.pos], na.rm=TRUE) - mean(Prx.split$CorMat[[1]][unsnr.dist.pos], na.rm=TRUE)))
		misfit.unsnrMat.MgCa[i.snr,i.timeuncertainty] <- abs((mean(splitvec.modified.modCorMat[[2]][,,i.snr,i.timeuncertainty][unsnr.dist.pos], na.rm=TRUE) - mean(Prx.split$CorMat[[2]][unsnr.dist.pos], na.rm=TRUE)))
		misfit.unsnrMat.TEX86[i.snr,i.timeuncertainty] <- abs((mean(splitvec.modified.modCorMat[[3]][,,i.snr,i.timeuncertainty][unsnr.dist.pos], na.rm=TRUE) - mean(Prx.split$CorMat[[3]][unsnr.dist.pos], na.rm=TRUE)))
		misfit.unsnrMat.bioind[i.snr,i.timeuncertainty] <- abs((mean(splitvec.modified.modCorMat[[4]][,,i.snr,i.timeuncertainty][unsnr.dist.pos], na.rm=TRUE) - mean(Prx.split$CorMat[[4]][unsnr.dist.pos], na.rm=TRUE)))
		misfit.unsnrMat.ice[i.snr,i.timeuncertainty] <- abs((mean(splitvec.modified.modCorMat[[5]][,,i.snr,i.timeuncertainty][unsnr.dist.pos], na.rm=TRUE) - mean(Prx.split$CorMat[[5]][unsnr.dist.pos], na.rm=TRUE)))
		misfit.unsnrMat.other[i.snr,i.timeuncertainty] <- abs((mean(splitvec.modified.modCorMat[[6]][,,i.snr,i.timeuncertainty][unsnr.dist.pos], na.rm=TRUE) - mean(Prx.split$CorMat[[6]][unsnr.dist.pos], na.rm=TRUE)))
		misfit.unsnrMat.all[i.snr,i.timeuncertainty] <- abs((mean(modified.modCorMat[,,i.snr,i.timeuncertainty][unsnr.dist.pos], na.rm=TRUE) - mean(Prx.split$CorMat[[7]][unsnr.dist.pos],na.rm=TRUE)))
	}
}

