
# additionalFunctions


# -1-

addNoise <- function(temperature, snr)
	#Input:
	#	temperature: temperature of the time series
	#	snr: signal-to-noise-ratio which is used to add defined noise
    {
    	if(snr==0){
    		noisy.temperature <- temperature
    	}else{
    		noisy.temperature <- (temperature + rnorm(length(temperature), mean=0, sd=sqrt(var(temperature, na.rm=TRUE)/snr)))
    	}
        #
        return(noisy.temperature)         
    }



# -2-

approx.nearest <- function(x, y, xi)
	#Input:
	#	x: numeric vector with time coordinates of the points to be interpolated
	#	y: corresponding observations to x
	#	xi: set of numeric values specifying where interpolation should take place
{
	result <- list()
	result$x <- xi
	result$y <- approx(c(x[1], x+c(diff(x)/2, 0)), c(y[1], y), xi, method="constant", f=1)$y
	output <- zoo(result$y, order.by=result$x)
	#
	return(output)
}



# -3-

corMatrix <- function(proxy_matrix)
	#Input:
	#	proxy_matrix: matrix with data whose correlation among each dataset pair should be calculated (time, numberOfTimeseries)
{
    t <- cor(proxy_matrix, use = "pairwise.complete")
    t[lower.tri(t)] <- NA
    diag(t) <- NA
    return(t)
}



# -4-

detr.timser <- function(time.series, detr=TRUE)
	#Input:
	#	time.series: zoo-object
	#	detr: linear detrending TRUE or FALSE
{
	if(detr==TRUE){
		result <- zoo(lm(coredata(time.series)~index(time.series), na.action=na.omit)$residuals, order.by=index(time.series))
	}else{
		result <- time.series
	}
	#
	return(result)
}



# -5-

DistanceOnSphere <- function(coordA, coordB, ellipsoid=TRUE){

    ## INPUT:
    ## coordA: vector giving latitude and longitude of point A: c(lat,lon)
    ## coordB: vector giving latitude and longitude of point B
    ## ellipsoid: if TRUE, use WGS84 reference ellipsoid for the Earth; if false,
    ##            assume simple spherical shape
    ## CONVENTION:
    ## * angle in decimal degree
    ## * latitude=[-90,90]: positive above equator, negative below equator
    ## * longitude=[-180,180]: positive to the East, negative to the West
    ##
    ## author Thomas Münch

    ## clean input variables
    coordA=as.vector(unlist(coordA))
    coordB=as.vector(unlist(coordB))

    ## same point has distance zero
    if (all(coordA==coordB)) return(0)
    
    ## conversion from degree to radian
    deg2rad=pi/180
    phiA=coordA[1]*deg2rad
    lambdaA=coordA[2]*deg2rad
    phiB=coordB[1]*deg2rad
    lambdaB=coordB[2]*deg2rad

    ## (equatorial) radius of the Earth [km]
    rE=6378.137
    
    ## distance for a simpe sphere

    if (ellipsoid==FALSE){
        
        aux=sin(phiA)*sin(phiB)+cos(phiA)*cos(phiB)*cos(lambdaA-lambdaB)
        zeta=acos(aux)

        res=rE*zeta

    ## distance for the WGS84 reference ellipsoid

    } else {

        ## oblateness factor for the Earth
        f=1/298.257223563

        F=(phiA+phiB)/2
        G=(phiA-phiB)/2
        l=(lambdaA-lambdaB)/2

        S=(sin(G))^2*(cos(l))^2+(cos(F))^2*(sin(l))^2
        C=(cos(G))^2*(cos(l))^2+(sin(F))^2*(sin(l))^2

        w=atan(sqrt(S/C))

        D=2*w*rE

        T=sqrt(S*C)/w
        H1=(3*T-1)/(2*C)
        H2=(3*T+1)/(2*S)

        res=D*(1+f*H1*(sin(F))^2*(cos(G))^2-f*H2*(cos(F))^2*(sin(G))^2)

    }

    
    return(res)
    
}



# -6-

distMat <- function(lat, long)
	#Input:
	#	lat,long: list or vector with coordinates of points whose distances between each other are illustrated in a distances matrix 
{
	distance <- array(NA, dim=c(length(lat), length(lat)))
	for(i in 1:(length(lat)-1)){
		for(n in (i+1):length(lat)){
			distance[i,n] <- DistanceOnSphere(c(lat[i], long[i]), c(lat[n], long[n]), ellipsoid=TRUE)
		}
	}
	return(distance)
}



# -7-

distMat.binned.corMat <- function(cor.Mat, dist.Mat, bin, binseq)
	#Input:
	#	cor.Mat: corraltion matrix
	#	dist.Mat: distance matrix belonging to the correlation matrix (distance between correlation pair position)
	#	bin: vector which names the boundaries of the bins
	#	binseq: centered values of the bins
	#Output:
	#	corFkt: correlation function over binned distances (zoo-object)
{
	corFkt <- numeric()
	for(i in 1:(length(bin)-1)){
		corFkt[i] <- mean(cor.Mat[which((dist.Mat >= bin[i]) & (dist.Mat < bin[i+1]))], na.rm=TRUE)
		if(i == (length(bin)-1)){
			corFkt[i] <- mean(cor.Mat[which((dist.Mat >= bin[i]) & (dist.Mat <= bin[i+1]))], na.rm=TRUE)
		}
	}
	return(zoo(corFkt, order.by=binseq))
}



# -8-

distortIrregAge <- function(age, time.uncertainty, points=c(2000,4000,6000,8000,10000))
	#Input:
	#	age: time points which should be distorted
	#	time.uncertainty: time uncertainty which should be applied
	#	points: ages of the original timeseries on which the shift due to time uncertainty should be applied
	#Output:
	#	distorted.age: distorted new time points
	{
		distort.point <- numeric()
		for(i in 1:length(points)){
			appl.time.un <- numeric()
			for(z in 1:(length(points)*1000)){
				appl.time.un[z] <- rnorm(1, mean=0, sd=time.uncertainty)
				if((appl.time.un[z] <= (mean(diff(points))/2)) & (appl.time.un[z] >= (-mean(diff(points))/2))){break}
			}
			distort.point[i] <- points[i] + appl.time.un[z]
		}			
        distorted.points <- c(1, distort.point) 
        distorted.age <- approx(c(1, points), distorted.points, age)$y
        #
        return(distorted.age)
    }



# -9-

gauss.filt.weights <- function(X, fc, dt, timser.length, int.method=c("linear","nearest"))
	#Input:
	#	X: time series (zoo-object)
	#	fc: cut-off frequency
	#	dt: regular inter-observation time step of the interpolation
	#	timser.length: oldest age of the time series
	#	int.method: kind of interpolation (linear, nearest neighbor)
	#Output:
	#	Xsmooth: filtered time series (zoo-object)
{
	if (!is.zoo(X)){stop("time series must be zoo-object")}
	tx <- index(X)
	x <- coredata(X)
	int.out <- seq(from=dt, to=max(tx), by=dt)
	Xsmooth <- rep.int(NA, length(seq(from=dt, to=timser.length, by=dt)))
		#
		if(int.method=="linear"){ X.int <- zoo(approx(tx, x, int.out, method="linear")$y, order.by=int.out) }
		if(int.method=="nearest"){ X.int <- approx.nearest(tx, x, int.out)}	
		fw <- gauss.weights(fc)
		filt.pos <- seq(from=-round(length(fw)/(2*dt))*dt, to=round(length(fw)/(2*dt))*dt, by=dt) + 1 + (length(fw)-1)/2
		if(max(filt.pos) > length(fw)){ filt.pos <- seq(from=-(round(length(fw)/(2*dt))-1)*dt, to=(round(length(fw)/(2*dt))-1)*dt, by=dt) + 1 + (length(fw)-1)/2 }
		#
	Xsmooth[1:length(int.out)] <- filter(X.int, fw[filt.pos]/sum(fw[filt.pos]))
	#
	return(zoo(Xsmooth, order.by=seq(from=dt, to=timser.length, by=dt)))
}



# -10-

gauss.weights <- function(fc)
	#Input:
	#	fc: cut-off frequency
{
	h <- 1/ (fc*pi*sqrt(2/log(2)))
	filt.position <- seq(from=round(-3*h), to=round(3*h), by=1)
	filter <- 1/sqrt(2 * pi * h^2) * (exp(-((filt.position^2))/(2 * h^2)))
	#
	return(filter/sum(filter))
}



# -11-

modifyIrregData <- function(irreg.temperature, irreg.age, time.un, snr, points=c(2000,4000,6000,8000,10000), tsc.in=200, detrending=FALSE, filt.position=10, int.method="linear", startTimser=10, endTimser=11700)
	#Input:	
	# 	irreg.temperature: temperature of the irregular sampled time series (no zoo-object!)
	#	irreg.age: ages of the irregular sampled time series
	#	time.un: time uncertainty which should be applied
	#	snr: signal-to-noise-ratio which should be applied (add white noise)
	#	points: ages of the original time series on which the shift due to time uncertainty should be applied
	#	tsc.in: period used for filtering in gauss.filt.weights 
	#	detrending: FALSE<-using the original data; TRUE<-subtracting a possible linear tendency within the original data
	#	filt.position: regular inter-observation time step of the interpolation
	#	int.method: interpolation method applied while using gauss.filt.weights ("linear" or "nearest")
	#	startTimser: first age of the filtered time series
	#	endTimser: oldest age of the time series
	#Output:
	#	res: gaussfiltered time series with/without noise, time uncertainty, detrending (no zoo-object!)
{
	if(snr==0){
		tmp <- zoo(irreg.temperature, order.by=distortIrregAge(age=irreg.age, time.uncertainty=time.un, points=points))
	}else{
		tmp <- zoo(addNoise(temperature=irreg.temperature, snr=snr), order.by=distortIrregAge(age=irreg.age, time.uncertainty=time.un, points=points))
	}
	tmp.cut <- tmp[which(index(tmp) <= endTimser)]
	#
	res <- coredata(gauss.filt.weights(detr.timser(tmp.cut, detr=detrending), 1/tsc.in, filt.position, endTimser, int.method))	
	#
	return(res)
}



# -12-

##' Simulate a random timeseries with a powerlaw spectrum
##'
##' Method: FFT white noise, rescale, FFT back, the result is scaled to variance 1
##' @title Simulate a random timeseries with a powerlaw spectrum
##' @param beta slope
##' @param N length of timeseries to be generated
##' @return vector containing the timeseries
##' @author Thomas Laepple
##' @export
SimPowerlaw <- function(beta, N)
{
  N2 <- (3^ceiling(log(N, base = 3)))
  df  <- 1 / N2
  f <- seq(from = df, to = 1/2, by = df)
  Filter <- sqrt(1/(f^beta))
  Filter <- c(max(Filter), Filter, rev(Filter))
  x   <- scale(rnorm(N2, 1))
  fx  <- fft(x)
  ffx <- fx * Filter
  result <- Re(fft(ffx, inverse = TRUE))[1:N]
  return(scale(result)[1:N])
}



# -13-

#' Resample a equidistant timeseries (e.g. model result) at the "timepoints" using block averaging
#' The blocks are divided at 1/2 time between the requested output points
#' For the first (and last) timepoint, the interval starting mean(diff(timepoints)) before (ending after)
#' are used.
#' Example usage is to downsample a model timeseries to mimick an integrating proxy
#' (e.g. water isotopes that are measured by melting pieces of ice)
#' @title Subsample (downsample) using block averaging
#' @param ts  ts object or vector containing the equidistant timeseries
#' @param timepoints  vector with the points in time
#' @return values at timepoints
#' @examples
#' input <- ts(SimPowerlaw(0.5,1000))
#' timepoints <- seq(from=50,to=950,by=50)
#' result <- SubsampleTimeseriesBlock(input,timepoints)
#' plot(input,main="Comparison of block avg. vs. simple interpolation",ylab="unitless")
#' points(timepoints,result,pch=19,col="red",lwd=3)
#' points(approx(time(input),c(input),timepoints),col="green",pch=10,lwd=3)
#' legend("bottom",col=c("black","red","green"),lwd=2,bty="n",c("High-resolution timeseries (input)","
#' Block Avg","interpolated values"))
#' @author Thomas Laepple
#' @export
SubsampleTimeseriesBlock <- function(ts, timepoints)
{
    result <- list()
    dt <- mean(diff(timepoints))
    timepoints.bound <- c(head(timepoints, 1) - dt/2, timepoints[-length(timepoints)] +
        diff(timepoints)/2, tail(timepoints + dt/2, 1))
    breaks <- .bincode(c(time(ts)), breaks = timepoints.bound,
        TRUE, TRUE)
    temp <- tapply(c(ts), breaks, mean, na.rm = TRUE)
    return(temp)
}



# -14-

splitMat <- function(cor.Mat, splitlist)
	#Input:
	#	cor.Mat: correlation matrix which should be split in a kind of 'submatrices'
	#	splitlist: list which is used to split the correlation matrix
{
	Split <- list()
	for(d in 1:length(splitlist)){
		Split[[d]] <- array(NA, dim=dim(cor.Mat))
		Split[[d]][splitlist[[d]], splitlist[[d]]] <- cor.Mat[splitlist[[d]], splitlist[[d]]]
	}
	return(Split)
}



# -15-

spatial.cor.analysis.splitvec <- function(data, metadat, detrending=FALSE, tsc.in=200, filt.position=10, int.method="linear", startTimser=10, endTimser=6000, null.hyp=FALSE, quant=c(0.05,0.95), powlaw.beta=1, rep=100, split=FALSE, splitlist=NA, binning=FALSE, bin=NA, binseq=NA, na=FALSE, nalist=naList)
	#Input:
	#	data: list with time series (zoo-objects)
	#	metadat: list with metadata (1-name, 2-latitude, 3-longitude, 4-mean.resolution, 5-proxy.type)
	#	detrending: FALSE <- using the original data; TRUE <- subtracting a possible linear tendency within the original data
	#	tsc.in: period used for filtering in gauss.filt 
	#	filt.position: regular inter-observation time step of the interpolation
	#	int.method: interpolation method applied while using gauss.filt.weights ("linear" or "nearest")
	#	startTimser: first age of the filtered time series
	#	endTimser: oldest age of the time series
	#	null.hyp: TRUE <- calculating quantiles; FALSE <- not calculating quantiles
	#	quant: (only used if null.hyp = TRUE) vector defining quantiles
	#	powlaw.beta: (only used if null.hyp = TRUE) exponent of the power law used to create surrogate time series
	#	rep: (only used if null.hyp = TRUE) number of surrogate time series per original time series
	#	split: TRUE <- splitting the correlation matrix (e.g. for distinguishing proxy types); FALSE <- not splitting the correlation matrix
	#	splitlist: (only used if split = TRUE) list with numbers of time series which is used for splitting the correlation matrix
	#	binning: TRUE <- calculating the spatial correlation function; FALSE <- not calculating the spatial correlation function
	#	bin: (only used if binning = TRUE) sequence used for defining the bins
	#	binseq: (only used if binning = TRUE) sequence containing the midpoint of the bins
	#	na: TRUE <- sets all calculated correlations of the corMatrix NA, if the mean.resolution of the time series is too worse
	#	naList: vector which names the time series whose mean.resolution is worse than the time scale applied
	#
	#Output:
	#	result ($dist,$CorMat,$qCorMat,$CorFkt,$qCorFkt,$QuantCorFkt,$numCorFkt)
{
	Distance <- distMat(metadat[[2]], metadat[[3]])
	#
	Data <- array(NA, dim=c(length(seq(from=startTimser, to=endTimser, by=filt.position)), length(data)))
	for(i in 1:length(data)){
		Data[,i] <- coredata(gauss.filt.weights(detr.timser(data[[i]], detr=detrending), 1/tsc.in, filt.position, endTimser, int.method))	
	}
	#
	print("Working on null.hyp")
	if(null.hyp==FALSE){
		result.CorMat <- corMatrix(Data)
			if(na==TRUE){
				for(u in 1:length(naList)){
					result.CorMat[,naList[u]] <- NA
					result.CorMat[naList[u],] <- NA
				}
			}else{
				result.CorMat <- result.CorMat	
			}
	}else{
		result.CorMat <- corMatrix(Data)
			if(na==TRUE){
				for(u in 1:length(naList)){
					result.CorMat[,naList[u]] <- NA
					result.CorMat[naList[u],] <- NA
				}
			}else{
				result.CorMat <- result.CorMat
			}			
		#
		RandData <- list()
		for(c in 1:rep){
			RandData[[c]] <- array(NA, dim=c(length(seq(from=startTimser, to=endTimser, by=filt.position)), length(data)))
		}
		for(c in 1:rep){
			print(paste("Working on rep=",c,sep=""))
			for(b in 1:length(data)){
				RandData[[c]][,b] <- coredata(gauss.filt.weights(detr.timser(zoo(SubsampleTimeseriesBlock(SimPowerlaw(powlaw.beta, endTimser), index(data[[b]])), order.by=index(data[[b]])), detr=detrending), 1/tsc.in, filt.position, endTimser, int.method))
			}
		}
		result.qCorMat <- list()
		for(c in 1:rep){
			result.qCorMat[[c]] <- corMatrix(RandData[[c]])
				if(na==TRUE){
					for(u in 1:length(naList)){
						result.qCorMat[[c]][,naList[u]] <- NA
						result.qCorMat[[c]][naList[u],] <- NA
					}
				}else{
					result.qCorMat <- result.qCorMat			
				}
		}
	}
	#
	print("Working on split")
	if(split==FALSE){
		if(null.hyp==TRUE){
			CorMat <- result.CorMat
			qCorMat <- result.qCorMat
		}else{
			CorMat <- result.CorMat
		}
	}else{
		if(null.hyp==TRUE){
			CorMat <- splitMat(result.CorMat, splitlist)
			CorMat[[(length(splitlist)+1)]] <- result.CorMat
			qCorMat <- list()
			for(c in 1:rep){
				qCorMat[[c]] <- splitMat(result.qCorMat[[c]], splitlist)
				qCorMat[[c]][[(length(splitlist)+1)]] <- result.qCorMat[[c]]
			}
		}else{
			CorMat <- splitMat(result.CorMat, splitlist)
			CorMat[[(length(splitlist)+1)]] <- result.CorMat
		}		
	}
	#
	print("Working on binning")
	if(binning==FALSE){
		if(null.hyp==TRUE){
			CorMat <- CorMat
			qCorMat <- qCorMat
		}else{
			CorMat <- CorMat
		}
	}else{
		if(split==FALSE){
			if(null.hyp==TRUE){
				CorFkt <- coredata(distMat.binned.corMat(CorMat, Distance, bin, binseq))
				qCorFkt <- array(NA, dim=c(length(binseq), rep))
				for(c in 1:rep){
					qCorFkt[,c] <- coredata(distMat.binned.corMat(qCorMat[[c]], Distance, bin, binseq))
				}
				QuantCorFkt <- array(NA, dim=c(length(quant), length(binseq)))
				for(i in 1:length(binseq)){
					QuantCorFkt[,i] <- quantile(qCorFkt[i,], quant, na.rm=TRUE)
				}	
			}else{
				CorFkt <- coredata(distMat.binned.corMat(CorMat, Distance, bin, binseq))
			}
		}else{
			if(null.hyp==TRUE){
				CorFkt <- list()
				for(i in 1:length(CorMat)){
					CorFkt[[i]] <- coredata(distMat.binned.corMat(CorMat[[i]], Distance, bin, binseq))
				}
				qCorFkt <- list()
				for(n in 1:(length(splitlist)+1)){
					qCorFkt[[n]] <- array(NA, dim=c(length(binseq), rep))
					for(c in 1:rep){
						qCorFkt[[n]][,c] <- coredata(distMat.binned.corMat(qCorMat[[c]][[n]], Distance, bin, binseq))
					}
				}
				QuantCorFkt <- list()
				for(n in 1:(length(splitlist)+1)){
					QuantCorFkt[[n]] <- array(NA, dim=c(length(quant), length(binseq)))
					for(i in 1:length(binseq)){
						QuantCorFkt[[n]][,i] <- quantile(qCorFkt[[n]][i,], quant, na.rm=TRUE)
					}
				}
			}else{
				CorFkt <- list()
				for(i in 1:length(CorMat)){
					CorFkt[[i]] <- coredata(distMat.binned.corMat(CorMat[[i]], Distance, bin, binseq))
				}
			}
		}
	}
	#
	print("Working on numCorFkt")
	if(binning==TRUE){
		if(split==TRUE){
			numCorFkt <- list()
			for(n in 1:(length(splitlist)+1)){
				numCorFkt[[n]] <- numeric()
				for(i in 1:(length(bin)-1)){
					numCorFkt[[n]][i] <- length(which(!is.na(match(which((Distance >= bin[i]) & (Distance < bin[i+1])), which(!is.na(CorMat[[n]]))))))
				}
			}
		}else{
			numCorFkt <- numeric()
			for(i in 1:(length(bin)-1)){
				numCorFkt[i] <- length(which(!is.na(match(which((Distance >= bin[i]) & (Distance < bin[i+1])), which(!is.na(CorMat))))))
			}
		}
	}else{
		numCorFkt <- NA
	}
	#
	result <- list("dist"=list(),"CorMat"=list(),"qCorMat"=list(),"CorFkt"=list(),"qCorFkt"=list(),"QuantCorFkt"=list(),"numCorFkt"=list())
	result$dist <- Distance
	result$CorMat <- CorMat
	if(null.hyp==TRUE){
		result$qCorMat <- qCorMat
	}else{
		result$qCorMat <- NA
	}
	if(binning==TRUE){
		result$CorFkt <- CorFkt
		result$numCorFkt <- numCorFkt
		if(null.hyp==TRUE){
			result$qCorFkt <- qCorFkt
			result$QuantCorFkt <- QuantCorFkt
		}else{
			result$qCorFkt <- NA
			result$QuantCorFkt <- NA
		}
	}else{
		result$CorFkt <- NA
		result$qCorFkt <- NA
		result$QuantCorFkt <- NA
		result$numCorFkt <- NA
	}
	#
	return(result)
}				


