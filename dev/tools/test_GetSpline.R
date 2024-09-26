library(sp)
library(raster)
library(terra)
library(sf)

library(rjson)
library(geojsonR)

########################################
numSite <- 1; cc <- 50; bb <- 1


########################################
params <- fromJSON(file='/usr3/graduate/mkmoon/GitHub/biomass/dev/PBM_Parameters.json')
source(params$setup$rFunctions)

## Get site name and image directory
geojsonDir <- params$setup$geojsonDir

strSite <- list.dirs(params$setup$outDir,full.names=F,recursive=F)[numSite]
print(strSite)

ckDir <- paste0(params$setup$outDir,strSite,'/chunk')
print(ckDir)


########################################
imgBase <- raster(paste0(params$setup$outDir,strSite,'/base_image.tif'))
numPix <- length(imgBase)
imgNum <- setValues(imgBase,1:numPix)
numChunks <- params$setup$numChunks
chunk <- numPix%/%numChunks

# ##
# ptShp <- shapefile(paste0(params$setup$workDir,'shp/mg_pts_',numSite,'.shp'))
# ptShp <- spTransform(ptShp,crs(imgBase))
# pixNums <- extract(imgNum,ptShp)

pixNums <- 400000

##
ppp <- 1
pixNum  <- pixNums[ppp]
  
ckNum <- sprintf('%03d',(pixNum%/%chunk+1))
file <- list.files(path=ckDir,pattern=glob2rx(paste0('*',ckNum,'.rda')),full.names=T)

load(file)
  

##
blue  <- band1[pixNum%%chunk,]
green <- band2[pixNum%%chunk,]
red   <- band3[pixNum%%chunk,]
nir   <- band4[pixNum%%chunk,]
phenYrs <- params$setup$phenStartYr:params$setup$phenEndYr


########################################
GetSpline <- function(blue, green, red, nir, dates, phenYrs, params, bb, yy){
  
  # Despike, calculate dormant value, fill negative VI values with dormant value
  log <- try({    
    
    pheno_pars <- params$phenology_parameters
    qa_pars    <- params$qa_parameters
    
    b2 <- blue/10000; b3 <- green/10000; b4 <- red/10000; b5 <- nir/10000
    
    i1   <- (b5 - b4) / (b5 + b4) # NDVI
    i2   <- (b5 - b3) / (b5 + b3) # GNDVI
    i3   <- b5 * i1               # NIRv
    vi   <- 2.5*(b5 - b4) / (b5 + 2.4*b4 + 1) # EVI2
    
    
    # Spikes check, and remove
    spikes     <- CheckSpike_MultiBand(b2, b4, vi, dates, pheno_pars)
    b2[spikes] <- NA; b3[spikes] <- NA; b4[spikes] <- NA; b5[spikes] <- NA
    i1[spikes] <- NA; i2[spikes] <- NA; i3[spikes] <- NA; vi[spikes] <- NA
    
    # Replace negative VIs with dormant value
    dormIms <- dates >= pheno_pars$dormStart & dates <= pheno_pars$dormEnd
    vi_dorm <- quantile(vi[dormIms & vi>0],probs=pheno_pars$dormantQuantile,na.rm=T)   # Calc vi dormant value using non-negative VIs
    
    
    #now calculate dormancy values and fill individual bands
    dormObs <- dormIms & vi < vi_dorm    #Defining dormant observations for bands as median on dates when vi < vi_dorm
    b2_dorm <- median(b2[dormObs], na.rm=T); b2[dormObs] <- b2_dorm
    b3_dorm <- median(b3[dormObs], na.rm=T); b3[dormObs] <- b3_dorm
    b4_dorm <- median(b4[dormObs], na.rm=T); b4[dormObs] <- b4_dorm
    b5_dorm <- median(b5[dormObs], na.rm=T); b5[dormObs] <- b5_dorm
    i1_dorm <- median(i1[dormObs], na.rm=T); i1[dormObs] <- i1_dorm
    i2_dorm <- median(i1[dormObs], na.rm=T); i2[dormObs] <- i2_dorm
    i3_dorm <- median(i1[dormObs], na.rm=T); i3[dormObs] <- i3_dorm
    vi[vi < vi_dorm] <- vi_dorm
    
    
    #
    splineStart <- as.Date(as.Date(paste0(phenYrs,'-01-01')) - pheno_pars$splineBuffer) 
    numDaysFit  <- 365 + (pheno_pars$splineBuffer * 2)    
    splineEnd   <- splineStart+(numDaysFit-1)
    all_dates   <- seq(min(splineStart), max(splineEnd), by="day")
    
    numYrs <- length(phenYrs)
    vecLength <- numDaysFit*numYrs
    
    #
    smoothMat <- matrix(NA, numDaysFit, numYrs)
    b2Mat <- smoothMat; b3Mat <- smoothMat; b4Mat <- smoothMat; b5Mat <- smoothMat
    i1Mat <- smoothMat; i2Mat <- smoothMat; i3Mat <- smoothMat
    
  },silent=TRUE)
  #If there is an error despiking or other initial steps, return NAs
  if(inherits(log, "try-error")){return(matrix(NA,58*length(phenYrs)))}   
  
  outMat <- matrix(NA,1,365)
  log <- try({    
    
    # Fit spline
    if(bb==1){
      smoothed_ts <- Smooth_VI(vi, dates, all_dates, pheno_pars, vi_dorm)
    }else if(bb==2){
      smoothed_ts <- Smooth_Bands(i1, dates, all_dates, pheno_pars)  
    }else if(bb==3){
      smoothed_ts <- Smooth_Bands(i2, dates, all_dates, pheno_pars)  
    }else if(bb==4){
      smoothed_ts <- Smooth_Bands(i3, dates, all_dates, pheno_pars)  
    }else if(bb==5){
      smoothed_ts <- Smooth_Bands(b2, dates, all_dates, pheno_pars)  
    }else if(bb==6){
      smoothed_ts <- Smooth_Bands(b3, dates, all_dates, pheno_pars)  
    }else if(bb==7){
      smoothed_ts <- Smooth_Bands(b4, dates, all_dates, pheno_pars)  
    }else if(bb==8){
      smoothed_ts <- Smooth_Bands(b5, dates, all_dates, pheno_pars)  
    }
    
  },silent=TRUE) #End of the try block
  
  if(inherits(log, "try-error")){
    outAll <- matrix(NA,1,365)
  }else{
    outAll <- smoothed_ts[all_dates >= as.Date(paste0(yy, "-01-01")) & all_dates < as.Date(paste0(yy+1, "-01-01"))]
    outAll <- outAll[1:365]
  }
  
  return(outAll)
  
}



