library(sp)
library(raster)
library(terra)
library(sf)

library(rjson)
library(geojsonR)

library(doMC)
library(doParallel)

###############################
args <- commandArgs()
print(args)

numSite <- as.numeric(substr(args[3],1,3))
cc      <- as.numeric(substr(args[3],4,6))
# numSite <- 1; cc <- 100


###############################
params <- fromJSON(file='/usr3/graduate/mkmoon/GitHub/biomass/dev/PBM_Parameters.json')
source(params$setup$rFunctions)



########################################
## Get site name and image directory
strSite <- list.dirs(params$setup$outDir,full.names=F,recursive=F)[numSite]
print(strSite)

imgDir <- paste0(params$setup$outDir,strSite,'/mosaic')
print(imgDir)



###############################
dfiles <- list.files(path=imgDir,pattern=glob2rx('*mosaic.tif'))
files  <- list.files(path=imgDir,pattern=glob2rx('*mosaic.tif'),full.names=T)

# Dates
yy <- substr(dfiles,3,4)
mm <- substr(dfiles,5,6)
dd <- substr(dfiles,7,8)
dates_all <- as.Date(paste(mm,'/',dd,'/',yy,sep=''),'%m/%d/%y')
dates <- unique(dates_all)

print(length(dates))

# Divide into chunks
imgBase <- raster(paste0(params$setup$outDir,strSite,'/base_image.tif'))

numCk <- params$setup$numChunks
chunk <- length(imgBase)%/%numCk

ckDir <- paste0(params$setup$outDir,strSite,'/chunk')
if (!dir.exists(ckDir)) {dir.create(ckDir)}


###############################
if(cc==numCk){
  chunks <- c((chunk*(cc-1)+1):length(imgBase))
}else{
  chunks <- c((chunk*(cc-1)+1):(chunk*cc))
}

band1 <- matrix(NA,length(chunks),length(dates))
band2 <- matrix(NA,length(chunks),length(dates))
band3 <- matrix(NA,length(chunks),length(dates))
band4 <- matrix(NA,length(chunks),length(dates))
for(i in 1:length(dates)){
  band1[,i] <- values(raster(files[i],1))[chunks]
  band2[,i] <- values(raster(files[i],2))[chunks]
  band3[,i] <- values(raster(files[i],3))[chunks]
  band4[,i] <- values(raster(files[i],4))[chunks]
  
  if(i%%20==0) print(i)
}
  

# Save
ckNum <- sprintf('%03d',cc)
save(band1,band2,band3,band4,dates,
     file=paste0(ckDir,'/chunk_',ckNum,'.rda'))



