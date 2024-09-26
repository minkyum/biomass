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

numSite <- as.numeric(substr(args[3],1, 3)) # site number
cc      <- as.numeric(substr(args[3],4, 6)) # chunk number
bb      <- as.numeric(substr(args[3],7, 8)) # vi or bands
yy      <- as.numeric(substr(args[3],9,12)) # year
# numSite <- 1; cc <- 50; bb <- 1; yy <- 2022



###############################
params <- fromJSON(file='/usr3/graduate/mkmoon/GitHub/biomass/dev/PBM_Parameters.json')
source(params$setup$rFunctions)



########################################
## Get site name, image directory and coordinate
strSite <- list.dirs(params$setup$outDir,full.names=F,recursive=F)[numSite]
print(strSite)

ckDir <- paste0(params$setup$outDir,strSite,'/chunk')
print(ckDir)

ckNum <- sprintf('%03d',cc)
file <- list.files(path=ckDir,pattern=glob2rx(paste0('*',ckNum,'.rda')),full.names=T)

load(file)


##########################################
numPix <- dim(band1)[1]
phenYrs <- params$setup$phenStartYr:params$setup$phenEndYr

s_mat <- matrix(NA,numPix,365)

for (i in 1:numPix){

  s_mat[i,] <- GetSpline(band1[i,],band2[i,],band3[i,],band4[i,],dates,phenYrs,params,bb,yy)
  
  if(i%%10000==0) print(i)
}


# Save
ckDir <- paste0(params$setup$outDir,strSite,'/chunk_spl')
if (!dir.exists(ckDir)) {dir.create(ckDir)}
ckDir <- paste0(ckDir,'/',sprintf('%02d',bb))
if (!dir.exists(ckDir)) {dir.create(ckDir)}

save(s_mat,file=paste0(ckDir,'/chunk_spl_',yy,'_',ckNum,'.rda'))




