# Richness and range size rarity from distribution predictions
# BINARY PRESENCE VERSION

library(tidyverse)
library(raster)
library(rgdal)
library(RColorBrewer)
library(stringr)
library(parallel)
library(doParallel)
library(foreach)
library(rasterVis)
library(sf)
ver = "VE8"

# Tables with species, model, area
cols = c("sciname","temperature","season","range_area")
splist = read_csv(paste0("/gpfs/gibbs/pi/jetz/from_loomis/data/results_30x30/birds/sdm_birds_",ver,"_0.csv"))[,cols]

# Plotting
bp <- rev(brewer.pal(9, "Spectral"))
colors = c("#0F70DE" ,bp[1:2], bp[3], bp[5:9], "#B30707")

# load country polygons
countries = st_read("/vast/palmer/home.mccleary/jc3893/30x30/shapefiles/", "countries") %>%
  st_transform(crs="+proj=cea +lat_ts=30 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs ")
canada = st_read("/vast/palmer/home.mccleary/jc3893/30x30/shapefiles/", "Canada") %>%
  st_transform(crs="+proj=cea +lat_ts=30 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs ")
# define coords system
plotcrs = "+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83"

# version with just countries
countries_outline <- countries[c(9,38),] %>%
  st_transform(crs=plotcrs, allow_ballpark=T) %>%
  st_crop(extent(-5771203, 2991315, -1690735, 4518085))
# select countries
countries = countries[c(9,44:54,56:80,82:87,89:97),] 

# version with all states/provs
countries_proj <- countries %>%
  st_union(canada) %>%
     st_transform(crs=plotcrs)

# wd
setwd(paste0("/gpfs/gibbs/pi/jetz/from_loomis/data/results_30x30/birds/",ver,"_30x30/"))
dir.create("agg/")

# parallel computing
n.cores <- 2
my.cluster <- parallel::makeCluster(
  n.cores,
  type = "PSOCK"
)
doParallel::registerDoParallel(cl = my.cluster)

# Richness estimation
# Season loop 
foreach(season = c("winter","summer")) %dopar% {
  library(tidyverse)
  library(raster)
  library(rgdal)
  library(RColorBrewer)
  library(stringr)
  library(parallel)
  library(doParallel)
  library(foreach)
  # Across temperature models
  for (temperature in c("means","variability","extremes")){
    # Subset
    splist_sub = splist[splist$temperature==temperature &
                          splist$season==season,]
    # Files
    files = list.files(pattern=paste0("_",season,"_",temperature,"_00.tif"),
                       recursive=T, full.names=T)
    files = str_subset(files, "PA")
    # remove sp not included
    files = files[grepl(paste0(splist_sub[,1][[1]],collapse = "|"), files)]
    # Richness- add species prediction rasters
    for (file in 1:length(files)){
      try({
      tif = files[file]
      sciname = gsub("./","",dirname(dirname(dirname(tif))))
      print(sciname)
        r_tif = raster(tif)
        if(file==1){rich = r_tif
        }else{
          rich = sum(rich, r_tif, na.rm=T)
        }
      })
    }
    # raster processing (crop and mask to US/Canada, project to new coords system)
      rich = rich %>%
        crop(countries) %>%
        mask(mask = countries) %>%
        projectRaster(crs=plotcrs) %>%
        trim()
      zrange=c(0,max(getValues(rich), na.rm=T))
      # plot, save
    writeRaster(rich, paste0("agg/",season,"_",temperature,"_rich.tif"),
                overwrite=T)
    jpeg(paste0("agg/",season,"_",temperature,"_rich.jpeg"),
         width=1000, height=1000, units="px")
    plot(rich, col=colors)
    plot(countries_outline, col="transparent", border="black",
         usePolypath=F, add=T)
    dev.off()
  }
}


# Example species (fig. 2)
sp_name = c("Ixoreus_naevius", "Sayornis_phoebe", "Icterus_cucullatus")
col = c("royalblue2", "yellow2", "orange")
for (sp in c(2)){ # 1:3
  if(sp==2){season="winter"}else{season="summer"}
jpeg(paste0("agg/",sp_name[sp],"_temp_comp.jpeg"),
     width=600, height=250, units="px")
par(mfrow=c(1,3), mar=c(1.1,1.1,1.1,0))
# loop models
for (temperature in c("means","variability","extremes")){
  # get raster, crop, mask, project
r <- raster(paste0(sp_name[sp],"/predictions/PA/Thresholded_predictions_",sp_name[sp],
                   "_",season,"_",temperature,"_00.tif")) %>%
  crop(countries) %>%
  mask(mask = countries) %>%
  projectRaster(crs=plotcrs)
# if temp is means, trim, otherwise disaggregate and process
if(temperature=="means"){r1 = r %>% trim()
r1test <- reclassify(r1, cbind(-Inf, 1.5, 2))
r1test[is.na(r1test[])] <- 0
plot(r1, col=col[sp], legend=F, xaxt='n', yaxt='n', box=F, axes=F)
}else{
rda <- r %>%
  projectRaster(r1)
# create levels to define range loss, no difference, or gain between layers
rda[is.na(rda[])] <- 0
rda <- reclassify(rda, cbind(.5, 1.5, 1))
r <- r1test-rda %>%
  trim()
# plotting
cpal <- c('gray80','transparent', col[sp], 'black')
r2 <- ratify(r)
# # for legend
# levelplot(r2,col.regions=cpal,att='ID',
#           scales=list(x=list(at=NULL), y=list(at=NULL)),
#           par.settings = list(axis.line = list(col = "transparent")), xaxt='n', yaxt='n', box=F, axes=F)
plot(r, col=cpal, legend=F, xaxt='n', yaxt='n', box=F, axes=F)}
plot(countries_proj, col="transparent", border="black",
     usePolypath=F, add=T)
}
dev.off()
}


# Delta biodiversity maps (fig. 5)
bp <- brewer.pal(9, "BrBG")
colors = rev(c(rep(bp[1],1), bp[2:4],  bp[6:8], rep(bp[9],1)))

setwd(paste0("/gpfs/gibbs/pi/jetz/from_loomis/data/results_30x30/birds/",ver,"_30x30/"))
  for(season in c("summer","winter")){
    # load richness raster
    r1 <- raster(paste0("agg/",season,"_means_rich.tif"))
    for (temperature in c("variability","extremes")){
      # load raster, project, calculate delta
      r <- raster(paste0("agg/",season,"_",temperature,"_rich.tif"))
      r1 <- projectRaster(r1, r)
      rd <- r-r1
      writeRaster(rd, paste0("agg/",season,"_",temperature,"_rich_delta.tif"),
                  overwrite=T)
      # aggregate for faster plotting
      rd <- aggregate(rd, 4)
      jpeg(paste0("agg/",season,"_",temperature,"_rich_delta.jpeg"),
           width=500, height=400, units="px")
      plot(rd, breaks=seq(-60,60,15), col=colors)
      plot(countries_proj, col="transparent", border="black", usePolypath=F, add=T)
      dev.off()
  }
}


# better visualizing delta richness - fig. 5
# define colors
bp <- rev(brewer.pal(9, "Spectral"))
colors = c("#0F70DE" ,bp[1:2], bp[3], bp[5:9], "#B30707")

bp2 <- brewer.pal(9, "YlOrBr")
bp2 <- brewer.pal(9, "Purples")
colors2 = rev(c(rep("lightgoldenrod3",3),"lightgoldenrod2","lightgoldenrod1",
                rep("white",2),bp2[2:8],rep(bp2[9],2)))
{
jpeg("agg/fig1.jpeg", width=750, height=400, units="px")
par(mfrow=c(2,3), mar=c(1.1,.6,1.1,0))
# loop seasons and plot rasters
for(season in c("summer","winter")){
      r <- raster(paste0("agg/",season,"_means_rich.tif"))
      r[r == 0] <- NA
      r = trim(r)
      plot(r, col=colors, zlim=c(0,184), legend=F,
           xaxt='n', yaxt='n', box=F, axes=F)
      plot(countries_outline, col="transparent", border="black",
           usePolypath=F, add=T, box=F)
      # loop variability/extremes predictions
      for (temperature in c("variability","extremes")){
       d <- raster(paste0("agg/",season,"_",temperature,"_rich_delta.tif"))
       d[d == 0] <- NA
       d = trim(d)
       plot(d, col=colors2, breaks=seq(-50,30,5), legend=F, # plot one with legend on for 
            xaxt='n', yaxt='n', box=F, axes=F)
       plot(countries_outline, col="transparent", border="black",
            usePolypath=F, add=T, box=F)
    }
  }
  dev.off()
}

# Absolute richness predictions - figure s10
jpeg("agg/figs1.jpeg", width=750, height=400, units="px")
par(mfrow=c(2,3), mar=c(1.1,.6,1.1,0))
for(season in c("summer","winter")){
  for (temperature in c("means","variability","extremes")){
    r <- raster(paste0("agg/",season,"_",temperature,"_rich.tif"))
    r[r == 0] <- NA
    r = trim(r)
    plot(r, col=colors, zlim=c(0,184), legend=F,
         xaxt='n', yaxt='n', box=F, axes=F)
    plot(countries_outline, col="transparent", border="black",
         usePolypath=F, add=T, box=F)
  }
}
dev.off()

