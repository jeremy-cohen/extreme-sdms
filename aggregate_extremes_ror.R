# Richness and range size rarity from distribution predictions
# RELATIVE OCCURRENCE RATE VERSION

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
    files = str_subset(files, "ROR")
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
    writeRaster(rich, paste0("agg/ror_",season,"_",temperature,"_rich.tif"),
                overwrite=T)
    # plot, save
    jpeg(paste0("agg/ror_",season,"_",temperature,"_rich.jpeg"),
         width=1000, height=1000, units="px")
    plot(rich, col=colors)
    plot(countries_outline, col="transparent", border="black",
         usePolypath=F, add=T)
    dev.off()
  }
}

# Example species (fig. s8)
sp_name = c("Ixoreus_naevius", "Sayornis_phoebe", "Icterus_cucullatus")
# define colors
rorcolors = c(brewer.pal(9, 'YlOrRd'), "#800026")
# loop species
for (sp in c(1:3)){ # 1:3
  if(sp==2){season="winter"}else{season="summer"}
jpeg(paste0("agg/ror_",sp_name[sp],"_temp_comp_ror.jpeg"),
     width=600, height=250, units="px")
par(mfrow=c(1,3), mar=c(1.1,1.1,1.1,0))
# loop models
for (temperature in c("means","variability","extremes")){
  # get raster, crop, mask, project
r <- raster(paste0(sp_name[sp],"/predictions/ROR/PO_predictions_",sp_name[sp],
                   "_",season,"_",temperature,"_00.tif")) %>%
  crop(countries) %>%
  mask(mask = countries) %>%
  projectRaster(crs=plotcrs) %>%
  trim()
# plotting
plot(r, col=rorcolors, breaks=c(0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1),
     xaxt='n', yaxt='n', legend=F, box=F, axes=F)
plot(countries_proj, col="transparent", border="black", usePolypath=F, add=T)
}
dev.off()
}


# Delta biodiversity maps (fig. s11)
bp <- brewer.pal(9, "BrBG")
colors = rev(c(rep(bp[1],1), bp[2:4],  bp[6:8], rep(bp[9],1)))
setwd(paste0("/gpfs/gibbs/pi/jetz/from_loomis/data/results_30x30/birds/",ver,"_30x30/"))
  for(season in c("summer","winter")){
    # load richness raster
    r1 <- raster(paste0("agg/ror_",season,"_means_rich.tif"))
    for (temperature in c("variability","extremes")){
      # load raster, project, calculate delta
      r <- raster(paste0("agg/ror_",season,"_",temperature,"_rich.tif"))
      r1 <- projectRaster(r1, r)
      rd <- r-r1
      writeRaster(rd, paste0("agg/ror_",season,"_",temperature,"_rich_delta.tif"),
                  overwrite=T)
      # aggregate for faster plotting
      rd <- aggregate(rd, 4)
      jpeg(paste0("agg/ror_",season,"_",temperature,"_rich_delta.jpeg"),
           width=500, height=400, units="px")
      plot(rd, breaks=seq(-60,60,15), col=colors, legend=F)
      plot(countries_proj, col="transparent", border="black", usePolypath=F, add=T)
      dev.off()
  }
}


# Better visualizing delta richness - fig. s11
# define colors
bp <- rev(brewer.pal(9, "Spectral"))
colors = c("#0F70DE" ,bp[1:2], bp[3], bp[5:9], "#B30707")

bp2 <- brewer.pal(9, "YlOrBr")
bp2 <- brewer.pal(9, "Purples")
colors2 = rev(c("lightgoldenrod3","lightgoldenrod2","lightgoldenrod1",
                "white","white",bp2[3:9]))
{
jpeg("agg/ror_fig1.jpeg", width=750, height=400, units="px")
par(mfrow=c(2,3), mar=c(1.1,.6,1.1,0))
# loop seasons and plot rasters
for(season in c("summer","winter")){
      r <- raster(paste0("agg/ror_",season,"_means_rich.tif"))
      r[r == 0] <- NA
      r = trim(r)
      plot(r, col=colors, zlim=c(0,184), legend=F,
           xaxt='n', yaxt='n', box=F, axes=F)
      plot(countries_outline, col="transparent", border="black",
           usePolypath=F, add=T, box=F)
      # loop variability/extremes predictions
      for (temperature in c("variability","extremes")){
       d <- raster(paste0("agg/ror_",season,"_",temperature,"_rich_delta.tif"))
       d[d == 0] <- NA
       d = trim(d)
       plot(d, col=colors2, breaks=seq(-40,20,5), legend=T,
            xaxt='n', yaxt='n', box=F, axes=F)
       plot(countries_outline, col="transparent", border="black",
            usePolypath=F, add=T, box=F)
    }
  }
  dev.off()
}

# Absolute richness predictions - figure s12
jpeg("agg/ror_figs1.jpeg", width=750, height=400, units="px")
par(mfrow=c(2,3), mar=c(1.1,.6,1.1,0))
for(season in c("summer","winter")){
  for (temperature in c("means","variability","extremes")){
    r <- raster(paste0("agg/ror_",season,"_",temperature,"_rich.tif"))
    r[r == 0] <- NA
    r = trim(r)
    plot(r, col=colors, zlim=c(0,184), legend=F,
         xaxt='n', yaxt='n', box=F, axes=F)
    plot(countries_outline, col="transparent", border="black",
         usePolypath=F, add=T, box=F)
  }
}
dev.off()
