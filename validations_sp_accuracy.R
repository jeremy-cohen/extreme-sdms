# Validating models against test site lists
library(lubridate)
library(sf)
library(gridExtra)
library(tidyverse)
library(raster)
library(RColorBrewer)
library(MetBrewer)
library(ggpubr)
library(rnaturalearth)
library(rnaturalearthdata)
library(Metrics)
library(data.table)
ver = "VE8"
setwd("/vast/palmer/home.mccleary/jc3893/30x30/")

# country polygons and grid to select sites
countries = st_read("shapefiles/", "countries") %>%
  st_transform(crs="+proj=cea +lat_ts=30 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs ")
bbox <- st_read("shapefiles/", "30x30BBox_NAmerica") %>%
  st_transform(crs="+proj=cea +lat_ts=30 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs ")
grid200 <- raster("grids/Global_reference_raster_200km_CEA.tif") %>%
  crop(bbox)
grid200 <- setValues(grid200, 1:ncell(grid200))

setwd("/vast/palmer/home.mccleary/jc3893/30x30/extremes/")

# observed richness and species at the top birded site per 200km cell
for (season in c("summer", "winter")){
  ebird <- read_csv(paste0("/gpfs/gibbs/pi/jetz/data/temporary_datasets/ebird_extremes_jeremy/ebird_2024_",
                           season,"_extracted.csv"))
  for (cell in 1:ncell(grid200)){
    if(cell==1 & season=="summer"){rm(out,out.0,out.02,out.05,out.1)}
    print(cell)
    rcell <- rasterFromCells(grid200, cell)
    ecell <- ebird[ebird$x > xmin(rcell) &
                     ebird$x < xmax(rcell) &
                     ebird$y > ymin(rcell) &
                     ebird$y < ymax(rcell),]
    if (nrow(ecell) >= 100){
      ecell$locid = paste(ecell$x, ecell$y)
      top <- modal(as.factor(ecell$locid))
      esite <- ecell[ecell$locid==top,]
      if (nrow(esite) >= 100){
        sp_pres <- ifelse(colMeans(esite[,c(7:695)]) > .02, 1, 0)
        siterich <- sum(as.numeric(sp_pres))
        total <- unlist(c(cell, season, esite[1,"x"], esite[1,"y"],
                            esite$locid[[1]], nrow(esite), siterich, sp_pres))
        assign(paste0("total.02"), total)
         }
        out <- rbind(out.02,total.0.02)}else{out.02 <- total.0.02}
      
    }
  }
colnames(out) <- c("cell","season","x","y","locality_id","visits","rich",
                   gsub(" ","_",names(sp_pres)))
rownames(out) <- NULL
write_csv(out, paste0("site_richness_validation_withsp_out.02.csv"))


# Add estimated richness across model types from models
out <- read_csv("site_richness_validation_withsp_out.02.csv") 
coordinates(out) <- c("x","y")
setwd(paste0("/gpfs/gibbs/pi/jetz/from_loomis/data/results_30x30/birds/",
             ver,"_30x30/"))
# loop seasons
for (season in c("summer", "winter")){
  out.s = out[out$season==season,]
  # stack projected richness estimates across 3 model types
  stack <- stack()
  for (temperature in c("means","variability","extremes")){
    r <- raster(paste0("agg/",season,"_",temperature,"_rich.tif")) %>%
      projectRaster(crs="+proj=cea +lat_ts=30 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs ")
    if(temperature!="means"){r <- projectRaster(r, stack)}
    stack <- stack(stack, r)
  }
  ests <- raster::extract(stack, out.s)
  colnames(ests) <- paste0("est.",c("means","variability","extremes"))
  assign(paste0("ests.",season), ests)
}
# attach richness estimates and delta richness between estimates/observed values
setwd("/vast/palmer/home.mccleary/jc3893/30x30/extremes/")
val = read_csv(paste0("site_richness_validation_withsp_out.02.csv"))
val <- cbind(data.frame(val), rbind(ests.summer, ests.winter))
val$optional <- NULL
val <- val[complete.cases(val),]
val = val[val$est.means>0,]
val[,697:699] <- round(val[,697:699],0)
val$delta.means <- val$est.means - val$rich
val$delta.variability <- val$est.variability - val$rich
val$delta.extremes <- val$est.extremes - val$rich
write_csv(val, paste0("site_richness_validation_withsp_out.02_",ver,"_comp.csv"))


#### get species presence/absence estimates for all model types
# Tables with species, model, area
cols = c("sciname","temperature","season","range_area")
# species list
splist = read_csv(paste0("/gpfs/gibbs/pi/jetz/from_loomis/data/results_30x30/birds/sdm_birds_",
                         ver,"_0.csv"))[,cols]
# working directory
setwd(paste0("/gpfs/gibbs/pi/jetz/from_loomis/data/results_30x30/birds/",
             ver,"_30x30/"))
# What species are in each cell? 
# loop seasons
for (season in c("winter","summer")){
  print(season)
  # loop model types
  for (temperature in c("means","variability","extremes")){ 
    print(temperature)
    # Subset
    splist_sub = splist[splist$temperature==temperature &
                          (splist$season==season |
                             splist$season=="resident"),]
    # Files
    files = c(list.files(pattern=paste0("_",season,"_",temperature,"_00.tif"),
                         recursive=T, full.names=T),
              list.files(pattern=paste0("_resident_",temperature,"_00.tif"),
                         recursive=T, full.names=T))
    files = str_subset(files, "PA")
    # remove sp not included
    files = files[grepl(paste0(splist_sub[,1][[1]],collapse = "|"), files)]
    # Get sp per cell
    stack <- stack()
    for (file in 1:length(files)){
      r <- raster(files[file])
      stack <- stack(stack, r)
    }
    # bring in site-level values
    setwd("/vast/palmer/home.mccleary/jc3893/30x30/extremes/")
    vals <- read_csv(paste0("site_richness_validation_withsp_out.0_",ver,"_comp.csv")) 
    vals <- vals[vals$season==season,]
    coordinates(vals) <- c("x","y")
    
    # use full raster to isolate cells with the sites
    setwd(paste0("/gpfs/gibbs/pi/jetz/from_loomis/data/results_30x30/birds/",
    ver,"_30x30/"))
    r <- raster(paste0("agg/",season,"_",temperature,"_rich.tif")) %>%
      projectRaster(crs="+proj=cea +lat_ts=30 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs ")
    r <- setValues(r, 1:ncell(r))
    ext <- raster::extract(r, vals)
    # get points for only the relevant cells
    r <- setValues(r, NA)
    values(r)[ext] =  1
    pts <- data.frame(rasterToPoints(r))
    coordinates(pts) <- c("x","y")
    # get species per point
    ex <- raster::extract(stack, pts)
    colnames(ex) <- sub("Thresholded_predictions_", "", colnames(ex))
    colnames(ex) <- sub("_resident.*", "", colnames(ex))
    colnames(ex) <- sub(paste0("_",season,".*"), "", colnames(ex))
    pts <- data.frame(pts)
    # save
    write_csv(cbind(pts[,c("x","y")], data.frame(ex)),
              paste0("/vast/palmer/home.mccleary/jc3893/30x30/extremes/sp_in_cells_",
                     season,"_",temperature,"_",ver,".csv"))
  }
}

# Compare observed site richness to estimated richness
# and observed/estimated species presence/absence
setwd("/vast/palmer/home.mccleary/jc3893/30x30/extremes/")
# load richness table
val <- read_csv(paste0("site_richness_validation_withsp_out.02_",ver,"_comp.csv"))
# define species columns
spnames <- names(val[8:696])
# loop seasons, subset data seasonally
for (season in c("summer", "winter")){
  vals <- val[val$season==season,]
  vals = vals[vals$est.means>0,]
  # loop model types
  for (temperature in c("means", "variability", "extremes")){
    # get richness and sp list per cell
    vals$est <- vals[,paste0("est.",temperature)][[1]] 
    # lists of sp observed 
    sp_est <- read_csv(paste0("sp_in_cells_",season,"_",temperature,".csv"))
    est_points <- sp_est
    coordinates(est_points) <- c("x","y")
    for (i in 1:nrow(vals)){
      cell <- vals[i,"cell"]
      #find closest cell to point, get sp
      val_coords <- vals[i, c("x","y")]
      coordinates(val_coords) <- c("x","y")
      crs(val_coords) <- crs(est_points) <- 
        "+proj=cea +lat_ts=30 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs "
      wm <- which.min(pointDistance(val_coords, est_points))
      spr <- sort(names(sp_est)[(which(sp_est[wm,3:ncol(sp_est)]==1)) +2])
      # site level species list
      sppt <- sort(names(vals)[(which(vals[i,spnames]==1)) +7])
      # same for absent species
      spr_abs <- setdiff(names(vals[spnames]), spr)
      sppt_abs <- setdiff(names(vals[spnames]), sppt)
      # get true pres, false pres, true abs, false abs
      tpres <- length(which(spr %in% sppt))
      fpres <- length(setdiff(spr, sppt))
      tabs <- length(which(spr_abs %in% sppt_abs))
      fabs <- length(setdiff(spr_abs, sppt_abs))
      pres_rate <- tpres/(tpres+fpres)
      abs_rate <- tabs/(tabs+fabs)
      fpres_rate <- fpres/(tpres+fpres)
      fabs_rate <- fabs/(tabs+fabs)
      ### species-level accuracy by location
      # first find which species are correctly or incorrectly pres/abs
      sp_corr_pres <- spr[which(spr %in% sppt)]
      sp_corr_abs <- spr_abs[which(spr_abs %in% sppt_abs)]
      sp_false_pres <- spr[which(spr %in% sppt_abs)]
      sp_false_abs <- spr_abs[which(spr_abs %in% sppt)]
      # then assign species level values
      loc_tpres <- which(spnames %in% sp_corr_pres)
      loc_tabs <- which(spnames %in% sp_corr_abs)
      loc_fpres <- which(spnames %in% sp_false_pres)
      loc_fabs <- which(spnames %in% sp_false_abs)
      pred <- rep(NA, length(spnames))
      pred[loc_tpres] <- "tpres"
      pred[loc_tabs] <- "tabs"
      pred[loc_fpres] <- "fpres"
      pred[loc_fabs] <- "fabs"
      # collect location level info
      pa <- c(cell[[1]], season, temperature, tpres, fpres, tabs, fabs, 
              pres_rate, abs_rate, fpres_rate, fabs_rate)
      if(i==1 & temperature=="means"){coll <- pa}else{coll <- rbind(coll, pa)}
      # collect species-level info
      spl <- c(cell[[1]], season, temperature, pred)
      if(i==1 & temperature=="means"){splev <- spl}else{splev <- rbind(splev, spl)}
    }
    # species-level metrics
    splev <- data.frame(splev)
    colnames(splev) <- c("cell","season","temperature",spnames)
    sub_splev <- splev[,4:ncol(splev)]
    foo = function(x){
      len = length(x)
      t_pres = sum(x == "tpres")/len
      t_abs = sum(x == "tabs")/len
      f_pres = sum(x == "fpres")/len
      f_abs = sum(x == "fabs")/len
      return(c(t_pres, t_abs, f_pres, f_abs))
    }
    sp_sum <- t(sapply(sub_splev, foo))
    sp_sum <- cbind(rownames(sp_sum),season, temperature, sp_sum)
    sp_sum <- data.frame(sp_sum)
    rownames(sp_sum) <- NULL
    colnames(sp_sum) <- c("spname","season","temperature","tpres","tabs","fpres","fabs")
    # collect info
    if(temperature=="means"){sptab <- sp_sum}else{sptab <- rbind(sptab, sp_sum)}
    
    # richness scatter plot - fig. s2
    if(season=="summer" & temperature=="means"){col="forestgreen"}
    if(season=="summer" & temperature=="variability"){col="orange"}
    if(season=="summer" & temperature=="extremes"){col="purple"}
    if(season=="winter" & temperature=="means"){col="darkolivegreen"}
    if(season=="winter" & temperature=="variability"){col="darkorange3"}
    if(season=="winter" & temperature=="extremes"){col="mediumpurple3"}
    if(temperature=="means"){ylab="Site Richness"}else{ylab=NULL}
    if(season=="summer"){xlab=NULL}else{xlab="Estimated richness\nfrom stacked SDMs"}
    gg <-  ggplot(vals, aes(est, rich)) +
      geom_point(size=2, pch=21, color=col) +
      geom_smooth(color="black", method="lm") +
      geom_abline(intercept = 0, slope = 1, lty=2, color="gray50") +
      theme(text = element_text(size=12),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            axis.line = element_line(color = "black"),
            axis.text.y = element_text(color = "black", size=12),
            axis.text.x = element_text(color = "black", size=12)) +
      xlim(0,max(c(vals$est.means,vals$est.variability,vals$est.extremes))) +
      xlab(xlab) +
      ylab(ylab)
    assign(paste0("gg.",season,".",temperature), gg)
    
    # define axis labels
    if(temperature=="means"){ylab="True negatives rate"}else{ylab=NULL}
    if(season=="summer"){xlab=NULL}else{xlab="True positives rate"}
    # compile data, name columns
    coll.df <- data.frame(coll)
    colnames(coll.df) <- c("cell","season","temperature","true_pres","false_pres",
                           "true_abs","false_abs","pres_rate","abs_rate",
                           "fpres_rate","fabs_rate")
    # fix column coding
    coll.df <- mutate_at(coll.df, c('pres_rate', 'abs_rate',
                                    'fpres_rate', 'fabs_rate'), as.numeric)
    coll.df <- mutate_at(coll.df, 'temperature', as.factor)
    coll.df$temperature <- factor(coll.df$temperature, 
                                  levels=c('means', 'variability', 'extremes'))
    # calculate deltas between means and other two models
    coll1 <- coll.df[coll.df$temperature=="means",]
    collc <- coll.df[coll.df$temperature!="means",]
    collc$delta_pres <- collc$pres_rate - coll1$pres_rate
    collc$delta_abs <- collc$abs_rate - coll1$abs_rate
    collc$delta_fpres <- collc$fpres_rate - coll1$fpres_rate
    collc$delta_fabs <- collc$fabs_rate - coll1$fabs_rate
    
    coll_sub <- coll.df[coll.df$temperature==temperature,] 
    # scatter plot pres vs abs rates (fig. s4)
    pos <-  ggplot(coll_sub, aes(pres_rate, abs_rate)) +
      geom_point(size=2) +
      geom_smooth(color="black", method="lm") +
      theme(text = element_text(size=12),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            axis.line = element_line(color = "black"),
            axis.text.y = element_text(color = "black", size=12),
            axis.text.x = element_text(color = "black", size=12)) +
      ylim(0.9,1) +
      xlim(0,1) +
      xlab(xlab) +
      ylab(ylab) +
      ggtitle(paste(temperature,season))
    assign(paste0("pos.",season,".",temperature), pos)
    
    # plot trends during the last run
    if(temperature=="extremes"){
      
      # delta comparisons
      rows <- nrow(vals)
      valp <- cbind(rep("points",rows),vals$rich)
      colnames(valp) <- c("temperature","value")
      valp <- data.frame(valp)
      valp$value <- as.numeric(valp$value)
      # compile, code columns, column names
      val2 <- cbind(c(rep("means",rows),rep("variability",rows),rep("extremes",rows)),
                    c(vals$delta.means,vals$delta.variability,vals$delta.extremes))
      colnames(val2) <- c("temperature","value")
      val2 <- data.frame(val2)
      val2$value <- as.numeric(val2$value)
      val2$temperature <- as.factor(val2$temperature)
      
      # fix column coding - species level data
      sptab <- data.frame(sptab)         
      sptab <- mutate_at(sptab, c('tpres', 'tabs', 'fpres', 'fabs'), as.numeric)   
      # rates (denominator set at minimum 1 to avoid dividing by 0)
      sptab$pres_rate <- sptab$tpres/ifelse((sptab$tpres+sptab$fpres)==0,
                                            1,(sptab$tpres+sptab$fpres))
      sptab$abs_rate <- sptab$tabs/ifelse((sptab$tabs+sptab$fabs)==0,
                                          1,(sptab$tabs+sptab$fabs))
      sptab$fpres_rate <- sptab$fpres/ifelse((sptab$tpres+sptab$fpres)==0,
                                             1,(sptab$tpres+sptab$fpres))
      sptab$fabs_rate <- sptab$fabs/ifelse((sptab$tabs+sptab$fabs)==0,
                                           1,(sptab$tabs+sptab$fabs))
      sptab <- mutate_at(sptab, c('pres_rate', 'abs_rate',
                                  'fpres_rate', 'fabs_rate'), as.numeric)
      sptab <- mutate_at(sptab, 'temperature', as.factor)
      sptab$temperature <- factor(sptab$temperature, 
                                  levels=c('means', 'variability', 'extremes'))
      # deltas to climate means models
      sp1 <- sptab[sptab$temperature=="means",]
      spc <- sptab[sptab$temperature!="means",]
      spc$delta_pres <- spc$pres_rate - sp1$pres_rate
      spc$delta_abs <- spc$abs_rate - sp1$abs_rate
      spc$delta_fpres <- spc$fpres_rate - sp1$fpres_rate
      spc$delta_fabs <- spc$fabs_rate - sp1$fabs_rate
      
      # violin plot - change in percent correctly identified as pres/abs (fig. 1)
      hsize = 1.5; tsize = 15
      if(season=="summer"){col="forestgreen"}else{col="darkolivegreen"}
      # presences
      # raw means
      pv <- ggplot(data=coll1, aes(x=temperature, y=pres_rate)) + 
        geom_hline(yintercept = mean(coll1$pres_rate), size = hsize, color = "gray") +
        geom_violin(trim=F, alpha = .75, fill=col) +
        scale_fill_manual(values = col)  +
        stat_summary(fun.data="mean_sdl", geom="crossbar", width=0.06, fill="white") +
        xlab(NULL) +
        ylab('True positives rate') +
        theme_light() +
        theme(legend.position = "none",
              text = element_text(size = tsize),
              axis.title.y = element_text(size = tsize),
              axis.text.x = element_text(size = tsize))
      # delta variability/extremes
      if(season=="summer"){cols=c("orange","purple")}else{
        cols=c("darkorange3","mediumpurple3")}
      pv2 <- ggplot(data=collc, aes(x=temperature, y=delta_pres, fill=temperature)) + 
        geom_hline(yintercept = 0, size = hsize, color = "gray") +
        geom_violin(trim=F, alpha = .75) +
        scale_fill_manual(values = cols) +
        stat_summary(fun.data="mean_sdl", geom="crossbar", width=0.06, fill="white") +
        ylab(expr(paste(Delta, " True positives rate"))) +
        xlab(NULL) +
        theme_light() +
        theme(legend.position = "none",
              text = element_text(size = tsize),
              axis.title.y = element_text(size = tsize),
              axis.text.x = element_text(size = tsize))
      assign(paste0("pv.",season), pv)
      assign(paste0("pv2.",season), pv2)
      
      # absences
      # raw means
      if(season=="summer"){col="forestgreen"}else{col="darkolivegreen"}
      nv <- ggplot(data=coll1, aes(x=temperature, y=abs_rate)) + 
        geom_hline(yintercept = mean(coll1$abs_rate), size = hsize, color = "gray") +
        geom_violin(trim=F, alpha = .75, fill=col) +
        scale_fill_manual(values = col)  +
        stat_summary(fun.data="mean_sdl", geom="crossbar", width=0.06, fill="white") +
        xlab(NULL) +
        ylab('True negatives rate') +
        theme_light() +
        theme(legend.position = "none",
              text = element_text(size = tsize),
              axis.title.y = element_text(size = tsize),
              axis.text.x = element_text(size = tsize))
      # delta variability/extremes
      if(season=="summer"){cols=c("orange","purple")}else{
        cols=c("darkorange3","mediumpurple3")}
      nv2 <- ggplot(data=collc, aes(x=temperature, y=delta_abs, fill=temperature)) + 
        geom_hline(yintercept = 0, size = hsize, color = "gray") +
        geom_violin(trim=F, alpha = .75) +
        scale_fill_manual(values = cols) +
        stat_summary(fun.data="mean_sdl", geom="crossbar", width=0.06, fill="white") +
        ylab(expr(paste(Delta, " True negatives rate"))) +
        xlab(NULL) +
        theme_light() +
        theme(legend.position = "none",
              text = element_text(size = tsize),
              axis.title.y = element_text(size = tsize),
              axis.text.x = element_text(size = tsize))
      assign(paste0("nv.",season), nv)
      assign(paste0("nv2.",season), nv2)
      
      # violin plot - change in percent positive or negative of false pres/abs (fig. s3)
      # presences
      # raw means
      if(season=="summer"){col="forestgreen"}else{col="darkolivegreen"}
      fpv <- ggplot(data=coll1, aes(x=temperature, y=fpres_rate)) + 
        geom_hline(yintercept = mean(coll1$fpres_rate), size = hsize, color = "gray") +
        geom_violin(trim=F, alpha = .75, fill=col) +
        scale_fill_manual(values = col)  +
        stat_summary(fun.data="mean_sdl", geom="crossbar", width=0.06, fill="white") +
        xlab(NULL) +
        ylab('False positives rate') +
        theme_light() +
        theme(legend.position = "none",
              text = element_text(size = tsize),
              axis.title.y = element_text(size = tsize),
              axis.text.x = element_text(size = tsize))
      # delta variability/extremes
      if(season=="summer"){cols=c("orange","purple")}else{
        cols=c("darkorange3","mediumpurple3")}
      fpv2 <- ggplot(data=collc, aes(x=temperature, y=delta_fpres, fill=temperature)) + 
        geom_hline(yintercept = 0, size = hsize, color = "gray") +
        geom_violin(trim=F, alpha = .75) +
        scale_fill_manual(values = cols) +
        stat_summary(fun.data="mean_sdl", geom="crossbar", width=0.06, fill="white") +
        ylab(expr(paste(Delta, " False positives rate"))) +
        xlab(NULL) +
        theme_light() +
        theme(legend.position = "none",
              text = element_text(size = tsize),
              axis.title.y = element_text(size = tsize),
              axis.text.x = element_text(size = tsize))
      assign(paste0("fpv.",season), fpv)
      assign(paste0("fpv2.",season), fpv2)
      
      # absences
      # raw means
      if(season=="summer"){col="forestgreen"}else{col="darkolivegreen"}
      fnv <- ggplot(data=coll1, aes(x=temperature, y=fabs_rate)) + 
        geom_hline(yintercept = mean(coll1$fabs_rate), size = hsize, color = "gray") +
        geom_violin(trim=F, alpha = .75, fill=col) +
        scale_fill_manual(values = col)  +
        stat_summary(fun.data="mean_sdl", geom="crossbar", width=0.06, fill="white") +
        xlab(NULL) +
        ylab('False negatives rate') +
        theme_light() +
        theme(legend.position = "none",
              text = element_text(size = tsize),
              axis.title.y = element_text(size = tsize),
              axis.text.x = element_text(size = tsize))
      # delta variability/extremes
      if(season=="summer"){cols=c("orange","purple")}else{
        cols=c("darkorange3","mediumpurple3")}
      fnv2 <- ggplot(data=collc, aes(x=temperature, y=delta_fabs, fill=temperature)) + 
        geom_hline(yintercept = 0, size = hsize, color = "gray") +
        geom_violin(trim=F, alpha = .75) +
        scale_fill_manual(values = cols) +
        stat_summary(fun.data="mean_sdl", geom="crossbar", width=0.06, fill="white") +
        ylab(expr(paste(Delta, " False negatives rate"))) +
        xlab(NULL) +
        theme_light() +
        theme(legend.position = "none",
              text = element_text(size = tsize),
              axis.title.y = element_text(size = tsize),
              axis.text.x = element_text(size = tsize))
      assign(paste0("fnv.",season), fnv)
      assign(paste0("fnv2.",season), fnv2)
      
      # means vs extremes pos/neg rates (fig. s3)
      coll.comp = data.frame(means.pres=coll.df[coll.df$temperature=="means",]$pres_rate,
                             variability.pres=coll.df[coll.df$temperature=="variability",]$pres_rate,
                             extremes.pres=coll.df[coll.df$temperature=="extremes",]$pres_rate,
                             means.abs=coll.df[coll.df$temperature=="means",]$abs_rate,
                             variability.abs=coll.df[coll.df$temperature=="variability",]$abs_rate,
                             extremes.abs=coll.df[coll.df$temperature=="extremes",]$abs_rate)
      
      ttype = c("variability","extremes")
      ttypename = c("Variability","Extremes")
      for (tt in 1:2){
        if(season=="summer" & tt==1){col="orange"}
        if(season=="summer" & tt==2){col="purple"}
        if(season=="winter" & tt==1){col="darkorange3"}
        if(season=="winter" & tt==2){col="mediumpurple3"}
        coll.comp$var.pres = coll.comp[,paste0(ttype[tt],".pres")]
        coll.comp$var.abs = coll.comp[,paste0(ttype[tt],".abs")]
        
        # presences
      ggpres <-  ggplot(coll.comp, aes(means.pres, var.pres)) +
        geom_point(size=2, color=col, pch=21) +
        geom_smooth(color="black", method="lm") +
        geom_abline(intercept = 0, slope = 1, lty=2, color="gray50") +
        theme(text = element_text(size=12),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              axis.line = element_line(color = "black"),
              axis.text.y = element_text(color = "black", size=12),
              axis.text.x = element_text(color = "black", size=12)) +
        xlab("Means - True presence rate") +
        ylab(paste(ttypename[tt],"- \nTrue presence rate"))
      assign(paste0("gg.",season,".",ttype[tt],".pres"), ggpres)
      # absences
      ggabs <-  ggplot(coll.comp, aes(means.abs, var.abs)) +
        geom_point(size=2, color=col, pch=21) +
        geom_smooth(color="black", method="lm") +
        geom_abline(intercept = 0, slope = 1, lty=2, color="gray50") +
        theme(text = element_text(size=12),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              axis.line = element_line(color = "black"),
              axis.text.y = element_text(color = "black", size=12),
              axis.text.x = element_text(color = "black", size=12)) +
        xlab("Means - True absence rate") +
        ylab(paste(ttypename[tt],"- \nTrue absence rate"))
      assign(paste0("gg.",season,".",ttype[tt],".abs"), ggabs)
      }
      
    }
    
    # stats for validation
    beta <- coef(summary(lm(est~rich, data=vals)))[2,1]
    r2 <- summary(lm(est~rich, data=vals))$r.squared
    rmse <- sqrt(mean((vals$rich - vals$est)^2))
    smape <- smape(vals$rich, vals$est)
    row <- c(season, temperature, beta, r2, rmse, smape)
    if(temperature=="means" & season=="summer"){
      table <- row}else{table <- rbind(table, row)}
  }
  # column naming, save
  colnames(coll) <- c("cell","season","temperature","true_pres","false_pres",
                      "true_abs","false_abs","pres_rate","abs_rate",
                      "fpres_rate","fabs_rate")
  write_csv(data.frame(coll), paste0("point_performance_",season,"_.02_",ver,".csv"))
}

{ # saving multipanel figures
gga <- ggarrange(gg.summer.means, gg.summer.variability, gg.summer.extremes,
                 gg.winter.means, gg.winter.variability, gg.winter.extremes,
                 nrow=2, ncol=3, widths = c(1.15,1,1))
ggsave(paste0("validation_scatter_.02_",ver,".jpeg"), gga, "jpeg",
       height=6, width=10, units="in")

posa <- ggarrange(pos.summer.means, pos.summer.variability, pos.summer.extremes,
                  pos.winter.means, pos.winter.variability, pos.winter.extremes,
                  nrow=2, ncol=3, widths = c(1.15,1,1))
ggsave(paste0("validation_pos_neg_rate_scatter.02_",ver,".jpeg"), posa, "jpeg",
       height=6, width=7, units="in")

pngg <- ggarrange(pv.summer, pv2.summer, nv.summer, nv2.summer, 
                  pv.winter, pv2.winter, nv.winter, nv2.winter, 
                  nrow=2, ncol=4, widths = c(1.5,3.5), 
                  labels=c("a)","","b)","","c)","","d)",""))
ggsave(paste0("validation_posneg_rate_violin_.02_",ver,".jpeg"), pngg, "jpeg",
       height=6, width=12, units="in")

fgg <- ggarrange(fpv.summer, fpv2.summer, fnv.summer, fnv2.summer, 
                 fpv.winter, fpv2.winter, fnv.winter, fnv2.winter, 
                 nrow=2, ncol=4, widths = c(1.5,3.5), 
                 labels=c("a)","","b)","","c)","","d)",""))
ggsave(paste0("validation_fposneg_rate_violin_.02_",ver,".jpeg"), fgg, "jpeg",
       height=6, width=12, units="in")

pagg <- ggarrange(gg.summer.variability.pres, gg.summer.extremes.pres,
                  gg.summer.variability.abs, gg.summer.extremes.abs,
                  gg.winter.variability.pres, gg.winter.extremes.pres,
                  gg.winter.variability.abs, gg.winter.extremes.abs,
                 nrow=2, ncol=4, widths = c(1,1,1,1), 
                 labels=c("a)","b)","c)","d)","e)","f)","g)","h)"))
ggsave(paste0("validation_posneg_rate_scatter_comp_.02_",ver,".jpeg"), pagg, "jpeg",
       height=6, width=12, units="in")

colnames(table) <- c("season","temperature","coef","r2","rmse","smape")
write_csv(data.frame(table), paste0("validation_table_.02_",ver,".csv"))
}



### maps of validation sites
  # transform for plotting
  vals_sp <- st_as_sf(vals, coords=c("x","y")) %>%
  st_set_crs("+proj=cea +lat_ts=30 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs") %>%
  st_transform(crs="+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs") %>%
  sf:::as_Spatial() %>%
  data.frame()
  # load country polygons
world <- ne_countries(scale='medium', returnclass = 'sf') %>%
  st_transform(crs="+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")

# points (fig s1)
map <- ggplot() +
  geom_sf(data=world) + 
  coord_sf(xlim=c(min(vals_sp$coords.x1)-50000, max(vals_sp$coords.x1)+50000),
           ylim=c(min(vals_sp$coords.x2)-50000, max(vals_sp$coords.x2)+50000), expand=F) +
  theme_bw() +
  theme(text = element_text(size=13),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"),
        axis.text.x = element_text(color = "black", size=8),
        axis.text.y = element_text(color = "black", size=8),
        plot.title = element_text(size = 22, face = "bold", hjust=.5)) +
  geom_point(vals_sp, mapping=aes(x=coords.x1, y=coords.x2)) +
  xlab("") +
  ylab("")
ggsave(paste0("validation_map.jpeg"), map, "jpeg",
       height=6, width=6, units="in")


# points colored by presence or absence rate (figs s5-s6)
for (season in c("summer","winter")){
# merge with pres/abs rates
coll = read_csv(paste0("point_performance_",season,"_.02_",ver,".csv"))
vals_sp_s = left_join(vals_sp, coll, by="cell")
# colors
if (season=="summer"){cols = c("forestgreen", "purple")}else{
  cols = c("darkolivegreen", "mediumpurple3")}
# subset
vals_sp_s_means = vals_sp_s[vals_sp_s$temperature=="means",]
# means, presence
mapmp <- ggplot() +
  geom_sf(data=world) + 
  coord_sf(xlim=c(min(vals_sp$coords.x1)-50000, max(vals_sp$coords.x1)+50000),
           ylim=c(min(vals_sp$coords.x2)-50000, max(vals_sp$coords.x2)+50000), expand=F) +
  theme_bw() +
  theme(text = element_text(size=13),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"),
        axis.text.x = element_text(color = "black", size=8),
        axis.text.y = element_text(color = "black", size=8)) +
  geom_point(vals_sp_s_means, mapping=aes(x=coords.x1, y=coords.x2, col=pres_rate)) +
  scale_colour_gradient(low="white", high=cols[1], limits=c(0,1), name="True Presence rate") +
  xlab("") +
  ylab("")

# means, absence
mapma <- ggplot() +
  geom_sf(data=world) + 
  coord_sf(xlim=c(min(vals_sp$coords.x1)-50000, max(vals_sp$coords.x1)+50000),
           ylim=c(min(vals_sp$coords.x2)-50000, max(vals_sp$coords.x2)+50000), expand=F) +
  theme_bw() +
  theme(text = element_text(size=13),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"),
        axis.text.x = element_text(color = "black", size=8),
        axis.text.y = element_text(color = "black", size=8)) +
  geom_point(vals_sp_s_means, mapping=aes(x=coords.x1, y=coords.x2, col=abs_rate)) +
  scale_colour_gradient(low="white", high=cols[1], limits=c(0.8,1), name="True Absence rate") +
  xlab("") +
  ylab("")

# extremes, presence
vals_sp_s_extremes = vals_sp_s[vals_sp_s$temperature=="extremes",]
mapep <- ggplot() +
  geom_sf(data=world) + 
  coord_sf(xlim=c(min(vals_sp$coords.x1)-50000, max(vals_sp$coords.x1)+50000),
           ylim=c(min(vals_sp$coords.x2)-50000, max(vals_sp$coords.x2)+50000), expand=F) +
  theme_bw() +
  theme(text = element_text(size=13),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"),
        axis.text.x = element_text(color = "black", size=8),
        axis.text.y = element_text(color = "black", size=8)) +
  geom_point(vals_sp_s_extremes, mapping=aes(x=coords.x1, y=coords.x2, col=pres_rate)) +
  scale_colour_gradient(low="white", high=cols[2], limits=c(0,1), name="True Presence rate") +
  xlab("") +
  ylab("")

# extremes, absence
mapea <- ggplot() +
  geom_sf(data=world) + 
  coord_sf(xlim=c(min(vals_sp$coords.x1)-50000, max(vals_sp$coords.x1)+50000),
           ylim=c(min(vals_sp$coords.x2)-50000, max(vals_sp$coords.x2)+50000), expand=F) +
  theme_bw() +
  theme(text = element_text(size=13),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"),
        axis.text.x = element_text(color = "black", size=8),
        axis.text.y = element_text(color = "black", size=8)) +
  geom_point(vals_sp_s_extremes, mapping=aes(x=coords.x1, y=coords.x2, col=abs_rate), cex=2) +
  scale_colour_gradient(low="white", high=cols[2], limits=c(0.8,1), name="True Absence rate") +
  xlab("") +
  ylab("")

gga <- ggarrange(mapmp, mapma, mapep, mapea, nrow=2, ncol=2,
                  labels=c("a)","b)","c)","d)"))
ggsave(paste0("validation_map_",season,".jpeg"), gga, "jpeg",
       height=7, width=15, units="in", bg="white")
}



