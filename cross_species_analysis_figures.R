# cross-scale patterns

{ library(tidyverse)
  library(RColorBrewer)
  library(MetBrewer)
  library(gridExtra)
  library(sf)
  library(ggpubr)
  library(car)
  library(raster)
  library(nlme)
  library(ape)
  library(phytools)
  library(car)
  ver <- "VE8"}

######### MODEL METRICS & AREA OF OCCURRENCE & PREDICTOR IMPORTANCES
data <- read_csv(paste0("~/extremes/sdm_birds_",ver,"_0.csv")) # load model outputs
# calculate % importance of climate, topographic, landcover suites
data$clim <- rowSums(data[,c(31:32,71:76)], na.rm=T)/rowSums(data[,31:76], na.rm=T)
data$top <- rowSums(data[,33:34], na.rm=T)/rowSums(data[,31:76], na.rm=T)
data$lc <- rowSums(data[,35:70], na.rm=T)/rowSums(data[,31:76], na.rm=T)
# remove hawaiian and seabirds (won't be needed)
list <- read_csv("~/extremes/species_list_2024.csv")
data = left_join(data, list[,c(1,8:9)], by=c("common_name"="Common"))
data = data[data$hawaii=="no" & data$seabird=="no",]
# Trait data
tpdata <- read_csv("~/extremes/species_list_2024_BirdTree-AVONET.csv")
tpdata$Scientific = gsub(" ","_",tpdata$Scientific)
# Combine to species level data
data = left_join(data, tpdata, join_by("sciname"=="Scientific"))
# save
write_csv(data, paste0("~/extremes/sdm_birds_",ver,"_0.csv"))

### Deltas (change between Means and Variability/Extremes models)
data <- read_csv(paste0("~/extremes/sdm_birds_",ver,"_0.csv")) %>%
  distinct(sciname, season, temperature, .keep_all=T)
# define model metrics
metrics <- c("AUC","SPSthreshold",
             "TSS","Sensitivity","Specificity")
# Get change between models
# loop seasons
rm(out)
for (season in c("summer", "winter")){
  # loop species
  for (i in 1:length(unique(data$sciname))){
    # define species, subset by species
    sci_name <- unique(data$sciname)[i]
    common <- unique(data$common_name)[i]
    sub <- data[data$sciname==sci_name & data$season==season,
                c("temperature","range_area",metrics,"clim","top","lc","prop_in_range")]
    traits = data[data$sciname==sci_name & data$season==season, c(73:94)][1,]
    # calculate differences for model metrics, area of range, proportion in range, compile
    if(nrow(sub)>0){
      for (temp in c("variability","extremes")){
        difs <- round(sub[sub$temperature==temp,3:10]-sub[1,3:10], 4)
        aocdif <- round(((sub[sub$temperature==temp,"range_area"]/
                            round(sub[1,"range_area"], 4) - 1) * 100),4)
        pirdif <- round(((sub[sub$temperature==temp,"prop_in_range"]/
                            round(sub[1,"prop_in_range"], 4) - 1) * 100),4)
        collect <- c(common, sci_name, season, temp, as.numeric(aocdif),
                     as.numeric(pirdif), as.numeric(difs), unlist(traits))
        if(temp!="variability"){spcollect <- rbind(spcollect, collect)}else{spcollect <- collect}
      }
      if(!exists("out")){out <- spcollect}else{out <- rbind(out, spcollect)}
    }
  }
}
# define columns and save
colnames(out) <- c("common_name", "sciname", "season", "temp",
                   "delta_range_area", "delta_prop_in_range",
                   paste0("delta_",metrics),
                   "delta_clim", "delta_top","delta_lc",names(traits))
rownames(out) <- NULL
out <- data.frame(out)
write_csv(out, paste0("~/extremes/birds_",ver,"_comp.csv"))


# Range area by model plot (fig 3)
out <- read_csv(paste0("~/extremes/sdm_birds_",ver,"_0.csv"))
out <- out[out$temperature=="means",]
out$range_area <- log(out$range_area) # log range area
# plot absolute ranges from means models (panel a)
aoc.means <- ggplot(out, aes(x=range_area, fill=season)) +
  geom_density(alpha=.3) +
  scale_fill_manual(values=c("forestgreen","darkolivegreen")) +
  geom_vline(xintercept=0, color="black", linetype="dashed", linewidth=1) +
  theme(legend.position = c(.23,.8),
        legend.title = element_blank(),
        legend.text=element_text(size=12),
        axis.text.y=element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(color = "black"),
        axis.text.x = element_text(color = "black", size=15),
        plot.title = element_text(size = 15, face = "bold", hjust=.5)) +
  xlim(5,18) +
  theme(text = element_text(size=13)) +
  xlab(expression(Range ~ size ~ (log(km^2)))) +
  ylab("Density") +
  ggtitle("means")

# plot delta ranges from other models (panels b-c) 
out2 <- read_csv(paste0("~/extremes/birds_",ver,"_comp.csv"))
for (temp in c("variability","extremes")){
  if(temp=="variability"){cols=c("orange","darkorange3")
  }else{cols=c("purple","mediumpurple3")}
  out3 = out2[out2$temp==temp,]
  out3$season = factor(out3$season)
  # plot
  aocgg <- ggplot(out3, aes(x=delta_range_area, fill=season)) +
    geom_density(alpha=.3, bounds=c(-100,100)) +
    scale_fill_manual(values=cols) +
    geom_vline(xintercept=0, color="gray50", linetype="dashed", cex=.7) +
    geom_vline(xintercept=median(out3[out3$season=="summer",]$delta_range_area),
               color="black", linetype="dashed", cex=1) +
    geom_vline(xintercept=median(out3[out3$season=="winter",]$delta_range_area),
               color="black", linetype="dashed", cex=1) +
    xlim(-40,40) +
    theme(legend.position = c(.78,.8),
          legend.title = element_blank(),
          legend.text=element_text(size=12),
          axis.text.y=element_blank(),
          axis.ticks.y = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(color = "black"),
          axis.text.x = element_text(color = "black", size=15),
          plot.title = element_text(size = 15, face = "bold", hjust=.5)) +
    theme(text = element_text(size=13)) +
    xlab(expression("% Difference in Range Size")) +
    ylab("") +
    ggtitle(temp)
  assign(paste0("aoc.",temp), aocgg)
}
gg.aoc <- ggarrange(aoc.means, aoc.variability, aoc.extremes, nrow=1, ncol=3)
ggsave(paste0("~/extremes/extremes_aoc_density_",ver,".jpeg"),
       gg.aoc, "jpeg", height=3, width=9, units="in")

# Violin plot for model performance metrics (fig. s7)
data <- read_csv(paste0("~/extremes/sdm_birds_",ver,"_0.csv"))
# absolute values for means
for(season in c("summer","winter")){
  if(season=="summer"){col="forestgreen"}else{col="darkolivegreen"}
  datameans <- data[data$temperature=="means" & data$season==season,]
  gg <- ggplot(data=datameans, aes(x=temperature, y=AUC)) +
    geom_hline(yintercept = median(datameans$AUC), size = 1.5, color = "gray") +
    geom_violin(trim=F, alpha = .6, fill=col) +
    scale_fill_manual(values = col)  +
    stat_summary(fun="median", geom="point", shape=21, size=4,
                 stroke=1.5, fill="white", col="black") +
    xlab(NULL) +
    ylab("AUC") +
    theme_light() +
    theme(legend.position = "none",
          text = element_text(size = 30))
  
  out <- read_csv(paste0("~/extremes/birds_",ver,"_comp.csv"))
  out2 <- out[out$season==season,]
  out2 <- mutate(out2, temp=factor(temp,levels=c("limits","variability", "extremes")))
  
  # remove extreme outliers (<2%, >98%) for prettier plotting - they are still included in statistical analyses
  quantiles <- quantile(out2$delta_AUC, probs=c(.02, .98), na.rm = T)
  out3 <- subset(out2, out2$delta_AUC > quantiles[1]  & out2$delta_AUC < quantiles[2] )
  
  # deltas for variability/extremes
  if(season=="summer"){cols=c("orange","purple")}else{
    cols=c("darkorange3","mediumpurple3")}
  gg2 <- ggplot(data=out3, aes(x=temp, y=delta_AUC, fill=temp)) +
    geom_hline(yintercept = 0, size = 1.5, color = "gray") +
    geom_violin(trim=F, alpha = .6) +
    scale_fill_manual(values = cols) +
    stat_summary(fun="median", geom="point", shape=21, size=4,
                 stroke=1.5, fill="white", col="black") +
    ylab(expr(paste(Delta, " AUC"))) +
    xlab(NULL) +
    theme_light() +
    theme(legend.position = "none",
          text = element_text(size = 30))
  gg2
  assign(paste0("gg.AUC.",season),gg)
  assign(paste0("gg2.AUC.",season),gg2)
}
# save
gg.ppm <- ggarrange(gg.AUC.summer, gg2.AUC.summer, gg.AUC.winter, gg2.AUC.winter,
                    nrow=2, ncol=2, widths = c(2,4), labels=c("a)","","b)",""),
                    font.label=list(color="black",size=25))
ggsave(paste0("~/extremes/extremes_auc_violin_",ver,"2.jpeg"),
       gg.ppm, "jpeg", height=10, width=9, units="in")

# Predictor importances (fig. s9)
# Violin plot - raw means
for(season in c("summer", "winter")){
  if(season=="summer"){col="forestgreen"}else{col="darkolivegreen"}
  for(metric in c("clim","top","lc")){
    if(metric=="clim"){metricname="Climate vars"}
    if(metric=="top"){metricname="Topographic vars"}
    if(metric=="lc"){metricname="Landcover vars"}
    datameans <- data[data$temperature=="means" & data$season==season,]
    datameans$metric <- datameans[,metric][[1]]
    gg <- ggplot(data=datameans, aes(x=temperature, y=metric)) +
      geom_hline(yintercept = mean(datameans$metric), size = 1.5, color = "gray") +
      geom_violin(trim=F, alpha = .6, fill=col) +
      scale_fill_manual(values = col)  +
      stat_summary(fun="median", geom="point", shape=21, size=4,
                   stroke=1.5, fill="white", col="black") +
      xlab(NULL) +
      ylab(metricname) +
      theme_light() +
      theme(legend.position = "none", text = element_text(size = 27))
    # variability/extremes deltas
    out2 <- out[out$season==season,]
    out2$metric <- out2[,paste0("delta_",metric)][[1]]
    if(metric=="range_area"){out2$metric <- ifelse(out2$metric>0,out2$metric^(1/3),(-(-out2$metric)^(1/3)))}
    if(metric=="range_area"){out2$metric <- log(out2$metric)}
    out2 <- mutate(out2, temp=factor(temp,levels=c("variability", "extremes")) )
    if(season=="summer"){cols=c("orange","purple")}else{
      cols=c("darkorange3","mediumpurple3")}
    gg2 <- ggplot(data=out2, aes(x=temp, y=metric, fill=temp)) +
      geom_hline(yintercept = 0, size = 1.5, color = "gray") +
      geom_violin(trim=F, alpha = .6) +
      scale_fill_manual(values = cols) +
      stat_summary(fun="median", geom="point", shape=21, size=4,
                   stroke=1.5, fill="white", col="black") +
      ylab(expr(paste(Delta, " ", !!metricname))) +
      xlab(NULL) +
      theme_light() +
      theme(legend.position = "none", text = element_text(size = 27))
    gg2
    assign(paste0("gg.",metric,".",season),gg)
    assign(paste0("gg2.",metric,".",season),gg2)
    
  }
}
gg.ppm <- ggarrange(gg.clim.summer, gg2.clim.summer, gg.top.summer,
                    gg2.top.summer, gg.lc.summer, gg2.lc.summer,
                    gg.clim.winter, gg2.clim.winter, gg.top.winter,
                    gg2.top.winter, gg.lc.winter, gg2.lc.winter,
                    nrow=2, ncol=6, widths = c(2,3.5,2,3.5,2,3.5),
                    labels=c("a)","","b)","","c)","","d)","","e)","","f)",""),
                    font.label=list(color="black",size=25))
ggsave(paste0("~/extremes/extremes_predimps_violin_",ver,".jpeg"),
       gg.ppm, "jpeg", height=10, width=24, units="in")


# How much of the prediction falls within expert range? (fig. 4)
# raw means
for(season in c("summer","winter")){
  if(season=="summer"){col="forestgreen"}else{col="darkolivegreen"}
  datameans <- data[data$temperature=="means" & data$season==season,]
  gg <- ggplot(data=datameans, aes(x=temperature, y=prop_in_range*100)) +
    geom_hline(yintercept = mean(datameans$prop_in_range)*100, size = 1.5, color = "gray") +
    geom_violin(trim=F, alpha = .6, fill=col) +
    scale_fill_manual(values = col)  +
    stat_summary(fun="median", geom="point", shape=21, size=4,
                 stroke=1.5, fill="white", col="black") +
    xlab(NULL) +
    ylab("Percent prediction\nwithin range") +
    theme_light() +
    theme(legend.position = "none",
          text = element_text(size = 23))
  # variability/extreme deltas
  out2 <- out[out$season==season,]
  out2 <- mutate(out2, temp=factor(temp,levels=c("variability", "extremes")) )
  
  # average delta proportion in range for variability/extremes vs means
  print(median(out2[out2$temp=="variability",]$delta_prop_in_range, na.rm=T))
  print(median(out2[out2$temp=="extremes",]$delta_prop_in_range, na.rm=T))
  
  # remove extreme values for plotting, original data are extremely long tailed
  # all data is included for statistical analysis
  quantiles <- quantile(out2$delta_prop_in_range, probs=c(.02, .98), na.rm = T)
  out2 <- subset(out2, out2$delta_prop_in_range > quantiles[1]  & out2$delta_prop_in_range < quantiles[2] )
  out2$delta_prop_in_range = as.numeric(out2$delta_prop_in_range)
  out2 = out2[out2$temp %in% c("variability","extremes"),]
  # plot
  if(season=="summer"){cols=c("orange","purple")}else{
    cols=c("darkorange3","mediumpurple3")}
  gg2 <- ggplot(data=out2, aes(x=temp, y=delta_prop_in_range, fill=temp)) +
    geom_hline(yintercept = 0, size = 1.5, color = "gray") +
    geom_violin(trim=F, alpha = .6) +
    scale_fill_manual(values = cols) +
    stat_summary(fun="median", geom="point", shape=21, size=4,
                 stroke=1.5, fill="white", col="black") +
    ylab(expr(paste(Delta, " Percent prediction"))) +
    xlab(NULL) +
    theme_light() +
    theme(legend.position = "none", text = element_text(size = 23))
  gg2
  
  assign(paste0("gg.",season),gg)
  assign(paste0("gg2.",season),gg2)
}
gg.ppm <- ggarrange(gg.summer, gg2.summer, gg.winter, gg2.winter,
                    nrow=2, ncol=2, widths = c(2,3.5),
                    labels=c("a)","","b)",""),
                    font.label=list(color="black",size=20))
ggsave(paste0("~/extremes/extremes_propinrange_violin_",ver,".jpeg"),
       gg.ppm, "jpeg", height=8, width=8, units="in")


# Autocorrelation values for table s1
for(season in c("summer","winter")){
  datas = data[data$season==season,]
  print(mean(datas$train.ac.mean))
  print(mean(datas$train.sc.sd))
  print(mean(datas$test.ac.mean))
  print(mean(datas$test.ac.sd))
}

# species counts for methods
length(unique(data$common_name))
nrow(data[data$season=="summer",])/3
nrow(data[data$season=="winter",])/3


