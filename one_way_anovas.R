## One-way anovas to explore model-level differences across species

library(WRS2)
setwd("~/extremes/")

# load species level model outputs and deltas
dat <- read.csv("birds_VE8_comp.csv", as.is=TRUE)
sdm <- read.csv("sdm_birds_VE8_0.csv", as.is=TRUE)

## Do range sizes differ by model type for summer/winter?
# loop seasons and species, calculate percent of maximum range across 3 models
rm(out)
for (season in c("summer", "winter")){
  for (i in 1:length(unique(data$sciname))){
    sci_name <- unique(data$sciname)[i]
    common <- unique(data$common_name)[i]
    sub <- data[data$sciname==sci_name & data$season==season,
                c("sciname","common_name","season","temperature","range_area")]
    # calculate differences
    if(nrow(sub)==3){
        sub$range_area_pct = sub$range_area/max(sub$range_area)
      if(!exists("out")){out <- sub}else{out <- rbind(out, sub)}
    }
  }
}
# models across species (one way anovas plus post-hoc tests by season)
t1s.range = t1waybt(range_area_pct ~ temperature, data=out[out$season=="summer",], tr = 0.2, nboot = 9999)
mcs.range = mcppb20(range_area_pct ~ temperature, data=out[out$season=="summer",], tr = 0.2, nboot = 9999)
write.csv(mcs.range$comp, "mcs_range.csv")
t1w.range = t1waybt(range_area_pct ~ temperature, data=out[out$season=="winter",], tr = 0.2, nboot = 9999)
mcw.range = mcppb20(range_area_pct ~ temperature, data=out[out$season=="winter",], tr = 0.2, nboot = 9999)
write.csv(mcw.range$comp, "mcw_range.csv")
# compile
col.range = rbind(t1s.range[1:2], t1w.range[1:2])

## Does proportion in range differ by model type for summer/winter?
# (one way anovas plus post-hoc tests by season)
t1s.pir = t1waybt(prop_in_range ~ temperature, data=sdm[sdm$season=="summer",], tr = 0.2, nboot = 9999)
mcs.pir = mcppb20(prop_in_range ~ temperature, data=sdm[sdm$season=="summer",], tr = 0.2, nboot = 9999)
write.csv(mcs.pir$comp, "mcs_pir.csv")
t1w.pir = t1waybt(prop_in_range ~ temperature, data=sdm[sdm$season=="winter",], tr = 0.2, nboot = 9999)
mcw.pir = mcppb20(prop_in_range ~ temperature, data=sdm[sdm$season=="winter",], tr = 0.2, nboot = 9999)
write.csv(mcw.pir$comp, "mcw_pir.csv")
# compile
col.pir = rbind(t1s.pir[1:2], t1w.pir[1:2])

## Does AUC/TSS differ by model type for summer/winter?
# (one way anovas plus post-hoc tests by season)
t1s.auc = t1waybt(AUC ~ temperature, data=sdm[sdm$season=="summer",], tr = 0.2, nboot = 9999)
mcs.auc = mcppb20(AUC ~ temperature, data=sdm[sdm$season=="summer",], tr = 0.2, nboot = 9999)
write.csv(mcs.auc$comp, "mcs_AUC.csv")
t1s.tss = t1waybt(TSS ~ temperature, data=sdm[sdm$season=="summer",], tr = 0.2, nboot = 9999)
mcs.tss = mcppb20(TSS ~ temperature, data=sdm[sdm$season=="summer",], tr = 0.2, nboot = 9999)
write.csv(mcs.tss$comp, "mcs_TSS.csv")
t1w.auc = t1waybt(AUC ~ temperature, data=sdm[sdm$season=="winter",], tr = 0.2, nboot = 9999)
mcw.auc = mcppb20(AUC ~ temperature, data=sdm[sdm$season=="winter",], tr = 0.2, nboot = 9999)
write.csv(mcw.auc$comp, "mcw_AUC.csv")
t1w.tss = t1waybt(TSS ~ temperature, data=sdm[sdm$season=="winter",], tr = 0.2, nboot = 9999)
mcw.tss = mcppb20(TSS ~ temperature, data=sdm[sdm$season=="winter",], tr = 0.2, nboot = 9999)
write.csv(mcw.tss$comp, "mcw_TSS.csv")
# compile
col.auc = rbind(t1s.auc[1:2], t1w.auc[1:2])
col.tss = rbind(t1s.tss[1:2], t1w.tss[1:2])

## Does contribution of variable class differ by model type for summer/winter?
# (one way anovas plus post-hoc tests by season)
t1s.clim = t1waybt(clim ~ temperature, data=sdm[sdm$season=="summer",], tr = 0.2, nboot = 9999)
mcs.clim = mcppb20(clim ~ temperature, data=sdm[sdm$season=="summer",], tr = 0.2, nboot = 9999)
write.csv(mcs.clim$comp, "mcs_clim.csv")
t1s.top = t1waybt(top ~ temperature, data=sdm[sdm$season=="summer",], tr = 0.2, nboot = 9999)
mcs.top = mcppb20(top ~ temperature, data=sdm[sdm$season=="summer",], tr = 0.2, nboot = 9999)
write.csv(mcs.top$comp, "mcs_top.csv")
t1s.lc = t1waybt(lc ~ temperature, data=sdm[sdm$season=="summer",], tr = 0.2, nboot = 9999)
mcs.lc = mcppb20(lc ~ temperature, data=sdm[sdm$season=="summer",], tr = 0.2, nboot = 9999)
write.csv(mcs.lc$comp, "mcs_lc.csv")
t1w.clim = t1waybt(clim ~ temperature, data=sdm[sdm$season=="winter",], tr = 0.2, nboot = 9999)
mcw.clim = mcppb20(clim ~ temperature, data=sdm[sdm$season=="winter",], tr = 0.2, nboot = 9999)
write.csv(mcw.clim$comp, "mcw_clim.csv")
t1w.top = t1waybt(top ~ temperature, data=sdm[sdm$season=="winter",], tr = 0.2, nboot = 9999)
mcw.top = mcppb20(top ~ temperature, data=sdm[sdm$season=="winter",], tr = 0.2, nboot = 9999)
write.csv(mcw.top$comp, "mcw_top.csv")
t1w.lc = t1waybt(lc ~ temperature, data=sdm[sdm$season=="winter",], tr = 0.2, nboot = 9999)
mcw.lc = mcppb20(lc ~ temperature, data=sdm[sdm$season=="winter",], tr = 0.2, nboot = 9999)
write.csv(mcw.lc$comp, "mcw_lc.csv")
# compile
col.clim = rbind(t1s.clim[1:2], t1w.clim[1:2])
col.top = rbind(t1s.top[1:2], t1w.top[1:2])
col.lc = rbind(t1s.lc[1:2], t1w.lc[1:2])





