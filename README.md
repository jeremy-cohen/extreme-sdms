## DATA FILES

### sdm_birds_v8.csv 
Species-level distribution model outputs. Each row is a species and columns contain model type information, including species modeled, modeling season, model type, model version, number of observations, estimated presence cells, model metrics (AUC, SPS threshold, TSS, sensitivity, specificity), and Random forest importance score for all model variables (any variable with _imp suffix). Cells containing NA have no associated model or path or information for a given species.

#### Column IDs
common_name: common name of species
sciname: scientific name of species
season: season at which model was fit
temperature: model type (means, variability, extremes)
version: model version
modName: model name
noPts: number of sampling points in model
range_area: estimated range size from model prediction
AUC: model area under the curve
SPSThresh: model threshold value
TSS: model true skill score
Sensitivity: model sensitivity
Specificity: model specificity
prop_in_range: Proportion of thresholded model prediction falling within expert range polygon
train.ac.mean: mean moran's I of training data across simulations
train.ac.sd: standard deviation of moran's I of training data across simulations
test.ac.mean: mean moran's I of testing data across simulations
test.ac.sd: standard deviation of moran's I of testing data across simulations
slope_pres_bio1_ehe: slope of relationship between bio1 (mean temperature) and extreme heat index among presence observations
slope_pres_bio1_ece: slope of relationship between bio1 (mean temperature) and extreme cold index among presence observations
slope_pres_bio12_spei: slope of relationship between bio12 (mean precipitation) and extreme drought index among presence observations
slope_abs_bio1_ehe: slope of relationship between bio1 (mean temperature) and extreme heat index among absence observations
slope_abs_bio1_ece: slope of relationship between bio1 (mean temperature) and extreme cold index among absence observations
slope_abs_bio12_spei: slope of relationship between bio12 (mean precipitation) and extreme drought index among absence observations
year_imp: importance score for year 
day_of_year_imp: importance score for julian date 
time_observations_started_imp: importance score for time observations started 
duration_minutes_imp: importance score for duration of survey 
effort_distance_km_imp: importance score for distance traveled 
number_observers_imp: importance score for number of observers 
bio1_imp: importance score for mean temperature in climate 
bio12_imp: importance score for total precipitation in climate          
tri_imp: importance score for topographic roughness index 
elev_imp: importance score for elevation 
pland_10_cropland_rainfed_imp: importance score for cropland (rainfed) landcover (type 1) 
pland_100_mosaic_tree_shrub_imp: importance score for mosaic tree/shrub landcover 
pland_11_cropland_rainfed_imp: importance score for cropland (rainfed) landcover (type 2) 
pland_110_mosaic_herbacious_imp: importance score for mosaic-herbaceous landcover 
pland_12_cropland_rainfed_imp: importance score for cropland (rainfed) landcover (type 3) 
pland_120_shrubland_imp: importance score for shrubland landcover (type 1) 
pland_121_shrubland_imp: importance score for shrubland landcover (type 2)
pland_122_shrubland_imp: importance score for shrubland landcover (type 3)
pland_130_grassland_imp: importance score for grassland landcover 
pland_140_lichens_mosses_imp: importance score for lichens/mosses landcover 
pland_150_sparse_imp: importance score for sparse landcover (1)
pland_152_sparse_imp: importance score for sparse landcover (2)
pland_153_sparse_imp: importance score for sparse landcover (3)
pland_160_flooded_freshwater_imp: importance score for flooded freshwater landcover 
pland_170_flooded_saltwater_imp: importance score for flooded saltwater landcover 
pland_180_flooded_shrub_imp: importance score for shrub landcover 
pland_190_urban_imp: importance score for urban landcover 
pland_20_cropland_irrigated_imp: importance score for irrigated cropland landcover 
pland_200_barren_imp: importance score for barren landcover (1)
pland_201_barren_imp: importance score for barren landcover (2)
pland_202_barren_imp: importance score for barren landcover (3)
pland_210_water_imp: importance score for water landcover 
pland_220_ice_imp: importance score for ice landcover 
pland_30_mosaic_cropland_imp: importance score for mosaic cropland landcover 
pland_40_mosaic_natural_veg_imp: importance score for natural vegetation landcover 
pland_50_evergreen_broadleaf_imp: importance score for evergreen broadleaf landcover (1)
pland_60_deciduous_broadleaf_imp: importance score for deciduous broadleaf landcover (1)
pland_61_deciduous_broadleaf_imp: importance score for deciduous broadleaf landcover (2)
pland_62_deciduous_broadleaf_imp: importance score for deciduous broadleaf landcover (3)
pland_70_evergreen_needleleaf_imp: importance score for evergreen needleleaf landcover (1)
pland_71_evergreen_needleleaf_imp: importance score for evergreen needleleaf landcover (2)
pland_72_evergreen_needleleaf_imp: importance score for evergreen needleleaf landcover (3)
pland_80_deciduous_needleleaf_imp: importance score for deciduous needleleaf landcover (1)
pland_81_deciduous_needleleaf_imp: importance score for deciduous needleleaf landcover (2)
pland_82_deciduous_needleleaf_imp: importance score for deciduous needleleaf landcover (3)
pland_90_mixed_forest_imp: importance score for mixed forest landcover
evi_imp: importance score for seasonal enhanced vegetation index 
bio4_imp: importance score for temperature seasonality
bio15_imp: importance score for precipitation seasonality
ehe_imp: importance score for extreme heat index
ece_imp: importance score for extreme cold index
spei_imp: importance score for drought index
clim: proportion of total importance score from climate variables
top: proportion of total importance score from topographic variables
lc: proportion of total importance score from landcover variables
Hand.Wing.Index: Hand-wing index value
Mass: body mass (g)
Habitat: habitat preference
Trophic.Level: trophic level
Trophic.Niche: trophic niche
Scientific.AVONET: AVONET species name
Scientific.BirdTree: BirdTree species name


### birds_VE8_comp.csv 
Differences between Variability/Extremes model and Means model by species/season.

#### Column IDs
common_name: common name of species
sciname: scientific name of species
season: season at which model was fit
temp: model type (means, variability, extremes)
delta_range_area: change in estimated range size from Means model
delta_prop_in_range: change in estimated proportion of range within expert range from Means model
delta_AUC: change in area under the curve from Means model
delta_SPSthreshold: Change in estimated threshold from Means model
delta_TSS: Change in TSS from Means model
delta_Sensitivity: Change in sensitivity from Means model
delta_Specificity: Change in specificity from Means model
delta_clim: change in proportion of total importance score from climate variables from Means model
delta_top: change in proportion of total importance score from topographic variables from Means model
delta_lc: change in proportion of total importance score from landcover variables from Means model


### species_list_2024.csv 
Index of all candidate species. Cells containing NA have no information about extinction for a given species.

#### Column IDs
Common: common name
Scientific: scientific name
Abbr: four letter species abbreviation
Code: species rarity code- 1 (common) or 2 (uncommon)
ebird_code: species six letter code
order: taxonomic order
Hawaii: Hawaiian species? 1-yes 0-no
seabird: Seabird? 1-yes 0-no

## Sharing/access info
Ebird data can be downloaded via www.ebird.org
AVONET trait data is available at https://onlinelibrary.wiley.com/doi/full/10.1111/ele.13898


## CODE FILES

1. data_update_ebird.R - Building and annotating the custom seasonal ebird datasets.
2. modeling_ranger_extremes_breeding.R - Species distribution modeling workflow including creation of species-level figures (summer).
3. modeling_ranger_extremes_nonbreeding.R - Species distribution modeling workflow including creation of species-level figures (winter).
4. fusing_outputs.R - combining all SDM outputs to create cross-species data tables.
5. aggregate_extremes.R - combining thresholded SDM predictions spatially to generate richness.
6. aggregate_extremes_ror.R - combining relative occurrence ratio SDM predictions spatially to generate richness.
7. validations_sp_accuracy.R - Site-level validations to compared observed vs. expected presences, absences, and richness.
8. one_way_anovas.R - one-way anova tests to explore cross-species trends across model types.
9. one_way_anovas.R - figures and tables to explore cross-species trends across model types.
10. trait_models.R - models with trait covariates accounting for phylogenetic structure to explore cross-species trends across model types.
