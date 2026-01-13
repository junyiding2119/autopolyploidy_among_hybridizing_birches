###############
####ENM
library(rgbif)
library(geosphere)
library(sdm)
library(raster)
library(usdm)
library(sf)
library(factoextra)
library(ggplot2)
library(gridExtra)

#######
#####Downloading GBIF data(accessed 16 April 2024)

sp <- shapefile("E://Experiment Data/02ENMstudyArea/data/50m_physical_NaturalEarth/ne_50m_admin_0_countries.shp")
sp <- sp[sp$ADMIN=="China",]

clim_lst <- list.files("E://Experiment Data/02ENMstudyArea/data/world_bioclim/Current_1981-2010_CHELSA",pattern=".tif$",full.names = T) 
clim <- raster::stack(clim_lst)

lst_soil <- list.files("E://Experiment Data/02ENMstudyArea/data/soil_world",pattern='.tif$',full.names = T)
soil <- raster::stack(lst_soil)

species <- "Betula costata Trautv."

occ_raw <- as.data.frame(occ_data(scientificName = species,
                                  hasCoordinate=TRUE, basisOfRecord = "PRESERVED_SPECIMEN",
						                      limit = 50000)[["data"]])
# remove erroneous coordinates, where either the latitude or longitude is missing
occ_clean <- subset(as.data.frame(occ_raw),
                   (!is.na(decimalLatitude))&(!is.na(decimalLongitude))&(year >= 1960)) 
cat(nrow(occ_raw)-nrow(occ_clean), "records are removed")
occ_clean <- occ_clean[,4:3]

# read additional populations from our records
add <- read.table("location_costata.txt", header = T)
occ_clean <- rbind(add, occ_clean)

# remove duplicated data based on latitude and longitude
dups <- duplicated(occ_clean[c("decimalLatitude","decimalLongitude")])
occ_unique <- occ_clean[!dups,]
cat(nrow(occ_clean)-nrow(occ_unique), "records are removed")

# thin occ data (keep one occurrence point per cell)
cells <- cellFromXY(clim[[1]],occ_unique)
dups <- duplicated(cells)
occ_single <- occ_unique[!dups,]
occ_single <- cbind(occ_single,rep("1",time=length(occ_single[,1])))
colnames(occ_single) <- c("Longitude", "Latitude", "occ")
coordinates(occ_single) <- ~ Longitude + Latitude
proj4string(occ_single) <- CRS(proj4string(clim[[1]]))
# remove aboard occ data 
over_result <- over(occ_single, sp)
in_china <- subset(occ_single, !is.na(over_result$scalerank))

### fine erroneous locations 
#notice that some locations are erroneous, such as some costata in the SW China
#we removed them manually
plot(sp)
plot(in_china,add=T)

write.table(in_china, file=paste("filter_occurences_",species,".txt", sep = ""),
            quote=FALSE, sep="\t", na="",
            row.names=FALSE, col.names=TRUE)
			
##we remove some utilis and albosinensis locations in the contact zone
## utilis individuals
utilis_polygon <- data.frame(c(77.479,75.443,102.067,106.681,103.213,99.188,77.479),
							 c(35.348,30.08,20.338,24.611,28.73,31.085,35.348))
Poly <- sp::Polygons(list(sp::Polygon(utilis_polygon)), ID = "A")
SpatialPoly <- sp::SpatialPolygons(list(Poly))

species <- "Betula utilis D.Don"

occ_raw <- as.data.frame(occ_data(scientificName = species,
                                  hasCoordinate=TRUE, basisOfRecord = "PRESERVED_SPECIMEN",
						                      limit = 50000)[["data"]])
# remove erroneous coordinates, where either the latitude or longitude is missing
occ_clean <- subset(as.data.frame(occ_raw),
                   (!is.na(decimalLatitude))&(!is.na(decimalLongitude))&(year >= 1960)) 
cat(nrow(occ_raw)-nrow(occ_clean), "records are removed")
occ_clean <- occ_clean[,4:3]

# read additional populations from our records
add <- read.table("location_utilis.txt", header = T)
occ_clean <- rbind(add, occ_clean)

# remove duplicated data based on latitude and longitude
dups <- duplicated(occ_clean[c("decimalLatitude","decimalLongitude")])
occ_unique <- occ_clean[!dups,]
cat(nrow(occ_clean)-nrow(occ_unique), "records are removed")

# thin occ data (keep one occurrence point per cell)
cells <- cellFromXY(clim[[1]],occ_unique)
dups <- duplicated(cells)
occ_single <- occ_unique[!dups,]
occ_single <- cbind(occ_single,rep("1",time=length(occ_single[,1])))
colnames(occ_single) <- c("Longitude", "Latitude", "occ")
coordinates(occ_single) <- ~ Longitude + Latitude
proj4string(occ_single) <- CRS(proj4string(clim[[1]]))

# remove aboard occ data 
over_result <- over(occ_single, sp)
in_china <- subset(occ_single, !is.na(over_result$scalerank))

tem <- as.data.frame(in_china@coords)
coordinates(tem) <- ~ Longitude + Latitude
over_result <- over(tem, SpatialPoly)

write.table(tem[which(over_result != "NA"),], 
			file=paste("filter_occurences_",species,".txt", sep = ""),
            quote=FALSE, sep="\t", na="",
            row.names=FALSE, col.names=TRUE)
			
## albosinensis individuals
albo_polygon <- data.frame(c(98.238,106.423,111.334,113.899,113.586,111.51,105.805,98.238),
						   c(36.323,30.211,30.949,34.733,38.993,39.901,40.799,36.323))
Poly <- sp::Polygons(list(sp::Polygon(albo_polygon)), ID = "A")
SpatialPoly <- sp::SpatialPolygons(list(Poly))

species <- "Betula albosinensis Burkill"

occ_raw <- as.data.frame(occ_data(scientificName = species,
                                  hasCoordinate=TRUE, basisOfRecord = "PRESERVED_SPECIMEN",
						          ,limit = 50000)[["data"]])
# remove erroneous coordinates, where either the latitude or longitude is missing
occ_clean <- subset(as.data.frame(occ_raw),
                   (!is.na(decimalLatitude))&(!is.na(decimalLongitude))&(year >= 1960)) 
cat(nrow(occ_raw)-nrow(occ_clean), "records are removed")
occ_clean <- occ_clean[,4:3]

# read additional populations from our records
add <- read.table("location_albosinensis.txt", header = T)
occ_clean <- rbind(add, occ_clean)

# remove duplicated data based on latitude and longitude
dups <- duplicated(occ_clean[c("decimalLatitude","decimalLongitude")])
occ_unique <- occ_clean[!dups,]
cat(nrow(occ_clean)-nrow(occ_unique), "records are removed")

# thin occ data (keep one occurrence point per cell)
cells <- cellFromXY(clim[[1]],occ_unique)
dups <- duplicated(cells)
occ_single <- occ_unique[!dups,]
occ_single <- cbind(occ_single,rep("1",time=length(occ_single[,1])))
colnames(occ_single) <- c("Longitude", "Latitude", "occ")
coordinates(occ_single) <- ~ Longitude + Latitude
proj4string(occ_single) <- CRS(proj4string(clim[[1]]))

tem <- as.data.frame(occ_single@coords)
coordinates(tem) <- ~ Longitude + Latitude
over_result <- over(tem, SpatialPoly)

write.table(tem[which(over_result != "NA"),], 
			file=paste("filter_occurences_",species,".txt", sep = ""),
            quote=FALSE, sep="\t", na="",
            row.names=FALSE, col.names=TRUE)

#######
##Environmental variables PCA 

# excluding the collinear variables
studyArea_clim <- crop(clim,extent(sp)) 
studyArea_clim <- mask(studyArea_clim,sp)

studyArea_soil <- crop(soil,extent(sp)) 
studyArea_soil <- mask(studyArea_soil,sp)

#studyArea_soil_resample <- resample(studyArea_soil,studyArea_clim)

climate_stack <- raster::stack(studyArea_clim,studyArea_soil)

# Define variable groups based on their characteristics
climate_annual_vars <- c("bio01","bio02", "bio03", "bio04", "bio07", "bio12", "bio15")
climate_monthly_vars <- c("bio05", "bio06", "bio13", "bio14")
climate_quarterly_vars <- c("bio08", "bio09", "bio10", "bio11", "bio16", "bio17", "bio18", "bio19")
soil_vars <- c("cfvo_5.15cm", "clay_5.15cm", "nitrogen_5.15cm", "ocd_5.15cm", "phh2o_5.15cm", "sand_5.15cm", "silt_5.15cm", "soc_5.15cm")

# Verify all required variables exist in the RasterStack
all_vars <- c(climate_annual_vars, climate_monthly_vars, climate_quarterly_vars, soil_vars)
var_names <- names(climate_stack)

if (!all(all_vars %in% var_names)) {
  missing_vars <- setdiff(all_vars, var_names)
  stop("Missing variables in RasterStack: ", paste(missing_vars, collapse = ", "))
}

# Sample data from RasterStack (50000 random pixels)
set.seed(1234)
sample_data <- sampleRandom(climate_stack, size = 100000, na.rm = TRUE, df = TRUE)
remaining_vars <- all_vars
current_data <- as.data.frame(sample_data[, remaining_vars, drop = FALSE])
threshold <- 10  # VIF threshold

# Stepwise VIF elimination with priority: 
#  1. Annual climate > 2. Monthly climate > 3. Soil > 4. Quarterly climate
#  Franklin J. 2010. Mapping Species Distributions: Spatial Inference and Prediction. Cambridge: Cambridge University Press.
repeat {
  vif_scores <- vif(current_data)
  
  # Exit loop if all VIF values are acceptable
  if (all(vif_scores$VIF <= threshold)) break
  
  high_vif <- vif_scores[vif_scores$VIF > threshold, ]
  if (nrow(high_vif) == 0) break  # Safety check
  
  # Priority 1: Remove annual climate variables
  high_vif_annual <- high_vif[high_vif$Variables %in% climate_annual_vars, ]
  if (nrow(high_vif_annual) > 0) {
    var_to_remove <- high_vif_annual$Variables[which.max(high_vif_annual$VIF)]
  } 
  else {
    # Priority 2: Remove monthly climate variables
    high_vif_monthly <- high_vif[high_vif$Variables %in% climate_monthly_vars, ]
    if (nrow(high_vif_monthly) > 0) {
      var_to_remove <- high_vif_monthly$Variables[which.max(high_vif_monthly$VIF)]
    } 
    else {
      # Priority 3: Remove soil variables
      high_vif_soil <- high_vif[high_vif$Variables %in% soil_vars, ]
      if (nrow(high_vif_soil) > 0) {
        var_to_remove <- high_vif_soil$Variables[which.max(high_vif_soil$VIF)]
      } 
      else {
        # Priority 4: Remove quarterly climate variables
        var_to_remove <- high_vif$Variables[which.max(high_vif$VIF)]
      }
    }
  }
  
  remaining_vars <- setdiff(remaining_vars, var_to_remove)
  current_data <- current_data[, remaining_vars, drop = FALSE]
  
  cat("Removed:", var_to_remove, " | Remaining variables:", length(remaining_vars), "\n")
  
  if (length(remaining_vars) <= 1) break
}

print(remaining_vars)

# Create final RasterStack with selected variables
selected_stack <- climate_stack[[remaining_vars]]
cat("Final retained variables:", paste(remaining_vars, collapse = ", "), "\n")

#Final retained variables: 
#bio03, bio08, bio15, bio18, bio19, 
#cfvo_5.15cm, clay_5.15cm, nitrogen_5.15cm, ocd_5.15cm, 
#phh2o_5.15cm, silt_5.15cm, soc_5.15cm 

writeRaster(studyArea_clim,
            # a series of names for output files
            filename=paste0("../environmental_data/",names(studyArea_clim),".tif"), 
            format="GTiff", ## the output format
            bylayer=TRUE, ## this will save a series of layers
            overwrite=T)

writeRaster(studyArea_soil,
            # a series of names for output files
            filename=paste0("../environmental_data/",names(studyArea_soil),".tif"), 
            format="GTiff", ## the output format
            bylayer=TRUE, ## this will save a series of layers
            overwrite=T)

######
##SDM - ENMTML

##### Clear all states
if(!is.null(dev.list())) dev.off()
par(mfrow=c(1,1))
rm(list=ls(all=TRUE))

library(raster) 
library(ENMTML) #version 1.0.0

ENMTML(
 pred_dir = "../environmental_data",
 proj_dir = NULL,
 result_dir = "../Result",
 occ_file = "../occurrences/filter_occurences_all_species.txt",
 sp = 'species',
 x = 'x',
 y = 'y',
 min_occ = 10,
 thin_occ = NULL,
 eval_occ = NULL,
 colin_var = NULL,
 imp_var = FALSE,
 sp_accessible_area = NULL,
 pseudoabs_method = c(method='GEO_ENV_KM_CONST', width='200'),
 pres_abs_ratio = 1,
 part=c(method='BOOT', replicates='10', proportion='0.7'),
 save_part = FALSE,
 save_final = TRUE,
 algorithm = c("BIO","DOM","ENF" ,"MXD","SVM","GLM","GAM","RDF","MLK","GAU"),
 thr = c(type='MAX_TSS'),
 msdm = NULL,
 ensemble=c(method='W_MEAN', metric='TSS'),
 extrapolation = FALSE,
 cores = 4
)

ENMTML(
 pred_dir = "../environmental_data",
 proj_dir = NULL,
 result_dir = "../Result",
 occ_file = "../occurrences/filter_occurences_all_species.txt",
 sp = 'species',
 x = 'x',
 y = 'y',
 min_occ = 10,
 thin_occ = NULL,
 eval_occ = NULL,
 colin_var = NULL,
 imp_var = FALSE,
 sp_accessible_area = NULL,
 pseudoabs_method = c(method='GEO_ENV_KM_CONST', width='200'),
 pres_abs_ratio = 1,
 part=c(method='BOOT', replicates='10', proportion='0.7'),
 save_part = FALSE,
 save_final = TRUE,
 algorithm = c("DOM", "GAU", "MXD", "RDF", "SVM"),
 thr = c(type='MAX_TSS'),
 msdm = NULL,
 ensemble=c(method='W_MEAN', metric='TSS'),
 extrapolation = FALSE,
 cores = 4
)


#plot
chinese_sp <- shapefile("E://Experiment Data/02ENMstudyArea/data/china_boundaries/provincial_boundaries/bou2_4p.shp") 
ens_lst <- list.files("E://paper/07 costatae evolution/data_analysis/ENM 20250724/Result/Ensemble/W_MEAN",pattern=".tif$",full.names = T)
ens <- raster::stack(ens_lst[2],ens_lst[3],ens_lst[4],ens_lst[6],ens_lst[1],ens_lst[5])
spplot(ens,cuts = 5,col.regions = rev(terrain.colors(7)),sp.layout = chinese_sp,names.attr=c('B. ashburneri', 'B. buggsii',"B. costata","B. utilis","B. albosinensis","B. ermanii"))

#binary maps
ens_binary_lst <- list.files("E://paper/07 costatae evolution/data_analysis/ENM 20250724/Result/Ensemble/W_MEAN/MAX_TSS",pattern=".tif$",full.names = T)

ens_binary <- raster::stack(ens_binary_lst[2],ens_binary_lst[3],ens_binary_lst[4],ens_binary_lst[6],ens_binary_lst[1],ens_binary_lst[5])

spplot(ens_binary,cuts = 1,col.regions = c("#F2F2F2", "#7992a6"),sp.layout = chinese_sp,names.attr=c('B. ashburneri', 'B. buggsii',"B. costata","B. utilis","B. albosinensis","B. ermanii"))
