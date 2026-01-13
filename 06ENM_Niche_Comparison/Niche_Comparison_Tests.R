###########
##Niche Comparison Tests (package ECOSPAT)
library(dismo)
library(rgeos)
library(sf)
library(ade4)
library(ecospat)
library(terra)
library(factoextra)
library(gridExtra)
library(stringr)
set.seed(1234)

setwd("E://paper/07 costatae evolution/data_analysis/ENM 20250724/00Main")

all_var_list <- list.files("../environmental_data/",pattern=".asc$",full.names = T)
all_var <- raster::stack(all_var_list)

sp <- shapefile("E://Experiment Data/02ENMstudyArea/data/50m_physical_NaturalEarth/ne_50m_admin_0_countries.shp")
sp<- sp[sp$ADMIN=="China",]

#extract the occurence climate data 
occ_all_species <- read.delim("../occurrences/filter_occurences_all_species.txt",header = T)

coordinates(occ_all_species) <- ~ x + y
proj4string(occ_all_species) <- CRS(proj4string(all_var[[1]]))

occ_data <- data.frame(extract(all_var, occ_all_species))
occ_data <- cbind(occ_data,occ_all_species@data)

#Creating the background climate data for each species
occ_all_species <- read.delim("../occurrences/filter_occurences_all_species.txt",header = T)
colnames(occ_all_species) <- c("Species","Longitude", "Latitude")
occ.albo <- occ_all_species[which(occ_all_species[,1] == "albosinensis"),]
occ.ash  <- occ_all_species[which(occ_all_species[,1] == "ash"),]
occ.cos  <- occ_all_species[which(occ_all_species[,1] == "costata"),]
occ.erma <- occ_all_species[which(occ_all_species[,1] == "ermanii"),]
occ.uti  <- occ_all_species[which(occ_all_species[,1] == "utilis"),]
occ.bugg <- occ_all_species[which(occ_all_species[,1] == "buggsii"),]


generate.background <- function(occ,sp,all_var,name,bg_points=1000,species) {

    coordinates(occ) <- ~ Longitude + Latitude
    proj4string(occ) <- CRS(proj4string(all_var[[1]]))
	
    # using binary maps to generate background
    r <- rast(paste("../Result/Ensemble/W_MEAN/MAX_TSS/",species,".tif", sep = ""))
    r[r == 0] <- NA
    
    # detect the biggest patches (groups of cells that are surrounded by cells that are NA)
    patched <- patches(r, directions = 8)
    polygons <- as.polygons(patched)
    polygons$area <- expanse(polygons)
    max_polygon <- polygons[which.max(polygons$area), ]
    
	# make a buffer to merge adjacent small polygons
    buffer_distance <- 800000
    buffered_largest <- st_buffer(st_as_sf(max_polygon), dist = buffer_distance)
    polys_sf <- st_as_sf(polygons)

    intersecting_indices <- st_intersects(buffered_largest, polys_sf)
    intersecting_shapes <- polys_sf[intersecting_indices[[1]], ]
    # merge those small buffer
    merged_shape_buffer <- sf:::as_Spatial(st_union(intersecting_shapes))
	
    #generate background points
    samp <- spsample(merged_shape_buffer, bg_points, type='random', iter=25)
    #thinning points
    cells <- cellFromXY(all_var[[1]],samp)
    cells <- unique(cells)
    
    background_xy <- xyFromCell(all_var[[1]], cells)
    background_xy <- stats::na.omit(as.data.frame(background_xy))
    colnames(background_xy) <- c("Longitude", "Latitude")
    
    #data extraction of background points and modification of the dataset
    coordinates(background_xy) <- ~ Longitude + Latitude
    proj4string(background_xy) <- CRS(proj4string(all_var[[1]]))
    
	# check
    plot(sp)
	plot(r, col = "gray30", legend = FALSE,add=T)
    plot(background_xy,col = "black",add=T)
    plot(occ,add=T,col="red")
    plot(merged_shape_buffer, border = "blue", lwd = 0.5,add=T)
    
    background_data <- stats::na.omit(data.frame(extract(all_var, background_xy)))
    background_data <- cbind(background_data,rep(name,time=length(background_data[,1])))
    
    colnames(background_data) <- colnames(occ_data)
    return(background_data)
}

#More points were generated for the larger polygons, whereas 500 points were generated for smaller regions (buggsii) 
#to ensure that the points accurately represented the background range without contributing 
#overlapping points from the polygon.

bg.albo  <- generate.background(occ=occ.albo,sp=sp,all_var=all_var,
                                bg_points=1000,name="bg_alb",species="albosinensis")
bg.ash   <- generate.background(occ=occ.ash ,sp=sp,all_var=all_var,
                                bg_points=500,name="bg_ash",species="ash")
bg.cos   <- generate.background(occ=occ.cos ,sp=sp,all_var=all_var,
                                bg_points=1000,name="bg_cos",species="costata")
bg.erma  <- generate.background(occ=occ.erma,sp=sp,all_var=all_var,
                                bg_points=1000,name="bg_erm",species="ermanii")
bg.uti   <- generate.background(occ=occ.uti ,sp=sp,all_var=all_var,
                                bg_points=1000,name="bg_uti",species="utilis")
bg.bugg  <- generate.background(occ=occ.bugg,sp=sp,all_var=all_var,
                                bg_points=300 ,name="bg_bug",species="buggsii")


####PCA-env analysis####
###Dataframe for PCA with presence and background points
PCA_input <- rbind(occ_data, bg.albo, bg.ash, bg.cos, bg.erma, bg.uti, bg.bugg)
write.table(PCA_input, 
			file="NicheComparisonTests_occ_back_input.txt",
            quote=FALSE, sep="\t", na="",
            row.names=FALSE, col.names=TRUE)

#Calibrating the PCA in the whole studay area (PCAenv)
pcaenv <- dudi.pca(PCA_input[,1:11], center = TRUE, scale = TRUE, scannf = FALSE, nf = 3)
#Only occurrence data
pca.clim <- prcomp(occ_data[,1:11], center = TRUE, scale = TRUE)

explained_variance <- pcaenv$eig / sum(pcaenv$eig) * 100
pca.clim$sdev^2 / sum(pca.clim$sdev^2) *100
col <- as.factor(PCA_input$species)
col <- as.factor(str_replace_all(PCA_input$species, c("bg_alb" = "albosinensis", "bg_ash" = "ash", "bg_cos" = "costata", "bg_erm" = "ermanii", "bg_uti" = "utilis", "bg_bug" = "buggsii")))

#plot PCA for occurrence data
scree.pc12 <- fviz_pca_biplot(pcaenv,axes = c(1,2),title = "(A) PC1 vs. PC2", 
                # Individuals
                geom = "point",repel = FALSE,mean.point=FALSE,
                col.ind = col,
                alpha.ind = 0.5,
                pointshape = 16, pointsize = 2,labelsize=4,
                palette = c("#E78AC3", "#FC8D62","#d8ba8e", "#8DA0CB", "#FFD92F", "#A6D854"),
                # Variables
                col.var = "black",arrowsize=0.3,
                legend.title = list(fill = "Species")) + theme_bw() + theme(panel.grid.major=element_blank(), panel.grid.minor = element_blank()) + labs(fill="Species") +
    theme(plot.title = element_text(color="black", size=30, face="bold"),axis.title = element_text(size = 30),
          axis.text = element_text(size = 20),axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"))+
    xlab("PC1 (34.7%)") + ylab("PC2 (26.2%)")
  
scree.pc23<-fviz_pca_biplot(pcaenv,axes = c(2,3),title = "(B) PC2 vs. PC3", 
                            # Individuals
                            geom = "point",repel = FALSE,mean.point=FALSE,
                            col.ind = col,
                            alpha.ind = 0.5,
                            pointshape = 16, pointsize = 2,labelsize=4,
                            palette = c("#E78AC3", "#FC8D62","#d8ba8e", "#8DA0CB", "#FFD92F", "#A6D854"),
                            # Variables
                            col.var = "black",arrowsize = 0.3,
                            legend.title = list(fill = "Species")) + theme_bw() + theme(panel.grid.major=element_blank(), panel.grid.minor = element_blank()) + labs(fill="Species") +
  theme(plot.title = element_text(color="black", size=30, face="bold"),axis.title = element_text(size = 30),
        axis.text = element_text(size = 20),axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"))+
  xlab("PC2 (26.2%)") + ylab("PC3 (14.6%)")
 
pdf("fviz_pca_climvars_12_23.pdf", width = 8, height = 12, bg = "white")
grid.arrange(scree.pc12, scree.pc23, ncol = 1)
dev.off()

scree.pc13<-fviz_pca_biplot(pcaenv,axes = c(1,3),title = "(A) PC1 vs. PC3", 
                            # Individuals
                            col.ind = col,
                            alpha.ind = 0.5,
                            pointshape = 16, pointsize = 2,labelsize=4,
                            palette = c("#E78AC3", "#FC8D62","#d8ba8e", "#8DA0CB", "#FFD92F", "#A6D854"),
                            # Variables
                            col.var = "black",arrowsize = 0.3,
                            legend.title = list(fill = "Species")) + theme_bw() + theme(panel.grid.major=element_blank(), panel.grid.minor = element_blank()) + labs(fill="Species") +
  theme(plot.title = element_text(color="black", size=30, face="bold"),axis.title = element_text(size = 30),
        axis.text = element_text(size = 20),axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"))+
  xlab("PC1 (32.2%)") + ylab("PC3 (14.6%)")
#save plots
pdf("fviz_pca_climvars_13.pdf", width = 8.5, height = 6, bg = "white")
print(scree.pc13)
dev.off()
 
#analyses of variance of PC scores withing and among species
#boxplot
vardat <- data.frame(cbind(pcaenv$li[,1], pcaenv$li[,2], pcaenv$li[,3]))
#vardat <- cbind(PCA_input$species,vardat)
vardat <- cbind(str_replace_all(PCA_input$species, c("bg_alb" = "albosinensis", "bg_ash" = "ash", "bg_cos" = "costata", "bg_erm" = "ermanii", "bg_uti" = "utilis", "bg_bug" = "buggsii")),vardat)
colnames(vardat) <- c("species","PC1","PC2","PC3")
spnames<-c(expression(italic("B. ashburneri")),expression(italic("B. buggsii")),expression(italic("B. costata")),expression(italic("B. utilis")),expression(italic("B. albosinensis")),expression(italic("B. ermanii")))

vardat$species <- factor(vardat$species, levels=c("ash","buggsii","costata","utilis","albosinensis","ermanii"))
vardat <- vardat[order(vardat$species),]

pdf("pcscores_perspecies_boxplots2.pdf", width=10, height=12, bg="white")
par(mfrow=c(3,1))
par(mar=c(1,10,10,2))
boxplot(vardat$PC1 ~ vardat$species, col=c("#FC8D62","#d8ba8e","#8DA0CB","#A6D854","#E78AC3","#FFD92F"), names=spnames, las=2, xlab=" ", ylab="", main ="", cex.axis=1, xaxt='n')
title(main = "(C) PC1", cex.main=1, adj=0)
mtext(text = "PC1 scores",
      side = 2,#side 1 = bottom
      line = 5, cex = 1)

par(mar=c(1,10,10,2))
boxplot(vardat$PC2 ~ vardat$species, col=c("#FC8D62","#d8ba8e","#8DA0CB","#A6D854","#E78AC3","#FFD92F"), names=spnames, las=2, xlab=" ", ylab="", main ="",cex.axis=1, xaxt='n')
title(main = "(D) PC2", cex.main=1, adj=0)
mtext(text = "PC2 scores",
      side = 2,#side 1 = bottom
      line = 5, cex = 1)

par(mar=c(1,10,10,2))
boxplot(vardat$PC3 ~ vardat$species, col=c("#FC8D62","#d8ba8e","#8DA0CB","#A6D854","#E78AC3","#FFD92F"), names=spnames, las=2, ylab="", xlab="", main ="",cex.axis=1)
title(main = "(E) PC3", cex.main=1, adj=0)
mtext(text = "PC3 scores",
      side = 2,#side 1 = bottom
      line = 5, cex = 1)

dev.off()

#anova
library(multcompView)

anova <- aov(vardat$PC1 ~ vardat$species)
tukey <- TukeyHSD(anova)

multcompLetters4(anova, tukey)

anova <- aov(vardat$PC2 ~ vardat$species)
tukey <- TukeyHSD(anova)

multcompLetters4(anova, tukey)

anova <- aov(vardat$PC3 ~ vardat$species)
tukey <- TukeyHSD(anova)

multcompLetters4(anova, tukey)

#range of median PC scores values across the six species
library(dplyr)
vardat %>%
  group_by(species) %>%
  summarise(median_pc1 = median(PC1)) %>%
  pull(median_pc1) %>%
  range()

vardat %>%
  group_by(species) %>%
  summarise(median_pc2 = median(PC2)) %>%
  pull(median_pc2) %>%
  range()

vardat %>%
  group_by(species) %>%
  summarise(median_pc3 = median(PC3)) %>%
  pull(median_pc3) %>%
  range()

#Calibrating the PCA in the whole studay area (PCAenv)
pcaenv <- dudi.pca(PCA_input[,1:11], center = TRUE, scale = TRUE, scannf = FALSE, nf = 2)

row.albo <- which(PCA_input[,12] == "albosinensis")
row.ash  <- which(PCA_input[,12] == "ash")
row.cos  <- which(PCA_input[,12] == "costata")
row.erma <- which(PCA_input[,12] == "ermanii")
row.uti  <- which(PCA_input[,12] == "utilis")
row.bugg <- which(PCA_input[,12] == "buggsii")

row.bg.albo <- which(PCA_input[,12] == "bg_alb")
row.bg.ash  <- which(PCA_input[,12] == "bg_ash")
row.bg.cos  <- which(PCA_input[,12] == "bg_cos")
row.bg.erma <- which(PCA_input[,12] == "bg_erm")
row.bg.uti  <- which(PCA_input[,12] == "bg_uti")
row.bg.bugg <- which(PCA_input[,12] == "bg_bug")


# PCA scores for the whole study area
scores.globclim <- pcaenv$li
# PCA scores for the species distribution
scores.albo <- scores.globclim[row.albo,]
scores.ash  <- scores.globclim[row.ash,]
scores.cos  <- scores.globclim[row.cos,]
scores.erma <- scores.globclim[row.erma,]
scores.uti  <- scores.globclim[row.uti,]
scores.bugg <- scores.globclim[row.bugg,]

#PCA scores for the whole study area of each species
scores.bg_albo <- rbind(scores.albo, scores.globclim[row.bg.albo,])
scores.bg_ash  <- rbind(scores.ash,  scores.globclim[row.bg.ash,])
scores.bg_cos  <- rbind(scores.cos,  scores.globclim[row.bg.cos,])
scores.bg_erma <- rbind(scores.erma, scores.globclim[row.bg.erma,])
scores.bg_uti  <- rbind(scores.uti,  scores.globclim[row.bg.uti,])
scores.bg_bugg <- rbind(scores.bugg, scores.globclim[row.bg.bugg,])

zalbo <- ecospat.grid.clim.dyn(glob=scores.globclim, glob1=scores.bg_albo, scores.albo, R=100)
zash  <- ecospat.grid.clim.dyn(glob=scores.globclim, glob1=scores.bg_ash , scores.ash, R=100)
zcos  <- ecospat.grid.clim.dyn(glob=scores.globclim, glob1=scores.bg_cos , scores.cos, R=100)
zerma <- ecospat.grid.clim.dyn(glob=scores.globclim, glob1=scores.bg_erma, scores.erma, R=100)
zuti  <- ecospat.grid.clim.dyn(glob=scores.globclim, glob1=scores.bg_uti , scores.uti, R=100)
zbugg <- ecospat.grid.clim.dyn(glob=scores.globclim, glob1=scores.bg_bugg, scores.bugg, R=100)

##############
## plot diploid parent and tetraploid niche dynamics
par(mfrow = c(1, 3)) 
par(mar = c(4, 4, 2, 1))

ecospat.plot.niche.dyn(zash, zalbo, quant=0.05, interest=FALSE,title= "Niche Overlap for ash and albo", 
                       name.axis1="PC1",name.axis2="PC2", col.abn= "#87D180", col.unf = "#02734A",
                       col.stab = "#2C69B0", col.pio = "#F4737A", col.exp = "#C52E19",transparenc=40)
ecospat.shift.centroids(scores.ash,scores.albo,scores.bg_ash,scores.bg_albo,col = "gray15")

ecospat.plot.niche.dyn(zash, zuti, quant=0.05, interest=FALSE,title= "Niche Overlap for ash and uti", 
                       name.axis1="PC1",name.axis2="PC2", col.abn= "#87D180", col.unf = "#02734A",
                       col.stab = "#2C69B0", col.pio = "#F4737A", col.exp = "#C52E19",transparenc=40)
ecospat.shift.centroids(scores.ash,scores.uti,scores.bg_ash,scores.bg_uti,col = "gray15")

ecospat.plot.niche.dyn(zcos, zerma, quant=0.05, interest=FALSE,title= "Niche Overlap for cos and erman", 
                       name.axis1="PC1",name.axis2="PC2", col.abn= "#87D180", col.unf = "#02734A",
                       col.stab = "#2C69B0", col.pio = "#F4737A", col.exp = "#C52E19",transparenc=40)
ecospat.shift.centroids(scores.cos,scores.erma,scores.bg_cos,scores.bg_erma,col = "gray15")


##############
## Niche Pioneering + Niche Abandonment

calculate.five.niche.dyn.index <- function(scores.globclim ,scores.clim.nat,scores.sp.nat,scores.clim.inv,scores.sp.inv) {
  # https://github.com/ecospat/ecospat/issues/65
  # gridding the native niche
  grid.clim.nat <- ecospat.grid.clim.dyn(glob=scores.globclim, 
                                         glob1=scores.clim.nat,
                                         sp=scores.sp.nat, R=100,
                                         th.sp=0)
  
  # gridding the invasive niche
  grid.clim.inv <- ecospat.grid.clim.dyn(glob=scores.globclim,
                                         glob1=scores.clim.inv,
                                         sp=scores.sp.inv, R=100,
                                         th.sp=0) 
  
  niche.dyn <- ecospat.niche.dyn.index(grid.clim.nat, grid.clim.inv, intersection = 0)
  
  
  # Code to derive abandonment, pionneering, stability, expansion, unfilling for the pooled niche
  
  #Grid the pooled distribution
  grid.clim.pooled<- ecospat.grid.clim.dyn(glob=scores.globclim,
                                           glob1=scores.globclim,
                                           sp=rbind(scores.sp.nat,scores.sp.inv), 
                                           R=100,
                                           th.sp=0) 
  
  # extract rasters per category
  abandonment_r<- (((grid.clim.nat$Z > 0) - (grid.clim.inv$Z > 0)) == 1) * # crop to environment only in the native range
    (grid.clim.nat$z.uncor>0) * # crop to native niche only
    grid.clim.pooled$z.uncor# retrieve niche density
  pioneering_r <- (((grid.clim.inv$Z > 0) - (grid.clim.nat$Z > 0)) == 1) * # environment only in the invaded range
    (grid.clim.inv$z.uncor >0)* # crop niche in the invaded range
    grid.clim.pooled$z.uncor 
  stability_r <- (((grid.clim.inv$Z > 0) * (grid.clim.nat$Z > 0)) == 1) * # environment in both native and invaded range
    ((grid.clim.inv$z.uncor * grid.clim.nat$z.uncor)>0)* # crop niche in the native and invaded range
    grid.clim.pooled$z.uncor
  expansion_r <- (((grid.clim.inv$Z > 0) * (grid.clim.nat$Z > 0)) == 1) * # environment in both native and invaded range
    ((grid.clim.inv$z.uncor >0) * (grid.clim.nat$z.uncor==0))* # crop niche in the invaded range
    grid.clim.pooled$z.uncor 
  unfilling_r<-(((grid.clim.inv$Z > 0) * (grid.clim.nat$Z > 0)) == 1) * # environment in both native and invaded range
    ((grid.clim.nat$z.uncor >0) * (grid.clim.inv$z.uncor==0))* # crop niche in the native range
    grid.clim.pooled$z.uncor 
  
  # compute proportions
  abandonment <- cellStats(raster(abandonment_r),'sum') / cellStats(raster(grid.clim.pooled$z.uncor),'sum')
  pioneering <- cellStats(raster(pioneering_r),'sum') / cellStats(raster(grid.clim.pooled$z.uncor),'sum')
  stability <- cellStats(raster(stability_r),'sum') / cellStats(raster(grid.clim.pooled$z.uncor),'sum')
  expansion <- cellStats(raster(expansion_r),'sum') / cellStats(raster(grid.clim.pooled$z.uncor),'sum')
  unfilling <- cellStats(raster(unfilling_r),'sum') / cellStats(raster(grid.clim.pooled$z.uncor),'sum')
  
  return(c(abandonment, pioneering, stability,expansion, unfilling))
}

calculate.five.niche.dyn.index(scores.globclim=scores.globclim, scores.clim.nat=scores.bg_ash, 
                               scores.sp.nat=scores.ash,scores.clim.inv=scores.bg_albo,
                               scores.sp.inv=scores.albo)

calculate.five.niche.dyn.index(scores.globclim=scores.globclim, scores.clim.nat=scores.bg_ash, 
                               scores.sp.nat=scores.ash,scores.clim.inv=scores.bg_uti,
                               scores.sp.inv=scores.uti)

calculate.five.niche.dyn.index(scores.globclim=scores.globclim, scores.clim.nat=scores.bg_cos, 
                               scores.sp.nat=scores.cos,scores.clim.inv=scores.bg_erma,
                               scores.sp.inv=scores.erma)

par(mfrow = c(2, 3)) 
par(mar = c(4, 4, 2, 1))

ecospat.plot.niche(zash,  title = expression(italic("B. ashburneri")), name.axis1 = "PC1", name.axis2 = "PC2", cor=F)
ecospat.plot.niche(zbugg, title = expression(italic("B. buggsii")), name.axis1 = "PC1", name.axis2 = "PC2", cor=F)
ecospat.plot.niche(zcos,  title = expression(italic("B. costata")), name.axis1 = "PC1", name.axis2 = "PC2", cor=F)
ecospat.plot.niche(zuti,  title = expression(italic("B. utilis")), name.axis1 = "PC1", name.axis2 = "PC2", cor=F)
ecospat.plot.niche(zalbo, title = expression(italic("B. albosinensis")), name.axis1 = "PC1", name.axis2 = "PC2", cor=F)
ecospat.plot.niche(zerma, title = expression(italic("B. ermanii")), name.axis1 = "PC1", name.axis2 = "PC2", cor=F)


###Overlap tests
ecospat.niche.overlap(zalbo, zash, cor=TRUE)$D
ecospat.niche.overlap(zalbo, zcos, cor=TRUE)$D
ecospat.niche.overlap(zalbo, zerma,cor=TRUE)$D
ecospat.niche.overlap(zalbo, zuti, cor=TRUE)$D
ecospat.niche.overlap(zalbo, zbugg,cor=TRUE)$D

ecospat.niche.overlap(zash, zcos,  cor=TRUE)$D
ecospat.niche.overlap(zash, zerma, cor=TRUE)$D
ecospat.niche.overlap(zash, zuti,  cor=TRUE)$D
ecospat.niche.overlap(zash, zbugg, cor=TRUE)$D

ecospat.niche.overlap(zcos, zerma, cor=TRUE)$D
ecospat.niche.overlap(zcos, zuti,  cor=TRUE)$D
ecospat.niche.overlap(zcos, zbugg, cor=TRUE)$D

ecospat.niche.overlap(zerma, zuti, cor=TRUE)$D
ecospat.niche.overlap(zerma, zbugg,cor=TRUE)$D

ecospat.niche.overlap(zuti, zbugg, cor=TRUE)$D

###Similarity tests
ecospat.niche.similarity.test(zalbo,zash, 1000, overlap.alternative = "higher",ncores= 4)$p.D
ecospat.niche.similarity.test(zash,zalbo, 1000, overlap.alternative = "higher",ncores= 4)$p.D
ecospat.niche.similarity.test(zalbo,zcos, 1000, overlap.alternative = "higher",ncores= 4)$p.D
ecospat.niche.similarity.test(zcos,zalbo, 1000, overlap.alternative = "higher",ncores= 4)$p.D
ecospat.niche.similarity.test(zalbo,zerma, 1000, overlap.alternative = "higher",ncores= 4)$p.D
ecospat.niche.similarity.test(zerma,zalbo, 1000, overlap.alternative = "higher",ncores= 4)$p.D
ecospat.niche.similarity.test(zalbo,zuti, 1000, overlap.alternative = "higher",ncores= 4)$p.D
ecospat.niche.similarity.test(zuti,zalbo, 1000, overlap.alternative = "higher",ncores= 4)$p.D
ecospat.niche.similarity.test(zalbo,zbugg, 1000, overlap.alternative = "higher",ncores= 4)$p.D
ecospat.niche.similarity.test(zbugg,zalbo, 1000, overlap.alternative = "higher",ncores= 4)$p.D

ecospat.niche.similarity.test(zash, zcos, 1000, overlap.alternative = "higher",ncores= 4)$p.D
ecospat.niche.similarity.test(zcos, zash, 1000, overlap.alternative = "higher",ncores= 4)$p.D
ecospat.niche.similarity.test(zash, zerma, 1000, overlap.alternative = "higher",ncores= 4)$p.D
ecospat.niche.similarity.test(zerma, zash, 1000, overlap.alternative = "higher",ncores= 4)$p.D
ecospat.niche.similarity.test(zash, zuti, 1000, overlap.alternative = "higher",ncores= 4)$p.D
ecospat.niche.similarity.test(zuti, zash, 1000, overlap.alternative = "higher",ncores= 4)$p.D
ecospat.niche.similarity.test(zash, zbugg, 1000, overlap.alternative = "higher",ncores= 4)$p.D
ecospat.niche.similarity.test(zbugg, zash, 1000, overlap.alternative = "higher",ncores= 4)$p.D

ecospat.niche.similarity.test(zcos, zerma, 1000, overlap.alternative = "higher",ncores= 4)$p.D
ecospat.niche.similarity.test(zerma, zcos, 1000, overlap.alternative = "higher",ncores= 4)$p.D
ecospat.niche.similarity.test(zcos, zuti, 1000, overlap.alternative = "higher",ncores= 4)$p.D
ecospat.niche.similarity.test(zuti, zcos, 1000, overlap.alternative = "higher",ncores= 4)$p.D
ecospat.niche.similarity.test(zcos, zbugg, 1000, overlap.alternative = "higher",ncores= 4)$p.D
ecospat.niche.similarity.test(zbugg, zcos, 1000, overlap.alternative = "higher",ncores= 4)$p.D

ecospat.niche.similarity.test(zerma, zuti, 1000, overlap.alternative = "higher",ncores= 4)$p.D
ecospat.niche.similarity.test(zuti, zerma, 1000, overlap.alternative = "higher",ncores= 4)$p.D
ecospat.niche.similarity.test(zerma, zbugg, 1000, overlap.alternative = "higher",ncores= 4)$p.D
ecospat.niche.similarity.test(zbugg, zerma, 1000, overlap.alternative = "higher",ncores= 4)$p.D

ecospat.niche.similarity.test(zuti, zbugg, 1000, overlap.alternative = "higher",ncores= 4)$p.D
ecospat.niche.similarity.test(zbugg, zuti, 1000, overlap.alternative = "higher",ncores= 4)$p.D

###Equivalency tests,Divergence
ecospat.niche.equivalency.test(zalbo,zash ,rep=1000, overlap.alternative = "higher",ncores= 4)$p.D
ecospat.niche.equivalency.test(zalbo,zcos ,rep=1000, overlap.alternative = "higher",ncores= 4)$p.D
ecospat.niche.equivalency.test(zalbo,zerma,rep=1000, overlap.alternative = "higher",ncores= 4)$p.D
ecospat.niche.equivalency.test(zalbo,zuti ,rep=1000, overlap.alternative = "higher",ncores= 4)$p.D
ecospat.niche.equivalency.test(zalbo,zbugg,rep=1000, overlap.alternative = "higher",ncores= 4)$p.D

ecospat.niche.equivalency.test(zash, zcos ,rep=1000, overlap.alternative = "higher",ncores= 4)$p.D
ecospat.niche.equivalency.test(zash, zerma,rep=1000, overlap.alternative = "higher",ncores= 4)$p.D
ecospat.niche.equivalency.test(zash, zuti ,rep=1000, overlap.alternative = "higher",ncores= 4)$p.D
ecospat.niche.equivalency.test(zash, zbugg,rep=1000, overlap.alternative = "higher",ncores= 4)$p.D

ecospat.niche.equivalency.test(zcos, zerma,rep=1000, overlap.alternative = "higher",ncores= 4)$p.D
ecospat.niche.equivalency.test(zcos, zuti ,rep=1000, overlap.alternative = "higher",ncores= 4)$p.D
ecospat.niche.equivalency.test(zcos, zbugg,rep=1000, overlap.alternative = "higher",ncores= 4)$p.D

ecospat.niche.equivalency.test(zerma,zuti ,rep=1000, overlap.alternative = "higher",ncores= 4)$p.D
ecospat.niche.equivalency.test(zerma,zbugg,rep=1000, overlap.alternative = "higher",ncores= 4)$p.D

ecospat.niche.equivalency.test(zuti, zbugg,rep=1000, overlap.alternative = "higher",ncores= 4)$p.D


higher
 ecospat.plot.overlap.test(ecospat.niche.similarity.test(zcos, zerma, 1000, overlap.alternative = "higher",ncores= 4), "D", "Similarity")
 
####PCA-env analysis - Analyze climate and soil data separately####
###Dataframe for Worldclim PCA with presence and background points
PCA_input <- read.delim("./NicheComparisonTests_occ_back_PCA_clim_input.txt",header=T)
pcaenv <- dudi.pca(PCA_input[,1:5], center = TRUE, scale = TRUE, scannf = FALSE, nf = 2)

row.albo <- which(PCA_input[,6] == "albosinensis")
row.ash  <- which(PCA_input[,6] == "ash")
row.cos  <- which(PCA_input[,6] == "costata")
row.erma <- which(PCA_input[,6] == "ermanii")
row.uti  <- which(PCA_input[,6] == "utilis")
row.bugg <- which(PCA_input[,6] == "buggsii")

row.bg.albo <- which(PCA_input[,6] == "bg_alb")
row.bg.ash  <- which(PCA_input[,6] == "bg_ash")
row.bg.cos  <- which(PCA_input[,6] == "bg_cos")
row.bg.erma <- which(PCA_input[,6] == "bg_erm")
row.bg.uti  <- which(PCA_input[,6] == "bg_uti")
row.bg.bugg <- which(PCA_input[,6] == "bg_bug")


# PCA scores for the whole study area
scores.globclim <- pcaenv$li
# PCA scores for the species distribution
scores.albo <- pcaenv$li[row.albo,]
scores.ash  <- pcaenv$li[row.ash,]
scores.cos  <- pcaenv$li[row.cos,]
scores.erma <- pcaenv$li[row.erma,]
scores.uti  <- pcaenv$li[row.uti,]
scores.bugg <- pcaenv$li[row.bugg,]

#PCA scores for the whole study area of each species
scores.bg_albo <- rbind(scores.albo, pcaenv$li[row.bg.albo,])
scores.bg_ash  <- rbind(scores.ash,  pcaenv$li[row.bg.ash,])
scores.bg_cos  <- rbind(scores.cos,  pcaenv$li[row.bg.cos,])
scores.bg_erma <- rbind(scores.erma, pcaenv$li[row.bg.erma,])
scores.bg_uti  <- rbind(scores.uti,  pcaenv$li[row.bg.uti,])
scores.bg_bugg <- rbind(scores.bugg, pcaenv$li[row.bg.bugg,])

zalbo <- ecospat.grid.clim.dyn(glob=scores.globclim, glob1=scores.bg_albo, scores.albo, R=100)
zash  <- ecospat.grid.clim.dyn(glob=scores.globclim, glob1=scores.bg_ash , scores.ash, R=100)
zcos  <- ecospat.grid.clim.dyn(glob=scores.globclim, glob1=scores.bg_cos , scores.cos, R=100)
zerma <- ecospat.grid.clim.dyn(glob=scores.globclim, glob1=scores.bg_erma, scores.erma, R=100)
zuti  <- ecospat.grid.clim.dyn(glob=scores.globclim, glob1=scores.bg_uti , scores.uti, R=100)
zbugg <- ecospat.grid.clim.dyn(glob=scores.globclim, glob1=scores.bg_bugg, scores.bugg, R=100)

par(mfrow = c(2, 3)) 
par(mar = c(4, 4, 2, 1))

ecospat.plot.niche(zash,  title = expression(italic("B. ashburneri")), name.axis1 = "PC1", name.axis2 = "PC2", cor=TRUE)
ecospat.plot.niche(zbugg, title = expression(italic("B. buggsii")), name.axis1 = "PC1", name.axis2 = "PC2", cor=TRUE)
ecospat.plot.niche(zcos,  title = expression(italic("B. costata")), name.axis1 = "PC1", name.axis2 = "PC2", cor=TRUE)
ecospat.plot.niche(zuti,  title = expression(italic("B. utilis")), name.axis1 = "PC1", name.axis2 = "PC2", cor=TRUE)
ecospat.plot.niche(zalbo, title = expression(italic("B. albosinensis")), name.axis1 = "PC1", name.axis2 = "PC2", cor=TRUE)
ecospat.plot.niche(zerma, title = expression(italic("B. ermanii")), name.axis1 = "PC1", name.axis2 = "PC2", cor=TRUE)
###
###Dataframe for Soil PCA
###
PCA_input <- read.delim("./NicheComparisonTests_occ_back_PCA_soil_input.txt",header=T)
pcaenv <- dudi.pca(PCA_input[,1:6], center = TRUE, scale = TRUE, scannf = FALSE, nf = 2)

row.albo <- which(PCA_input[,7] == "albosinensis")
row.ash  <- which(PCA_input[,7] == "ash")
row.cos  <- which(PCA_input[,7] == "costata")
row.erma <- which(PCA_input[,7] == "ermanii")
row.uti  <- which(PCA_input[,7] == "utilis")
row.bugg <- which(PCA_input[,7] == "buggsii")

row.bg.albo <- which(PCA_input[,7] == "bg_alb")
row.bg.ash  <- which(PCA_input[,7] == "bg_ash")
row.bg.cos  <- which(PCA_input[,7] == "bg_cos")
row.bg.erma <- which(PCA_input[,7] == "bg_erm")
row.bg.uti  <- which(PCA_input[,7] == "bg_uti")
row.bg.bugg <- which(PCA_input[,7] == "bg_bug")


# PCA scores for the whole study area
scores.globclim <- pcaenv$li
# PCA scores for the species distribution
scores.albo <- pcaenv$li[row.albo,]
scores.ash  <- pcaenv$li[row.ash,]
scores.cos  <- pcaenv$li[row.cos,]
scores.erma <- pcaenv$li[row.erma,]
scores.uti  <- pcaenv$li[row.uti,]
scores.bugg <- pcaenv$li[row.bugg,]

#PCA scores for the whole study area of each species
scores.bg_albo <- rbind(scores.albo, pcaenv$li[row.bg.albo,])
scores.bg_ash  <- rbind(scores.ash,  pcaenv$li[row.bg.ash,])
scores.bg_cos  <- rbind(scores.cos,  pcaenv$li[row.bg.cos,])
scores.bg_erma <- rbind(scores.erma, pcaenv$li[row.bg.erma,])
scores.bg_uti  <- rbind(scores.uti,  pcaenv$li[row.bg.uti,])
scores.bg_bugg <- rbind(scores.bugg, pcaenv$li[row.bg.bugg,])

zalbo <- ecospat.grid.clim.dyn(glob=scores.globclim, glob1=scores.bg_albo, scores.albo, R=100)
zash  <- ecospat.grid.clim.dyn(glob=scores.globclim, glob1=scores.bg_ash , scores.ash, R=100)
zcos  <- ecospat.grid.clim.dyn(glob=scores.globclim, glob1=scores.bg_cos , scores.cos, R=100)
zerma <- ecospat.grid.clim.dyn(glob=scores.globclim, glob1=scores.bg_erma, scores.erma, R=100)
zuti  <- ecospat.grid.clim.dyn(glob=scores.globclim, glob1=scores.bg_uti , scores.uti, R=100)
zbugg <- ecospat.grid.clim.dyn(glob=scores.globclim, glob1=scores.bg_bugg, scores.bugg, R=100)

par(mfrow = c(2, 3)) 
par(mar = c(4, 4, 2, 1))

ecospat.plot.niche(zash,  title = expression(italic("B. ashburneri")), name.axis1 = "PC1", name.axis2 = "PC2", cor=TRUE)
ecospat.plot.niche(zbugg, title = expression(italic("B. buggsii")), name.axis1 = "PC1", name.axis2 = "PC2", cor=TRUE)
ecospat.plot.niche(zcos,  title = expression(italic("B. costata")), name.axis1 = "PC1", name.axis2 = "PC2", cor=TRUE)
ecospat.plot.niche(zuti,  title = expression(italic("B. utilis")), name.axis1 = "PC1", name.axis2 = "PC2", cor=TRUE)
ecospat.plot.niche(zalbo, title = expression(italic("B. albosinensis")), name.axis1 = "PC1", name.axis2 = "PC2", cor=TRUE)
ecospat.plot.niche(zerma, title = expression(italic("B. ermanii")), name.axis1 = "PC1", name.axis2 = "PC2", cor=TRUE)