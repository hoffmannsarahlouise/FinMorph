#'
#' Supplementary R script to accompany:
#' 
#' Comparative Morphology of Shark Pectoral Fins
#' 
#' Hoffmann, S.L. Buser, T.J., and Porter, M.E. 
#' 
#' Submitted to Journal of Morphology
#' 
#' Please direct any questions or comments regarding the contents of this script
#' to: Thaddaeus Buser (thaddaeus.buser@gmail.com)
#' 
# Assign working drive
setwd("//Analysis_in_R_200705")
#
# load libraries and custom functions
library(vegan)
library(phytools)
library(geomorph)
library(geiger)
library(cluster)
library(ape)
library(pvclust)
library(MASS)
library(nlme)
source("my.custom.functions.R") #Be sure keep the custom functions script in your working drive 
#
#############################################
#
# Read in primary data
fin.measures <- read.csv(file = "finmorph_all_07022020_trim.csv", header = T, row.names = 1)
#
shark.tree <- read.nexus(file = "mtDNA_tree_run3_MCC.tree")
plotTree(shark.tree)
#
#
fin.landmarks <- readland.tps(file = "Fin_landmarks.TPS", specID = "imageID",readcurves = T)
#'
####'
####'  SECTION I: PECTORAL FIN SHAPE
####'  
#' 
#' We gathered outline data for the pectoral fins as well but found that
#' the outlines were not necessary to capture the shape of the fins well,
#' so, to save ourselves an extraneous number of dependant variables, we'll
#' trip our dataset to just include the landmarks
fin.landmarks <- fin.landmarks[1:5,,]
#
landmark.colorz2 <-c("green", "red", "blue", "orange", "black") #for plotting purposes below
# Generalized Procrustes analysis of landmark data
sharks.gpa_test <- gpagen(A=fin.landmarks) 
# Plot the aligned landmark arrays color-coded as above
dev.off()
plotAllSpecimens(sharks.gpa_test$coords, label=T,
                 plot.param = list(pt.bg = landmark.colorz2, mean.cex=.5))
#
#
dev.off()
# Principal component analysis of aligned landmark data
shark.morphopace <- plotTangentSpace(sharks.gpa_test$coords, axis1 = 1, axis2 = 2)
#read in species classifier for each specimen
classifier_sharkTPS <- read.csv(file = "shark.landmark.classifier.csv")
# Average landmark positions for each species
species.aves.shark_TPS <-my.landmark.species.average(classifier = classifier_sharkTPS, aligned.coords = sharks.gpa_test$coords)
#plot again
dev.off()
plotAllSpecimens(species.aves.shark_TPS, label=T,
                 plot.param = list(pt.bg = landmark.colorz2, mean.cex=.5))
#General Procrustes alignment for spp. averages
sharks.sppave.gpa <- gpagen(A=species.aves.shark_TPS) 
#plot again
plotAllSpecimens(sharks.sppave.gpa$coords, label=T,
                 plot.param = list(pt.bg = landmark.colorz2, mean.cex=.5))
#pca
shark.sppave.morphopace.full <- plotTangentSpace(sharks.sppave.gpa$coords, axis1 = 1, axis2 = 2)
#                          PC1     PC2     PC3     PC4     PC5      PC6
#Standard deviation     0.1342 0.04971 0.04190 0.03316 0.02107 0.007983
#Proportion of Variance 0.7553 0.10366 0.07366 0.04612 0.01863 0.002670
#Cumulative Proportion  0.7553 0.85892 0.93258 0.97870 0.99733 1.000000
plotOutliers(sharks.sppave.gpa$coords)
#Major outlier: Squatina_dumeril
#Minor outliers: Alopias_vulpinus, Prionace_glauca
#
#######################################
#
#####permutations of the dataset BEGIN ##
#no Squatina
shark.sppave.no.squatina <- sharks.sppave.gpa$coords[,,c(dimnames(sharks.sppave.gpa$coords)[[3]][(dimnames(sharks.sppave.gpa$coords)[[3]]!="Squatina_dumeril")])]
length(sharks.sppave.gpa$coords)
length(shark.sppave.no.squatina)
#no Squatina, Alopias, or Prionace  
shark.sppave.no.outliers <- shark.sppave.no.squatina[,,c(dimnames(shark.sppave.no.squatina)[[3]][(dimnames(shark.sppave.no.squatina)[[3]]!="Alopias_vulpinus")])]
length(shark.sppave.no.outliers)
shark.sppave.no.outliers <- shark.sppave.no.outliers[,,c(dimnames(shark.sppave.no.outliers)[[3]][(dimnames(shark.sppave.no.outliers)[[3]]!="Prionace_glauca")])]
length(shark.sppave.no.outliers)
#
#no squatina dataset
shark.sppave.no.squatina.gpa <- gpagen(A=shark.sppave.no.squatina) 
  #plot again
plotAllSpecimens(shark.sppave.no.squatina.gpa$coords, label=T,
                 plot.param = list(pt.bg = landmark.colorz2, mean.cex=.5))
  #pca
shark.sppave.no.squatina.gpa.morphospace <- plotTangentSpace(shark.sppave.no.squatina.gpa$coords, axis1 = 1, axis2 = 2)
#                           PC1     PC2     PC3     PC4     PC5      PC6
#Standard deviation     0.08101 0.04277 0.03392 0.02290 0.02029 0.008174
#Proportion of Variance 0.62240 0.17343 0.10908 0.04972 0.03904 0.006340
#Cumulative Proportion  0.62240 0.79583 0.90491 0.95463 0.99366 1.000000
plotOutliers(shark.sppave.no.squatina.gpa$coords) # looks good
#
#
#no outliers dataset
shark.sppave.no.outliers.gpa <- gpagen(A=shark.sppave.no.outliers) 
#plot again
plotAllSpecimens(shark.sppave.no.outliers.gpa$coords, label=T,
                 plot.param = list(pt.bg = landmark.colorz2, mean.cex=.5))
#pca
shark.sppave.no.outliers.gpa.morphospace <- plotTangentSpace(shark.sppave.no.outliers.gpa$coords, axis1 = 1, axis2 = 2)
#                           PC1     PC2     PC3     PC4     PC5      PC6
#Standard deviation     0.06655 0.04216 0.02996 0.02246 0.01664 0.008547
#Proportion of Variance 0.55650 0.22334 0.11278 0.06341 0.03479 0.009180
#Cumulative Proportion  0.55650 0.77985 0.89263 0.95603 0.99082 1.000000
plotOutliers(shark.sppave.no.outliers.gpa$coords) # looks good
#
####### Permutations END ##
#'
#' We'll look at the phylomorphospace of pectoral fin shape now. To
#' make the process easier, we'll pull out the first to PC axes and
#' assign them to an object.
shark.sppave.full.pc1pc2 <- shark.sppave.morphopace.full$pc.scores[,1:2]
name.check(shark.tree, shark.sppave.full.pc1pc2) #name need to be in the same order
# Rearrange the taxa 
shark.sppave.full.pc1pc2 <- shark.sppave.full.pc1pc2[shark.tree$tip.label,]
row.names(shark.sppave.full.pc1pc2) == shark.tree$tip.label #checks out!
#  
#### Trim phylogeny to match reduced datasets
# No Squatina
# Do some prunin' a la the continuous_tutorial
TreeOnly.squatina <- setdiff(shark.tree$tip.label ,dimnames(shark.sppave.no.squatina.gpa$coords)[[3]])
TreeOnly.squatina # Enter the name of the object we just created to see what’s in it.
#Trimmed tree for no Squatina
shark.tree.nosquatina <- drop.tip(shark.tree,TreeOnly.squatina)
plot(shark.tree.nosquatina)
#Pull out PC1, PC2 and organize row order
shark.sppave.nosquatina.pc1pc2 <- shark.sppave.no.squatina.gpa.morphospace$pc.scores[,1:2]
name.check(shark.tree.nosquatina, shark.sppave.nosquatina.pc1pc2)
shark.sppave.nosquatina.pc1pc2 <- shark.sppave.nosquatina.pc1pc2[shark.tree.nosquatina$tip.label,]
row.names(shark.sppave.nosquatina.pc1pc2) == shark.tree.nosquatina$tip.label
#
#
# No outliers
# Do some prunin' a la the continuous_tutorial
TreeOnly.outliers <- setdiff(shark.tree$tip.label ,dimnames(shark.sppave.no.outliers.gpa$coords)[[3]])
TreeOnly.outliers # Enter the name of the object we just created to see what’s in it.
#Trimmed tree for no Squatina
shark.tree.nooutliers <- drop.tip(shark.tree,TreeOnly.outliers)
plot(shark.tree.nooutliers)
#Pull out PC1, PC2 and organize row order
shark.sppave.nooutliers.pc1pc2 <- shark.sppave.no.outliers.gpa.morphospace$pc.scores[,1:2]
name.check(shark.tree.nooutliers, shark.sppave.nooutliers.pc1pc2)
shark.sppave.nooutliers.pc1pc2 <- shark.sppave.nooutliers.pc1pc2[shark.tree.nooutliers$tip.label,]
row.names(shark.sppave.nooutliers.pc1pc2) == shark.tree.nooutliers$tip.label
#
#
############################### Phylomorphospace and stats #######
##########Full dataset
#read in new color info for spp aves
shark.sppave.colorz <- read.csv( file = "Fin_data_Sppave_200705_shape.csv", row.names = 1, header = T)
#'
#'Check the name order as above and resort to align
name.check(shark.tree, shark.sppave.colorz)
shark.sppave.colorz <- shark.sppave.colorz[shark.tree$tip.label,]
row.names(shark.sppave.colorz) == shark.tree$tip.label
#
plotGMPhyloMorphoSpace(shark.tree, shark.sppave.full.pc1pc2, 
                       tip.labels = F, yaxis = 2, node.labels = F,
                       plot.param = list(t.bg=as.character(shark.sppave.colorz$Ecotype.new.col),
                                         txt.cex=0.5)
)
legend(x= "topleft", legend = c(unique(as.character(shark.sppave.colorz$Ecotype.new))) , 
       pch = 21, pt.cex = 1.5, col = "black", pt.bg = c(unique(as.character(shark.sppave.colorz$Ecotype.new.col)))   )
title(main = "Sharks Ecology - Full Dataset")
#
#color by taxonomic order
plotGMPhyloMorphoSpace(shark.tree, shark.sppave.full.pc1pc2, 
                       tip.labels = T, yaxis = 2, node.labels = F,
                       plot.param = list(t.bg=as.character(shark.sppave.colorz$Order.col),
                                         txt.cex=0.5)
)
legend(x= 0.0, y = 0.35, legend = c(unique(as.character(shark.sppave.colorz$Order))) , 
       pch = 21, pt.cex = 1.5, col = "black", pt.bg = c(unique(as.character(shark.sppave.colorz$Order.col)))   )
title(main = "Sharks Taxonomic Order - Full Dataset")
#
#color by taxonomic family
plotGMPhyloMorphoSpace(shark.tree, shark.sppave.full.pc1pc2, 
                       tip.labels = F, yaxis = 2, node.labels = F,
                       plot.param = list(t.bg=as.character(shark.sppave.colorz$Family.col),
                                         txt.cex=0.5)
)
legend(x= 0.0, y = 0.35, legend = c(unique(as.character(shark.sppave.colorz$Family))) , 
       pch = 21, pt.cex = 1.5, col = "black", pt.bg = c(unique(as.character(shark.sppave.colorz$Family.col)))   )
title(main = "Sharks Taxonomic Family - Full Dataset")
#
##### shape transformation
##### Shape change associated with each PC axis ##### 
ave.shark.fin.shape.full <- mshape(sharks.sppave.gpa$coords) #average shape of all specimens
#'
#' Now we'll import the landmark coordinate data of our ouline
shark.fin.ref.specimen <- readland.tps(file = "Porbeagle_3_NAST_2016.TPS", readcurves = F, specID = "imageID")
dev.off() #reset the display
# Transform the outline to take on the average shape
shark.fin.shape.full.ave.outline <- warpRefOutline(file ="Porbeagle_3_NAST_2016_coords.TPS", coord =shark.fin.ref.specimen[,,1], ref=ave.shark.fin.shape.full)
#' Now we can use this outline of the average shape to visualize the shape change captured by each
#' of the PC axes. We'll start with just an un-transformed average shape and TPS grid to give ourselves a point
#' of reference
plotRefToTarget(ave.shark.fin.shape.full, ave.shark.fin.shape.full, outline = shark.fin.shape.full.ave.outline$outline, method = "TPS")
title(main = "Untransformed- Full Dataset")
#'
#' Next, we'll visualize the shape change associated with the most negative observed value of PC1
plotRefToTarget(ave.shark.fin.shape.full, shark.sppave.morphopace.full$pc.shapes$PC1min, outline = shark.fin.shape.full.ave.outline$outline, method = "TPS")
title(main = "PC1 negative- Full Dataset")
#'
#' Now the most positive value of PC1
plotRefToTarget(ave.shark.fin.shape.full, shark.sppave.morphopace.full$pc.shapes$PC1max, outline = shark.fin.shape.full.ave.outline$outline, method = "TPS")
title(main = "PC1 positive- Full Dataset")
#'
#' Most negative value of PC2
plotRefToTarget(ave.shark.fin.shape.full, shark.sppave.morphopace.full$pc.shapes$PC2min, outline = shark.fin.shape.full.ave.outline$outline, method = "TPS")
title(main = "PC2 negative- Full Dataset")
#' Most positive value of PC2
plotRefToTarget(ave.shark.fin.shape.full, shark.sppave.morphopace.full$pc.shapes$PC2max, outline = shark.fin.shape.full.ave.outline$outline, method = "TPS")
title(main = "PC2 positive- Full Dataset")
#' Most negative value of PC3
plotRefToTarget(ave.shark.fin.shape.full, shark.sppave.morphopace.full$pc.shapes$PC3min, outline = shark.fin.shape.full.ave.outline$outline, method = "TPS")
title(main = "PC3 negative- Full Dataset")
#' Most positive value of PC3
plotRefToTarget(ave.shark.fin.shape.full, shark.sppave.morphopace.full$pc.shapes$PC3max, outline = shark.fin.shape.full.ave.outline$outline, method = "TPS")
title(main = "PC3 positive- Full Dataset")
#' Most negative value of PC4
plotRefToTarget(ave.shark.fin.shape.full, shark.sppave.morphopace.full$pc.shapes$PC4min, outline = shark.fin.shape.full.ave.outline$outline, method = "TPS")
title(main = "PC4 negative- Full Dataset")
#' Most positive value of PC4
plotRefToTarget(ave.shark.fin.shape.full, shark.sppave.morphopace.full$pc.shapes$PC4max, outline = shark.fin.shape.full.ave.outline$outline, method = "TPS")
title(main = "PC4 positive- Full Dataset")
#
# Regular MANOVA
gdf_sharkTPS_habitat.full.reg <- geomorph.data.frame(coords =shark.sppave.morphopace.full$pc.scores,
                                                habitat = shark.sppave.colorz$Ecology)
#
procD.lm(coords~habitat, data = gdf_sharkTPS_habitat.full.reg)
#
### Phylogenetic MANOVA
#' First we'll make a geomorph-formatted dataframe to hold our shape variables, 
#' phylogenetic data, and habitat type for each species
gdf_sharkTPS_habitat.full <- geomorph.data.frame(coords =shark.sppave.morphopace.full$pc.scores, phy=shark.tree,
                                            habitat = shark.sppave.colorz$Ecology)
#
attributes(gdf_sharkTPS_habitat.full) #looks like they're all there
#
# Now we can perform our phylogenetic ANOVA
procD.pgls(coords ~ habitat, phy = phy, data = gdf_sharkTPS_habitat.full)
#' It looks like there is NOT a significant difference in average fin shape
#' across the different habitats.
#
########### No Squatina dataset
# Remove Squatina from colorz
shark.sppave.colorz.nosquatina <- shark.sppave.colorz[shark.tree.nosquatina$tip.label,]
row.names(shark.sppave.colorz.nosquatina) == shark.tree.nosquatina$tip.label
#
shark.sppave.colorz.nosquatina$Ecotype.new <- as.factor(as.character(shark.sppave.colorz.nosquatina$Ecotype.new))
shark.sppave.colorz.nosquatina$Ecotype.new.col <- as.factor(as.character(shark.sppave.colorz.nosquatina$Ecotype.new.col))
#
plotGMPhyloMorphoSpace(shark.tree.nosquatina, shark.sppave.nosquatina.pc1pc2, 
                       tip.labels = T, yaxis = 2, node.labels = F,
                       plot.param = list(t.bg=as.character(shark.sppave.colorz.nosquatina$Ecotype.new.col),
                                         txt.cex=0.5)
)
legend(x= "topleft", legend = c(unique(as.character(shark.sppave.colorz.nosquatina$Ecotype.new))) , 
       pch = 21, pt.cex = 1.5, col = "black", pt.bg = c(unique(as.character(shark.sppave.colorz.nosquatina$Ecotype.new.col)))   )
title(main = "Sharks Ecology - No Squatina")
#
#color by taxonomic order
plotGMPhyloMorphoSpace(shark.tree.nosquatina, shark.sppave.nosquatina.pc1pc2, 
                       tip.labels = F, yaxis = 2, node.labels = F,
                       plot.param = list(t.bg=as.character(shark.sppave.colorz.nosquatina$Order.col),
                                         txt.cex=0.5)
)
legend(x= 0.0, y = 0.2, legend = c(unique(as.character(shark.sppave.colorz.nosquatina$Order))) , 
       pch = 21, pt.cex = 1.5, col = "black", pt.bg = c(unique(as.character(shark.sppave.colorz.nosquatina$Order.col)))   )
title(main = "Sharks Taxonomic Order - No Squatina Dataset")
#
#color by taxonomic family
plotGMPhyloMorphoSpace(shark.tree.nosquatina, shark.sppave.nosquatina.pc1pc2, 
                       tip.labels = T, yaxis = 2, node.labels = F,
                       plot.param = list(t.bg=as.character(shark.sppave.colorz.nosquatina$Family.col),
                                         txt.cex=0.5)
)
legend(x= 0.05, y = 0.2, legend = c(unique(as.character(shark.sppave.colorz.nosquatina$Family))) , 
       pch = 21, pt.cex = 1.5, col = "black", pt.bg = c(unique(as.character(shark.sppave.colorz.nosquatina$Family.col)))   )
title(main = "Sharks Taxonomic Family - No Squatina Dataset")
#
##### shape transformation
##### Shape change associated with each PC axis ##### 
ave.shark.fin.shape.nosquatina <- mshape(shark.sppave.no.squatina.gpa$coords) #average shape of all specimens
#'
shark.fin.shape.nosquatina.ave.outline <- warpRefOutline(file ="Porbeagle_3_NAST_2016_coords.TPS", coord =shark.fin.ref.specimen[,,1], ref=ave.shark.fin.shape.nosquatina)
#' Now we can use this outline of the average shape to visualize the shape change captured by each
#' of the PC axes. We'll start with just an un-transformed average shape and TPS grid to give ourselves a point
#' of reference
plotRefToTarget(ave.shark.fin.shape.nosquatina, ave.shark.fin.shape.nosquatina, outline = shark.fin.shape.nosquatina.ave.outline$outline, method = "TPS")
title(main = "Untransformed- No Squatina Dataset")
#'
#' Next, we'll visualize the shape change associated with the most negative observed value of PC1
plotRefToTarget(ave.shark.fin.shape.nosquatina, shark.sppave.no.squatina.gpa.morphospace$pc.shapes$PC1min, outline = shark.fin.shape.nosquatina.ave.outline$outline, method = "TPS")
title(main = "PC1 negative- No Squatina Dataset")
#'
#' Now the most positive value of PC1
plotRefToTarget(ave.shark.fin.shape.nosquatina, shark.sppave.no.squatina.gpa.morphospace$pc.shapes$PC1max, outline = shark.fin.shape.nosquatina.ave.outline$outline, method = "TPS")
title(main = "PC1 positive- No Squatina Dataset")
#'
#' Most negative value of PC2
plotRefToTarget(ave.shark.fin.shape.nosquatina, shark.sppave.no.squatina.gpa.morphospace$pc.shapes$PC2min, outline = shark.fin.shape.nosquatina.ave.outline$outline, method = "TPS")
title(main = "PC2 negative- No Squatina Dataset")
#' Most positive value of PC2
plotRefToTarget(ave.shark.fin.shape.nosquatina, shark.sppave.no.squatina.gpa.morphospace$pc.shapes$PC2max, outline = shark.fin.shape.nosquatina.ave.outline$outline, method = "TPS")
title(main = "PC2 positive- No Squatina Dataset")
#' Most negative value of PC3
plotRefToTarget(ave.shark.fin.shape.nosquatina, shark.sppave.no.squatina.gpa.morphospace$pc.shapes$PC3min, outline = shark.fin.shape.nosquatina.ave.outline$outline, method = "TPS")
title(main = "PC3 negative- No Squatina Dataset")
#' Most positive value of PC3
plotRefToTarget(ave.shark.fin.shape.nosquatina, shark.sppave.no.squatina.gpa.morphospace$pc.shapes$PC3max, outline = shark.fin.shape.nosquatina.ave.outline$outline, method = "TPS")
title(main = "PC3 positive- No Squatina Dataset")
#' Most negative value of PC4
plotRefToTarget(ave.shark.fin.shape.nosquatina, shark.sppave.no.squatina.gpa.morphospace$pc.shapes$PC4min, outline = shark.fin.shape.nosquatina.ave.outline$outline, method = "TPS")
title(main = "PC4 negative- No Squatina Dataset")
#' Most positive value of PC4
plotRefToTarget(ave.shark.fin.shape.nosquatina, shark.sppave.no.squatina.gpa.morphospace$pc.shapes$PC4max, outline = shark.fin.shape.nosquatina.ave.outline$outline, method = "TPS")
title(main = "PC4 positive- No Squatina Dataset")
#
# Regular MANOVA
gdf_sharkTPS_habitat.nosquatina.reg <- geomorph.data.frame(coords =shark.sppave.no.squatina.gpa.morphospace$pc.scores,
                                                     habitat = shark.sppave.colorz.nosquatina$Ecotype.new)
#
summary(procD.lm(coords~habitat, data = gdf_sharkTPS_habitat.nosquatina.reg))
#
### Phylogenetic MANOVA
#' First we'll make a geomorph-formatted dataframe to hold our shape variables, 
#' phylogenetic data, and habitat type for each species
gdf_sharkTPS_habitat.nosquatina <- geomorph.data.frame(coords =shark.sppave.no.squatina.gpa.morphospace$pc.scores, phy=shark.tree.nosquatina,
                                                 habitat = shark.sppave.colorz.nosquatina$Ecotype.new)
#
attributes(gdf_sharkTPS_habitat.nosquatina) #looks like they're all there
#
# Now we can perform our phylogenetic ANOVA
summary(procD.pgls(coords ~ habitat, phy = phy, data = gdf_sharkTPS_habitat.nosquatina))
#' It looks like there is NOT a significant difference in average fin shape
#' across the different habitats.
#
########### No outliers dataset #######
# Remove all outliers from colorz
shark.sppave.colorz.nooutliers <- shark.sppave.colorz[shark.tree.nooutliers$tip.label,]
row.names(shark.sppave.colorz.nooutliers) == shark.tree.nooutliers$tip.label

plotGMPhyloMorphoSpace(shark.tree.nooutliers, shark.sppave.nooutliers.pc1pc2, 
                       tip.labels = F, yaxis = 2, node.labels = F,
                       plot.param = list(t.bg=as.character(shark.sppave.colorz.nooutliers$Ecotype.new.col),
                                         txt.cex=0.5)
)
legend(x= 0.1, y = 0.2, legend = c(unique(as.character(shark.sppave.colorz.nooutliers$Ecotype.new))) , 
       pch = 21, pt.cex = 1.5, col = "black", pt.bg = c(unique(as.character(shark.sppave.colorz.nooutliers$Ecotype.new.col)))   )
title(main = "Sharks Ecology - No Outliers")
#
#color by taxonomic order
plotGMPhyloMorphoSpace(shark.tree.nooutliers, shark.sppave.nooutliers.pc1pc2, 
                       tip.labels = T, yaxis = 2, node.labels = F,
                       plot.param = list(t.bg=as.character(shark.sppave.colorz.nooutliers$Order.col),
                                         txt.cex=0.5)
)
legend(x= -0.12, y = 0.12, legend = c(unique(as.character(shark.sppave.colorz.nooutliers$Order))) , 
       pch = 21, pt.cex = 1.5, col = "black", pt.bg = c(unique(as.character(shark.sppave.colorz.nooutliers$Order.col)))   )
title(main = "Sharks Taxonomic Order - No Outliers Dataset")
#
#color by taxonomic family
plotGMPhyloMorphoSpace(shark.tree.nosquatina, shark.sppave.nosquatina.pc1pc2, 
                       tip.labels = T, yaxis = 2, node.labels = F,
                       plot.param = list(t.bg=as.character(shark.sppave.colorz.nosquatina$Family.col),
                                         txt.cex=0.5)
)
legend(x= 0.05, y = 0.2, legend = c(unique(as.character(shark.sppave.colorz.nosquatina$Family))) , 
       pch = 21, pt.cex = 1.5, col = "black", pt.bg = c(unique(as.character(shark.sppave.colorz.nosquatina$Family.col)))   )
title(main = "Sharks Taxonomic Family - No Squatina Dataset")
#
##### shape transformation
##### Shape change associated with each PC axis ##### 
ave.shark.fin.shape.nooutliers <- mshape(shark.sppave.no.outliers.gpa$coords) #average shape of all specimens
#'
shark.fin.shape.nooutliers.ave.outline <- warpRefOutline(file ="Porbeagle_3_NAST_2016_coords.TPS", coord =shark.fin.ref.specimen[,,1], ref=ave.shark.fin.shape.nooutliers)
#' Now we can use this outline of the average shape to visualize the shape change captured by each
#' of the PC axes. We'll start with just an un-transformed average shape and TPS grid to give ourselves a point
#' of reference
plotRefToTarget(ave.shark.fin.shape.nooutliers, ave.shark.fin.shape.nooutliers, outline = shark.fin.shape.nooutliers.ave.outline$outline, method = "TPS")
title(main = "Untransformed- No Outliers Dataset")
#'
#' Next, we'll visualize the shape change associated with the most negative observed value of PC1
plotRefToTarget(ave.shark.fin.shape.nooutliers, shark.sppave.no.outliers.gpa.morphospace$pc.shapes$PC1min, outline = shark.fin.shape.nooutliers.ave.outline$outline, method = "TPS")
title(main = "PC1 negative- No Outliers Dataset")
#'
#' Now the most positive value of PC1
plotRefToTarget(ave.shark.fin.shape.nooutliers, shark.sppave.no.outliers.gpa.morphospace$pc.shapes$PC1max, outline = shark.fin.shape.nooutliers.ave.outline$outline, method = "TPS")
title(main = "PC1 positive- No Outliers Dataset")
#'
#' Most negative value of PC2
plotRefToTarget(ave.shark.fin.shape.nooutliers, shark.sppave.no.outliers.gpa.morphospace$pc.shapes$PC2min, outline = shark.fin.shape.nooutliers.ave.outline$outline, method = "TPS")
title(main = "PC2 negative- No Outliers Dataset")
#' Most positive value of PC2
plotRefToTarget(ave.shark.fin.shape.nooutliers, shark.sppave.no.outliers.gpa.morphospace$pc.shapes$PC2max, outline = shark.fin.shape.nooutliers.ave.outline$outline, method = "TPS")
title(main = "PC2 positive- No outliers Dataset")
#' Most negative value of PC3
plotRefToTarget(ave.shark.fin.shape.nooutliers, shark.sppave.no.outliers.gpa.morphospace$pc.shapes$PC3min, outline = shark.fin.shape.nooutliers.ave.outline$outline, method = "TPS")
title(main = "PC3 negative- No Outliers Dataset")
#' Most positive value of PC3
plotRefToTarget(ave.shark.fin.shape.nooutliers, shark.sppave.no.outliers.gpa.morphospace$pc.shapes$PC3max, outline = shark.fin.shape.nooutliers.ave.outline$outline, method = "TPS")
title(main = "PC3 positive- No Outliers Dataset")
#' Most negative value of PC4
plotRefToTarget(ave.shark.fin.shape.nooutliers, shark.sppave.no.outliers.gpa.morphospace$pc.shapes$PC4min, outline = shark.fin.shape.nooutliers.ave.outline$outline, method = "TPS")
title(main = "PC4 negative- No Outliers Dataset")
#' Most positive value of PC4
plotRefToTarget(ave.shark.fin.shape.nooutliers, shark.sppave.no.outliers.gpa.morphospace$pc.shapes$PC4max, outline = shark.fin.shape.nooutliers.ave.outline$outline, method = "TPS")
title(main = "PC4 positive- No Outliers Dataset")
#
# Regular MANOVA
gdf_sharkTPS_habitat.nooutliers.reg <- geomorph.data.frame(coords =shark.sppave.no.outliers.gpa.morphospace$pc.scores,
                                                           habitat = shark.sppave.colorz.nooutliers$Ecotype.new)
#
summary(procD.lm(coords~habitat, data = gdf_sharkTPS_habitat.nooutliers.reg))
# 
### Phylogenetic MANOVA
#' First we'll make a geomorph-formatted dataframe to hold our shape variables, 
#' phylogenetic data, and habitat type for each species
gdf_sharkTPS_habitat.nooutliers <- geomorph.data.frame(coords =shark.sppave.no.outliers.gpa.morphospace$pc.scores, phy=shark.tree.nooutliers,
                                                       habitat = shark.sppave.colorz.nooutliers$Ecotype.new)
#
attributes(gdf_sharkTPS_habitat.nooutliers) #looks like they're all there
#
# Now we can perform our phylogenetic ANOVA
summary(procD.pgls(coords ~ habitat, phy = phy, data = gdf_sharkTPS_habitat.nooutliers))
#' It looks like there is NOT a significant difference in average fin shape
#' across the different habitats.
#'
#'
####'
####'  SECTION II: PECTORAL FIN INTERNAL ANATOMY
####'  
######## Begin measurement data analysis
# We'll do a PCA of the internal pectoral fin anatomy variables across all specimens
test.pca <- prcomp(fin.measures[,18:37])
plot(test.pca$x[,"PC1"], test.pca$x[,"PC2"])
text(as.numeric(test.pca$x[,"PC1"]), test.pca$x[,"PC2"],
     as.character(fin.measures$Genus.species))
test.pca
#
#' Now we can use the "broken stick" method to estimate the number of "significant" PC axes.
#' See methods for discussion of "broken stick" method. 
screeplot(test.pca, bstick = T, type = "lines")
#' Looks like the first PC axis definitely explains more than would be expected by chance
#' and the second explains about as much as would be expected by chance, while the rest 
#' explain less.
#
# Let's calculate the % var captured by each PC axis
test.pca$sdev/(sum(test.pca$sdev))
# Let's color code each taxonomic family to see how they're distributed
# in morphospace
plot(test.pca$x[,"PC1"], test.pca$x[,"PC2"])
points(as.numeric(test.pca$x[,"PC1"]), test.pca$x[,"PC2"],
       pch =21, bg =as.character(fin.measures$Family.col))
legend(x= -40, y = 30, legend = c(unique(as.character(fin.measures$Family))) , 
       pch = 21, pt.cex = 1.5, col = "black", pt.bg = c(unique(as.character(fin.measures$Family.col)))   )
title(main = "Measurement PCA- Family Coloration")
#'
#' Okay, now we'll use a simple loop to take the average value of these 
#' traits for each species.
dataset.matrix <- data.frame(fin.measures[,6],fin.measures[,18:37])
colnames(dataset.matrix)[1]<-"Species"
#
new.matrix <- matrix(nrow = length(unique(dataset.matrix[,"Species"])), ncol = ncol(dataset.matrix))
row.names(new.matrix) <- as.character(unique(dataset.matrix[,"Species"]))
colnames(new.matrix) <- colnames(dataset.matrix)
#
for(j in 2:ncol(new.matrix)){
#
  for(i in 1:nrow(new.matrix)){
    new.matrix[i,j] <- mean(dataset.matrix[dataset.matrix[,"Species"] == row.names(new.matrix)[i], j])
#    
  }
}
new.matrix <- new.matrix[,2:ncol(new.matrix)]
#
shark.sppave.linear.data <- new.matrix
#
############################### SPECIES AVERAGES BEGIN ################
#' Before we can test for relationships, we want to ensure that the set of taxa in 
#' our morphological dataset perfectly matches the set of taxa in our phylogeny. 
#' To do so, we will remove any taxa that are not found in both datasets.
## NOTE: the following six lines of code were modified from the "Continuous Trait R 
## Tutorial" written by D. Luke Mahler (2014), available at: 
## http://treethinkers.org/tutorials/continuous-traits/
#'
TreeOnly <- setdiff(shark.tree$tip.label ,rownames(shark.sppave.linear.data))
TreeOnly # Enter the name of the object we just created to see what’s in it.
#'
#Trimmed tree
shark.tree.pruned <- drop.tip(shark.tree,TreeOnly)
plot(shark.tree.pruned) #Let's take a look at our trimmed phylogeny
#'
#Make sure the names are all in order
name.check(shark.tree.pruned, shark.sppave.linear.data) #good
#Double check order
shark.tree.pruned$tip.label== row.names(shark.sppave.linear.data) #not yet
#' Rearrange to match the order of taxa in the phylo
shark.sppave.linear.data <- shark.sppave.linear.data[shark.tree.pruned$tip.label,]
#Double check order
shark.tree.pruned$tip.label== row.names(shark.sppave.linear.data) #good!
#
# Read in species average dataset
shark.sppave.colorz.full <- read.csv( file = "Fin_data_Sppave_200705.csv", row.names = 1, header = T)
shark.sppave.colorz.full <- shark.sppave.colorz.full[rownames(shark.sppave.linear.data),]
#
shark.sppave.linear.data.pca <- prcomp(shark.sppave.linear.data)
plot(shark.sppave.linear.data.pca$x[,"PC1"], shark.sppave.linear.data.pca$x[,"PC2"])
text(as.numeric(shark.sppave.linear.data.pca$x[,"PC1"]), shark.sppave.linear.data.pca$x[,"PC2"],
     as.character(row.names(shark.sppave.linear.data)))
shark.sppave.linear.data.pca
#
#' Now we can use the "broken stick" method to estimate the number of "significant" PC axes.
#' See methods for discussion of "broken stick" method. 
screeplot(shark.sppave.linear.data.pca, bstick = T, type = "lines")
# Looks like only the first PC axis explains more variance than is expected by chance
#'
#' See the percent of total variation summarized by each PC axis
shark.sppave.linear.data.pca$sdev/(sum(shark.sppave.linear.data.pca$sdev))
#' Make an object to help us color-code our plot points 
measure.colorz <- shark.sppave.colorz.full[row.names(shark.sppave.linear.data.pca$x),]
#'
#' Plot the morphospace of the first two PC axes and color-code the points
#' by taxonomic family
plot(shark.sppave.linear.data.pca$x[,"PC1"], shark.sppave.linear.data.pca$x[,"PC2"])
points(as.numeric(shark.sppave.linear.data.pca$x[,"PC1"]), shark.sppave.linear.data.pca$x[,"PC2"],
       pch =21, bg =as.character(measure.colorz$Family.col))
legend(x= 15, y = -2, legend = c(unique(as.character(measure.colorz$Family))) , 
       pch = 21, pt.cex = 1.5, col = "black", pt.bg = c(unique(as.character(measure.colorz$Family.col)))   )
title(main = "Measurement PCA- Family Coloration")
#
# Let's make it a little more formal
linear.sppave.full.PC1PC2 <- shark.sppave.linear.data.pca$x[,1:2]
name.check(shark.tree.pruned, linear.sppave.full.PC1PC2)
#
shark.sppave.colorz.full <- shark.sppave.colorz.full[shark.tree.pruned$tip.label,]
row.names(shark.sppave.colorz.full) == shark.tree.pruned$tip.label
#
linear.sppave.full.PC1PC2 <- linear.sppave.full.PC1PC2[shark.tree.pruned$tip.label,]
row.names(linear.sppave.full.PC1PC2) == shark.tree.pruned$tip.label
#
# Now we'll plot the phylomorphospace and color-code by ecomorphotype of each species
plotGMPhyloMorphoSpace(shark.tree.pruned, linear.sppave.full.PC1PC2, 
                       tip.labels = F, yaxis = 2, node.labels = F,
                       plot.param = list(t.bg=as.character(shark.sppave.colorz.full$Ecotype.new.col),
                                         txt.cex=0.5)
)
legend(x = "topleft", legend = c(unique(as.character(shark.sppave.colorz.full$Ecotype.new))) , 
       pch = 21, pt.cex = 1.5, col = "black", pt.bg = c(unique(as.character(shark.sppave.colorz.full$Ecotype.new.col)))   )
title(main = "Sharks Ecology - Linear Data, Full Dataset")
#
####################### Phylogenetic MANOVA
# We'll make some objects to hold our data and make running the 
# phyMANOVA easier
shark.groups <- as.factor(measure.colorz$Ecotype.new)
names(shark.groups) <- row.names(shark.sppave.linear.data)
shark.sppave.linear.data
# Check the order
row.names(shark.sppave.linear.data) == names(shark.groups)
row.names(shark.sppave.linear.data) == shark.tree.pruned$tip.label
row.names(shark.sppave.linear.data.pca$x) == shark.tree.pruned$tip.label
# Start with a regular MANOVA
manova(shark.sppave.linear.data.pca$x~shark.groups)
# Now we'll do a phyMANOVA. The function that we'll use for this first attempt
# is limited in the number of dependant variables that can be used by the number
# of observations, so we'll use the maximum number of PC axes that the function 
# can handle. 
shark.pc.data <- shark.sppave.linear.data.pca$x[,1:9]
shark.MANOVA <- aov.phylo(shark.pc.data~shark.groups, phy = shark.tree.pruned, nsim = 10000, test = "Wilks")
print(attributes(shark.MANOVA)$summary) # summary table
#                                                         
# Since the above test is limited in the number of dependant variables that it can consider,
# let's try with the geomorph MANOVA and phyMANOVA, which are not so limited.
linear.shark.full.gdf <- geomorph.data.frame(phy = shark.tree.pruned, habitat = shark.groups , shape =shark.sppave.linear.data)
attributes(linear.shark.full.gdf)
#
#regular MANOVA
summary(procD.lm(shape~habitat, iter = 10000, seed = 12345, data = linear.shark.full.gdf))
#phyMANOVA
summary(procD.pgls(shape~habitat, phy=phy, iter = 10000, seed = 12345, data = linear.shark.full.gdf))
# Significant differences between ecomorphotypes on all accounts!
#