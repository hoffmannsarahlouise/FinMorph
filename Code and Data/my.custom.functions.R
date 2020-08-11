#' Running list of custom functions
#' 
#' Corresponding author: T.J. Buser    
#' email: BuserT@OregonState.edu
#' 
#' 
#' Custom functions:
#' 
#' 
######## Species Average for landmark data ###############
#'
#requires two variables:
# 1) "classifier" - must be an n by 2 matrix, with "ID" as the first column header
#' and "Species" as the second column header. Each row should contain the information
#' for a given specimen, such that "ID" is the entry in the array (see below) for that 
#' individual, and "Species" is the species of that inividual.
#' 2) "aligned.coords" - must be a p x k x n array of procrustes-aligned landmark data. 
#' Where "p" is the number of landmarks, "k" is the number of landmark dimensions 
#' (i.e., "2" or "3"), and "n" is the number of specimens. This format is the default 
#' output for the procrustes alignment functions found in the package "geomorph," 
#' specifically "gpagen()" and "bilat.symmetry." For gpagen, you'll need to specify the
#' "$coords" componenet of the functions output, for bilat.symmetry, you'll need to 
#' specify either "$symm.shape" or "$asymm.shape."
#'
#'
#'
my.landmark.species.average <- function(classifier, aligned.coords){
#'  
  require(geomorph)
  # Number of species
  species.n <- as.character(unique(classifier$Species))
  
  # Set up an array to hold the data for each species
  # Need the number of coordinates, dimensions, and species
  #pick the first entry to get coordinates and dimensions
  dummy.dim <- aligned.coords[,,1] 
  #set up array with appropriate attributes
  coords.array <- array(0, dim = c(nrow(dummy.dim), ncol(dummy.dim), length(species.n)))   #do this only the first time
  #label with species names
  dimnames(coords.array) <- list(c(), c(), c(species.n))  #do this only the first time
  
  
  #### loop section
  #i <-1
  
  for(i in 1:length(species.n)){
    # Species to be worked on for loop i
    fish.sp <-classifier[classifier$Species==species.n[i],]  ##
    if(length(fish.sp$ID)==1){
      dummy1 <- as.numeric(row.names(fish.sp))
      species.coords <- aligned.coords[,,c(dummy1)]
      coords.array[,,species.n[i]] <- species.coords
    }
    else{
      # Rows containing relevant species
      dummy1 <- as.numeric(row.names(fish.sp))                  ## 
      # Array of coordinate data for each individual of pertinent species
      species.coords <- aligned.coords[,,c(dummy1)]           ##
      # matrix of average coordinates data for pertinent species
      mean.species.coords <- mshape(species.coords)            ##
      # put the data in its appropriate place in the array (i.e., for the pert. species)
      coords.array[,,species.n[i]] <- mean.species.coords ##
    }
  }
  print(coords.array)
}
#'
#' Proper syntax for using the function:
#' my.landmark.species.average(classifier = , aligned.coords = )
#' 
#' 
#' 
#' 
########## End #####################
#'
#
##################### Convert MorphoJ landmark data (2D matrix) into 3D array format ########
#'
#' Function for converting MorphoJ landmark data output into an array format that is readable
#' by functions in the package "geomorph." 
#'
MorphoJ.Convert <- function(num.landmarks, num.dimensions, coordinate.data){
  
  test.array <- array(0, dim = c(num.landmarks, num.dimensions, nrow(coordinate.data)))   #do this only the first time
  dummy.names <- as.character(row.names(coordinate.data))
  dimnames(test.array) <- list(c(), c(), c(dummy.names))
  
  for(i in 1:length(dummy.names)){
    species.matrix <- matrix(as.numeric(coordinate.data[dummy.names[i],]), nrow = num.landmarks, ncol = num.dimensions, byrow = T)
    test.array[,,dummy.names[i]] <- species.matrix   
  }
  return(test.array)
}
#'
#'
#' Proper syntax for using the function:
#' MorphoJ.Convert(num.landmarks = , num.dimensions = , coordinate.data = )
#'
#'
#'
############### End #################
#'
#'
#'
########################################################
#
############ Pairwise Procrustes distance matrix generator function ##################
#
# "data.array"  must be a p x k x n array of procrustes-aligned landmark data. 
#' Where "p" is the number of landmarks, "k" is the number of landmark dimensions 
#' (i.e., "2" or "3"), and "n" is the number of specimens. This format is the default 
#' output for the procrustes alignment functions found in the package "geomorph," 
#' specifically "gpagen()" and "bilat.symmetry." For gpagen, you'll need to specify the
#' "$coords" componenet of the functions output, for bilat.symmetry, you'll need to 
#' specify either "$symm.shape" or "$asymm.shape."
#
########################################################
#
proc.dist.matrix <- function(data.array){
  
  name.dummy <- dimnames(data.array)
  length(name.dummy[[3]])
  name.vector <- name.dummy[[3]]
  comp.matrix <- matrix(data = 0, length(name.dummy[[3]]), length(name.dummy[[3]]))
  row.names(comp.matrix) <- name.dummy[[3]]
  colnames(comp.matrix) <- name.dummy[[3]]
  
  for (i in 1:length(name.vector)){
    #i = 1
    for (j in 1:length(name.vector)){
      #j = 6
      
      comp.matrix[i,j] <- procdist(data.array[,,name.vector[[i]]], data.array[,,name.vector[[j]]])
      
    }
  }
  return(comp.matrix)
}

#'
######## End ##############
#'
#'
#'
################# Procrustes distance from average shape ###############
#'
#'
#' Function for calculating a matrix of pairwise comparisons of procrustes
#' distance from one shape (in a group) to the mean shape (for said group).
#'
#' Function input is a 3D array of procrustes-aligned landmark coordinate
#' data. 
#'
proc.dist.from.ave <- function(data.array){
  
  ave.shape <- mshape(data.array)
  name.dummy <- dimnames(data.array)
  name.vector <- name.dummy[[3]]
  comp.matrix <- matrix(data = 0, nrow = 1, ncol = length(name.dummy[[3]]))
  row.names(comp.matrix) <- c("Average_shape")
  colnames(comp.matrix) <- name.dummy[[3]]
  
  
  for (j in 1:length(name.vector)){
    
    comp.matrix[,j] <- procdist(ave.shape, data.array[,,name.vector[[j]]])
    
  }
  
  return(comp.matrix)
}
#'
#'
#'
############ End ##############
#'
#'
#'
# Reading in, reversing, and converting to .tps format
#
#
my.landmark.array.reflector <- function(landmark.array.keep, landmark.array.reverse, species){
  landmark.array <- landmark.array.reverse
  landmark.array.dummy <- landmark.array
  for(i in 1:length(dimnames(landmark.array)[[3]])){
    landmark.array.dummy[,,i][,1] <- landmark.array[,,i][,1]*(-1)
  }
  writeland.tps(landmark.array.dummy, paste(species, ".landmark.array.reversed.side.tps", sep = ""))
  writeland.tps(landmark.array.keep, paste(species, ".landmark.array.kept.side.tps", sep = ""))
}
#'
#'
#'
############ End ##############
#'
#'
#'
# Curve resampling function

my.curve.resampler <- function(landmark.array, fixed.landmarks, curve.resample){
  
  require("geomorph")
  #my.curve.resampler(landmark.array = preoperc_lands, fixed.landmarks = 11, curve.resample = 38)
  
  
  #landmark.array <- preoperc_lands
  ########## landmark.array 
  #' must be a p x k x n array of un-alinged landmark data. 
  #' Where "p" is the number of landmarks, "k" is the number of landmark dimensions 
  #' (i.e., "2" or "3"), and "n" is the number of specimens. This format is the 
  #' standard for functions from the package "geomorph" that read landmark
  #' data into R (e.g., readland.nts).
  
  #fixed.landmarks = 11
  ########## fixed.landmarks
  #' The number of fixed landmarks in your dataset. These MUST be arranged such that
  #' they come before the semilandmarks. For example, if you have 10 landmarks and
  #' 50 semilandmarks, they must be arranged such that, for a given specimen, rows
  #' 1-10 are the landmark coordinates and 11-60 are the semilandmarks that make up
  #' the curve. At this time, the function can only handle semilandmarks that define
  #' a single curve (i.e., all semilandmarks must be part of the same curve).
  
  #curve.resample = 38
  ########## curve.resample
  #' The number of evenly-spaced points to be generated along the curve. This does not
  #' need to be the same number of points that were originally used to define the curve.
  
  
  
  # Number of specimens
  specimens.n <- as.character(dimnames(landmark.array)[[3]])
  
  # Set up an array to hold the data for each specimen
  # Need the number of coordinates, dimensions, and species
  #pick the first entry to get coordinates and dimensions
  dummy.dim <- landmark.array[,,1] 
  #set up array with appropriate attributes
  coords.array <- array(0, dim = c((sum(fixed.landmarks, curve.resample)+1), ncol(dummy.dim), length(specimens.n)))   #do this only the first time
  #label with specimen names
  dimnames(coords.array) <- list(c(), c(), c(specimens.n))  #do this only the first time
  
  #i <-1
  
  for(i in 1:length(specimens.n)){
    coords.dummy <- landmark.array[,,i]
    curvestart.dummy <- coords.dummy[as.numeric(fixed.landmarks),] 
    curvecoords.dummy <- coords.dummy[as.numeric(fixed.landmarks):as.numeric(nrow(coords.dummy)),]
    curve.resampled.coords.dummy <- digit.curves(curvestart.dummy, curvecoords.dummy, nPoints = as.numeric(curve.resample), closed = F)   
    
    new.spec.coords.dummy <- matrix(NA, nrow = (sum(fixed.landmarks, nrow(curve.resampled.coords.dummy)) - 1), ncol = ncol(coords.dummy))
    new.spec.coords.dummy[1:as.numeric(fixed.landmarks),] <-  coords.dummy[1:as.numeric(fixed.landmarks),]
    new.spec.coords.dummy[as.numeric(fixed.landmarks):nrow(coords.array[,,i]),] <- curve.resampled.coords.dummy
    
    coords.array[,,i] <- new.spec.coords.dummy
  }
  col1 <- (as.numeric(fixed.landmarks):(nrow(coords.array[,,1])-2))
  col2 <- (as.numeric(fixed.landmarks) + 1):(nrow(coords.array[,,1])-1)
  col3 <- (as.numeric(fixed.landmarks) + 2):nrow(coords.array[,,1])
  row.length <- length((as.numeric(fixed.landmarks) + 2):nrow(coords.array[,,1]))
  sliders.matrix <- matrix(c(col1, col2, col3), nrow = row.length, ncol = 3)     
  colnames(sliders.matrix) <- c("before", "slide", "after") 
  write.csv(sliders.matrix, file = "my.slider.matrix.csv", row.names = F)
  print(coords.array)
}




################## End ###################

##' Function to read .dta and write a .tps file
##' 
##' Specify the file name in quotes and the function will read in the .dts file
##' then save it as .TPS format in the same folder, with the same file name but
##' appended with ".TPS"
##' 
##' Usefull for concatenating multiple .dta landmark files
##' 
##' **** Requires package geomorph ******
##' 
##' Thaddaeus Buser (busert@oregonstate.edu) August 2018


DTA2TPS.converter <- function(dta.filename){
  
  require(geomorph)
  
  dummy1 <- readland.nts(paste(dta.filename))
  writeland.tps(dummy1,paste(dta.filename,".TPS", sep = ""))
}


################## End ###################


my.plotspec <- function (spec, digitspec, fixed = NULL, ptsize = 1, centered = FALSE, 
                         ...) 
{
  mesh <- NULL
  if (inherits(spec, "shape3d") == TRUE || inherits(spec, "mesh3d") == 
      TRUE) {
    if (centered == TRUE) {
      specimen <- scale(as.matrix(t(spec$vb)[, -4]), scale = FALSE)
      spec$vb <- rbind(t(specimen), 1)
    }
    if (centered == FALSE) {
      specimen <- as.matrix(t(spec$vb)[, -4])
    }
    mesh <- spec
    if (is.null(mesh$material)) {
      mesh$material <- "gray"
    }
  }
  else if (inherits(spec, "matrix") == FALSE) {
    stop("File is not a shape3d/mesh3d object or xyz matrix")
  }
  else if (inherits(spec, "matrix") == TRUE && dim(spec)[2] == 
           3) {
    if (centered == TRUE) {
      specimen <- scale(spec, scale = FALSE)
    }
    if (centered == FALSE) {
      specimen <- spec
    }
  }
  else {
    stop("File is not matrix in form: vertices by xyz")
  }
  if (is.null(dim(digitspec)) || dim(digitspec)[2] != 3) {
    stop("Digitized file is not xyz matrix in form: p x k")
  }
  plot3d(specimen[, 1], specimen[, 2], specimen[, 3], size = ptsize, 
         aspect = FALSE, ...)
  if (!is.null(mesh)) {
    shade3d(mesh, add = TRUE)
  }
  if (!is.null(fixed)) {
    points3d(digitspec[1:fixed, ], aspect = FALSE, size = 10, 
             col = "red")
    points3d(digitspec[(fixed + 1):nrow(digitspec), ], aspect = F, 
             size = 10, col = "green")
  }
  else {
    points3d(digitspec, aspect = F, size = 30, col = "red")
  }
}

###################################### END ########

