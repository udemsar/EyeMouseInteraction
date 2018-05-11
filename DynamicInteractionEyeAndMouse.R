# ------------------------------------------------
# Measure of Dynamic Interaction between gaze and mouse trajectories
# ------------------------------------------------
# Supplementary Info 1 
# ------------------------------------------------
# Code provided as is and can be used or modified freely. 
#
# We would however appreciate if you cite the paper:
#
# Demsar U and Coltekin A, "Quantifying Gaze and Mouse Interactions on 
# Spatial Visual Interfaces with a New Methodology from Geographic Information Science"
# ------------------------------------------------
# Code author: Urska Demsar
# University of St Andrews
# St Andrews, Scotland, UK
# http://udemsar.com
# Date: 13 July 2016
# ------------------------------------------------
# This is the main file that implements the calculation of space-time densities
# around mouse and eye trajectories and the measure of dynamic interaction
# between the two trajectories. For description of the algorithm, see the paper.
# ------------------------------------------------
# Supporting R files:
# DecayFunctions.R
# DensityAroundOnePoint.R
# DensityAroundTrajectory.R
# SquaredDistance2points.R
# WriteCSVtableFromRectangular3DArrays.R
# 
# For details on space-time density calculation in these files see the following paper:
# Demsar U, Buchin K, van Loon EE and Shamoun-Baranes J, 2015, 
# Stacked space-time densities: a geovisualisation approach to explore 
# dynamics of space use over time. GeoInformatica, 19(1):85-115.
# DOI 10.1007/s10707-014-0207-5
#
# For details on Brownian bridge decay in these files see the following paper:
# Horne JS, Garton EO, Krone SM, Lewis JS, 2007, 
# Analyzing animal movements using Brownian bridges, Ecology 88(9):2354-2363
# ------------------------------------------------
# Input: two data files with eye and mouse trajectories for one task and all participants
#
# Each file must be a csv file including the following five variables:
# TaskNo, PartNo, GazePointX/MousePointX,	GazePointY/MousePointY and Time, where
#
# - TaskNo - number of current task
# - PartNo - number of participants (this are used as trajectory IDs)
# - GazePointX/MousePointX - x coordinate of each eye or mouse trajectory point 
# in the screen coordinate system and in pixels
# - GazePointY/MousePointY - y coordinate of each eye or mouse trajectory point 
# in the screen coordinate system and in pixels
# - Time - time in ms since the start of experiment, where eye and mouse trajectories
# both start a 0 and where both movement records are synchronised with each other
#
# ------------------------------------------------
# The following also need to be manually specified: 
# - volume extent (min/max X, min/max Y, min/max T)
# - voxel size
# - kernel size 
# - kernel type
# ------------------------------------------------
# Outputs: 
# - csv files with eye and mouse densities per participant given as volumes, i.e. with 
# four columns: x, y, t, density value
# - csv file of the sum volume per participant (eye + mouse density)
# - csv file with measure of dynamic interaction for each participant
# ------------------------------------------------


# ------------------------------------------------
# Parameters to be set/changed by the user for own data

# Set working directory to the current directory where the data should be located
setwd(getwd())

# Set voxel size in pixels  
voxelSize <- 10 # for both foveal and parafoveal interaction

# Set kernel size in pixels
kernelSize <- 50 # for foveal interaction
#kernelSize <- 110 # for parafoveal interaction

# Threshold for density sum decision function, calculated across the entire voxel
a <- (3/(pi*kernelSize^2))*voxelSize^2;

# Set spatio-temporal extent of the density volumes:
# X & Y axis - stimulus extent on the screen:
minXcoord <- 510 
maxXcoord <- 1410 
minYcoord <- 90 
maxYcoord <- 990 

# Select kernel type: 1-linear, 2-bisquare, 3-Gaussian, 4-Brownian bridges
method <- 1
#method <- 2
#method <- 3
#method <- 4

# If Brownian bridges, the user sets sigmas 1 and 2 based on uncertainty:
# sigma11 = uncertainty parameter of movement between two points
# sigma12 = uncertainty parameter around the trajectory point
# See Horne et al. 2007 (as per above)
if (method == 4) {
  sigma11 <- 3 # to be set by the user - this is just for testing
  sigma12 <- 1 # to be set by the user - this is just for testing
} else {
  sigma11 <- 0
  sigma12 <- 0
}

# Name of input and output files
# Input files
eye_trajectories.file <- 'Task_1_EyeTrajectories_Example.csv'
mouse_trajectories.file <- 'Task_1_MouseTrajectories_Example.csv'

# Output file
output_interaction.file <- 'ValuesOfDynamicInteraction.csv'

# ------------------------------------------------
# General setup & read data

# read additional library
library(plot3D)

# Read fuctions in additional files
# Function to calculate density for each trajectory 
source("DensityAroundOnePoint.R") 
source("DensityAroundTrajectory.R") 
# Functions for saving density as as csv coordinate file
source("WriteCSVtableFromRectangular3DArrays.R") 
# Functions for geometric and kernel calculations
source("SquaredDistance2points.R")
source("DecayFunctions.R")

# Read data
dfEyeAll <- read.csv(eye_trajectories.file,stringsAsFactors=FALSE)
head(dfEyeAll)
dfMouseAll <- read.csv(mouse_trajectories.file,stringsAsFactors=FALSE)
head(dfMouseAll)

# -----------------------------------------------
# Read unique participant IDs and how many there are
listOfParticipants <- unique(dfEyeAll$PartNo)
noOfParticipants <- length(listOfParticipants)

# Initialise list of dynamic interaction values for all participants
dynamicInteractionIndices <- vector(mode="numeric",length=noOfParticipants)

# -----------------------------------------------
# Loop across participants
for (shk in 1:noOfParticipants) { 

  # --------------------------
  # Find the eye and mouse trajectories of this particular participant
  participantEyeIndices <- which(dfEyeAll$PartNo == listOfParticipants[shk])
  participantMouseIndices <- which(dfMouseAll$PartNo == listOfParticipants[shk])
    
  # Extract trajectories of this one participant
  trajEyeData <- dfEyeAll[participantEyeIndices,]
  trajMouseData <- dfMouseAll[participantMouseIndices,]
    
  # Check that all ok
  head(trajEyeData)
  head(trajMouseData)  
  tail(trajEyeData)
  tail(trajMouseData)  
    
  # --------------------------
  # Build three 3D arrays of x, y, z coordinates and initialise the density volumes
  
  # Read temporal extent from the eye&mouse trajectories
  # Z axis - temporal extent, in miliseconds or minutes
  minZcoord <- min(trajEyeData$Time,trajMouseData$Time) # start point in time
  maxZcoord <- max(trajEyeData$Time,trajMouseData$Time) # end point in time
    
  startX <- floor(minXcoord/voxelSize)*voxelSize
  startY <- floor(minYcoord/voxelSize)*voxelSize
  startZ <- floor(minZcoord/voxelSize)*voxelSize
  endX <- ceiling(maxXcoord/voxelSize)*voxelSize
  endY <- ceiling(maxYcoord/voxelSize)*voxelSize
  endZ <- ceiling(maxZcoord/voxelSize)*voxelSize
    
  xvoxels <- (endX-startX)/voxelSize
  yvoxels <- (endY-startY)/voxelSize
  zvoxels <- (endZ-startZ)/voxelSize
    
  # Build the volume
  x <- seq(startX,endX,by=voxelSize)
  y <- seq(startY,endY,by=voxelSize)
  z <- seq(startZ,endZ,by=voxelSize)
  M <- mesh(x,y,z)
  xcoord <- M$x
  ycoord <- M$y
  zcoord <- M$z

  # Initialise individual participant eye and mouse densities as zeros everywhere 
  # in a 3D array of the same size as the xcoord/ycoord/vcoord volumes
  eyeDensity <- array(data=0,dim=dim(xcoord))
  mouseDensity <- array(data=0,dim=dim(xcoord))
  
  #-------------------------------------------
  # Find ID of the current trajectories and respective lenghts
  currentTrajID <- shk
  currentEyeTrajLength <- length(trajEyeData$Time)
  currentMouseTrajLength <- length(trajMouseData$Time)
  
  # currentLines will be a submatrices in trajData starting with coordinate data only
  currentEyeLine <- trajEyeData[,3:5]
  currentMouseLine <- trajMouseData[,3:5]
  
  #-------------------------------------------
  # Calculate eye density for current participant
  print('*********')
  print(paste(sprintf('Participant no.: %d',shk), sprintf('out of %d', noOfParticipants)))
  print(paste(sprintf('Eye density for trajectory no.: %d',currentTrajID),sprintf('with this many points: %d', currentEyeTrajLength) ))
    
  # Sometimes there may be only one point in a trajectory, meaning that there won't be a segment needed for
  # stacked density calculation. If that is the case, current Line will only have 1 point, so we just calculate
  # density around this one point.
  
  if (currentEyeTrajLength > 1) {
     
    # if the trajectory has two points or more, we calculate stacked densities around each segment
    eyeDensity <- DensityAroundTrajectory(currentEyeLine,kernelSize,xcoord,ycoord,zcoord,method,voxelSize,sigma11,sigma12) 
    
    } else {
      # if however the trajectory has only one point, then we calculate the density around this point
    eyeDensity <- DensityAroundOnePoint(currentEyeLine,kernelSize,xcoord,ycoord,zcoord,method,voxelSize,sigma11,sigma12)

  } # if (currentTrajLength > 1)
    
  # Export eye density 
  eyeDensity.file <- paste(c(listOfParticipants[shk],"_Density_Eye.csv"),collapse="")
  WriteCSVtableFromFourRectangular3DArrays(xcoord,ycoord,zcoord,eyeDensity,eyeDensity.file)
  
  #-------------------------------------------
  # Calculate mouse density for current participant
  print('*********')
  print(paste(sprintf('Participant no.: %d',shk), sprintf('out of %d', noOfParticipants)))
  print(paste(sprintf('Mouse density for trajectory no.: %d',currentTrajID),sprintf('with this many points: %d', currentMouseTrajLength) ))
  
  # Sometimes there may be only one point in a trajectory, meaning that there won't be a segment needed for
  # stacked density calculation. If that is the case, current Line will only have 1 point, so we just calculate
  # density around this one point.
  
  if (currentMouseTrajLength > 1) {
    
    # if the trajectory has two points or more, we calculate stacked densities around each segment
    mouseDensity <- DensityAroundTrajectory(currentMouseLine,kernelSize,xcoord,ycoord,zcoord,method,voxelSize,sigma11,sigma12) 
    
  } else {
    # if however the trajectory has only one point, then we calculate the density around this point
    mouseDensity <- DensityAroundOnePoint(currentMouseLine,kernelSize,xcoord,ycoord,zcoord,method,voxelSize,sigma11,sigma12)
    
  } # if (currentTrajLength > 1)
  
  # Export mouse density 
  mouseDensity.file <- paste(c(listOfParticipants[shk],"_Density_Mouse.csv"),collapse="")
  WriteCSVtableFromFourRectangular3DArrays(xcoord,ycoord,zcoord,mouseDensity,mouseDensity.file)
  
  #-------------------------------------------
  # Calculate density sum for current participant
  
  sumDensity <- eyeDensity+mouseDensity
  
  # Export density sum
  sumDensity.file <- paste(c(listOfParticipants[shk],"_Density_Sum.csv"),collapse="")
  WriteCSVtableFromFourRectangular3DArrays(xcoord,ycoord,zcoord,sumDensity,sumDensity.file)
  
  #-------------------------------------------
  # Calculate level of dynamic interaction for current participant

  # Find voxels that have the sum value of more than threshold a and cound them  
  interactionVoxels <- which(sumDensity > a)
  noInteractionVoxels <- length(interactionVoxels)
  
  # Normalise the count with total no. of voxels in the sumDensity
  dynamicInteraction <- noInteractionVoxels/length(sumDensity)
  
  # Attach this to the dynamicInteractionIndices table for each participant
  dynamicInteractionIndices[shk] <- dynamicInteraction
  
} # for (shk in 1:noOfParticpants)

# Export table of interaction indices and participant numbers
outputTable <- cbind(listOfParticipants,dynamicInteractionIndices)
write.csv(outputTable,output_interaction.file,row.names=FALSE)

