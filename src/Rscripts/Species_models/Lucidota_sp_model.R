# Lucidota sp 2D model 

#### Setting working directory ####
whereami <- getwd()
where.pieces <- strsplit(whereami,"/")
if(where.pieces[[1]][2]=="home") {
  print("On the cluster!")
  setwd("/home/waldrop@chapman.edu/firefly")
}

#### Load functions ####
source("./src/Rscripts/model_fxns.R")

n <- 76     # set seed number
reps <- 50  # number of replicates

#### Parameters ####
parameters <- list(
  Species = "Lucidota_sp",
  # Measurements from morphology on antennomere 6
  width.olf.hair = 9.67e-6,         # mean width of an olfactory sensillum
  width.mech.hair.prox = 5.16e-6,   # mean width of an mechanosensory sensillum (proximal)
  width.mech.hair.med = 2.82e-6,    # mean width of an mechanosensory sensillum (medial)
  width.mech.hair.dis = 1.28e-6,    # mean width of an mechanosensory sensillum (distal)
  olf.hair.length = 13.6e-6,        # Mean length of an olfactory sensillum
  mech.hair.length = 77.0e-6,       # Mean length of a mechanosensory sensillum
  mech.hair.angle = (29/360)*2*pi,  # Hair angle of a mechanosensory sensillum
  mean.dist.all = 25.5e-6,          # Mean distance between any two sensilla in the array
  sd.dist.all= 7e-6,              # Standard deviation of distances between sensilla
  frac.olf = 0.59,                  # Fraction of all sensilla that are olfactory
  density.all = 0.0026/(1e-6)^2,    # Mean density of all sensilla 
  width.ant = 490e-6,               # width of antennomere 6
  # Parameters for IBM model: 
  np = 3,                           # Number of points per IBM smallest grid
  # Notes ~ Rule of thumb: 2 boundary points per fluid grid point. 
  #        vertex pts too far apart: flow thru boundary, too close: numerical weirdness
  L = 400e-6,                       # Length of computational domain (m)
  N = 1024,                         # Number of Cartesian grid meshwidths at the finest level of the AMR grid
  length.strip = 10e-6,              # small strip
  base.width = 10e-6
)

# Calculated Parameters for 2D model
parameters$area.strip <- parameters$length.strip * parameters$width.ant
parameters$dx <- (2.0 * parameters$L) / (parameters$N * parameters$np)            # Cartesian mesh width (m)

parameters$total.hairs <- ceiling(parameters$density.all * parameters$area.strip)
parameters$num.olf.hairs <- ceiling(parameters$frac.olf * parameters$total.hairs)
parameters$num.mech.hairs <- parameters$total.hairs - parameters$num.olf.hairs
parameters$overlap <- floor((parameters$mech.hair.length * cos(parameters$mech.hair.angle)) / parameters$mean.dist.all)

print(paste("Creating model for", parameters$Species, "which will have", parameters$num.olf.hairs, "olfactory hairs."))

dir.create(paste0("./data/vertex-files/",parameters$Species,"/"), recursive = TRUE, showWarnings = FALSE)
dir.create(paste0("./data/csv-files/",parameters$Species,"/"), recursive = TRUE, showWarnings = FALSE)

#### Make the model!####
set.seed(n)
make.model(parameters, reps, F)

