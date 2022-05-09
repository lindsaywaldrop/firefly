# Generates grid for use with Constraint IB Method IBAMR
# GtD = gap to diameter ratio between hairs
# dist = distance from antennule
# theta = angle of center hair with positive x-axis
#
# To run this R script on Bridges, enter the following commands: 
# - module add R
# - R   ## This will start R!##
# - source('generate_grid2d.R')  ## Follow prompts if installing packages
# - quit()
# - n

#### Loads required packages ####

plotit <- 1
# plot the hairs? yes = 1, no = 0
norows <- 1  
hair.density <- 1.0

mainDir1 <- "./data/vertex-files"
mainDir2 <- "./data/csv-files"
subDir1 <- paste(nohairs, "hair_files", sep = "")
dir.create(file.path(mainDir1, subDir1), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(mainDir2, subDir1), recursive = TRUE, showWarnings = FALSE)

#### Defines functions ####
circle <- function(center, radius, L, dx){
  require(pracma)
  x_grid <- seq(-(radius + 0.01), radius + 0.01, by = dx)
  y_grid <- seq(-(radius + 0.01), radius + 0.01, by = dx)
  whole_grid <- meshgrid(x_grid, y_grid)
  
  THETA <- c(seq(0, 2 * pi, length = 250), 0)
  RHO <- array(1, length(THETA)) * radius
  Z <- array(1, length(RHO)) * 0
  nap <- matrix(c(THETA, RHO), nrow = length(THETA), ncol = 2)
  points <- pracma::pol2cart(nap)
  points <- as.data.frame(points)
  
  In <- inpolygon(whole_grid$X, whole_grid$Y, points$x, points$y, boundary = FALSE)
  Xin <- whole_grid$X[In]
  Yin <- whole_grid$Y[In]
  
  X <- Xin + center[1]
  Y <- Yin + center[2]
  circ <- data.frame(X, Y)
  return(circ)
}

plotahair <- function(hairxCenterx, hairxCentery, hdia, L, dx, no, plotit){
  h1 <- circle(c(hairxCenterx, hairxCentery), 0.5 * hdia, L, dx)
  if(plotit == 1){
    points(hairxCenterx, hairxCentery, pch = 19, cex = 2.5)
    text(hairxCenterx, hairxCentery, labels = no, col = "red")
  }
  return(h1)
}
makehairs <- function( th, GtD, number, nohairs, plotit = 0){
  
  #th=0
  #nohairs <- 25
  np <- 3
  L <- 2.0         # length of computational domain (m)
  N <- 4096        # number of Cartesian grid meshwidths at the finest level of the AMR grid
  dx <- (2.0 * L) / (N * np)  # Cartesian mesh width (m)
  # Notes ~ Rule of thumb: 2 boundary points per fluid grid point. 
  #        vertex pts too far apart: flow thru boundary, too close: numerical weirdness
  NFINEST <- 5  # NFINEST = 4 corresponds to a uniform grid spacing of h=1/64
  kappa_target <- 1.0e-2        # target point penalty spring constant (Newton)
  
  #hdia <- (3/8)*0.01     # Diameter of hair
  #adia <- 10*hdia     # Diameter of flagellum
  hdia <- 0.002     # Diameter of hair
  adia <- 0.1     # Diameter of flagellum
  
  th2 <- (th / 180) * pi      # Angle off positive x-axis in radians
  #GtD <- 5.0      # Gap width to diameter ratio
  dist <- 2 * hdia     # Distance between antennule and hair 
  mindGap <- (0.5 * adia) + (0.5 * hdia) + dist  # Calculate distance between hair centers
  width <- (GtD * hdia) + hdia
  
  #### Calculates center positions (x,y) of each hair ####
  hair1Centerx <- mindGap * cos(th2)
  hair1Centery <- mindGap * sin(th2)
  
  hair2Centerx <- hair1Centerx - width * sin(th2)
  hair2Centery <- hair1Centery + width * cos(th2)
  
  hair3Centerx <- hair1Centerx + width * sin(th2)
  hair3Centery <- hair1Centery - width * cos(th2)
  
  hair4Centerx <- hair1Centerx + width * cos(th2 + (30/180) * pi)
  hair4Centery <- hair1Centery + width * sin(th2 + (30/180) * pi)
  
  hair5Centerx <- hair1Centerx + width * cos(th2 - (30/180) * pi)
  hair5Centery <- hair1Centery - width * sin(-th2 + (30/180) * pi)
  
  hair6Centerx <- hair2Centerx + width * cos(th2 + (30/180) * pi)
  hair6Centery <- hair2Centery + width * sin(th2 + (30/180) * pi)
  
  hair7Centerx <- hair3Centerx + width * cos(th2 - (30/180) * pi)
  hair7Centery <- hair3Centery - width * sin(-th2 + (30/180) * pi)
  
  hair8Centerx <- hair7Centerx + width * cos(th2 - (30/180) * pi)
  hair8Centery <- hair7Centery - width * sin(-th2 + (30/180) * pi)
  
  hair9Centerx <- hair7Centerx + width * cos(th2 + (30/180) * pi)
  hair9Centery <- hair7Centery + width * sin(th2 + (30/180) * pi)
  
  hair10Centerx <- hair5Centerx + width * cos(th2 + (30/180) * pi)
  hair10Centery <- hair5Centery + width * sin(th2 + (30/180) * pi)
  
  hair11Centerx <- hair4Centerx + width * cos(th2 + (30/180) * pi)
  hair11Centery <- hair4Centery + width * sin(th2 + (30/180) * pi)
  
  hair12Centerx <- hair6Centerx + width * cos(th2 + (30/180) * pi)
  hair12Centery <- hair6Centery + width * sin(th2 + (30/180) * pi)
  
  hair13Centerx <- hair8Centerx + width * cos(th2 - (30/180) * pi)
  hair13Centery <- hair8Centery - width * sin(-th2 + (30/180) * pi)
  
  hair14Centerx <- hair9Centerx + width * cos(th2 - (30/180) * pi)
  hair14Centery <- hair9Centery - width * sin(-th2 + (30/180) * pi)
  
  hair15Centerx <- hair10Centerx + width * cos(th2 - (30/180) * pi)
  hair15Centery <- hair10Centery - width * sin(-th2 + (30/180) * pi)
  
  hair16Centerx <- hair11Centerx + width * cos(th2 - (30/180) * pi)
  hair16Centery <- hair11Centery - width * sin(-th2 + (30/180) * pi)
  
  hair17Centerx <- hair12Centerx + width * cos(th2 - (30/180) * pi)
  hair17Centery <- hair12Centery - width * sin(-th2 + (30/180) * pi)
  
  hair18Centerx <- hair12Centerx + width * cos(th2 + (30/180) * pi)
  hair18Centery <- hair12Centery + width * sin(th2 + (30/180) * pi)
  
  hair19Centerx <- hair13Centerx + width * cos(th2 - (30/180) * pi)
  hair19Centery <- hair13Centery - width * sin(-th2 + (30/180) * pi)
  
  hair20Centerx <- hair14Centerx + width * cos(th2 - (30/180) * pi)
  hair20Centery <- hair14Centery - width * sin(-th2 + (30/180) * pi)
  
  hair21Centerx <- hair15Centerx + width * cos(th2 - (30/180) * pi)
  hair21Centery <- hair15Centery - width * sin(-th2 + (30/180) * pi)
  
  hair22Centerx <- hair16Centerx + width * cos(th2 - (30/180) * pi)
  hair22Centery <- hair16Centery - width * sin(-th2 + (30/180) * pi)
  
  hair23Centerx <- hair17Centerx + width * cos(th2 - (30/180) * pi)
  hair23Centery <- hair17Centery - width * sin(-th2 + (30/180) * pi)
  
  hair24Centerx <- hair18Centerx + width * cos(th2 - (30/180) * pi)
  hair24Centery <- hair18Centery - width * sin(-th2 + (30/180) * pi)
  
  hair25Centerx <- hair18Centerx + width * cos(th2 + (30/180) * pi)
  hair25Centery <- hair18Centery + width * sin(th2 + (30/180) * pi)
  
  #### Produces points within defined hairs ####
  
  # Antennule
  ant <- circle(c(0, 0), 0.5 * adia, L, dx);  # Produces points that define antennule
  aN <- size(ant$X, 2)                   # Records number of points inside antennule
  if(plotit == 1){
    plot(0, 0, xlim = c(-0.1, 0.1), ylim = c(-0.1, 0.1), 
         pch = 19, cex = 4.5) #Plots antennule
    text(0, 0, labels = "Ant", col = "red")
  }
  
  for (i in 1:nohairs){
    hairx <- eval(as.name(paste("hair", i, "Centerx", sep = "")))
    hairy <- eval(as.name(paste("hair", i, "Centery", sep = "")))
    h <- plotahair(hairx, hairy, hdia, L, dx, i, plotit)
    assign(paste("h", i, sep = ""), h)
  }
  hN <- size(h$X, 2)
  disp(paste("Array number: ", number, ", number of points per hair: ", hN, sep = ""))
  
  #### Write points to vertex file ####
  
  totalN <- aN + nohairs * hN  # Calculates total number of points (first line of vertex file)
  
  filename <- paste("./data/vertex-files/", nohairs, "hair_files/hairs", 
                    number, ".vertex", sep = "")   # Defines file name
  if(file.exists(filename)) file.remove(filename)  # Deletes file with that name if it exists
  cat(as.character(totalN), sep = "\n", file = filename, append = FALSE)
  # new code with fwrite in data.table package
  fwrite(data.frame(ant), file = filename, append = TRUE, sep = " ", nThread = 2)
  for (k in 1:nohairs){
    hair <- eval(as.name(paste("h", k, sep = "")))
    fwrite(data.frame(hair), file = filename, append = TRUE, sep = " ", nThread = 2)
  }
  
  allhairs <- data.frame("a" = c(aN, 0, 0))
  names(allhairs) <- "a"
  for (i in 1:nohairs){
    hairx <- eval(as.name(paste("hair", i, "Centerx", sep = "")))
    hairy <- eval(as.name(paste("hair", i, "Centery", sep = "")))
    allhairs <- cbind(allhairs, "h" = c(hN, hairx, hairy))
  }
  
  fwrite(allhairs, file = paste("./data/csv-files/", nohairs, 
                                "hair_files/hairs", number, ".csv", sep = ""), 
         row.names = FALSE)
  
}

##### Input parameter definitions ####

#parameters <- read.table(paste("./data/parameters/allpara_", endrun, ".txt", sep = ""))
#names(parameters) <- c("angle", "gap", "Re")

parameters <- data.frame(angle=90, gap=4, Re=10)

for (j in startrun:endrun){
  #ptm <- proc.time()
  makehairs(parameters$angle[j], parameters$gap[j], j, nohairs, plotit)
  #t <- (proc.time() - ptm)
  #message(t)
}
