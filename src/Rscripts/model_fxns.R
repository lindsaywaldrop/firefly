# Install required packages
packages <- c("pracma", "data.table")
package.check <- lapply(
  packages,
  FUN <- function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE, repos = "http://cran.us.r-project.org")
      library(x, character.only = TRUE)
    }
  }
)


#### Functions ####

make.peg <- function(dx, center.point, height, width){
  x_grid <- seq(center.point[1] - 0.5 * width - 0.5e-6, 
                center.point[1] + 0.5 * width + 0.5e-6, 
                by = dx)
  y_grid <- seq(center.point[2] - 1e-6, center.point[2] + height + 0.5 * width + 1e-6, 
                by = dx)
  whole_grid <- meshgrid(x_grid, y_grid)
  
  THETA <- seq(0, pi, length = 25)
  RHO <- array(1, length(THETA)) * 0.5 * width
  nap <- matrix(c(THETA, RHO), nrow = length(THETA), ncol = 2)
  pts <- pracma::pol2cart(nap)
  pts <- as.data.frame(pts)
  pts$y <- pts$y + (height - width) + (dx/2) + center.point[2]
  pts$x <- pts$x + center.point[1]
  pts2<- data.frame("x" = c(center.point[1] - 0.5 * width, 0.5 * width + center.point[1]), 
                    "y" = c(center.point[2] + (dx/2), center.point[2] + dx))
  pts.all <- rbind(pts, pts2)
  
  In <- inpolygon(whole_grid$X, whole_grid$Y, pts.all$x, pts.all$y, boundary = FALSE)
  Xin <- whole_grid$X[In]
  Yin <- whole_grid$Y[In]
  circ <- data.frame("x" = Xin, "y" = Yin)
  return(circ)
  
}

make.flying.peg <- function(dx, center.point, width){
  height <- width * 1.15
  x_grid <- seq(center.point[1] - 0.5 * width - 0.5e-6, 
                center.point[1] + 0.5 * width + 0.5e-6, 
                by = dx)
  y_grid <- seq(center.point[2] - 0.5 * height - 0.5e-6, 
                center.point[2] + 0.5 * height + 0.5e-6, 
                by = dx)
  whole_grid <- meshgrid(x_grid, y_grid)
  
  circ.seq <- seq(0, 2 * pi, length = 50)
  x.seq <- array(width * cos(circ.seq)) + center.point[1]
  y.seq <- array(height * sin(circ.seq)) + center.point[2]
  pts.all <- data.frame("x" = x.seq, "y" = y.seq)
  
  In <- inpolygon(whole_grid$X, whole_grid$Y, pts.all$x, pts.all$y, boundary = FALSE)
  Xin <- whole_grid$X[In]
  Yin <- whole_grid$Y[In]
  circ <- data.frame("x" = Xin, "y" = Yin)
  return(circ)
}

write.vertex<-function(Species, pegs, rep){
  # Write vertex file
  require(data.table)
  filename <- paste0("./data/vertex-files/", Species, "/", Species, "_", rep, ".vertex")
  if(file.exists(filename)) file.remove(filename)  # Deletes file with that name if it exists
  cat(as.character(length(pegs$x)), sep = "\n", file = filename, append = FALSE)
  fwrite(pegs, file = filename, append = TRUE, sep = " ", nThread = 2)
}

write.pegids <- function(Species, peg.ids, rep){
  filename <- paste0("./data/csv-files/", Species, "/", Species, "_", rep, ".csv")
  if(file.exists(filename)) file.remove(filename)  # Deletes file with that name if it exists
  write.csv(peg.ids, row.names = FALSE, file = filename)
}

#### Create the model ####

make.model <- function(parameters, rep, plotit = FALSE){
  ## Function that creates 2D model from parameter list. Reps are the repetitions 
  ## of slices through the antennomere. 
  
  # Decide on distribution of olfactory hairs
  hairs <- matrix(0, nrow = rep, ncol = parameters$total.hairs)
  while(sum(hairs) != ceil(parameters$frac.olf * rep * parameters$total.hairs)){
    for(pep in 1:rep) hairs[pep, ] <- rbinom(parameters$total.hairs, c(1, 0), 
                                             c(parameters$frac.olf))
  }
  
  # Decide on linear hair positions
  positions <- matrix(0, nrow = rep, ncol = parameters$total.hairs)
  positions[,1] <- runif(nrow(hairs), -parameters$width.ant, -(parameters$width.ant * (1 - 0.10)))
  for(i in 2:ncol(hairs)){
    positions[,i] <- positions[,i - 1] + rnorm(1, parameters$mean.dist.all, 
                                                       parameters$sd.dist.all)
  }
  
  # Recenter hair array on x=0
  mid.point.x <- sum(range(positions)) / 2
  positions <- positions - mid.point.x
  
  # Plot hairs on x line
  if(plotit == TRUE){
    plot(x = positions, y = rep(0, length(positions)))
    points(x = positions[hairs == 1], rep(0, length(positions[hairs == 1])), 
           pch = 19, col = "red")
    lines(c(-200e-6, 200e-6), c(0, 0))
  }
  
  for (kip in 1:rep){
    
  # Make fixed pegs
  length.peg <- ifelse(hairs[kip,] == 1, parameters$olf.hair.length, 
                       parameters$mech.hair.length * 0.10)
  width.peg <- ifelse(hairs[kip,] == 1,  parameters$width.olf.hair, parameters$width.mech.hair.prox)
  peg.ids <- data.frame("start_id" = rep(NA, sum(hairs[kip,])), 
                        "numPts" = rep(NA, sum(hairs[kip,])))
  pegs <- data.frame()
  count.hairs <- 0
  for(k in 1:2){
    for(j in 1:parameters$total.hairs){
      if(k == 1 & hairs[kip, j] == 1){
        count.hairs <- count.hairs + 1
        peg.ids$start_id[count.hairs] <- length(pegs$x) + 1
        peg.try <- make.peg(parameters$dx, c(positions[kip, j], parameters$base.width), 
                            length.peg[j], width.peg[j])
        pegs <- rbind(pegs, peg.try)
        print(length(pegs$x))
        peg.ids$numPts[count.hairs] <- length(peg.try$x)
      }
      else if(k == 2 & hairs[j] == 0){
        peg.try <- make.peg(parameters$dx, c(positions[kip, j], parameters$base.width), 
                            length.peg[j], width.peg[j])
        pegs <- rbind(pegs, peg.try)
      }
    }
  }
  
  # Make flying pegs
  if(parameters$overlap > 0){
    for (j in 1:parameters$overlap){
      fly.above <- rnorm(parameters$num.mech.hairs, 
                         mean = ((1/parameters$overlap) * 1.5 * j * parameters$mech.hair.length * 
                                   sin(parameters$mech.hair.angle) + parameters$base.width),
                         sd = 0.05 * ((1/parameters$overlap) * j * parameters$mech.hair.length * 
                                        sin(parameters$mech.hair.angle) + parameters$base.width))
      fly.across <- runif(1, -0.5*parameters$width.ant, 
                          -0.5*parameters$width.ant + 0.10 * parameters$width.ant)
      for(i in 2:parameters$num.mech.hairs){
        fly.across <- c(fly.across, fly.across[i - 1] + rnorm(1, 2.5 * parameters$mean.dist.all, parameters$sd.dist.all))
      }
      widths <- ifelse(j == 1, parameters$width.mech.hair.med, parameters$width.mech.hair.dis)
      for (k in 1:parameters$num.mech.hairs){
        center.point <- c(fly.across[k], fly.above[k])
        if(center.point[1] < 0.5 * parameters$width.ant){
          fly.peg <- make.flying.peg(parameters$dx, center.point, widths)
          pegs <- rbind(pegs, fly.peg) 
        }
      }
      rm(fly.across, fly.above)
    }
  }
  
  # Make base plate
  x_grid <- seq(-0.5*parameters$width.ant, 0.5*parameters$width.ant, by = parameters$dx)
  y_grid <- seq(0, parameters$base.width, by = parameters$dx)
  whole_grid <- meshgrid(x_grid, y_grid)
  pegs <- rbind(pegs, data.frame("x" = as.vector(whole_grid$X), "y" = as.vector(whole_grid$Y)))
  
  # Plot pegs
  if(plotit==TRUE){
    plot(pegs, xlim=c(-parameters$L, parameters$L), ylim = c(0, 0.5*parameters$L))
  }
  
  # Write out vertex and csv files with peg info
  write.vertex(parameters$Species, pegs, kip)
  write.pegids(parameters$Species, peg.ids, kip)
  rm(pegs, length.peg, width.peg)
  }
  return(rowSums(hairs))
}

