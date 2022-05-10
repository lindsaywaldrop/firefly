# Lucidota sp 2D model 

#### Load Required Libraries ####
library(pracma)
library(data.table)

#### Parameters ####
filename <- "Lucidota_sp"
# Measurements from morphology on antennomere 6
width.olf.hair <- 9.67e-6         # mean width of an olfactory sensillum
width.mech.hair.prox <- 5.16e-6   # mean width of an mechanosensory sensillum (proximal)
width.mech.hair.med <- 2.82e-6    # mean width of an mechanosensory sensillum (medial)
width.mech.hair.dis <- 1.28e-6    # mean width of an mechanosensory sensillum (distal)
olf.hair.length <- 13.6e-6        # Mean length of an olfactory sensillum
mech.hair.length <- 77.0e-6       # Mean length of a mechanosensory sensillum
mech.hair.angle <- (29/360)*2*pi  # Hair angle of a mechanosensory sensillum
mean.dist.all <- 25.5e-6          # Mean distance between any two sensilla in the array
sd.dist.all <- 7e-6              # Standard deviation of distances between sensilla
frac.olf <- 0.59                  # Fraction of all sensilla that are olfactory
density.all <- 0.0026/(1e-6)^2    # Mean density of all sensilla 
width.ant <- 490e-6               # width of antennomere 6

# Parameters for IBM model: 
np <- 2                           # Number of points per IBM smallest grid
L <- 400e-6                       # Length of computational domain (m)
N <- 2048                         # Number of Cartesian grid meshwidths at the finest level of the AMR grid
dx <- (2.0 * L) / (N * np)        # Cartesian mesh width (m)
# Notes ~ Rule of thumb: 2 boundary points per fluid grid point. 
#        vertex pts too far apart: flow thru boundary, too close: numerical weirdness
NFINEST <- 5                      # NFINEST = 4 corresponds to a uniform grid spacing of h=1/64
kappa_target <- 1.0e-2            # target point penalty spring constant (Newton)

# Parameters for 2D model
length.strip <- 10e-6 # um, small strip
area.strip <- length.strip*width.ant
total.hairs <- ceiling(density.all*area.strip)
num.olf.hairs <- ceiling(frac.olf*total.hairs)
num.mech.hairs <- total.hairs-num.olf.hairs
overlap <- floor((mech.hair.length*cos(mech.hair.angle))/mean.dist.all)
base.width <- 10e-6

#### Functions ####

make.peg <- function(dx, center.point, height, width){
  x_grid <- seq(center.point[1]-0.5*width-1e-6, center.point[1]+0.5*width+1e-6, by = dx)
  y_grid <- seq(center.point[2]-1e-6, center.point[2]+height+0.5*width+1e-6, by = dx)
  whole_grid <- meshgrid(x_grid, y_grid)
  
  THETA <- seq(0, pi, length = 25)
  RHO <- array(1, length(THETA)) * 0.5*width
  nap <- matrix(c(THETA, RHO), nrow = length(THETA), ncol = 2)
  pts <- pracma::pol2cart(nap)
  pts <- as.data.frame(pts)
  pts$y <- pts$y + (height-width) + (dx/2) + center.point[2]
  pts$x <- pts$x + center.point[1]
  pts2<- data.frame("x" = c(center.point[1]-0.5*width, 0.5*width+center.point[1]), 
                    "y" = c(center.point[2]+(dx/2), center.point[2]+dx))
  pts.all <- rbind(pts,pts2)
  
  In <- inpolygon(whole_grid$X, whole_grid$Y, pts.all$x, pts.all$y, boundary = FALSE)
  Xin <- whole_grid$X[In]
  Yin <- whole_grid$Y[In]
  
  circ <- data.frame("x" = Xin, "y" = Yin)
  return(circ)
  
}

make.flying.peg <- function(dx, center.point, width){
  height <- width*1.15
  x_grid <- seq(center.point[1]-0.5*width-2e-6, center.point[1]+0.5*width+2e-6, by = dx)
  y_grid <- seq(center.point[2]-0.5*height-2e-6, center.point[2]+0.5*height+2e-6, by = dx)
  whole_grid <- meshgrid(x_grid, y_grid)
  
  circ.seq <- seq(0, 2*pi, length = 50)
  x.seq <- array(width*cos(circ.seq)) + center.point[1]
  y.seq <- array(height*sin(circ.seq)) + center.point[2]
  pts.all <- data.frame("x"=x.seq, "y"=y.seq)
  
  In <- inpolygon(whole_grid$X, whole_grid$Y, pts.all$x, pts.all$y, boundary = FALSE)
  Xin <- whole_grid$X[In]
  Yin <- whole_grid$Y[In]
  
  circ <- data.frame("x" = Xin, "y" = Yin)
  return(circ)
}

#### Create the vertex file ####

hairs <- 1
while(sum(hairs)!=num.olf.hairs){
  hairs <- runif(total.hairs, 0, 1)
  hairs <- ifelse(hairs<0.5, 1, 0)
}

positions <- runif(1, 0, 50e-6)
for(i in 2:length(hairs)){
  positions <- c(positions, positions[i-1]+rnorm(1, mean.dist.all, sd.dist.all))
}

mid.point.x <- sum(range(positions))/2
positions <- positions-mid.point.x
plot(x=positions, y =rep(0,length(positions)))
points(x=positions[hairs==1], rep(0,length(positions[hairs==1])), pch=19, col="black")
lines(c(-200e-6, 200e-6), c(0, 0))

# Make base plate
x_grid <- seq(-0.5*L, 0.5*L, by = dx)
y_grid <- seq(0, base.width, by = dx)
whole_grid <- meshgrid(x_grid, y_grid)
pegs <- data.frame("x" = as.vector(whole_grid$X), "y" = as.vector(whole_grid$Y))

# Make fixed pegs
length.peg <- ifelse(hairs==1, olf.hair.length, mech.hair.length*0.10)
width.peg <- ifelse(hairs==1,  width.olf.hair, width.mech.hair.prox)
for(j in 1:total.hairs){
  peg.try <- make.peg(dx, c(positions[j],base.width), length.peg[j], width.peg[j])
  pegs <- rbind(pegs, peg.try)
}

# Make flying pegs
for (j in 1:overlap){
  fly.above <- rnorm(num.mech.hairs, mean = ((1/overlap)*1.5*j*mech.hair.length*sin(mech.hair.angle)+base.width),
                     sd = 0.05*((1/overlap)*j*mech.hair.length*sin(mech.hair.angle)+base.width))
  fly.across <- runif(1, -0.5*L, -0.5*L+100e-6)
  for(i in 2:num.mech.hairs){
    fly.across <- c(fly.across, fly.across[i-1]+rnorm(1, 2.5*mean.dist.all, sd.dist.all))
  }
  widths <- ifelse(j==1, width.mech.hair.med, width.mech.hair.dis)
  for (k in 1:num.mech.hairs){
    center.point <- c(fly.across[k], fly.above[k])
    fly.peg <- make.flying.peg(dx, center.point, widths)
    pegs <- rbind(pegs,fly.peg)
  }
  rm(fly.across, fly.above)
}

# Plot pegs
plot(pegs, xlim=c(-0.5*L, 0.5*L), ylim = c(0, 0.5*L))

# Write vertex file
filename <- paste0(filename,".vertex")
if(file.exists(filename)) file.remove(filename)  # Deletes file with that name if it exists
cat(as.character(length(pegs$x)), sep = "\n", file = filename, append = FALSE)
fwrite(pegs, file = filename, append = TRUE, sep = " ", nThread = 2)
