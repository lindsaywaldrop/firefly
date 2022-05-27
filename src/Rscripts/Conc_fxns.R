#### Required functions for running analysis on concentration data ####

magnitude<-function(u,v) sqrt(u^2 + v^2)

retrieve.conc <- function(Species, run){
  require(R.matlab)
  print.interval <- 10000
  init.dat <- readMat(paste0("./results/odorcapture/",Species,"/initdata_",sprintf("%04i",run),".mat"))
  dat <- readMat(paste0("./results/odorcapture/",Species,"/hairs_c_",sprintf("%04i",run),".mat"))
  cdat <- readMat(paste0("./results/odorcapture/",Species,"/c_",sprintf("%04i",run),".mat"))
  cmax <- max(cdat$c.1)
  steps.number <- length(dat) 
  hairs.number <- length(dat[[1]])
  conc.data <- matrix(NA, ncol = hairs.number, nrow = steps.number)
  colnames(conc.data) <- paste("hair", as.character(1:hairs.number), sep = "")
  rownames(conc.data) <- as.character(1:steps.number)
  final.total <- 0
  for (i in 1:steps.number){
    #csum <- sum(sum(cdat[[i]]/cmax))
    for (j in 1:hairs.number){
    a <- dat[[paste0("hairs.c.", i)]][[j]][[1]]
    if (length(a) == 0) { 
      conc.data[i, j] <- 0 
    } else {
      conc.data[i, j] <- sum(a)/cmax}
    }
  }
  conc.df <- as.data.frame(conc.data)
  #final.conc9 <- sum(conc.df[nrow(conc.df),])
  conc.df$time <- seq(from=0, by = print.interval*init.dat$dt[1], length.out=nrow(conc.data))
  return(conc.df)
}

diditrun <- function(Species, run){
  if(file.exists(paste0("./results/odorcapture/",Species,"/initdata_",sprintf("%04i",run),".mat"))){
    a <- TRUE
  } else {
    a <- FALSE
  }
  return(a)
}

find.final.sum <- function(df, rmna = TRUE){
  df<-subset(df, select= -time)
  if(ncol(df)==1){
    a <- df[nrow(df),]
  } else{
    a <- rowSums(df[nrow(df),], na.rm = rmna)
  }
  return(a)
}

plot.concvtime <- function(df){
  require(tidyr)
  require(ggplot2)
  require(viridis)
  conc.df.longer <- pivot_longer(df, cols=-c("time"), names_to = "hair", values_to = "Concentration")
  p <- ggplot(conc.df.longer, aes(time, Concentration, color = hair)) + 
        geom_point() + geom_line() +
        scale_color_viridis(discrete = TRUE) +
        theme_minimal()
  return(p)
}

plot.concmap <- function(Species, run, step.number, float=FALSE){
  require(R.matlab)
  require(tidyr)
  require(ggplot2)
  require(viridis)
  init.dat <- readMat(paste0("./results/odorcapture/",Species,"/initdata_",sprintf("%04i",run),".mat"))
  dots <- read.table(paste0("./data/vertex-files/",Species,"/",Species,"_",run,".vertex"), skip = 1)
  cdat <- readMat(paste0("./results/odorcapture/",Species,"/c_",sprintf("%04i",run),".mat"))
  cdat.origin <- as.data.frame(cdat[[step.number]])
  cdat.origin[cdat.origin<0] <-0
  if(float){
    cmax <- max(cdat.origin)
    cmin <- 0
  } else {
    cmax <- max(cdat$c.1)
    cmin <- 0
  }
  colnames(cdat.origin) <- init.dat$y
  cdat.origin$x <- init.dat$x
  cdat.origin.long <- pivot_longer(as.data.frame(cdat.origin), cols=-"x", names_prefix = "V", names_to="y")
  cdat.origin.long$y <- as.numeric(cdat.origin.long$y)
  p <- ggplot(cdat.origin.long, aes(x, y, fill = value)) + geom_tile() +
        geom_point(data = dots, mapping = aes(x = V1, y = V2), 
                   shape = ".", color = "white", fill = NA) +
        scale_fill_viridis(option = "A", name="Conc", limits=c(cmin,cmax)) + 
        theme_minimal() 
  return(p)
}

plot.velmap <- function(Species, run, component = "Um"){
  require(R.matlab)
  require(ggplot2)
  require(viridis)
  init.dat <- readMat(paste0("./results/odorcapture/",Species,"/initdata_",sprintf("%04i",run),".mat"))
  vel.dat <- readMat(paste0("./results/odorcapture/",Species,"/velocity_",sprintf("%04i",run),".mat"))
  dots <- read.table(paste0("./data/vertex-files/",Species,"/",Species,"_",run,".vertex"), skip = 1)
  uvel.origin <- as.data.frame(vel.dat$u.flick)
  vvel.origin <- as.data.frame(vel.dat$v.flick)
  if(component=="Um") {
    value <- magnitude(uvel.origin,vvel.origin)
  } else if (component=="Ux"){
    value <- uvel.origin
  } else if (component=="Uy"){
    value <- vvel.origin
  } else {
    stop("Velocity component unknown, please select Um, Ux, or Uy!")
  }
  colnames(value) <- init.dat$y
  value$x <- init.dat$x
  value.long <- pivot_longer(value, cols=-"x", names_to = "y")
  value.long$y <- as.numeric(value.long$y)
  p <- ggplot(value.long, aes(x, y, fill = value)) + 
        geom_tile() +
        scale_fill_viridis(name=component) + 
        geom_point(data = dots, mapping = aes(V1, V2), color = "white", fill = NA) +
        theme_minimal()
  return(p)
}







