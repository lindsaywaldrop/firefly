#### Required functions for running analysis on concentration data ####

#### Analysis Functions ####

load.specieslist <- function(){
  if(file.exists("../data/parameters/Specimen-list.csv")){
    specimen.list <- read.csv("../data/parameters/Specimen-list.csv") 
  } else {
    specimen.list <- read.csv("./data/parameters/Specimen-list.csv") 
  }
  drop.these <- which(specimen.list$Measured != "yes")
  specimen.list <- specimen.list[-drop.these,]
  specimen.list$Species <- sub(" ", "_", specimen.list$Species)
  specimen.list$Species <- sub("[.]","",specimen.list$Species)
  return(specimen.list)
}

retrieve.conc <- function(Species, set.name, run){
  require(R.matlab)
  print.interval <- 10000
  if(file.exists(paste0("./results/odorcapture/",set.name,"/",Species,
                        "/initdata_",sprintf("%04i",run),".mat"))){
    init.dat <- readMat(paste0("./results/odorcapture/",set.name,"/",Species,"/initdata_",sprintf("%04i",run),".mat"))
    dat <- readMat(paste0("./results/odorcapture/",set.name,"/",Species,"/hairs_c_",sprintf("%04i",run),".mat"))
    cdat <- readMat(paste0("./results/odorcapture/",set.name,"/",Species,"/c_",sprintf("%04i",run),".mat"))
  } else {
    init.dat <- readMat(paste0("../results/odorcapture/",set.name,"/",Species,"/initdata_",sprintf("%04i",run),".mat"))
    dat <- readMat(paste0("../results/odorcapture/",set.name,"/",Species,"/hairs_c_",sprintf("%04i",run),".mat"))
    cdat <- readMat(paste0("../results/odorcapture/",set.name,"/",Species,"/c_",sprintf("%04i",run),".mat"))
  }
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

diditrun <- function(Species, set.name, run){
  if(file.exists(paste0("./results/odorcapture/",set.name,"/",Species,"/initdata_",sprintf("%04i",run),".mat"))){
    a <- TRUE
  } else if(file.exists(paste0("../results/odorcapture/",set.name,"/",Species,"/initdata_",sprintf("%04i",run),".mat"))) {
    a <- TRUE
  }else {
    a <- FALSE
  }
  return(a)
}

find.final.sum <- function(df, rmna = TRUE){
  df<-subset(df, select= -time)
  df <-df[nrow(df),]
  a <- ifelse(df < 1e-4, 0, df)
  a <- as.data.frame(a)
  b <- sum(a, na.rm = rmna)
  return(b)
}

capture.analysis <- function(set.name, no.runs){
  
  specimen.list <- load.specieslist()
  final.conc <- matrix(NA, ncol = length(specimen.list$Species), nrow = no.runs)
  colnames(final.conc) <- specimen.list$Species
  numofhairs <- matrix(NA, ncol = length(specimen.list$Species), nrow = no.runs)
  colnames(numofhairs) <- specimen.list$Species
  
  for (pep in 1:length(specimen.list$Species)){
    print(specimen.list$Species[pep])
    for (j in 1:no.runs){
      if(diditrun(specimen.list$Species[pep], set.name, j)) {
        conc.df <- retrieve.conc(specimen.list$Species[pep], set.name, j)
        numofhairs[j, pep] <- ncol(conc.df) - 1
        final.conc[j, pep] <- find.final.sum(conc.df)
      } else {
        print(paste0("not run:", specimen.list$Species[pep], " # ", j))
        numofhairs[j, pep] <- 0
        final.conc[j, pep] <- 0
      }
    }
  }
  final.conc2 <- as.data.frame(final.conc)
  final.conc2$run <- seq(1, no.runs)
  numofhairs <- as.data.frame(numofhairs)
  num.hairs <- pivot_longer(numofhairs, cols=everything(), names_to = "Species", 
                            values_to = "NumofHairs")
  finalconc.long <- pivot_longer(final.conc2, cols = -"run", names_to = "Species", 
                                 values_to = "Concentration")
  finalconc.long$Signal <- rep(NA, nrow(finalconc.long))
  for (u in 1:length(specimen.list$Species)){
    for (k in 1:nrow(finalconc.long)){
      if(finalconc.long$Species[k] == specimen.list$Species[u]){
        finalconc.long$Signal[k] <- as.character(specimen.list$Signal[u])
      }
    }
  }
  finalconc.long$Signal <- factor(finalconc.long$Signal)
  finalconc.long$NumofHairs <- num.hairs$NumofHairs
  finalconc.long <- finalconc.long[ , c(1, 2, 4, 5, 3)]
  write.csv(finalconc.long, paste0("../results/r-csv-files/finalconc-long_", 
                                   set.name, "_", no.runs, ".csv"), 
            row.names = FALSE)
  
  perhair <- as.data.frame(final.conc / numofhairs)
  perhair$run <- seq(1, no.runs)
  perhair.long <- pivot_longer(perhair, cols = -"run", names_to = "Species", 
                               values_to = "Concentration")
  perhair.long$Signal <- rep(NA, nrow(perhair.long))
  for (u in 1:length(specimen.list$Species)){
    for (k in 1:nrow(perhair.long)){
      if(perhair.long$Species[k] == specimen.list$Species[u]){
        perhair.long$Signal[k] <- as.character(specimen.list$Signal[u])
      }
    }
  }
  perhair.long$Signal <- factor(perhair.long$Signal)
  perhair.long <- perhair.long[, c(1, 2, 4, 3)]
  
  write.csv(perhair.long, paste0("../results/r-csv-files/perhair-long_", set.name, 
                                 "_", no.runs, ".csv"), 
            row.names = FALSE)
  return(NULL)
}

make.sums <- function(df){
  df$Species <- factor(df$Species)
  x.sums <- rep(NA, nlevels(df$Species))
  x.signals <- rep(NA, nlevels(df$Species))
  for (i in 1:nlevels(df$Species)){
    x <- subset(df, df$Species==levels(df$Species)[i])
    x.sums[i] <- sum(x$Concentration)
    x.signals[i] <- x$Signal[1]
  }
  return(data.frame("Species"= levels(df$Species),
                    "Signal" = x.signals,
                    "Concentration" = x.sums))
}

prep.phyloanova <- function(df, Species.list){
  signals <- df$Signal
  signals[signals=="both"] <- "chemical"
  signals <- factor(signals)
  names(signals) <- Species.list
  sig.sums <- df$Concentration
  names(sig.sums) <- Species.list
  return(list(signals, sig.sums))
}

#### Plotting Functions ####

magnitude<-function(u, v) sqrt(u^2 + v^2)

plot.concvtime <- function(df, yrange = NULL){
  require(tidyr)
  require(ggplot2)
  require(viridis)
  
  conc.df.longer <- pivot_longer(df, cols = -c("time"), names_to = "hair", 
                                 values_to = "Concentration")
  if(is.null(yrange)) yrange <- range(conc.df.longer$Concentration)
  p <- ggplot(conc.df.longer, aes(time, Concentration, color = hair)) + 
        geom_point() + geom_line() +
        scale_color_viridis(discrete = TRUE) +
        ylim(yrange) +
        xlab("time (s)") + 
        theme_minimal()
  return(p)
}

plot.concmap <- function(Species, set.name, run, step.number, float = FALSE){
  require(R.matlab)
  require(tidyr)
  require(ggplot2)
  require(viridis)
  if(file.exists(paste0("./results/odorcapture/", set.name, "/", Species, 
                        "/initdata_", sprintf("%04i", run), ".mat"))){
    init.dat <- readMat(paste0("./results/odorcapture/", set.name, "/", Species,
                               "/initdata_", sprintf("%04i", run), ".mat"))
    dots <- read.table(paste0("./data/vertex-files/", Species, "/", Species, 
                              "_", run, ".vertex"), skip = 1)
    cdat <- readMat(paste0("./results/odorcapture/", set.name, "/", Species, 
                           "/c_", sprintf("%04i", run), ".mat"))
  } else if (file.exists(file.exists(paste0("../results/odorcapture/", set.name, 
                                            "/", Species, "/initdata_", sprintf("%04i", run), 
                                            ".mat")))){
    init.dat <- readMat(paste0("../results/odorcapture/", set.name, "/", Species, 
                               "/initdata_", sprintf("%04i", run), ".mat"))
    dots <- read.table(paste0("../data/vertex-files/", Species, "/", Species, 
                              "_", run, ".vertex"), skip = 1)
    cdat <- readMat(paste0("../results/odorcapture/", set.name, "/", Species, 
                           "/c_", sprintf("%04i", run), ".mat"))
  } else {
    stop("Correct files not found.")
  }
  cdat.origin <- as.data.frame(cdat[[step.number]])
  cdat.origin[cdat.origin < 0] <- 0
  if(float){
    cmax <- max(cdat.origin)
    cmin <- 1e-20
  } else {
    cmax <- max(cdat$c.1)
    cmin <- 1e-20
  }
  colnames(cdat.origin) <- init.dat$y
  cdat.origin$x <- init.dat$x
  cdat.origin.long <- pivot_longer(as.data.frame(cdat.origin), cols = -"x", 
                                   names_prefix = "V", names_to = "y")
  cdat.origin.long$y <- as.numeric(cdat.origin.long$y)
  cdat.origin.long$value[cdat.origin.long$value<1e-16] <- 1e-16
  p <- ggplot(cdat.origin.long, aes(x, y, fill = value)) + geom_raster() +
        geom_point(data = dots, mapping = aes(x = V1, y = V2), 
                   shape = 19, color = "white", fill = NA, size=1) +
        scale_fill_viridis(option = "A", name = "Conc", limits = c(cmin, cmax), 
                           trans="sqrt") + 
        ylab("y (m)") + xlab ("x (m)") + 
        theme_minimal() 
  return(p)
}

plot.velmap <- function(Species, set.name, run, component = "Um"){
  require(R.matlab)
  require(ggplot2)
  require(viridis)
  if(file.exists(paste0("./results/odorcapture/", set.name, "/", Species, 
                        "/initdata_", sprintf("%04i", run), ".mat"))){
    init.dat <- readMat(paste0("./results/odorcapture/", set.name, "/", Species, 
                               "/initdata_", sprintf("%04i", run), ".mat"))
    vel.dat <- readMat(paste0("./results/odorcapture/", set.name, "/", Species, 
                              "/velocity_", sprintf("%04i", run), ".mat"))
    dots <- read.table(paste0("./data/vertex-files/", Species, "/", Species, "_", 
                              run, ".vertex"), skip = 1)
  } else if (file.exists(paste0("../results/odorcapture/", set.name, "/", Species, 
                                "/initdata_", sprintf("%04i", run), ".mat"))){
    init.dat <- readMat(paste0("../results/odorcapture/", set.name, "/", Species, 
                               "/initdata_", sprintf("%04i", run), ".mat"))
    vel.dat <- readMat(paste0("../results/odorcapture/", set.name, "/", Species, 
                              "/velocity_", sprintf("%04i", run), ".mat"))
    dots <- read.table(paste0("../data/vertex-files/", Species, "/", Species, "_", 
                              run, ".vertex"), skip = 1)
  } else {
    stop("Correct files not found.")
  }
  uvel.origin <- as.data.frame(vel.dat$u.flick)
  vvel.origin <- as.data.frame(vel.dat$v.flick)
  if(component == "Um") {
    value <- magnitude(uvel.origin, vvel.origin)
  } else if (component == "Ux"){
    value <- uvel.origin
  } else if (component == "Uy"){
    value <- vvel.origin
  } else {
    stop("Velocity component unknown, please select Um, Ux, or Uy!")
  }
  colnames(value) <- init.dat$y
  value$x <- init.dat$x
  value.long <- pivot_longer(value, cols = -"x", names_to = "y")
  value.long$y <- as.numeric(value.long$y)
  p <- ggplot(value.long, aes(x, y, fill = value)) + 
        geom_tile() +
        scale_fill_viridis(name = component) + 
        geom_point(data = dots, mapping = aes(V1, V2), color = "white", fill = NA) +
        theme_minimal()
  return(p)
}







