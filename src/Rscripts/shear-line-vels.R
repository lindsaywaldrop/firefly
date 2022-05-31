#################################################################################################################
#################################################################################################################
###
### Velocity across horizontal line calculations
###
#################################################################################################################


#### Parameters ####
#species <- "Pyropyga_nigricans"
#n <- 10

shear.line.vels <- function(species, n){
  # Parameters that should not change
  output <- c("Um_profile", "Ux_profile", "Uy_profile")
  
  #### Main analysis ####
  
  # Checks for and makes new directory for time series data
  dir.create(file.path(paste("./results/r-csv-files/", species, "/", 
                             sep = "")), recursive = TRUE, showWarnings = FALSE)
  for (j in 1:length(output)){
    #print(j)
    for (k in 1:n){
      
      #print(paste("Simulation:", k))
      # Gets file list.
      f <- list.files(path=paste("./results/visit/", species, "/sim", k, "/", sep = ""),
                      pattern = output[j],
                      all.files=FALSE)
      data <- read.table(paste("./results/visit/", 
                               species, "/sim", k, "/", f[1], sep = "")) 	# Loads data for time step
      if(k == 1){
        all.data <- matrix(NA, nrow = nrow(data), ncol=(n+1))
        all.data[,1] <- data$V1
        all.data[,2] <- data$V2
      } else {
        all.data[,k+1] <- data$V2
      }
    }
    
    mean_profile <- rowMeans(all.data[,2:(n+1)])
    y <- all.data[,1]
    plot(mean_profile,y, type="l", col="blue")
    
    save.data <- data.frame("y" = y, "mean_profile" = mean_profile)
    
    #### Checking and Saving Data ####
    #complete<-as.numeric(sum(is.na(Qall2)))
    #message("~.*^*~Completeness check~*^*~.~\n",
    #        "Number of NAs: ",complete)
    #if (complete==0){
    #  message("Set complete. Saving now!")
    write.table(save.data, file = paste("./results/r-csv-files/", species, "/Shear_", 
                                        output[j],"_", n, "_", Sys.Date(), ".csv", sep = ""), 
                sep = ",", row.names = FALSE)
    rm(save.data, data, f, mean_profile, y, all.data)
    #} else {
    #  message("Set not complete, did not save")
    #}
  }
}

