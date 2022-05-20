#################################################################################################################
#################################################################################################################
###
### Velocity across horizontal line calculations
###
#################################################################################################################


#### Parameters ####
#species <- "Pyropyga_nigricans"
n <- 10

# Parameters that should not change
output <- c("Um_line20", "Um_line30", "Um_line50")

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
  x <- all.data[,1]
  plot(x,mean_profile, type="l", col="blue")
  
  save.data <- data.frame("x" = x, "mean_profile" = mean_profile)
  
  #### Checking and Saving Data ####
  #complete<-as.numeric(sum(is.na(Qall2)))
  #message("~.*^*~Completeness check~*^*~.~\n",
  #        "Number of NAs: ",complete)
  #if (complete==0){
  #  message("Set complete. Saving now!")
  write.table(save.data, file = paste("./results/r-csv-files/", species, "/Across_", 
                                      output[j],"_", n, "_", Sys.Date(), ".csv", sep = ""), 
              sep = ",", row.names = FALSE)
  rm(save.data, data, f, mean_profile, x, all.data)
  #} else {
  #  message("Set not complete, did not save")
  #}
}

