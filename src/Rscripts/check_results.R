rm(CL_data)
t_step <- 1e-9
step_interval <- 100
CL_data <- read.csv("src/matlab/CL_data.csv", header=F)
CL_data$sum.all <- rowSums(CL_data)
CL_data$time <- seq(0,step_interval*(nrow(CL_data)-1), 
                    by=step_interval)*(t_step/(step_interval/10))
plot(sum.all~time, data=CL_data, type="l",
     ylim = range(CL_data$sum.all))
points(CL_data$time, CL_data$sum.all, col="red", pch=19)

