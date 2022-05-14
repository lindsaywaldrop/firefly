library(R.matlab)
library(tidyr)
library(ggplot2)
library(viridis)

init.dat <- readMat("./results/odorcapture/Lucidota_sp/initdata_0001.mat")
dat <- readMat("./results/odorcapture/Lucidota_sp/hairs_c_0001.mat")
cdat <- readMat("./results/odorcapture/Lucidota_sp/c_0001.mat")

dt <- init.dat$dt
steps.number <- length(dat) 
hairs.number <- length(dat[[1]])
conc.data <- matrix(NA, ncol = hairs.number, nrow = steps.number)
colnames(conc.data) <- paste("hair", as.character(1:hairs.number), sep = "")
rownames(conc.data) <- as.character(1:steps.number)
for (i in 1:steps.number){
  csum <- sum(sum(cdat[[i]]))
  for (j in 1:hairs.number){
    a <- dat[[paste("hairs.c.", i, sep = "")]][[j]][[1]]
    if (length(a) == 0) { conc.data[i, j] <- 0 }
    else {conc.data[i, j] <- sum(a)/csum/cmax}
  }
}

conc.df <- as.data.frame(conc.data)
conc.df$time <- 0:(nrow(conc.data)-1)*dt*500
conc.df.longer <- pivot_longer(conc.df, cols=-c("time"), names_to = "hair", values_to = "Concentration")

ggplot(conc.df.longer, aes(time, Concentration, color = hair)) + 
  geom_point() + geom_line() +
  scale_color_viridis(discrete = TRUE) +
  theme_minimal() 



cdat.origin <- as.data.frame(cdat$c.20)
cmax <- max(cdat$c.1)
cdat.origin$x <- 1:nrow(cdat.origin)
cdat.origin.long <- pivot_longer(as.data.frame(cdat.origin), cols=-"x", names_prefix = "V", names_to="y")
cdat.origin.long$y <- as.numeric(cdat.origin.long$y)
ggplot(cdat.origin.long, aes(x,y, fill=value))+geom_tile() +
  scale_fill_viridis(option="A", limits=c(-2e-4,cmax)) +theme_minimal()
