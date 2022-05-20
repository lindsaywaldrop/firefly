library(R.matlab)
library(tidyr)
library(ggplot2)
library(viridis)

Species <- "Ellychnia_californica"
run <- 9
L<- 400e-6
steps <- 37

init.dat <- readMat(paste0("./results/odorcapture/",Species,"/initdata_000",run,".mat"))
dat <- readMat(paste0("./results/odorcapture/",Species,"/hairs_c_000",run,".mat"))
cdat <- readMat(paste0("./results/odorcapture/",Species,"/c_000",run,".mat"))
vel.dat <- readMat(paste0("./results/odorcapture/",Species,"/velocity_000",run,".mat"))
#vel.dat <- readMat(paste0("./results/ibamr/",Species,"/viz_IB2d",run,".mat"))
dots <- read.table(paste0("./data/vertex-files/",Species,"/",Species,"_",run,".vertex"), skip = 1)

dt <- init.dat$dt
steps.number <- length(dat) 
hairs.number <- length(dat[[1]])
conc.data <- matrix(NA, ncol = hairs.number, nrow = steps.number)
colnames(conc.data) <- paste("hair", as.character(1:hairs.number), sep = "")
rownames(conc.data) <- as.character(1:steps.number)
final.total <- 0
for (i in 1:steps.number){
  csum <- sum(sum(cdat[[i]]))
  for (j in 1:hairs.number){
    a <- dat[[paste0("hairs.c.", i)]][[j]][[1]]
    if(i==steps) final.total <- final.total + sum(dat[[paste0("hairs.c.", steps)]][[j]][[1]])
    if (length(a) == 0) { 
      conc.data[i, j] <- 0 
      } else {conc.data[i, j] <- sum(a)/csum/cmax}
  }
}



conc.df <- as.data.frame(conc.data)
final.conc9 <- sum(conc.df[nrow(conc.df),])
conc.df$time <- 0:(nrow(conc.data)-1)*dt*10000
conc.df.longer <- pivot_longer(conc.df, cols=-c("time"), names_to = "hair", values_to = "Concentration")

ggplot(conc.df.longer, aes(time, Concentration, color = hair)) + 
  geom_point() + geom_line() +
  scale_color_viridis(discrete = TRUE) +
  theme_minimal() 
ggsave("concbyhair.png")

for (bleb in 1:steps){
#  bleb = steps
  cdat.origin <- as.data.frame(cdat[[bleb]])
  cmax <- max(cdat$c.1)
  colnames(cdat.origin) <- init.dat$y
  cdat.origin$x <- init.dat$x
  cdat.origin.long <- pivot_longer(as.data.frame(cdat.origin), cols=-"x", names_prefix = "V", names_to="y")
  cdat.origin.long$y <- as.numeric(cdat.origin.long$y)
  ggplot(cdat.origin.long, aes(x,y, fill=value)) + geom_tile() +
    geom_point(data=dots, mapping=aes(x=V1,y=V2), shape = ".", color="white", fill=NA) +
    scale_fill_viridis(option="A") +theme_minimal() 
  ggsave(paste0("c_",bleb,".png"),width=9, height=5)
}
