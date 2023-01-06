source("./src/Rscripts/Conc_fxns.R")
library(patchwork)
library(tidyr)
set.name <- "water1en11"
Species <- "Photinus_macdermotti"
run <- 5
step.number <- 24

conc.df <- retrieve.conc(Species, set.name, run)
plot.concvtime(conc.df)
ggsave("Bwickrun4-conc.pdf", width=6, height=2)
find.final.sum(conc.df)

plot.concmap(Species, set.name, run, step.number, T)
conc.df.longer <- pivot_longer(conc.df, cols = -c("time"), names_to = "hair", 
                                         values_to = "Concentration")
yrange <- c(-0.0001,max(conc.df.longer$Concentration))

for (yup in 1:nrow(conc.df)){
  print(yup)
  conc.df.sub <- conc.df
  conc.df.sub$hair1[(yup+1):nrow(conc.df)] <- NA 
  conc.df.sub$hair2[(yup+1):nrow(conc.df)] <- NA 
  p1 <- plot.concmap(Species, set.name, run, yup, F)
  p2 <- plot.concvtime(conc.df.sub, yrange)
  p1/p2 + plot_layout(heights=c(2,1))
  ggsave(paste0("./doc/anim/", Species, "_", run, "_anim-", yup, ".png"), 
         last_plot(), width=5, height=5)
}

plot.velmap(Species, set.name, run)
ggsave("Bwickrun4-vel.pdf", width=5, height=4)
