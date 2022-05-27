source("./src/Rscripts/Conc_fxns.R")

setwd("/Volumes/Isengard/IBAMR/firefly")
Species <- "Bicellonycha_wickershamorum"
run <- 4
step.number <- 2

conc.df <- retrieve.conc(Species, run)
plot.concvtime(conc.df)
find.final.sum(conc.df)

plot.concmap(Species, run, step.number,F)

for (yup in 1:60){
  p <- plot.concmap(Species, run, yup)
  ggsave(paste0("./doc/anim/", Species, "_", run, "_anim-", yup, ".png"), 
         p)
}
