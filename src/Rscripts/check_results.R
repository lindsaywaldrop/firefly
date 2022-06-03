source("./src/Rscripts/Conc_fxns.R")

set.name <- "air1en11"
Species <- "Bicellonycha_wickershamorum"
run <- 11
step.number <- 24

conc.df <- retrieve.conc(Species, set.name, run)
plot.concvtime(conc.df)
find.final.sum(conc.df)

plot.concmap(Species, set.name, run, step.number,T)

for (yup in 1:60){
  p <- plot.concmap(Species, set.name, run, yup, T)
  ggsave(paste0("./doc/anim/", Species, "_", run, "_anim-", yup, ".png"), 
         p)
}

plot.velmap(Species, set.name, run)
