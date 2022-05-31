source("./src/Rscripts/Conc_fxns.R")

set.name <- "water1en11"
Species <- "Pyractomena_angulata"
run <- 20
step.number <- 10

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
