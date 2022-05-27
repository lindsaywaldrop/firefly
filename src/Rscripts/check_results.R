source("./src/Rscripts/Conc_fxns.R")

set.name <- "water1en15"
Species <- "Pyractomena_lucifera"
run <- 6
step.number <- 1

conc.df <- retrieve.conc(Species, set.name, run)
plot.concvtime(conc.df)
find.final.sum(conc.df)

plot.concmap(Species, set.name, run, step.number,T)

for (yup in 1:60){
  p <- plot.concmap(Species, run, yup)
  ggsave(paste0("./doc/anim/", Species, "_", run, "_anim-", yup, ".png"), 
         p)
}

plot.velmap(Species, set.name, run)
