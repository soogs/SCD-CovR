
# 1. conditions ####
conditions <- list(components = c(2, 3),
                   dimension = c("low", "high"),
                   vafx = c(0.9, 0.5),
                   vafy = c(0.9, 0.5),
                   weak = c("common", "distinctive"),
                   relevant = c("common", "distinctive"),
                   reps = c(1:2))

condition_df <- data.frame(components = NA,
                           dimension = NA,
                           vafx = NA,
                           vafy = NA,
                           weak = NA,
                           relevant = NA,
                           reps = NA)

counts <- 0
for (dimensionz in 1:2){
  for (componentsz in  1:2){
    for (vafxz in 1:2){
      for (vafyz in 1:2){
        for (weakz in 1:2){
          for (relevantz in 1:2){
            for (repsz in 1:2){
              
              counts <- counts + 1
              
              component_now <- conditions$components[componentsz]
              dimension_now <- conditions$dimension[dimensionz]
              vafx_now <- conditions$vafx[vafxz]
              vafy_now <- conditions$vafy[vafyz]
              weak_now <- conditions$weak[weakz]
              relevant_now <- conditions$relevant[relevantz]
              reps_now <- conditions$reps[repsz]
              
              condition_df[counts,] <- c(component_now,
                                         dimension_now,
                                         vafx_now,
                                         vafy_now,
                                         weak_now,
                                         relevant_now,
                                         reps_now)
              
              print(counts) 
            }
          }
        }
      }
    }
  }
}

# dimension (2) * error x (2) * error y(2) * weak component (2) * relevant component (2)
# = 32 conditions 
# analyzing it twice, leading to 64 conditions