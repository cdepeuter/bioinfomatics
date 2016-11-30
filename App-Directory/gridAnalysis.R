gridSearch = TRUE
grid.overlap <- seq(from=10, to=20, by=2)
grid.intervals <- seq(from=20, to=40, by=5)
grid.bins <- seq(from=10, to=20, by=2)
grid.results <- data.frame()
#build data frame
#loop through all grid values, set them to mapper values, source analysis, add results to dataframe

for(o in grid.overlap){
  overlap = o
  for(i in grid.intervals){
    intervals = i
    for(b in grid.bins){
      bins = b
      print(paste("Doing mapper analysis with (overlap, intervals, bins)", toString(o), toString(i), toString(b), sep="-"))
      source("./analysis.R")
      
      results <- c(o, i, b, m1$num_vertices, m1.bhi, hclusters.bhi, kclust.bhiDiffGenes)
      grid.results <- rbind(grid.results, results)
    }
  }
}
colnames(grid.results) <- c("overlap", "intervals", "bins", "num verticies", "mapperBHI", "hclustBHI", "kclustBHI")

gridSearch = FALSE