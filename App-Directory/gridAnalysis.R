gridSearch = TRUE
#grid.overlap <- seq(from=10, to=30, by=2)
grid.overlap <- c(5,10,15, 20)
#grid.intervals <- c(18, 20, 22, 24, 26, 28)
grid.intervals <- seq(from=10, to=30, by=2)
grid.bins <- seq(from=10, to=20, by=1)
grid.results <- data.frame()
#build data frame
#loop through all grid values, set them to mapper values, source analysis, add results to dataframe

for(ovlp in grid.overlap){
  m1.overlap = ovlp
  for(intv in grid.intervals){
    m1.intervals = intv
    for(mbins in grid.bins){
      m1.bins = mbins
      
      debug.print(paste("Doing mapper analysis with (overlap, intervals, bins)", ovlp, intv, mbins, sep="-"))
      
      #catch exceptions cuz mapper doesnt like some inputs were gonna give
      try({
        source("./analysis.R");

        results <- c(ovlp, intv, mbins, m1$num_vertices, length(diff_genes), m1.bhi, hclusters.bhi, kclust.bhi, m1.bhiDiffGenes, hclusters.bhiDiffGenes, kclust.bhiDiffGenes)
        grid.results <- rbind(grid.results, results)
        debug.print(results)
      })
     
    }
  }
}
colnames(grid.results) <- c("overlap", "intervals", "bins", "num verticies","numDiffexp", "mapperBHI", "hclustBHI", "kclustBHI", "mapperDiffBHI", "hclustBHIDiff", "kclustBHIDiff")

gridSearch = FALSE

print(mean(grid.results$mapperBHI))
print(mean(grid.results$hclustBHI))
print(mean(grid.results$kclustBHI))
print(mean(grid.results$mapperDiffBHI))
print(mean(grid.results$hclustBHIDiff))
print(mean(grid.results$kclustBHIDiff))