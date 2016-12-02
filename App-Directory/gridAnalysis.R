gridSearch = TRUE
grid.overlap <- seq(from=10, to=24, by=2)
#grid.overlap <- c(12, 14)
#grid.intervals <- c(20,40)
grid.intervals <- seq(from=20, to=50, by=5)
grid.bins <- seq(from=10, to=30, by=5)
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
      try(
        source("./analysis.R")
      )

      
      results <- c(ovlp, intv, mbins, m1$num_vertices, length(diff_genes), m1.bhi, hclusters.bhi, kclust.bhiDiffGenes, m1.bhiDiffGenes, hclusters.bhiDiffGenes, kclust.bhiDiffGenes)
      debug.print(results)
      grid.results <- rbind(grid.results, results)
    }
  }
}
colnames(grid.results) <- c("overlap", "intervals", "bins", "num verticies","numDiffexp", "mapperBHI", "hclustBHI", "kclustBHI", "mapperDiffBHI", "hclustBHIDiff", "kclustBHIDiff")

gridSearch = FALSE