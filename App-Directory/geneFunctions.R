geneFuncs <- featureData[, "GO:Function ID"]

getFunctionIdsForGene <- function(gene){
  theseFunctions <- as.character(featureData[gene, "GO:Function ID"])
  return(unlist(strsplit(theseFunctions, "///")))
}


getFunctionNamesForGene <- function(gene){
  theseFunctions <- as.character(featureData[gene, "GO:Function"])
  return(unlist(strsplit(theseFunctions, "///")))
}
getFunctionsFromGeneList <- function(genes, max = length(genes)){
  byGenes <- lapply(genes, getFunctionNamesForGene)
  return(sort(table(unlist(byGenes)),decreasing=TRUE)[1:max])
}

stringifyData <- function(rw){
  funcs <- names(rw)
  vals <- unname(rw)
  firstCollapse <- apply(t(rbind(funcs, vals)), 1, paste, collapse = ":")
  return(paste(firstCollapse, collapse = ", "))
}

getFunctionsInAnalysisFormat <- function(genes, max = length(genes)){
  byGenes <- lapply(genes, getFunctionNamesForGene)
  return(byGenes)
} 
genesandtheirfunctions = lapply(rownames(affy_fil),getFunctionsInAnalysisFormat)

genesandtheirfunctions = lapply(genesandtheirfunctions,unlist)
listofgenefunctions = unique(unlist(genesandtheirfunctions))


debug.print("GENE FUNCTIONS")

#what are the functions of the diff expressed genes
functionDiffexpressedFrequency <- getFunctionsFromGeneList(diff_genes)
#how many gene functions do we have
#numFunctions <- length(unique(unlist(funcByGene)))

#for mapper clusters, find their top function
m1.clusterFunctions <- lapply(lapply(cluster_list, getGeneIdsByMapperCluster), getFunctionsFromGeneList)
m1.topFunctionsByCluster <- lapply(lapply(cluster_list, getGeneIdsByMapperCluster), getFunctionsFromGeneList, 4)
m1.topDiffexpFunctionsByCluster <- lapply(lapply(cluster_list, getGeneIdsByMapperCluster, justDiffexp=TRUE), getFunctionsFromGeneList, 10)



#functions from all genes
allFunctions <- getFunctionsFromGeneList(geneIds)
allFunctions.total <- sum(allFunctions[complete.cases(allFunctions)])

hclusters.clusterFunctions <- lapply(hclusters.regularClusteredGenes, getFunctionsFromGeneList)
hclusters.topHClusterFunctions <-  lapply(hclusters.regularClusteredGenes, getFunctionsFromGeneList, 4)
hclusters.topHClusterDiffexpFunctions <- lapply(hclusters.regularClusteredGenesDiffexp, getFunctionsFromGeneList, 10)
#Gold standard genes, going to need to map these to affymetrix ids

#write table for UI
m1.topFunctionString <- lapply(m1.topFunctionsByCluster, stringifyData)
m1.topDiffexpFunctionString <- lapply(m1.topDiffexpFunctionsByCluster, stringifyData)
m1.tableData <- rbind(as.integer(cluster_list), pct_diffexp, num_diffexp, totalInMapperClusters, m1.topFunctionString, m1.topDiffexpFunctionString )
rownames(m1.tableData) <- c("cluster", "% diffexp", "num_diffexp", "total", "functions", "Diffexp functions");

hclusters.topHClustFunctionString <- lapply(hclusters.topHClusterFunctions, stringifyData)
hclusters.topHClustDiffexpFunctionString <- lapply(hclusters.topHClusterDiffexpFunctions, stringifyData)
hclusters.tableData <- rbind(as.integer(cluster_list), pct_diffexp_reg, num_diffexp_by_reg, hclusters.totalInCluster, hclusters.topHClustFunctionString, hclusters.topHClustDiffexpFunctionString)


#what are the functions of just the differentially expressed genes
debug.print("building mapper nodes")

MapperNodes <- mapperVertices(m1, geneIds)
MapperLinks <- mapperEdges(m1)
rnk <- round(pct_diffexp*100)
MapperNodes$pctdiffexp <- round(pct_diffexp*100)
unq <-  unique(rnk)
colorRampMap <- colorRampPalette(c('blue', 'red'))(max(unq))[rank(unq)]
jsColorString <- paste(paste("[\"", paste(colorRampMap, collapse="\",\"")), "\"]")

