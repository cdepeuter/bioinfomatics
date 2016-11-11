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


#what are the functions of the diff expressed genes
functionDiffexpressedFrequency <- getFunctionsFromGeneList(diff_genes)
#how many gene functions do we have
#numFunctions <- length(unique(unlist(funcByGene)))

#for mapper clusters, find their top function
topFunctionsByCluster <- lapply(lapply(cluster_list, getGeneIdsByMapperCluster), getFunctionsFromGeneList, 4)
topDiffexpFunctionsByCluster <- lapply(lapply(cluster_list, getGeneIdsByMapperCluster, justDiffexp=TRUE), getFunctionsFromGeneList, 10)



#functions from all genes
allFunctions <- getFunctionsFromGeneList(geneIds)

heirarchicalClusterFunctions <- lapply(regularClusteredGenes, getFunctionsFromGeneList)
topHClusterFunctions <-  lapply(regularClusteredGenes, getFunctionsFromGeneList, 4)
topHClusterDiffexpFunctions <- lapply(regularClusteredGenesDiffexp, getFunctionsFromGeneList, 10)
#Gold standard genes, going to need to map these to affymetrix ids

#write table for UI
topFunctionString <- lapply(topFunctionsByCluster, stringifyData)
topDiffexpFunctionString <- lapply(topDiffexpFunctionsByCluster, stringifyData)
tbldata <- rbind(as.integer(cluster_list), pct_diffexp, num_diffexp, totalInMapperClusters, topFunctionString, topDiffexpFunctionString )
rownames(tbldata) <- c("cluster", "% diffexp", "num_diffexp", "total", "functions", "Diffexp functions");

topHClustFunctionString <- lapply(topHClusterFunctions, stringifyData)
topHClustDiffexpFunctionString <- lapply(topHClusterDiffexpFunctions, stringifyData)
htbldata <- rbind(as.integer(cluster_list), pct_diffexp_reg, num_diffexp_by_reg, totalInCluster, topHClustFunctionString, topHClustDiffexpFunctionString)


#what are the functions of just the differentially expressed genes


MapperNodes <- mapperVertices(m1, geneIds)
MapperLinks <- mapperEdges(m1)
rnk <- round(pct_diffexp*100)
MapperNodes$pctdiffexp <- round(pct_diffexp*100)
unq <-  unique(rnk)
colorRampMap <- colorRampPalette(c('blue', 'red'))(max(unq))[rank(unq)]
jsColorString <- paste(paste("[\"", paste(colorRampMap, collapse="\",\"")), "\"]")

