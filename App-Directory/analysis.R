#get gene omnibus
#https://www.ncbi.nlm.nih.gov/gds/?term=metastatic
# need to make sure whatever datasets we get from here have features we can look for (highly/lowly metastatic)
#gds5437 https://www.ncbi.nlm.nih.gov/sites/GDSbrowser?acc=GDS5437
#setwd("~/Documents/columbia/bioinformatics/project/App-Directory")


#mapper inputs
if(!exists("gridSearch") || !gridSearch ){
  m1.overlap = 15
  m1.intervals = 16
  m1.bins = 22
}


debug.print("analysis run")
if(!exists("affy_exp")){
  debug.print("loading data")
  source("./loadData.R")
}


time1 <-proc.time()


m1 <- mapper1D(
  distance_matrix = corr,
  filter_values = filterforfil3,
  num_intervals = m1.intervals,
  percent_overlap = m1.overlap,
  num_bins_when_clustering = m1.bins)

g1 <- graph.adjacency(m1$adjacency, mode="undirected")

cluster_list <- seq(from = 1, to = m1$num_vertices)

time2 <- proc.time()

debug.print(paste("mapper done", time2-time1))
#so what genes end up where, get colors for plot
geneToMapperCluster = list()
for(i in 1:max(m1$level_of_vertex)){
  verts = m1$points_in_level[[i]]; 
  for(p in verts){
    geneToMapperCluster[p] = i;
  }
}

geneToMapperVertex = list()
for(i in 1:m1$num_vertices){
  verts = m1$points_in_vertex[[i]]; 
  for(p in verts){
    geneToMapperVertex[p] = i
  }
}

getVertexForGene <- function(gene){
  indx <- which(geneIds == gene); 
  return(unlist(geneToMapperVertex[indx]))
}
getClusterForGene <- function(gene){
  indx <- which(geneIds == gene);
  return(unlist(geneToMapperCluster[indx]))
}


diff_exp_clusters <- unlist(lapply(diff_genes, getClusterForGene))
diff_exp_vertices <- unlist(lapply(diff_genes, getVertexForGene))

#which genes are differentially expressed by index
diff_gene_nums <- lapply(diff_genes, function(gene){which(geneIds == gene)})

getPropDiffexp <- function(vert){
  diff_expressed <- getNumDiff(vert)
  return(diff_expressed/length(m1$points_in_vertex[[vert]]))
}

getNumDiff<- function(vert){
  return(sum(m1$points_in_vertex[[vert]] %in% diff_gene_nums))
}

pct_diffexp  <- unlist(lapply(cluster_list, getPropDiffexp))
num_diffexp  <- unlist(lapply(cluster_list, getNumDiff))
 

#get actual gene name by geneId, get that by cluster
dtt <- Table(gds)
refToId <- data.frame(as.character(dtt$IDENTIFIER))
rownames(refToId) <- as.character(dtt$ID_REF)
geneIdToName <- function(id){return(as.character(refToId[id, 1]))}
getGeneIdByIndex <- function(indx){return(geneIds[indx])}
getindexofgene<- function(geneinput){return(which(geneIds==geneinput))}
getGeneIdsByMapperCluster <- function(cluster, justDiffexp = FALSE){
  genes <- lapply(m1$points_in_vertex[[cluster]], getGeneIdByIndex);
  if(justDiffexp){
    genes <- genes[unlist(genes) %in% diff_genes]
  }
  return(unlist(genes))
}

getGeneNamesByCluster <-function(cluster){
  return(unlist(lapply(getGeneIdsByMapperCluster(cluster),geneIdToName)))
}

allClustersGenes <- lapply(cluster_list, getGeneNamesByCluster)
totalInMapperClusters <- unlist(lapply(allClustersGenes, length))# how does normal clustering do?
# if multiple clusters, whys that
# regardless of multiple or single, what genes were not significantly expressed
# but are in the same cluster as the significantly expressed ones any significance to these


#normal hierarchical clustering

time3<-proc.time()

debug.print(paste("staring h clust", time3-time2))
hclusters <- hclust(dsy)
dend <- as.dendrogram(hclusters)
#cut tree where clusters = num verticies for mapper
hclusters.clusts <- cutree(hclusters, k=m1$num_vertices)
geneToClusterReg <- list()
for(i in 1:length(hclusters.clusts)){
  geneToClusterReg[i] = as.numeric(hclusters.clusts[i])
}

#number of diffexp genes in regular cluster
getNum_reg <- function(clusternumber){
  inThisCluster <- which( geneToClusterReg == clusternumber );
  return(length(which(inThisCluster %in% diff_gene_nums )));
}

#proportion diffexp genes reg cluster
getPropDiffexp_reg <- function(clusternumber){
  inThisCluster <- which( geneToClusterReg == clusternumber );
  whatsDiffExp <- which( inThisCluster %in% diff_gene_nums );
  return(length(whatsDiffExp)/length(inThisCluster));
}


cluster_list <- seq(from = 1, to = m1$num_vertices)
pct_diffexp_reg  <- unlist(lapply(cluster_list, getPropDiffexp_reg))
num_diffexp_by_reg <-unlist(lapply(cluster_list, getNum_reg))

clusters <- lapply(cluster_list,
                   function(clusternumber){
                     return(which(geneToClusterReg %in% clusternumber))
                   }
            )
hclusters.totalInCluster <- unlist(lapply(clusters, length))
getGenesInCluster <- function(list, justDiffexp = FALSE){
  genes <- unlist(lapply(list, getGeneIdByIndex))
  if(justDiffexp){
    genes <- genes[unlist(genes) %in% diff_genes]
  }
  return(genes)
}

hclusters.regularClusteredGenes <- lapply(clusters, getGenesInCluster)
hclusters.regularClusteredGenesDiffexp <- lapply(clusters, getGenesInCluster, justDiffexp = TRUE)



Genebelongstocluster <- vector()
for(i in 1:m1$num_vertices){
  for(j in 1:length(m1$points_in_vertex[[i]])){
    Genebelongstocluster = c(Genebelongstocluster,i)
  }
}

time4 <- proc.time()
debug.print(paste("stariting k clust", time4-time3))

#do kmeans clustering
source("./kmeans.R")

#do gene function analysis
source("./geneFunctions.R")


#do BHI evaluation
source("./BHIcomparison.R")

#do pathway analysis
source(pathwayFile)

