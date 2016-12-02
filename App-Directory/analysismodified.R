#get gene omnibus
#https://www.ncbi.nlm.nih.gov/gds/?term=metastatic
# need to make sure whatever datasets we get from here have features we can look for (highly/lowly metastatic)
#gds5437 https://www.ncbi.nlm.nih.gov/sites/GDSbrowser?acc=GDS5437


#gds_set_name <- "GDS5437"
gds_set_name <- "GDS1439"


#set columns for healthy/sick data
if(gds_set_name == "GDS1439"){
  gstats.healthy <- 1:6
  gstats.sick <- 7:19
  annotation_file <- "hgu133a.db"
  biocLite("hgu133a.db")
  library(hgu133a.db)
  
}else if(gds_set_name == "GDS5437"){
  gstats.healthy <- 10:14
  gstats.sick <- 1:9
  annotation_file <- "moe430a.db"
  biocLite("moe430a.db")
  library(moe430a.db)
}

pval <- 0.02
max_diffexp <- 50
samplesize <-  1000
set.seed(21)

#TODO better way to decide this or just add as inputs on UI
#mapper inputs
overlap = 20
intervals = 10
bins = 40


if(!exists("affy_exp")){
  source("loadData.R")
}

#TODO check for redundant genes
#http://www3.stat.sinica.edu.tw/statistica/oldpdf/A12n112.pdf

#do regular gene expression analysis, find differentially expressed genes
disease_state<-as.character(pData(eset)[,2])
combn <- factor(disease_state);
design <- model.matrix(~combn);
fit <- lmFit(affy_exp,design);
efit <- eBayes(fit);

#get differentially expressed genes
diff_gene_table <- topTable(efit, number = max_diffexp ,p.value =  pval);
diff_genes <- rownames(diff_gene_table)


#grab some random genes to go with diff expressed, how many do we want
df <- tbl_df(affy_exp)
fil<-affy_exp[diff_genes,]
fil2<-affy_exp[sample(nrow(affy_exp),samplesize),]
fil3 <- unique(rbind(fil, fil2))
geneIds <- rownames(fil3)

#TODO shuffle here make sure everythings still good
#fil3 <- fil3[sample(nrow(fil3)),]

#get pairwise distance for mapper, 
#TODO different methods for correlation/distance
count=0
degpos = list()
filT = t(fil3)
corr = cor(filT,method="spearman")

#doesnt look like these are being used
#for(i in 1:length(corr[,1])){if(rownames(corr)[i] %in% diff_genes){degpos[count]=i;count=count+1}}
#degpos=unlist(degpos)

#do mapper, what clusters do the differentially expressed genes end up in

#filter function for mapper
#TODO generalize disease_state_benign, benignData and diseaseData indices
disease_state_benign = factor("benign",levels=c("benign","primary","metastatic"))

allpca = princomp(affy_exp)

#principal components for benign data
benigndata = affy_exp[,gstats.healthy]
benignpca = princomp(benigndata)
weights_beningn = loadings(benignpca)[,1]
reduced_dim_benign = benigndata %*% weights_beningn 
design_benign = model.matrix(~disease_state_benign)
fit_benign = lmFit(design_benign)
coeff_benign = coefficients(fit_benign)

#principal component analysis for disease data
diseasedata = affy_exp[,gstats.sick]
diseasepca = princomp(diseasedata)
weights_disease = loadings(diseasepca)[,1]
reduced_dim_disease = diseasedata %*% weights_disease
modeled_disease_filter = reduced_dim_disease %*% coeff_benign
filterforfil3 = modeled_disease_filter[geneIds,]

m1 <- mapper1D(
  distance_matrix = corr,
  filter_values = filterforfil3,
  num_intervals = intervals,
  percent_overlap = overlap,
  num_bins_when_clustering = bins)

library(igraph)
g1 <- graph.adjacency(m1$adjacency, mode="undirected")

cluster_list <- seq(from = 1, to = m1$num_vertices)

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


#now that we got vertex % plot
# plot(g1, layout = layout.auto(g1), 
#      vertex.color = colorRampPalette(c('blue', 'red'))(length(pct_diffexp))[rank(pct_diffexp)]
# )


#get actual gene name by geneId, get that by cluster
dtt <- Table(gds)
refToId <- data.frame(as.character(dtt$IDENTIFIER))
rownames(refToId) <- as.character(dtt$ID_REF)
geneIdToName <- function(id){return(as.character(refToId[id, 1]))}
getGeneIdByIndex <- function(indx){return(geneIds[indx])}
getindexofgene<- function(geneinput){return(which(geneIds==geneinput))}
getGeneIdsByCluster <- function(cluster){
  genes <- lapply(m1$points_in_vertex[[cluster]], getGeneIdByIndex);
  return(unlist(genes))
}

getGeneNamesByCluster <-function(cluster){
  return(unlist(lapply(getGeneIdsByCluster(cluster),geneIdToName)))
}

allClustersGenes <- lapply(cluster_list, getGeneNamesByCluster)
totalInMapperClusters <- unlist(lapply(allClustersGenes, length))

#getMode
#modeVert <-  names(sort(-table(diff_exp_vertices)))[1]
#modeClust <- names(sort(-table(diff_exp_clusters)))[1]
#howManyVertsInMode = length(which( unlist(diff_exp_vertices) == modeVert))


# how does normal clustering do?
# if multiple clusters, whys that
# regardless of multiple or single, what genes were not significantly expressed
# but are in the same cluster as the significantly expressed ones any significance to these


#normal hierarchical clustering
dsy <- daisy(fil3)
hclusters <- hclust(dsy)
dend <- as.dendrogram(hclusters)
#cut tree where clusters = num verticies for mapper
clusts <- cutree(hclusters, k=m1$num_vertices)
geneToClusterReg <- list()
for(i in 1:length(clusts)){
  geneToClusterReg[i] = as.numeric(clusts[i])
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
totalInCluster <- unlist(lapply(clusters, length))
clusters <- lapply(cluster_list,
                   function(clusternumber){
                     return(which(geneToClusterReg %in% clusternumber))
                   }
)

getGenesInCluster <- function(list){
  return(unlist(lapply(list, getGeneIdByIndex)))
}


source("../geneFunctions.R")


#write table for UI
topFunctionString <- lapply(topFunctionsByCluster, stringifyData)
tbldata <- rbind(as.integer(cluster_list), pct_diffexp, num_diffexp, totalInMapperClusters, topFunctionString )
rownames(tbldata) <- c("cluster", "% diffexp", "num_diffexp", "total", "functions");

topHClustFunctionString <- lapply(topHClusterFunctions, stringifyData)
regularClusteredGenes <- lapply(clusters, getGenesInCluster)
htbldata <- rbind(as.integer(cluster_list), pct_diffexp_reg, num_diffexp_by_reg, totalInCluster, topHClustFunctionString)



MapperNodes <- mapperVertices(m1, geneIds)
MapperLinks <- mapperEdges(m1)
rnk <- round(pct_diffexp*100)
MapperNodes$pctdiffexp <- round(pct_diffexp*100)
unq <-  unique(rnk)
colorRampMap <- colorRampPalette(c('blue', 'red'))(max(unq))[rank(unq)]
jsColorString <- paste(paste("[\"", paste(colorRampMap, collapse="\",\"")), "\"]")



# 
# map<- mapper(
#   corr,
#   filter_values = filterforfil3,
#   num_intervals = intervals,
#   percent_overlap = overlap,
#   num_bins_when_clustering = bins)
# 
# cluster_list_map <- seq(from = 1, to = map$num_vertices)
# pct_diffexp_map  <- unlist(lapply(cluster_list_map, getPropDiffexp))
# num_diffexp_map  <- unlist(lapply(cluster_list_map, getNumDiff))
# 
# MapperNodes <- mapperVertices(map, geneIds)
# MapperLinks <- mapperEdges(map)
# MapperNodes$pctdiffexp <- pct_diffexp_map
# MapperNodes$rankdiffexp <- rank(pct_diffexp_map)
# unq <-  unique(as.integer(rank(pct_diffexp_map)))



#USE HEAT MAP PLOT, WERE GONNA NEED LESS GENES


# 
# #this was shitty kmeans plot
# prclust = princomp(regClust$centers)
# weights = loadings(prclust)[,1:2]
# actual = regClust$centers %*% weights
# actual = data.frame(actual)
# actual$percent = pct_diffexp_reg
# actual$color[actual$percent <.3 & actual$percent > 0]='red1'
# actual$color[actual$percent <.6 & actual$percent >= .3]='red2'
# actual$color[actual$percent>=.6]='red3'
# actual$color[actual$percent==0]="blue"
# plot(actual[,1:2],col=actual$color,bg=actual$color,pch=21,cex=2)


#look at the pathways of the important genes
#can we make any insight from those

#use TOPAseq for some topological analysis of the pathways

#how do we evaluate these results?
#BHI starts here
source("https://bioconductor.org/biocLite.R")
biocLite("annotate")
biocLite("GO.db")
library(annotate)
library(GO.db)

Genebelongstocluster <-vector()
for(i in 1:m1$num_vertices){
  for(j in 1:length(m1$points_in_vertex[[i]])){
    Genebelongstocluster=c(Genebelongstocluster,i)
  }
}
names(Genebelongstocluster)=names(unlist(allClustersGenes))
if(require("Biobase") && require("annotate") && require("GO.db") &&
   require("hgu133a.db")) {
  bhiofmapper=BHI(Genebelongstocluster, annotation=annotation_file, names=names(Genebelongstocluster), category="all")
}
#bhiclusts= BHI(clusts,annotation="moe430a.db",names=names(clusts),category = "all")
