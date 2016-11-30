#get gene omnibus
#https://www.ncbi.nlm.nih.gov/gds/?term=metastatic
# need to make sure whatever datasets we get from here have features we can look for (highly/lowly metastatic)
#gds5437 https://www.ncbi.nlm.nih.gov/sites/GDSbrowser?acc=GDS5437
#setwd("~/Documents/columbia/bioinformatics/project/App-Directory")


gds_set_name <- "GDS5437"
#gds_set_name <- "GDS1439"


#set columns for healthy/sick data
if(gds_set_name == "GDS1439"){
  gstats.healthy <- 1:6
  gstats.sick <- 7:19
}else if(gds_set_name == "GDS5437"){
  gstats.healthy <- 10:14
  gstats.sick <- 1:9
}

pval <- 0.02
max_diffexp <- 100
samplesize <-  1500
set.seed(3)

#TODO better way to decide this or just add as inputs on UI
#mapper inputs
if(!gridSearch || !exists("gridSearch")){
  overlap = 12
  intervals = 40
  bins = 15
}



if(!exists("affy_exp")){
  source("./loadData.R")
}

#TODO filter out genes whos expression does not change across samples

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


#get pairwise distance for mapper, 

count=0

filT = t(fil3)

#corr = cor(filT,method="spearman")
#TODO test different methods for correlation/distance
corr = cor(filT,method="pearson")

#doesnt look like these are being used
#degpos = list()
#for(i in 1:length(corr[,1])){if(rownames(corr)[i] %in% diff_genes){degpos[count]=i;count=count+1}}
#degpos=unlist(degpos)

#do mapper, what clusters do the differentially expressed genes end up in

#filter function for mapper
benign.disease_state = factor("benign",levels=c("benign","primary","metastatic"))

allpca = princomp(affy_exp)

#principal components for benign data
benign.data = affy_exp[,gstats.healthy]
benign.pca = princomp(benign.data)
benign.weights = loadings(benign.pca)[,1]
reduced_dim_benign = benign.data %*% benign.weights
benign.design = model.matrix(~benign.disease_state)
benign.fit = lmFit(benign.design)
coeff_benign = coefficients(benign.fit)

#principal component analysis for disease data
diseasedata = affy_exp[,gstats.sick]
diseasepca = princomp(diseasedata)
weights_disease = loadings(diseasepca)[,1]
reduced_dim_disease = diseasedata %*% weights_disease
modeled_disease_filter = reduced_dim_disease %*% coeff_benign
filterforfil3 = modeled_disease_filter[geneIds,]


print(paste("From analysis.R, Doing mapper analysis with (soverlap, intervals, bins)", toString(o), toString(i), toString(b), sep="-"))
m1 <- mapper1D(
  distance_matrix = corr,
  filter_values = filterforfil3,
  num_intervals = intervals,
  percent_overlap = overlap,
  num_bins_when_clustering = bins)

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

getGenesInCluster <- function(list, justDiffexp = FALSE){
  genes <- unlist(lapply(list, getGeneIdByIndex))
  if(justDiffexp){
    genes <- genes[unlist(genes) %in% diff_genes]
  }
  return(genes)
}

hclusters.regularClusteredGenes <- lapply(clusters, getGenesInCluster)
hclusters.regularClusteredGenesDiffexp <- lapply(clusters, getGenesInCluster, justDiffexp = TRUE)
#do gene function analysis
source("./geneFunctions.R")


Genebelongstocluster <- vector()
for(i in 1:m1$num_vertices){
  for(j in 1:length(m1$points_in_vertex[[i]])){
    Genebelongstocluster = c(Genebelongstocluster,i)
  }
}

#do kmeans clustering
source("./kmeans.R")


#do BHI evaluation
source("./BHIcomparison.R")

