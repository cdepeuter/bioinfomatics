#get gene omnibus
#https://www.ncbi.nlm.nih.gov/gds/?term=metastatic
# need to make sure whatever datasets we get from here have features we can look for (highly/lowly metastatic)
#gds5437 https://www.ncbi.nlm.nih.gov/sites/GDSbrowser?acc=GDS5437
source("https://bioconductor.org/biocLite.R"); ## try http:// if https:// URLs are not supported
biocLite("affy");
biocLite("limma");
biocLite("GEOquery");
library(affy);
library(limma);
library(GEOquery);
library(plyr);
library(dplyr)

#is teh data normalized

gds_set_name <- "GDS1439"
gds <- getGEO(gds_set_name);
# Convert the downloaded dataset into an ExpressionSet object(functions you might need: "GDS2eSet")
eset <- GDS2eSet(gds,do.log2=TRUE);


# Convert ExpressionSet to a matrix(functions you might need: "as.matrix")
affy_exp <- as.matrix(eset);
# Obtain the affymetrix ID of each row, we need to map these ID to Uniprot gene ID in the next step 
affy_exp_names <- rownames(affy_exp);

#do regular gene expression analysis, find differentially expressed genes

disease_state<-as.character(pData(eset)[,2])
combn <- factor(disease_state);
#levels(combn)=c("metastatic (4T1) breast carcinoma","nonmetastatic (67NR) breast carcinoma","naÃ¯ve")
design <- model.matrix(~combn);
fit <- lmFit(affy_exp,design);
efit <- eBayes(fit);
diff_gene_table <- topTable(efit, number = nrow(affy_exp),p.value =  0.000001);
diff_genes <- rownames(diff_gene_table)
#TODO turn these into acutal gene ids

#filtering out for initial mapper, lets use 200 + diff expressed
df <- tbl_df(affy_exp)
fil<-affy_exp[diff_genes,]
fil2<-affy_exp[sample(nrow(affy_exp), 1500),]
fil3 <- unique(rbind(fil, fil2))
geneIds <- rownames(fil3)
count=0
degpos=list()
filT=t(fil3)
corr=cor(filT,method="spearman")
for(i in 1:length(corr[,1])){if(rownames(corr)[i] %in% diff_genes){degpos[count]=i;count=count+1}}
degpos=unlist(degpos)

#get distance matrix
#maybe use daisy, lab3

#do mapper, what clusters do the differentially expressed genes end up in

overlap = 10
intervals = 20
bins = 20


filter=list()
for(i in 1:nrow(fil3)){filter[[i]]=mean(corr[i,])}
filter=unlist(filter)
library(TDAmapper)
m1 <- mapper1D(
  distance_matrix = corr,
  filter_values = filter,
  num_intervals = intervals,
  percent_overlap = overlap,
  num_bins_when_clustering = bins)

library(igraph)

g1 <- graph.adjacency(m1$adjacency, mode="undirected")

#dnt plot yet
#plot(g1, layout = layout.auto(g1) )


#so what genes end up where
getGeneIdByIndex <- function(indx){return(geneIds[indx])}
getGenesByCluster <- function(cluster){genes<-lapply(m1$points_in_level[[cluster]], getGeneIdByIndex); return(unlist(genes))}
toCluster = list()
for(i in 1:max(m1$level_of_vertex)){verts = m1$points_in_level[[i]]; for(p in verts){toCluster[p] = i}}
toVertex = list()
for(i in 1:m1$num_vertices){verts = m1$points_in_vertex[[i]]; for(p in verts){toVertex[p] = i}}
getVertexForGene <- function(gene){indx <- which(geneIds == gene); return(unlist(toVertex[indx]))}
getClusterForGene <- function(gene){indx <- which(geneIds == gene); return(unlist(toCluster[indx]))}
diff_exp_clusters <- unlist(lapply(diff_genes, getClusterForGene))

diff_exp_vertices <- unlist(lapply(diff_genes, getVertexForGene))
print(unlist(diff_exp_clusters))

#for each vertex get proportion diff_exp_genes
diff_gene_nums <- lapply(diff_genes, function(gene){which(geneIds == gene)})
getProp <- function(vert){sum(m1$points_in_vertex[[vert]] %in% diff_gene_nums)/length(m1$points_in_vertex[[vert]])}
pct_diffexp  <- unlist(lapply(seq(from = 1, to = m1$num_vertices), getProp))


#plot
plot(g1, layout = layout.auto(g1), vertex.color=colorRampPalette(c('blue', 'red'))(length(pct_diffexp))[rank(pct_diffexp)])



#getMode
#modeVert <-  names(sort(-table(diff_exp_vertices)))[1]
#modeClust <- names(sort(-table(diff_exp_clusters)))[1]
#howManyVertsInMode = length(which( unlist(diff_exp_vertices) == modeVert))
#how does normal clustering do?


#if multiple clusters, whys that


#regardless of multiple or single, what genes were not significantly expressed
# but are in the same cluster as the significantly expressed ones any significance to these


#look at the pathways of the important genes
#can we make any insight from those

#use TOPAseq for some topological analysis of the pathways

#how do we evaluate these results?
