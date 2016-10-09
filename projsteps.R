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

gds_set_name <- "GDS5437";
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
diff_gene_table <- topTable(efit, number = nrow(affy_exp), p.value = 0.02);
diff_genes <- rownames(diff_gene_table)
#TODO turn these into acutal gene ids

#filtering out for initial mapper, lets use 200 + diff expressed
df <- tbl_df(affy_exp)
fil<-affy_exp[diff_genes,]
fil2<-affy_exp[sample(nrow(affy_exp), 500),]
fil3 <- rbind(fil, fil2)
count=0
degpos=list()
fil3=t(fil3)
corr=cor(fil3,method="spearman")
for(i in 1:length(corr[,1])){if(rownames(corr)[i] %in% diff_genes){degpos[count]=i;count=count+1}}
degpos=unlist(degpos)

#get distance matrix
#maybe use daisy, lab3

#do mapper, what clusters do the differentially expressed genes end up in

#par(mfrow=c(2,1))
overlap = 10
intervals = 20
bins = 20
#filter = 2*cos(0.5*(1:100)) #what is the filter and why, should capture somethign about data
filter=list()
for(i in 1:545){filter[[i]]=mean(corr[i,])}
filter=unlist(filter)

m1 <- mapper1D(
  distance_matrix = corr,
  filter_values = filter,
  num_intervals = intervals,
  percent_overlap = overlap,
  num_bins_when_clustering = bins)

library(igraph)

g1 <- graph.adjacency(m1$adjacency, mode="undirected")
plot(g1, layout = layout.auto(g1) )


#if multiple clusters, whys that

#regardless of multiple or single, what genes were not significantly expressed
# but are in the same cluster as the significantly expressed ones any significance to these


#look at the pathways of the important genes
#can we make any insight from those

#use TOPAseq for some topological analysis of the pathways

#how do we evaluate these results?
