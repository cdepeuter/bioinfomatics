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



#do regular gene expression analysis, find differentially expressed genes

#get distance matrix
#maybe use daisy, lab3

#do mapper, what clusters do the differentially expressed genes end up in
 
#if multiple clusters, whys that

#regardless of multiple or single, what genes were not significantly expressed
# but are in the same cluster as the significantly expressed ones any significance to these


#look at the pathways of the important genes
#can we make any insight from those

#use TOPAseq for some topological analysis of the pathways

#how do we evaluate these results?


