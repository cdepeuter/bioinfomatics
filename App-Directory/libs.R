source("https://bioconductor.org/biocLite.R"); ## try http:// if https:// URLs are not supported
biocLite("affy");
biocLite("limma");
biocLite("clValid")
biocLite("hopach")
biocLite("annotate")
biocLite("mygene")
biocLite("annotate")
biocLite("GO.db")
library(TDAmapper);
library(networkD3);
library(HSAUR);
library(shiny);
library(igraph);
library(scatterplot3d);
biocLite("genefilter")

#use unstable TDAMapper for up to date functions
#use proxy for ubuntu bug https://github.com/hadley/devtools/issues/877
#with_config(use_proxy("97.77.104.22", 3128), devtools::install_github("paultpearson/TDAmapper"))

library(affy);
library(limma);
library(GEOquery);
library(mygene);
library(annotate);
library(tidyverse);
library(cluster);
library(annotate)
library(GO.db)
library(Biobase)
library(genefilter)
library(clValid)
library(annotate)
library(GO.db)



