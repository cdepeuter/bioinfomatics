# source("https://bioconductor.org/biocLite.R"); ## try http:// if https:// URLs are not supported
# biocLite("affy");
# biocLite("limma");
# biocLite("clValid")
# biocLite("hopach")
# biocLite("annotate")
# biocLite("mygene")
# biocLite("annotate")
# biocLite("GO.db")
# biocLite("moe430a.db")
# biocLite("hgu133a.db")
# biocLite("genefilter")

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
library(TDAmapper);
library(networkD3);
library(HSAUR);
library(shiny);
library(igraph);
library(scatterplot3d);
library(annotate)
library(moe430a.db)
library(GO.db)
library(Biobase)
library(genefilter)
library(clValid)
library(annotate)
library(hgu133a.db)
library(GO.db)

