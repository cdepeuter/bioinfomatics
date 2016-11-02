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
topFunctionsByCluster <- lapply(lapply(cluster_list, getGeneIdsByCluster), getFunctionsFromGeneList, 4)

#functions from all genes
allFunctions <- getFunctionsFromGeneList(geneIds)

heirarchicalClusterFunctions <- lapply(regularClusteredGenes, getFunctionsFromGeneList)
topHClusterFunctions <-  lapply(regularClusteredGenes, getFunctionsFromGeneList, 4)

#Gold standard genes, going to need to map these to affymetrix ids
#TODO office hours, talk to UN
goldStandardBreastCancer  <-  c('BRCA1', 'BRCA2', 'ATM', 'BARD1', 'BRIP1', 'CASP8', 'CDH1', 'CHEK2', 'CTLA4', 'CYP19A1','FGFR2', 'H19', 'LSP1', 'MAP3K1',
                                'MRE11', 'NBN', 'PALB2', 'PTEN', 'RAD51', 'RAD51C', 'STK11', 'TERT', 'TOX3', 'TP53', 'XRCC2', 'XRCC3')
#goldStandardLungCancer <- c('ATK1', 'ALK', 'BRAF', 'DDR2', 'EGFR', 'ERBB2', 'KRAS', 'MAP2K1', 'NRAS', 'PIK3CA', 'PTEN', 'RET', 'RIT1', 'ROS1')
goldStandardProstateCancer <- c('AR', 'BRCA1', 'BRCA2', 'CD82', 'CDH1', 'CHEK2', 'EHBP1', 'ELAC2', 'EP300', 'EPHB2', 'EZH2', 'FGFR2', 'FGFR4', 'GNMT',
                                'HNF1B', 'HOXB13', 'HPCX', 'IGF2', 'ITGA6', 'KLF6', 'LRP2', 'MAD1L1', 'MED12', 'MSMB', 'MSR1', 'MXI1', 'NBN', 'PCAP', 
                                'PCNT', 'PLXNB1', 'PTEN', 'RNASEL', 'SRD5A2', 'STAT3', 'TGFBR1', 'WRN', 'WT1', 'ZFHX3')

#breastQuery <- queryMany(goldStandardBreastCancer, scopes="symbol", fields=c("uniprot", "ensembl.gene", "reporter"), species="human")
#prostateQuery <- queryMany(goldStandardProstateCancer, scopes="symbol", fields=c("uniprot", "ensembl.gene", "reporter"), species="human")'

#why do these queries cme back with some info #### when done in batches but not one at a time
# 
# thisGene <- featureData[1, "Gene ID"]
# geneSearch <- getGene(thisGene)