#how do we evaluate these results?
#BHI starts here

names(Genebelongstocluster)=names(unlist(allClustersGenes))

bhiofmapper = BHI(Genebelongstocluster, annotation="moe430a.db", names=names(Genebelongstocluster), category="all")

bhiclusts = BHI(clusts,annotation="moe430a.db",names=names(clusts),category = "all")


diffgenesinmapperclust=Genebelongstocluster[which(names(Genebelongstocluster) %in% diff_genes)]
diffgenesinhclust=clusts[which(names(clusts) %in% diff_genes)]
diffgenesinkclust=kclust$cluster[which(names(kclust$cluster) %in% diff_genes)]
bhiofdiffgenesinmapper=BHI(diffgenesinmapperclust, annotation="moe430a.db", names=names(diffgenesinmapperclust), category="all")
bhifofdiffgenesinhclust=BHI(diffgenesinhclust,annotation="moe430a.db",names=names(diffgenesinhclust),category="all")
bhiofdiffgenesinkclust=BHI(diffgenesinkclust,annotation="moe430a.db",names=names(diffgenesinkclust),category="all")


#write to table
bhis <- c(bhiofdiffgenesinmapper, bhifofdiffgenesinhclust, bhiofdiffgenesinkclust)
names(bhis) <- c("MAPPER", "H-Clust", "K-means")