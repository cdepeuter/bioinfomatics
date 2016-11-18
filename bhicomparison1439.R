diffgenesinmapperclust=Genebelongstocluster[which(names(Genebelongstocluster) %in% diff_genes)]
diffgenesinhclust=clusts[which(names(clusts) %in% diff_genes)]
diffgenesinkclust=kclust$cluster[which(names(kclust$cluster) %in% diff_genes)]
bhiofdiffgenesinmapper=BHI(diffgenesinmapperclust, annotation="hgu133a.db", names=names(diffgenesinmapperclust), category="all")
bhifofdiffgenesinhclust=BHI(diffgenesinhclust,annotation="hgu133a.db",names=names(diffgenesinhclust),category="all")
bhiofdiffgenesinkclust=BHI(diffgenesinkclust,annotation="hgu133a.db",names=names(diffgenesinkclust),category="all")
