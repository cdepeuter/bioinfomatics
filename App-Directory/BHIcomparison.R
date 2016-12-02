diffgenesinmapperclust=Genebelongstocluster[which(names(Genebelongstocluster) %in% diff_genes)]
diffgenesinhclust=clusts[which(names(clusts) %in% diff_genes)]
diffgenesinkclust=kclust$cluster[which(names(kclust$cluster) %in% diff_genes)]
bhiofdiffgenesinmapper=BHI(diffgenesinmapperclust, annotation=annotation_file, names=names(diffgenesinmapperclust), category="all")
bhifofdiffgenesinhclust=BHI(diffgenesinhclust,annotation=annotation_file,names=names(diffgenesinhclust),category="all")
bhiofdiffgenesinkclust=BHI(diffgenesinkclust,annotation=annotation_file,names=names(diffgenesinkclust),category="all")
