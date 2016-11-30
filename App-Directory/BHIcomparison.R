#how do we evaluate these results?
#BHI starts here

names(Genebelongstocluster)=names(unlist(allClustersGenes))

m1.bhi = BHI(Genebelongstocluster, annotation="moe430a.db", names=names(Genebelongstocluster), category="all")

hclusters.bhi = BHI(hclusters.clusts,annotation="moe430a.db",names=names(hclusters.clusts),category = "all")


m1.diffGenesInClust = Genebelongstocluster[which(names(Genebelongstocluster) %in% diff_genes)]
hclusters.diffGenesInClusts = hclusters.clusts[which(names(hclusters.clusts) %in% diff_genes)]
kclust.diffGenesInClusts = kclust$cluster[which(names(kclust$cluster) %in% diff_genes)]
m1.bhiDiffGenes = BHI(m1.diffGenesInClust, annotation="moe430a.db", names=names(m1.diffGenesInClust), category="all")
hclusters.bhiDiffGenes = BHI(hclusters.diffGenesInClusts,annotation="moe430a.db",names=names(hclusters.diffGenesInClusts),category="all")
kclust.bhiDiffGenes = BHI(kclust.diffGenesInClusts,annotation="moe430a.db",names=names(kclust.diffGenesInClusts),category="all")


#write to table
bhis <- c(m1.bhiDiffGenes, hclusters.bhiDiffGenes, kclust.bhiDiffGenes)
names(bhis) <- c("MAPPER", "H-Clust", "K-means")