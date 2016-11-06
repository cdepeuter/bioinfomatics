numcenters=m1$num_vertices
kclust=kmeans(fil3,centers=numcenters,iter.max=20)
bhikclust=BHI(kclust$cluster,annotation="moe430a.db",names=names(kclust$cluster),category="all")
