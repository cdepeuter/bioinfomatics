#fishersTest

#x=matrix(nrow=2,ncol=2)
#fisher=list()
#min=1
#for(i in 1:length(pct_diffexp)){
#x[1,1]=length(which(m1$points_in_vertex[[i]] %in% diff_gene_nums))
#print(x[1,1])
#x[1,2]=length(diff_genes)-x[1,1]
#print(x[1,2])
#x[2,1]=length(m1$points_in_vertex[[i]])-x[1,1]
#print(x[2,1])
#x[2,2]=numgenes-length(diff_genes)-x[2,1]
#print(x[2,2])
#fisher[[i]]=fisher.test(x,alternative = "greater")
#if(min>fisher[[i]]$p.value){min=fisher[[i]]$p.value}
#}
#fisherreg=list()
#minreg=1
#for(i in 1:length(pct_diffexp)){
#x[1,1]=length(which(which(toClusterReg==i) %in% diff_gene_nums ))
#print(x[1,1])
#x[1,2]=length(diff_genes)-x[1,1]
#print(x[1,2])
#x[2,1]=length(which(toClusterReg==i))-x[1,1]
#print(x[2,1])
#x[2,2]=numgenes-length(diff_genes)-x[2,1]
#print(x[2,2])
#fisherreg[[i]]=fisher.test(x,alternative = "greater")
#if(minreg>fisherreg[[i]]$p.value){minreg=fisherreg[[i]]$p.value}
#}