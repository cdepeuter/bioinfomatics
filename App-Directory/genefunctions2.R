names(genesandtheirfunctions) = rownames(affy_fil)


genefunctionindex=function(genefunction){
  l=list()
  for(i in 1:length(genesandtheirfunctions)){
    if(genefunction %in% genesandtheirfunctions[[i]]){
      l = c(l,i)
    }
  }
  return(unlist(l))
}

a12 = lapply(listofgenefunctions,genefunctionindex)

a22 =lapply(a12,function(x){return(rownames(affy_fil)[x])})
names(a22) = listofgenefunctions

numdiffgenesingf = function(gfgenes){return(length(which(gfgenes %in% diff_genes)))}
diffgenesingf = function(gfgenes){return(unlist(gfgenes[which(gfgenes %in% diff_genes)]))}
gfdiffgenes =lapply(a22,diffgenesingf)

gfdiffcount=list()#this list has the number of differentially expressed genes for the pathway at the same index
for(i in 1:length(listofgenefunctions)){
  gfdiffcount[[i]]=numdiffgenesingf(a22[[i]])
}

gfdiffcount = unlist(gfdiffcount)

contingency=matrix(nrow=2,ncol=2);
#function to execute fisher test for each gf
getpvaluefromfisher=function(i){
  contingency[1,1]=gfdiffcount[i];
  contingency[1,2]=length(a22[[i]])-length(contingency[1,1]);
  contingency[2,1]=length(diff_genes)-length(contingency[1,1]);
  contingency[2,2]=length(rownames(affy_fil))-length(diff_genes)-length(contingency[1,2]);
  return(fisher.test(as.matrix(contingency),alternative = "greater"))
}

pvaluefromfishergf=lapply(seq(from=1,to=length(listofgenefunctions)),getpvaluefromfisher)

deg_gf_fdr_func <-function(i){return(p.adjust(pvaluefromfishergf[[i]]$p.value,method="fdr"))}
deg_gf_fdr=unlist(lapply((seq(from=1,to=length(listofgenefunctions))),deg_gf_fdr_func))

#get the pvalue for each pathway
pvaluegf=vector()
for(i in 1:length(pvaluefromfishergf)){pvaluegf[i]=pvaluefromfishergf[[i]]$p.value}

#construct the resultant matrix which is of the form pathwayname,fdr,pvalue and order it by increasing p-values
result_matrix_gf=cbind(listofgenefunctions,pvaluegf,deg_gf_fdr)
colnames(result_matrix_gf)=c('Gene Function','P-value','FDR')
deg_order_gf=order(pvaluegf)
result_matrix_gf=result_matrix_gf[deg_order_gf,]
result_matrix_gf = tbl_df(result_matrix_gf)
result_matrix_gf = result_matrix_gf[order(result_matrix_gf$`P-value`),]


goodfunctions = unlist(result_matrix_gf[1:15,1])
indexofgoodfunctions = which(listofgenefunctions %in% goodfunctions)
goodfunctiongenes = gfdiffgenes[indexofgoodfunctions]


clustersofgoodgfgenes = list()
clustersofgoodgfgenesH = list()
clustersofgoodgfgenesK = list()

for(i in 1:length(goodfunctiongenes)){
  clustersofgoodgfgenes[[i]] = unlist(lapply((goodfunctiongenes[[i]]),getVertexForGene))
  clustersofgoodgfgenesH[[i]] = hclusters.clusts[which(names(hclusters.clusts) %in% goodfunctiongenes[[i]])]
  clustersofgoodgfgenesK[[i]] = kclust$cluster[which(names(kclust$cluster) %in% goodfunctiongenes[[i]])]
}
