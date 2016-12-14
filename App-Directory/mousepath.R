uni2 <- as.matrix(read.delim(file = "pathways_summary.txt",skip=1,header=T,sep ="\t"))
uni2 = tbl_df(uni2)
uni2 = uni2[,-1]
uni2 = uni2[,-2]
mappingmouse = AnnotationDbi::select(moe430a.db,keys=diff_genes,"SYMBOL")
genetosymbol = function(geneid){return(mappingmouse$SYMBOL[which(mappingmouse$PROBEID %in% geneid)[1]])}
#filtering out the pathways with unacceptable occurences
uni3 = count(uni2,PW.NAME)
uni3 = uni3[((uni3$n>5) & (uni3$n<500)),] 
#function to map gene names to uniprot names
l3 = list() #this list shall be of the form where the first element in each sublist is the pathway,
#and all the following elements are the genes in that pathway
colnames(uni3)=c("V2","n")
colnames(uni2)=c("V1","V2")
for(i in 1:length(uni3$V2)){
  l3[[i]]=c(as.character(uni3$V2[i]),unlist(list(as.character(uni2$V1[uni2$V2 %in% uni3$V2[i]]))))
}

diff_genes_mouse = unlist(lapply(diff_genes,genetosymbol))
diff_genes_mouse = diff_genes_mouse[-which(is.na(diff_genes_mouse))]
#for BIOLOGICAL FUNCTIONS , just needs to be in the same form as l3

#function to calculate number of differentially expressed genes in each pathway
numdiffgenesinpathway = function(pathwaygenes){return(length(which(pathwaygenes %in% diff_genes_mouse)))}

diffgenesinpathway = function(pathwaygenes){return(unlist(pathwaygenes[which(pathwaygenes %in% diff_genes_mouse)]))}


pathwaydiffcount=list()#this list has the number of differentially expressed genes for the pathway at the same index
for(i in 1:length(l3)){
  pathwaydiffcount[[i]]=numdiffgenesinpathway(l3[[i]][-1])
}

pathwaydiffgenes=list()#this list has the number of differentially expressed genes for the pathway at the same index
for(i in 1:length(l3)){
  pathwaydiffgenes[[i]]=diffgenesinpathway(l3[[i]][-1])
}

pathwaydiffcount=unlist(pathwaydiffcount)
contingency=matrix(nrow=2,ncol=2);
#function to execute fisher test for each pathway
getpvaluefromfisher=function(i){
  contingency[1,1]=pathwaydiffcount[i];
  contingency[1,2]=length(l3[[i]][-1])-length(contingency[1,1]);
  contingency[2,1]=length(diff_genes)-length(contingency[1,1]);
  contingency[2,2]=length(rownames(affy_fil))-length(diff_genes)-length(contingency[1,2]);
  return(fisher.test(as.matrix(contingency),alternative = "greater"))
}

pvaluefromfisher=lapply(seq(from=1,to=144),getpvaluefromfisher)

#function to find false discovery rates
deg_pathway_fdr_func <-function(i){return(p.adjust(pvaluefromfisher[[i]]$p.value,method="fdr"))}
deg_pathway_fdr=unlist(lapply((seq(from=1,to=144)),deg_pathway_fdr_func))

#get the pvalue for each pathway
pvalue=vector()
for(i in 1:length(pvaluefromfisher)){pvalue[i]=pvaluefromfisher[[i]]$p.value}

#construct the resultant matrix which is of the form pathwayname,fdr,pvalue and order it by increasing p-values
result_matrix=cbind(uni3$V2,pvalue,deg_pathway_fdr)
colnames(result_matrix)=c('Pathway','P-value','FDR')
deg_order=order(pvalue)
result_matrix=result_matrix[deg_order,]
