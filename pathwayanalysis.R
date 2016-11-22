uni2 <- as.matrix(read.delim(file = "UniProt2Reactome.txt", header = F, sep = "\t")); 
bgGenes <- unique(human_reactome[,1]);
uni2=tbl_df(uni2)
uni3=count(uni2,V2)
uni3=uni3[((uni3$n>5) & (uni3$n<500)),]
l3=list()
for(i in 1:length(uni3$V2)){
  l3[[i]]=c(as.character(uni3$V2[i]),unlist(list(as.character(uni2$V1[uni2$V2 %in% uni3$V2[i]]))))
}
u133a <- as.matrix(read.delim(file = "uniprot2U133A.txt", header = T, sep = "\t"));
u133a=tbl_df(u133a)
genetouniprot <- function(geneid){return (u133a$Uniprot_ID[which(u133a$U133A_ID %in% geneid)])}
diff_genes_uniprot<-unlist(lapply(diff_genes,genetouniprot))
numdiffgenesinpathway=function(pathwaygenes){return(length(which(pathwaygenes %in% diff_genes_uniprot)))}
pathwaydiffcount=list()
for(i in 1:length(l3)){
  pathwaydiffcount[[i]]=numdiffgenesinpathway(l3[[i]][-1])
}
pathwaydiffcount=unlist(pathwaydiffcount)
contingency=matrix(nrow=2,ncol=2);
getpvaluefromfisher=function(i){
contingency[1,1]=pathwaydiffcount[i];
contingency[1,2]=length(l3[[i]][-1])-length(contingency[1,1]);
contingency[2,1]=length(diff_genes)-length(contingency[1,1]);
contingency[2,2]=length(fil3[,1])-length(diff_genes)-length(contingency[1,2]);
return(fisher.test(as.matrix(contingency),alternative = "greater"))
}
pvaluefromfisher=lapply(seq(from=1,to=2223),getpvaluefromfisher)
deg_pathway_fdr_func <-function(i){return(p.adjust(pvaluefromfisher[[i]]$p.value,method="fdr"))}
deg_pathway_fdr=unlist(lapply((seq(from=1,to=2223)),deg_pathway_fdr_func))
pvalue=vector()
for(i in 1:length(pvaluefromfisher)){pvalue[i]=pvaluefromfisher[[i]]$p.value}
result_matrix=cbind(uni3$V2,pvalue,deg_pathway_fdr)
colnames(result_matrix)=c('Pathway','P-value','FDR')
deg_order=order(pvalue)
result_matrix=result_matrix[deg_order,]
