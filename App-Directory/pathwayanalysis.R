#gonna explain what each variable is, feel free to change names as you please
uni2 <- as.matrix(read.delim(file = "UniProt2Reactome.txt", header = F, sep = "\t")); #reading in pathway data
uni2 = tbl_df(uni2)
uni3 = count(uni2,V2)
uni3 = uni3[((uni3$n>5) & (uni3$n<500)),] #filtering out the pathways with unacceptable occurences
l3 = list() #this list shall be of the form where the first element in each sublist is the pathway,
#and all the following elements are the genes in that pathway

for(i in 1:length(uni3$V2)){
  l3[[i]]=c(as.character(uni3$V2[i]),unlist(list(as.character(uni2$V1[uni2$V2 %in% uni3$V2[i]]))))
}
#reading in uniprot data
u133a <- as.matrix(read.delim(file = "uniprot2U133A.txt", header = T, sep = "\t"));
u133a = tbl_df(u133a)

#function to map gene names to uniprot names
genetouniprot <- function(geneid){return (u133a$Uniprot_ID[which(u133a$U133A_ID %in% geneid)])}

#changing diff gene names to their respective uniprot names
diff_genes_uniprot <- unlist(lapply(diff_genes,genetouniprot))

#for BIOLOGICAL FUNCTIONS , just needs to be in the same form as l3

#function to calculate number of differentially expressed genes in each pathway
numdiffgenesinpathway = function(pathwaygenes){return(length(which(pathwaygenes %in% diff_genes_uniprot)))}

pathwaydiffcount=list()#this list has the number of differentially expressed genes for the pathway at the same index
for(i in 1:length(l3)){
  pathwaydiffcount[[i]]=numdiffgenesinpathway(l3[[i]][-1])
}

pathwaydiffcount=unlist(pathwaydiffcount)
contingency=matrix(nrow=2,ncol=2);
#function to execute fisher test for each pathway
getpvaluefromfisher=function(i){
contingency[1,1]=pathwaydiffcount[i];
contingency[1,2]=length(l3[[i]][-1])-length(contingency[1,1]);
contingency[2,1]=length(diff_genes)-length(contingency[1,1]);
contingency[2,2]=length(fil3[,1])-length(diff_genes)-length(contingency[1,2]);
return(fisher.test(as.matrix(contingency),alternative = "greater"))
}

pvaluefromfisher=lapply(seq(from=1,to=2223),getpvaluefromfisher)

#function to find false discovery rates
deg_pathway_fdr_func <-function(i){return(p.adjust(pvaluefromfisher[[i]]$p.value,method="fdr"))}
deg_pathway_fdr=unlist(lapply((seq(from=1,to=2223)),deg_pathway_fdr_func))

#get the pvalue for each pathway
pvalue=vector()
for(i in 1:length(pvaluefromfisher)){pvalue[i]=pvaluefromfisher[[i]]$p.value}

#construct the resultant matrix which is of the form pathwayname,fdr,pvalue and order it by increasing p-values
result_matrix=cbind(uni3$V2,pvalue,deg_pathway_fdr)
colnames(result_matrix)=c('Pathway','P-value','FDR')
deg_order=order(pvalue)
result_matrix=result_matrix[deg_order,]
