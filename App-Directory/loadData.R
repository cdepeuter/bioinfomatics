debug.debug=TRUE
debug.print <- function(x){if(debug.debug){print(x)}}


gds_set_name <- "GDS5437"
#gds_set_name <- "GDS1439"


#set columns for healthy/sick data
if(gds_set_name == "GDS1439"){
  gstats.healthy <- 1:6
  gstats.sick <- 7:19
  annotation_file <- "hgu133a.db"
  gstats.poapct <- .6
  gstats.poanum <- 8
  gstats.iqr <- .75
  pval <- .01
}else if(gds_set_name == "GDS5437"){
  gstats.healthy <- 10:14
  gstats.sick <- 1:9
  annotation_file <- "moe430a.db"
  gstats.poapct <- .2
  gstats.poanum <- 7.5
  gstats.iqr <- .25
  pval <- .05
}

pval <- 0.01
max_diffexp <- 1000


gds <- getGEO(gds_set_name);
eset <- GDS2eSet(gds,do.log2=TRUE);
affy_exp <- as.matrix(eset);
affy_exp_names <- rownames(affy_exp);
featureData<-fData(eset)

#do pre-filtering
#poverA -> gene expression value above pct, in at least num samples
f1 <- pOverA(gstats.poapct, gstats.poanum)
filterFunction <- filterfun(f1)
mask <- genefilter(affy_exp, filterFunction)
f2 <- function(x){
  IQR(x) > gstats.iqr
}
filterFunction2 <- filterfun(f2)

mask2 <- genefilter(affy_exp, filterFunction2)

affy_fil <- affy_exp[mask &mask2,]
debug.print("Dimension of filtered data")
debug.print(dim(affy_fil))


#do regular gene expression analysis, find differentially expressed genes
disease_state<-as.character(pData(eset)[,2])
combn <- factor(disease_state);
design <- model.matrix(~combn);
fit <- lmFit(affy_fil,design);
efit <- eBayes(fit);

#get differentially expressed genes
diff_gene_table <- topTable(efit, number = max_diffexp ,p.value =  pval);
diff_genes <- rownames(diff_gene_table)
geneIds <- rownames(affy_fil)

debug.print(paste("number diffexp genes ", length(diff_genes)))


#get pairwise distance for mapper, 

count=0

filT = t(affy_fil)

#corr = cor(filT,method="spearman")
#TODO test different methods for correlation/distance
debug.print("Getting distance/correlation")
corr = cor(filT,method="pearson")
#corr = cor(filT,method="spearman")
dsy <- daisy(affy_fil)

#do mapper, what clusters do the differentially expressed genes end up in
#filter function for mapper
benign.disease_state = factor("benign",levels=c("benign","primary","metastatic"))

allpca = princomp(affy_fil)

#principal components for benign data
benign.data = affy_fil[,gstats.healthy]
benign.pca = princomp(benign.data)
benign.weights = loadings(benign.pca)[,1]
reduced_dim_benign = benign.data %*% benign.weights
benign.design = model.matrix(~benign.disease_state)
benign.fit = lmFit(benign.design)
coeff_benign = coefficients(benign.fit)

#principal component analysis for disease data
diseasedata = affy_fil[,gstats.sick]
diseasepca = princomp(diseasedata)
weights_disease = loadings(diseasepca)[,1]
reduced_dim_disease = diseasedata %*% weights_disease
modeled_disease_filter = reduced_dim_disease %*% coeff_benign
filterforfil3 = modeled_disease_filter[geneIds,]

