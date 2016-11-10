gds <- getGEO(gds_set_name);
eset <- GDS2eSet(gds,do.log2=TRUE);
affy_exp <- as.matrix(eset);
affy_exp_names <- rownames(affy_exp);
featureData<-fData(eset)