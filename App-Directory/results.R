t1 <- read.table("gridResults1439Daisy_poa_4_10p5-pearsonSmallerParams.csv")
t2 <- read.table("gridResults1439Daisy_poa_4_10p5_no_iqr_pearson.csv")
t3 <-read.table("gridResults1439Daisy_poa_8_p6_iqr_75_pearson.csv")
pearsonTable1439 <- rbind(t1,t2,t3)
print(c(mean(pearsonTable1439$mapperBHI), mean(pearsonTable1439$hclustBHI), mean(pearsonTable1439$kclustBHI)))
print(c(sd(pearsonTable1439$mapperBHI), sd(pearsonTable1439$hclustBHI), sd(pearsonTable1439$kclustBHI)))

t4 <-read.table("gridResults1439Daisy_poa_45_9p5-spearman.csv", sep=",")
t5 <- read.table("gridResults1439Daisy_poa_4_10p5_noiqr_spearman.csv")
spearmanTable1439 <- rbind(t4,t5)

allValues<- rbind(pearsonTable1439, pearsonTable1439)
print(c(mean(spearmanTable1439$mapperBHI), mean(spearmanTable1439$hclustBHI), mean(spearmanTable1439$kclustBHI)))
print(c(sd(spearmanTable1439$mapperBHI), sd(spearmanTable1439$hclustBHI), sd(spearmanTable1439$kclustBHI)))

spearmanTable5437 <- read.table("5437_small_params_spearman_poa_8_p25_daisy.csv")
print(c(mean(spearmanTable5437$mapperBHI), mean(spearmanTable5437$hclustBHI), mean(spearmanTable5437$kclustBHI)))
print(c(sd(spearmanTable5437$mapperBHI), sd(spearmanTable5437$hclustBHI), sd(spearmanTable5437$kclustBHI)))


colors<-rainbow(3)
xrange <- range(allValues$num.verticies)
yrange<-range(.2, max(max(allValues$mapperBHI), max(allValues$hclustBHI), max(allValues$kclustBHI)))
sortedVerts <- sort(unique(allValues$num.verticies))
flattenMapperBHI <- unlist(lapply(sortedVerts, function(x){mean(allValues[allValues$num.verticies == x, "mapperBHI"])}))
flattenHclustBHI <- unlist(lapply(sortedVerts, function(x){mean(allValues[allValues$num.verticies == x, "hclustBHI"])}))
flattenKclustBHI <- unlist(lapply(sortedVerts, function(x){mean(allValues[allValues$num.verticies == x, "kclustBHI"])}))
jump <- 5


plot(xrange, yrange, type="n", xlab="Num Vertices", ylab="BHI")
colors<-rainbow(3)



lines(sortedVerts, flattenMapperBHI, col = colors[1])
lines(sortedVerts, flattenHclustBHI, col = colors[2])
lines(sortedVerts, flattenKclustBHI, col = colors[3])

legend('bottomright', c("Mapper", "HClust", "K Clust"), col=colors, lty=c(1,1))



# # lines(smooth, smoothHclustBHI, col = colors[2])
# # lines(smooth, smoothKclustBHI, col = colors[3])
# 
# smooth <- sortedVerts[seq(from = 1, to=length(sortedVerts), by=jump)]
# smoothMapperBHI <- unlist(lapply(seq(from=1, to=length(smooth)), function(x){mean(flattenMapperBHI[x:x+jump])}))
# smoothHclustBHI <- unlist(lapply(seq(from=1, to=length(smooth)), function(x){mean(flattenHclustBHI[x:x+jump])}))
# smoothKclustBHI <- unlist(lapply(seq(from=1, to=length(smooth)),function(x){mean(flattenKclustBHI[x:x+jump])}))