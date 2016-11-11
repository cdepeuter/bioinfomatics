numcenters <- m1$num_vertices
kclust <- kmeans(fil3, centers=numcenters, iter.max=50)
bhikclust <- BHI(kclust$cluster, annotation="moe430a.db", names=names(kclust$cluster), category="all")


#do pca for top 2

# tbpca <- tbl_df(pca[[2]])
# top2 <- tbpca[, 1:2]
# pcaseach <- as.matrix(features[,4:length(features)]) %*% as.matrix(top2)
# 
# 
# pcaMeans <- kmeans(pcaseach, input$clusters)
# 
# output$plot1 <- renderPlot({
#   par(mar = c(5.1, 4.1, 0, 1))
#   plot(pcaseach,
#        col = clusters()$cluster,
#        pch = 20, cex = 3)
#   points(clusters()$centers, pch = 4, cex = 4, lwd = 4)
# })