## Topological Data analysis
The purpose of this project is to compare the clustering of the TDA clustering method Mapper with K-Means and Hierarchical Clustering


### Run shiny app
1. Set working directory to App-Directory of this repo.
2. Source loadData.R
3. Source analysis.R
4. `shiny::runApp()`

### To work with a new GDS file
We have built this code with the idea of being able to insert other gene expression data with as few changes necessary as possible; however, if using a new GDS file a new condition will need to be added after those starting at line 10 in loadData.R, with the necessary variables filled in.

![alt text](https://github.com/cdepeuter/bioinfomatics/blob/master/App-Directory/mapperimgs/gds1437-15-17-17.png "Sample Mapper output")
