library(TDAmapper)

par(mfrow=c(2,1))
overlap = 50
intervals = 10
bins = 10
filter = 2*cos(0.5*(1:100))
n=100

#triangle, 
#x=c(c(1:100),c(1:100), rep(1,100))
#y=c(c(1:100),rep(100, 100), c(1:100) )


#x+10sin(x)
#x=c(c(0:n))
#y=c(c(10*sin(0:n)+0:n))


#weird random
#x <- c(0:n, n:0)
#y <- c(c(0, cumsum(stats::rnorm(n))), rev(c(0, cumsum(stats::rnorm(n)))))


#fig 8
x=2*cos(0.5*(1:100))
y=sin(1:100)


data <-data.frame(x=x, y=y)
plot(data)

m1 <- mapper1D(
  distance_matrix = dist(data),
  filter_values = filter,
  num_intervals = intervals,
  percent_overlap = overlap,
  num_bins_when_clustering = bins)

library(igraph)

g1 <- graph.adjacency(m1$adjacency, mode="undirected")
plot(g1, layout = layout.auto(g1) )
