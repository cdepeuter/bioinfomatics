
forceNetwork(Nodes = MapperNodes, Links = MapperLinks, 
             Source = "Linksource", Target = "Linktarget",
             Value = "Linkvalue", NodeID = "Nodename",
             Group = "Nodegroup", opacity = 1, 
             linkDistance = 10, charge = -400) 

ColourScale <- 'd3.scale.ordinal()
            .rank(pct_diffexp)/max(rank(pct_diffexp))
.range(["#FF6900", "#694489"]);'
JS(ColourScale)
