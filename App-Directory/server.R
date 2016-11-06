function(input, output) {
  source("../analysis.R")
  
  # data <- eventReactive(input$go, {
  #   rnorm(input$num)
  # })
  
  
  
  output$mapperPlot <- renderPlot({
    par(mar = c(1, 1, 1, 1))
    plot(g1, layout = layout.auto(g1), 
         main = gds_set_name,
         vertex.color = colorRampPalette(c('blue', 'red'))(length(pct_diffexp))[rank(pct_diffexp)]
    )
  })
  output$hTable <- renderTable(t(htbldata))
  output$diffExpTable <- renderTable(t(tbldata))
  
  ColourScale <- paste(cbind('d3.scale.ordinal().range(',jsColorString,');'))
  
    output$force <- renderForceNetwork({
      forceNetwork(Nodes = MapperNodes, Links = MapperLinks, 
                   Source = "Linksource", Target = "Linktarget",
                   Value = "Linkvalue", NodeID = "Nodename", colourScale = JS(ColourScale),
                   Group = "pctdiffexp", opacity = 1, Nodesize = "Nodesize",
                   linkDistance = 10, charge = -150, legend = TRUE)  
    })
  
  output$hclustPlot <-renderPlot({ plot(dend)})

}