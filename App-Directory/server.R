function(input, output) {
  source("../analysis.R")
  
  # data <- eventReactive(input$go, {
  #   rnorm(input$num)
  # })
  
  
  
  output$mapperPlot <- renderPlot({
    par(mar = c(5.1, 4.1, 0, 1))
    plot(g1, layout = layout.auto(g1), 
         vertex.color = colorRampPalette(c('blue', 'red'))(length(pct_diffexp))[rank(pct_diffexp)]
    )
  })
  output$table <- renderTable(tbldata)
  
  output$hclustPlot <-renderPlot({ plot(dend)})

}