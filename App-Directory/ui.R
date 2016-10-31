fluidPage(
          tags$head(
                 tags$link(
                   rel = "stylesheet",
                   type = "text/css",
                   href = "style.css"
                 )
               ),
  # navlistPanel(
  #   tabPanel("Mapper", 
  #            tags$head(tags$script(src="https://ajax.googleapis.com/ajax/libs/jquery/1.12.4/jquery.min.js")),
  #            tags$head(tags$script(src="cdp.js")),
  #            #includeCSS("www/style.css"),
  #            tags$head(
  #              tags$link(
  #                rel = "stylesheet",
  #                type = "text/css",
  #                href = "style.css"
  #              )
  #            ),
  #            # headerPanel('Store Clusters'),
  #            # sidebarPanel(
  #            #   numericInput('clusters', 'Cluster count', 3,
  #            #                min = 1, max = 9)
  #            # ),
  #            mainPanel(
  #              plotOutput('mapperPlot'),
  #              tags$div(class="clusterInfoTable", 
  #                   tableOutput('table')
  #              )
  # 
  #              #tags$img(src = "/rplt.png")
  #            )
  #   ),
  #   tabPanel("Force Network",
  #      
  #       #tags$img(src = "/rplt.png")
  #   ),
  #   tabPanel("Hierarchical Clusters", 
  #            mainPanel(
  #              plotOutput('hclustPlot')
  #            )
  #   )
  # )
  forceNetworkOutput("force")
  
)