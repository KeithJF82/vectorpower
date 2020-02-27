library(shiny)

# Define UI ----
ui <- fluidPage(
  titlePanel("Example 1"),
  sidebarLayout(
    sidebarPanel(
      p("Viewing available data for input into main population calculations"),
      radioButtons("xdata","X-axis",c("Data line","Mosquito density parameter")),
      radioButtons("ydata","Y-axis",c("Mosquito density parameter",
                                      "Annual EIR",
                                      "Year-round average slide prevalence",
                                      "Year-round average PCR prevalence",
                                      "Clinical incidence/year")),
      submitButton(text="Update Display"),
      width=4
    ),
    mainPanel(
      plotOutput("graph"),
      tableOutput("contents"),
      width=8
    )
  )
)

# Define server logic ----
server <- function(input, output, session) {
  library(vectorpower)
  
  data <- reactive({
    xvalues=switch(input$xdata,"Data line" = "N_M","Mosquito density parameter" = "M")
    benchmark=switch(input$ydata,"Mosquito density parameter" = "M",
                                      "Annual EIR" = "EIR",
                                      "Year-round average slide prevalence"="slide_prev",
                                      "Year-round average PCR prevalence"="pcr_prev",
                                      "Clinical incidence/year"="clin_inc")
    output_location = system.file("extdata",package="vectorpower")
    dataset_location = "DemoFolder1"
    dataset_folder=paste(output_location,dataset_location,sep="/")
    data <- get_folder_data(input_folder=dataset_folder,xvalues=xvalues,yvalues = benchmark, plot_flag = FALSE)
    return(data)
  })
  
  output$graph <- renderPlot({
    
    mv_values=data()[[1]]
    n_mv_values=length(mv_values)
    benchmark_values=data()[[2]]
    ylab=input$ydata
    
    text_size=1.0
    par(mar = c(4,4,1,1))
    matplot(mv_values,benchmark_values,type="p+l",pch=2,lty=1,col=1,lwd=1.5,
            xlab="Data line",ylab=ylab,cex.lab=text_size,cex.axis=text_size)
    
    return(plot)
    
    }, width=600, height=600,res=72)
}

# Run the app ----
shinyApp(ui = ui, server = server)
