library(shiny)

# Define UI ----
ui <- fluidPage(
  titlePanel("Example X"),
  sidebarLayout(
    sidebarPanel(
      p("Loading results of main population computations and viewing benchmark vs time graphs for multiple baseline data
        sets and intervention parameter values"),
      radioButtons("benchmark","Benchmark",c("EIR/day","Slide prevalence","PCR prevalence","Clinical incidence/year")),
      sliderInput("n_mv_range","Baseline data set(s)",1,min=1,max=1,value=c(1,1)),
      sliderInput("n_int_range","Intervention parameter value(s)",1,min=1,max=1,value=c(1,1)),
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
  library(vectorpower) # TODO: switch to library() before pushing to repository
  # library(devtools)
  # setwd("C:/Users/kjfras16/Documents/0 - Git repositories")
  # load_all("vectorpower")
  
  input_list <- reactive({
    input_list <- readRDS(file=
                        url("https://github.com/KeithJF82/vectorpower/raw/beta/inst/extdata/examples/mainpop_data_example.Rds"))
    # input_list <- 
    #   readRDS("C:/Users/kjfras16/Documents/0 - Git repositories/vectorpower/inst/extdata/examples/mainpop_data_example.Rds")
    return(input_list)
    })
  
  observe({
    updateSliderInput(session,"n_mv_range",max=input_list()$n_mv_values)
    updateSliderInput(session,"n_int_range",max=input_list()$n_int_values)
    })
  
  set_n_mv <- reactive({
    if(length(input$n_mv_range)==1){ set_n_mv=input$n_mv_range[1] } else {
      set_n_mv=c(input$n_mv_range[1]:input$n_mv_range[2]) }
    return(set_n_mv)
  })
  
  set_n_int <- reactive({
    if(length(input$n_int_range)==1){ set_n_int=input$n_int_range[1] } else {
      set_n_int=c(input$n_int_range[1]:input$n_int_range[2]) }
    return(set_n_int)
  })
  
  data <- reactive({
    benchmark=switch(input$benchmark,"EIR/day" = "EIR",
                                      "Slide prevalence" = "slide_prev",
                                      "PCR prevalence"="pcr_prev",
                                      "Clinical incidence/year"="clin_inc")
    data <- get_mainpop_data(input_list=input_list(),set_n_mv=set_n_mv(),set_n_int=set_n_int(),benchmark=benchmark)
    return(data)
  })
  
  output$graph <- renderPlot({
    
    benchmark=names(data()$data[2])
    set_n_mv = data()$set_n_mv
    set_n_int = data()$set_n_int
    n_mv_values=length(set_n_mv)
    n_int_values=length(set_n_int)
    n_curves=n_mv_values*n_int_values
    time_values=data()$data[[1]]
    benchmark_values=data()$data[[2]]
    xlim=c(min(time_values),max(time_values)*1.4)
    ylim=c(min(benchmark_values),max(benchmark_values))
    ylab=input$benchmark
    
    titles=rep(" ",n_curves)
    titles=rep(" ",n_curves)
    colours=1+c(1:n_mv_values)
    ltypes=rep(0,n_curves)
    for(i in 1:n_int_values){
      for(j in 1:n_mv_values){
        nt=((i-1)*n_mv_values)+j
        titles[nt]=paste("n_int=",set_n_int[i],", n_mv=",set_n_mv[j],sep="")
        ltypes[nt]=i
      }
    }
    
    text_size=1.0
    par(mar = c(4,4,1,1))
    if(is.na(dim(benchmark_values)[3]) || n_curves==1){
      matplot(time_values,benchmark_values,type="l",lty=set_n_int,col=colours,lwd=1.5,
              xlab="time (days)",ylab=ylab,xlim=xlim,ylim=ylim,cex.lab=text_size,cex.axis=text_size)
    } else {
      matplot(time_values,benchmark_values[,1,],type="l",lty=1,col=colours,lwd=1.5,
              xlab="time (days)",ylab=ylab,xlim=xlim,ylim=ylim,cex.lab=text_size,cex.axis=text_size)
      for(i in 2:n_int_values){
        matplot(time_values,benchmark_values[,i,],type="l",lty=i,col=colours,lwd=1.5,add=TRUE,xaxt="n",yaxt="n")
      }
    }
    legend("bottomright", inset=0.01, legend=titles, lty=ltypes,col=rep(colours,n_int_values), 
           horiz=FALSE,bg='white',cex=text_size)
    
    return(plot)
    
    }, width=600, height=600,res=72)
}

# Run the app ----
shinyApp(ui = ui, server = server)
