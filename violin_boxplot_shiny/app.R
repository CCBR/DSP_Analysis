library(shiny)
library(bslib)
library(GeomxTools)

# Load the plot script
source("violin_boxplot.R")

# Load the input data
count.data <- readRDS("data/Tosato_q3_counts_vessel.RDS")
annotation.data <- readRDS("data/Tosato_filtered_annotation_vessel.RDS")

# Define the column names to offer
annotation.cols <- colnames(annotation.data %>% 
                              select_if(is.character))

# Grab an example gene
#example.gene <- count.data$gene[1]

# Define UI
ui <- page_sidebar(
  # App title ----
  title = "Violin/Boxplot for a Gene of Interest",
  # Sidebar panel for inputs ----
  sidebar = sidebar(
    
    selectizeInput("annotation", 
                   "Annotation to display", 
                   choices = NULL, 
                   multiple = FALSE), 
    
    selectizeInput("gene", 
                   "Gene to display", 
                   choices = NULL, 
                   multiple = FALSE), 
    
    checkboxInput("summary_stats", 
                  "Display Summary Stats", 
                  value = FALSE), 
    
    checkboxInput("compare", 
                  "Wilcoxon Rank Sum Test", 
                  value = FALSE), 
    
    selectizeInput("compare.groups", 
                   "Groups to Compare", 
                   choices = NULL, 
                   multiple = FALSE), 
    
  ),
  
  plotOutput(outputId = "Violin_Boxplot")
)


# Define server logic required to create the plot
server <- function(input, output, session) {
  
  # Reactive gene name selection
  observe({
    updateSelectizeInput(session, 
                         "gene", 
                         choices = rownames(count.data), 
                         server = TRUE)
  })
  
  # Reactive annotation selection
  observe({
    updateSelectizeInput(session, 
                         "annotation", 
                         choices = annotation.cols, 
                         server = TRUE)
    
    updateSelectizeInput(session, 
                         "compare.groups", 
                         choices = unique(annotation.data[[input$annotation]]), 
                         server = TRUE)
  })
  
  # Reactive summary stats toggle
  summary_stats <- reactive({input$summary_stats})
  
  # Reactive compare toggle
  compare <- reactive({input$compare})
  
  # Groups to compare
  compare.list <- list(input$compare.groups)
  
  output$Violin_Boxplot <- renderPlot({
    
    violin_boxplot(counts = count.data,  
                   annotation.df = annotation.data, 
                   gene.list = input$gene, 
                   annotation.field = input$annotation, 
                   display.summary.stat = input$summary_stats, 
                   compare = input$compare, 
                   compare.groups = compare.list)
    
  })
  
}

shinyApp(ui = ui, server = server)
