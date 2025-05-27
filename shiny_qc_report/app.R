library(shiny)
library(bslib)
library(GeomxTools)

# Load the plot scripts
source("PCA_plot.R")

# Load the input data
data.folder <- "/Users/cauleyes/CPTR/DSP_Analysis/test_datasets/Human_Kidney/qc_report/"
counts.data <- read.csv(paste0(data.folder, "human_kidney_test_quantile_norm_counts.csv"))
annotation.data <- read.csv(paste0(data.folder, "human_kidney_test_annotation.csv"))

# Define the column names to offer
annotation.cols <- colnames(annotation.data %>% 
                              select_if(is.character))

# Reformat the counts and annotation data for mapping and transforming
rownames(counts.data) <- counts.data$gene
counts.data <- counts.data %>% select(-gene)
colnames(counts.data) <- sub("\\.dcc", "", colnames(counts.data))
colnames(counts.data) <- gsub("\\.", "-", colnames(counts.data))

rownames(annotation.data) <- sub("\\.dcc", "", annotation.data$sample_ID)

# Define UI
ui <- page_sidebar(
  # App title ----
  title = "DSP Analysis - QC Plots",
  # Sidebar panel for inputs ----
  sidebar = sidebar(
    
    selectizeInput("annotation", 
                   "Annotation to display", 
                   choices = NULL, 
                   multiple = FALSE)
    
  ),
  
  navset_card_underline(
    nav_panel("PCA", 
              plotOutput("PCA_plot")))
  
)


# Define server logic required to create the plot
server <- function(input, output, session) {
  
  # Reactive annotation selection
  observe({
    updateSelectizeInput(session, 
                         "annotation", 
                         choices = annotation.cols, 
                         server = TRUE, 
                         selected = "slide_name")
    
  })
  
  output$PCA_plot <- renderPlot({
    
    PCA_plot(counts = counts.data, 
             annotation.df = annotation.data, 
             annotation.field = input$annotation)
    
  })
  
}

shinyApp(ui = ui, server = server)
