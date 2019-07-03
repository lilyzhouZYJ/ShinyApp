library(shiny)
library(plotly)
library(stringr)

# Define UI ----
ui <- fluidPage(
  titlePanel("Shiny App"),
  
  sidebarLayout(
    sidebarPanel(
      checkboxInput("use_exon",
                    label="Check to use exon; uncheck to use gene name",
                    value=TRUE),
      textInput("exon",
                label="Insert exon info",
                value="CNV_1_865282_865854"),
      textInput("gene",
                label="Insert gene name",
                value="LINC01128"),
      checkboxInput("svd",
                    label="SVD deletion: 10",
                    value=FALSE),
      textInput("sample",
                label="Sample to highlight",
                value="M_090252"),
      checkboxInput("parent",
                    label="Display both parents",
                    value=FALSE)
    ),
    mainPanel(
      plotOutput("plot", width=1000, height=600)
    )
  )
)

# Define server logic ----
server <- function(input, output) {
  dataInput <- reactive({
    if(input$use_exon==TRUE){
      ex <- str_replace(input$exon, "X", "23")
      ex <- str_replace(ex, "Y", "24")
      ex <- str_split(ex, "_")
      input_start <- as.numeric(ex[[1]][3])
      input_stop <- as.numeric(ex[[1]][4])
      input_chr <- ex[[1]][2]
      if(input$svd==TRUE){
        data <- get(paste('dataset_10_chr',input_chr, sep=''))
        start_ind <- data$index[which(data$start>=input_start)[1]]
        end_ind <- data$index[max(which(data$stop<=input_stop))]
        data <- filter(data, index>=start_ind-5, index<=end_ind+5)
      }
      else{
        data <- get(paste('dataset_chr', input_chr, sep = ''))
        start_ind <- data$index[which(data$start>=input_start)[1]]
        end_ind <- data$index[tail(which(data$stop<=input_stop), n=1)]
        data <- filter(data, index>=start_ind-5, index<=end_ind+5)
      }
      data
    }
    else{
      row <- gene_list[which(gene_list$hgnc_symbol==input$gene),]
      input_chr <- row$chromosome_name
      if(input_chr=="X") {input_chr="23"}
      else if(input_chr=="Y") {input_chr="24"}
      input_start <- row$start_position
      input_stop <- row$end_position
      if(input$svd==TRUE){
        data <- get(paste('dataset_10_chr',input_chr, sep=''))
        start_ind <- data$index[which(data$start>=input_start)[1]]
        end_ind <- data$index[max(which(data$stop<=input_stop))]
        data <- filter(data, index>=start_ind-5, index<=end_ind+5)
      }
      else{
        data <- get(paste('dataset_chr', input_chr, sep = ''))
        start_ind <- data$index[which(data$start>=input_start)[1]]
        end_ind <- data$index[tail(which(data$stop<=input_stop), n=1)]
        data <- filter(data, index>=start_ind-5, index<=end_ind+5)
      }
      data
    }
  })
  
  output$plot <- renderPlot({
    data <- dataInput()
    plot <- ggplot(data, aes(x=exon, y=zrpkm, group=sample)) + 
      geom_smooth(method="loess", formula="y ~ x", span=0.7, se=FALSE, color="lightgray", size=1) +
      theme_minimal() +
      xlab("Exons") +
      ylab("ZRPKM") +
      theme(
        axis.text.x = element_text(size=15, angle=90),
        axis.text.y = element_text(size=15),
        axis.title = element_text(size=18)
      )
    if(input$parent==TRUE){
      plot <- plot %+% filter(data, !str_detect(sample, input$sample))
      plot + geom_line(data=filter(data,str_detect(sample,input$sample)), aes(color=sample),size=2) +
        geom_point(data=filter(data,str_detect(sample,input$sample)), aes(color=sample),size=4)
    }
    else{
      plot <- plot %+% filter(data, sample!=input$sample)
      plot + geom_line(data=filter(data,sample==input$sample), aes(color=sample),size=2) +
        geom_point(data=filter(data,sample==input$sample),aes(color=sample),size=4)
    }
  })
}

# Run the app ----
shinyApp(ui = ui, server = server)