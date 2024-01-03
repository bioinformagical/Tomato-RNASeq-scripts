#
# This is a Shiny web application that draws a barplot showing
# scaled gene expression counts.
# You can run the application by clicking the 'Run App' button above.
# 
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
source("Common.R")
sample_sheet = getSampleExcelSpreadsheet()  
scaled_counts=getScaledCounts(assembly = "SL4")
scaled_counts$assembly="SL4"
more_scaled_counts = getScaledCounts(assembly = "SL5")
more_scaled_counts$assembly = "SL5"
scaled_counts = rbind(scaled_counts,more_scaled_counts)

# Define UI for application that draws a barplot
ui <- fluidPage(

    # Application title
    titlePanel("RNA-Seq expression profile, scaled counts"),

    # Sidebar with text box to enter a gene name
    sidebarLayout(
        sidebarPanel(
            textInput("gene_name","Enter gene (from SL4 or SL5 assembly)",
                      value="Solyc10G001086")
        ),

        # Show a plot of the generated distribution
        mainPanel(
           plotOutput("barplot")
        )
    )
)

# Define server logic required to draw a barplot
server <- function(input, output) {

    output$barplot <- renderPlot({
      gene_name=input$gene_name
      assembly = scaled_counts[gene_name,"assembly"]
      description = scaled_counts[gene_name,"description"]
      main = gene_name
      ylab="counts per million (cpm)"
      if (is.na(description) | description == "NA") {
        main = paste0(gene_name," (from assembly ",
                     scaled_counts[gene_name,"assembly"],
                     ", no functional description available)")
      }
      else {
        main = paste0(gene_name," (",description,")")
      }
      indexes = which(names(scaled_counts)%in%c("description","assembly"))
      dat = scaled_counts[,-indexes]
      sample_colors = getSampleColors(names(dat),
                                      sample_sheet = sample_sheet)
      rowdata=dat[gene_name,names(dat)]
      rowdata=as.numeric(rowdata)
      names(rowdata)=names(sample_colors)
      barplot(rowdata,
              main=main,
              col=sample_colors,
              las=2,
              ylab=ylab)
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
