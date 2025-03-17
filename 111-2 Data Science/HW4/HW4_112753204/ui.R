library(shiny)
library(ggplot2)
library(DT)

fluidPage(
  
  titlePanel("Data Science HW4 - Iris Dataset Analysis"),
  h4("資科碩一 112753204 陳品伃"),
  #Page1
  tabsetPanel(
    tabPanel("PCA",
             h3("PCA Plot - Iris"),
             
             sidebarPanel(
               sliderInput('input_size', 'Input data size', min=6, max=nrow(iris),
                           value=min(150, nrow(iris)), round=0),
               selectInput('x', 'X',choices = c("PC1", "PC2", "PC3", "PC4"), selected = "PC1"),
               selectInput('y', 'Y',choices = c("PC1", "PC2", "PC3", "PC4"), selected = "PC2")),
             
             mainPanel(
               plotOutput("pca_plot"))
    ),
    #Page2
    tabPanel("PCA Table",
             fluidRow(
               column(width = 8,
                      h3("1. Summary of PCA"),
                      verbatimTextOutput("pca_summary"))
             ),
             fluidRow(
               column(width = 10,
                      h3("2. Rotation Table"),
                      tableOutput(outputId = "rotation"))
             ),
             fluidRow(
               column(width = 10,
                      h3("3. Center Table"),
                      tableOutput(outputId = "center"))
             ),
             fluidRow(
               column(width = 10,
                      h3("4. X Table"),
                      DTOutput(outputId = "xtable"))
             )
    ),
    #Page3
    tabPanel("CA",
             h3("CA Plot - Iris"),
             
             fluidRow(
              sidebarPanel(
                sliderInput('ca_input_size', 'Input data size', min=6, max=nrow(iris),
                             value=min(150, nrow(iris)), round=0)),
              mainPanel(
                plotOutput("ca_plot"))
             ),
             fluidRow(
               column(width = 12,
                      h3("CA Summary"),
                      verbatimTextOutput(outputId = "ca_summary"))
             )
      ),
    #Page4
    tabPanel("CA Kmeans",
             h3("CA Kmeans Plot - Iris"),
             
               sidebarPanel(
                 sliderInput('data_size', 'Input data size', min=11, max=nrow(iris),
                             value=min(150, nrow(iris)), round=0),
                 sliderInput('k_size', 'Center(k)', min=3, max=10,
                             value=min(10, 3), round=0)),
             
               mainPanel(
                 plotOutput("kmeans_plot"))
    )
  )
)

