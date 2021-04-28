#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinythemes)
library(tidyverse)

file <- read.table('sif_cbioportal_brca.tsv', header = TRUE)

# Define UI for application that draws a histogram
ui <- shinyUI(fluidPage(#shinythemes::themeSelector(),
                theme = shinytheme('superhero'),

    # Application title
    titlePanel("App Breast Cancer"),
    # Sidebar with a slider input for number of bins 
    sidebarLayout( position = 'left' ,# posso indicare al posizone dove mettere la sidebar
        sidebarPanel('Options', # posso anche mettere un sottotitolo nella sidebar
            #fileInput(
             #   inputId = 'otherfile',
              #  label = 'choose file'),
            checkboxGroupInput(
                inputId = 'Resources',
                label = 'Data resources',
                choices = c('brca_metabric','brca_igr_2015','brca_mbcproject_wagle_2017','breast_msk_2018','brca_tcga_pan_can_atlas_2018'),
                selected = c('brca_metabric','brca_igr_2015','brca_mbcproject_wagle_2017','breast_msk_2018','brca_tcga_pan_can_atlas_2018')),
            checkboxGroupInput(
                inputId = 'Types',
                label = 'Breast cancer subtypes',
                choices = c('HER2+','HR+','TNBC'),
                selected = c('HER2+','HR+','TNBC')),
            checkboxGroupInput(
                inputId = 'Class',
                label = 'Tumor classification',
                choices = c('Primary','Metastasis'),
                selected = c('Primary','Metastasis'))
        ),
    
        # Show a plot of the generated distribution
        mainPanel('main panel',
              tabsetPanel(
                tabPanel('main', 
                  fluidRow( # fluidrow serve come comnado per mettere dove vooglio i vari plot all'interno degli output
                    column(6,plotOutput(outputId = 'myplot')),
                    column(6,plotOutput(outputId = 'classplot' ))
                    ),
                  fluidRow(
                    column(3,  tableOutput(outputId = 'classtable'))
                   )),
                tabPanel('sec',
                           # downloadButton(
                            #  outputId = 'downloadPlot',
                             # label = 'Download Plot'),
                           )
          )
        )
)
))

# Define server logic required to draw a histogram
server <- function(input, output) {
  
  #data manilupation
  
  rawdata <- file %>% # questo e' il file riadattato per il primo grafico iterattivo dove posso variare tutte le variabili
    group_by(data,class,type) %>%
    summarise(n.samples=n_distinct(sample.id),n.patients=n_distinct(patient.id))
  
  #in questo modo sono risucito a rendere reattivo il file rawdata per tutte e tre i chechboxinput, in questo modo posson modificare i plot in base agli input dei checkbox della ui
  newData <- reactive({
    data1 <- rawdata
    data1 <- subset(data1, data %in% input$Resources)
    data1 <- subset(data1, type %in% input$Types)
    data1 <- subset(data1, class %in% input$Class)})
  
  rawdata2 <- file %>% # file per grafico e tabella tra metastatici e primari 
    group_by(class,type) %>%
    summarise(n.samples=n_distinct(sample.id),n.patients=n_distinct(patient.id))
  
  # rendo iterattivo solo per le resources e i types, le classi dovrebbero rimanere costanti 
  newData2 <- reactive({
    data2 <- rawdata2
    data2 <- subset(data2, type %in% input$Types)
    data2 <- subset(data2, class %in% input$Class)})
  
  all <- file %>% # in questo caso %>% indica al file sif di applicare l aseguente funzione che segue la 'pipe'
    group_by(class,type) %>% # in questo caso il comando group_by singifica che il file vien raggruppato per classe e tipo, poi una volta applicato usa la seconda funzione che segue 
    summarise(n.samples=n_distinct(sample.id),n.patients=n_distinct(patient.id)) %>% # summarise sommato creand onumer ocolonne pazienti e sample facendo la conta per i subtypes
    add_column(data='all',.before = 'class') # aggiunto colonna all prima di classe nella data frame
  
  rawdata_table <- rawdata %>%
                   rbind(all)
  
  #unicamente per la dataTable
  newData_table <- reactive({
  data3 <- rawdata_table
  # data3 <- subset(data3, data %in% input$Resources) # problema con all perche' non e' nell'input ( quindi o non lo rendo reattivo, o inserisco nuovo input oopure boh?)
  data3 <- subset(data3, type %in% input$Types)
  data3 <- subset(data3, class %in% input$Class)
  })
  
# non lo lascia fare, devo vedere come unire due reactive object o seno pensare altra soluzione 
  # rawdata_table <- reactive({
  #   data3 <- rawdata
  #   data3 <- subset(data3, data %in% input$Resources)
  #   data3 <- subset(data3, type %in% input$Types)
  #   data3 <- subset(data3, class %in% input$Class)})
  # 
  # all2 <- reactive({
  #   data4<- all
  #   data4<- subset(data4, type %in% input$Types)
  #   data4<- subset(data4, class %in% input$Class)
  # })
  
  #non se po fa 
  #Data_talbe <- rawdata_table %>%
   #             rbind(all2)

       #plotting
    output$myplot <- renderPlot({
        # generate bins based on input$bins from ui.R
        ggplot(newData(), aes(x=type, y=n.samples, fill=class)) +
            geom_bar(stat="identity") + 
            scale_fill_manual(values=c('#999999','#E69F00'))})
    output$classplot <- renderPlot({
        ggplot(newData2(), aes(x=type,y=n.samples,fill=class)) +
           geom_bar(stat="identity") + theme(aspect.ratio = 1,legend.position = "none") +
           scale_fill_manual(values=c('#999999','#E69F00')) +
           geom_text(aes(label=n.samples)) + # questo comando serve per aggiungere la numerazione delle quantita', in qunato non si riesce a capire bene il nuemro reale senza indicazione
           facet_wrap(~class)})# serve per fare la divisione per classi --> in questo modo passo da un istogramma a 2 gfafici ad istogramma divisi tra metastatici e primari 
    
    #DataTable
    output$classtable <- renderTable(newData_table()) # non posso mettere due oggetti di fila?
    
}

# Run the application 
shinyApp(ui = ui, server = server)
