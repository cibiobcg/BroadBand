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
library(shinyjs) 

# cercare modo per quelle cazzo di booleane

file <- read.table('sif_cbioportal_brca.tsv', header = TRUE)

# cna_false <- file %>%
#             filter(snv.data == FALSE)
  
# Define UI for application that draws a histogram

ui <- shinyUI(fluidPage(#shinythemes::themeSelector(),
                theme = shinytheme('superhero'),
                shinyjs::useShinyjs(), # questo mi serve per 'controllare' i vari click per gli input riguradanti snv e cna  
    # Application title
    titlePanel("App Breast Cancer"),
    # Sidebar with a slider input for number of bins 
    sidebarLayout( position = 'left' ,# posso indicare al posizone dove mettere la sidebar
        sidebarPanel('Options', # posso anche mettere un sottotitolo nella sidebar
            #fileInput(
             #   inputId = 'otherfile',
              #  label = 'choose file'),
            checkboxInput(
              inputId = 'ALL',
              label = 'all',
              value = TRUE),
            checkboxInput(
              inputId = 'CNA',
              label = 'Somatic copy number alterations',
              value = FALSE),
            checkboxInput(
              inputId = 'SNV',
              label = 'Somatic single nucleotide variants',
              value = FALSE),
            checkboxGroupInput(
                inputId = 'Resources',
                label = 'Data resources',
                choices = c('brca_metabric','brca_igr_2015','brca_mbcproject_wagle_2017','breast_msk_2018','brca_tcga_pan_can_atlas_2018'),
                selected = c('brca_metabric','brca_igr_2015','brca_mbcproject_wagle_2017','breast_msk_2018','brca_tcga_pan_can_atlas_2018')),
            # selectInput(
            #   inputId = 'Data',
            #   label = 'Data',
            #   choices = list('All'= 1, 'Somatic copy number alterations' = 2, 'Somatic single nucleotides variants'= 3),
            #   selected = 1
            # ),
            checkboxGroupInput(
                inputId = 'Types',
                label = 'Breast cancer subtypes',
                choices = c('HER2+','HR+','TNBC'),
                selected = c('HER2+','HR+','TNBC')),
            checkboxGroupInput(
                inputId = 'Class',
                label = 'Tumor classification',
                choices = c('Primary','Metastasis'),
                selected = c('Primary','Metastasis')),
    conditionalPanel( # questo per creare un checbox aggiuntivo quando si passa alla secondo pannello 
              condition = 'input.tabs== 2',
                checkboxGroupInput(inputId = 'mama',label = 'mama',
                  choices = c('mama','papa','fiol','fioi'))
        )),
        # Show a plot of the generated distribution
        mainPanel(
              tabsetPanel(type= 'tabs', id = 'tabs',
                tabPanel(id='main', 'First screen',
                  fluidRow( # fluidrow serve come comnado per mettere dove vooglio i vari plot all'interno degli output
                    column(6,plotOutput(outputId = 'myplot')),
                    column(6,plotOutput(outputId = 'classplot' ))
                    ),
                  fluidRow(
                    column(3,  tableOutput(outputId = 'classtable'))
                   )),
                tabPanel(id ='sec', value = 2,
                         'Second screen',
                         
                  # fluidRow(
                    
                  # )
                  )
                           # downloadButton(
                            #  outputId = 'downloadPlot',
                             # label = 'Download Plot'),
                           )
                         )
                      )
                    )
                 )

# Define server logic required to draw a histogram
server <- function(input, output) {
  
  # primo filtro su snv e cna
  
  # observeEvent(input$ALL,{
  #   if(input$ALL == TRUE){file %>%
  #       filter(cna.data == TRUE | snv.data == TRUE)
  #   } else {NULL}
  # })
  # observeEvent(input$CNA, {
  #   if(input$CNA == TRUE){file %>%
  #       filter(cna.data == TRUE)
  #   } else {NULL}
  # })
  # observeEvent(input$SNV,{
  #   if(input$SNV == TRUE){file %>%
  #       filter(snv.data == TRUE)
  #     }else{NULL}
  # })
  #in questo modo sono risucito a rendere reattivo il file rawdata per tutte e tre i chechboxinput, in questo modo posson modificare i plot in base agli input dei checkbox della ui
  newData <- reactive({
    rawdata <- file %>%       
      group_by(data,class,type,cna.data,snv.data) %>%
      summarise(n.samples=n_distinct(sample.id),n.patients=n_distinct(patient.id))
    data1 <- rawdata
    # if (input$ALL == TRUE){
    #   data1 <- data1 %>%
    #     filter(cna.data == TRUE | snv.data == TRUE)
    # } else { data1 <- data1}
    data1 <- subset(data1, data %in% input$Resources)
    data1 <- subset(data1, type %in% input$Types)
    data1 <- subset(data1, class %in% input$Class)
    data1 <- subset(data1, cna.data %in% input$CNA)
    data1 <- subset(data1, snv.data %in% input$SNV)
    })
  #rendo interattivo solo per tutti e tre i gli input ed ? resa interattiva!!
  
  newData2 <- reactive({
    data2 <- file
    data2 <- subset(data2, type %in% input$Types)
    data2 <- subset(data2, data %in% input$Resources)
    data2 <- subset(data2, class %in% input$Class)
    data2 <- subset(data2, cna.data %in% input$CNA)
    data2 <- subset(data2, snv.data %in% input$SNV)
    data2 <- data2 %>%
        group_by(class,type) %>%
        summarise(n.samples= n_distinct(sample.id),n.patients=n_distinct(patient.id)) 
    })
  
  # devo rendere reattivo all in modo che anche quest cambi con il cambiare della soruces, quindi cambino i valori man mano
  # all <- file %>% # in questo caso %>% indica al file sif di applicare l aseguente funzione che segue la 'pipe'
  #  group_by(class,type) %>% # in questo caso il comando group_by singifica che il file vien raggruppato per classe e tipo, poi una volta applicato usa la seconda funzione che segue 
  #  summarise(n.samples=n_distinct(sample.id),n.patients=n_distinct(patient.id)) %>% # summarise sommato creand onumer ocolonne pazienti e sample facendo la conta per i subtypes
  #  add_column(data='all',.before = 'class')
  
  #unicamente per la dataTable
  
  # !!! eliminare le due colonne di cnv e snv come output?
  newData_table <- reactive({
    rawdata <- file %>%       
      group_by(data,class,type,cna.data,snv.data) %>%
      summarise(n.samples=n_distinct(sample.id),n.patients=n_distinct(patient.id))
    data3 <- rawdata 
    data3 <- subset(data3, data %in% input$Resources) # problema con all perche' non e' nell'input ( quindi o non lo rendo reattivo, o inserisco nuovo input oopure boh?)
    data3 <- subset(data3, type %in% input$Types)
    data3 <- subset(data3, class %in% input$Class)
    data3 <- subset(data3, cna.data %in% input$CNA)
    data3 <- subset(data3, snv.data %in% input$SNV)
    all2 <- data3 %>%
          group_by(class,type) %>%
          summarise(n.samples_sum=sum(n.samples),n.patients_sum=sum(n.patients))%>% # sum somma, che e' diverso da n_dsitincr il quale e' un equivalnte della funzionelenght
          add_column(data='all',.before = 'class')
    all2 <- all2 %>%
          rename(n.samples =n.samples_sum, n.patients = n.patients_sum)
    data3 <- data3 %>%
           rbind(all2)
  #data3 <- data3$cna.data == NULL
    })
  
# gestione combinazione tasti
  observeEvent(input$ALL,{
    if(input$ALL == TRUE){
               shinyjs::disable('CNA')
               shinyjs::disable('SNV')
    } else{
              shinyjs::enable('CNA')
              shinyjs::enable('SNV')
    }})
  observeEvent(input$CNA,{
    if(input$CNA == TRUE){
                shinyjs::disable('ALL')
    } else{
                shinyjs::enable('ALL')
    }})
  observeEvent(input$SNV,{
    if(input$SNV == TRUE){
              shinyjs::disable('ALL')
    } else {
              shinyjs::enable('ALL')
    }})
  
  
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
           geom_text(aes(label=n.samples)) + # questo comando serve per aggiungere la numerazione delle quantita', in qunato non si riesce a capire bene il numero reale senza indicazione
           facet_wrap(~class)})# serve per fare la divisione per classi --> in questo modo passo da un istogramma a 2 gfafici ad istogramma divisi tra metastatici e primari 
    
    #DataTable
    output$classtable <- renderTable(newData_table()) # non posso mettere due oggetti di fila?
    
}

# Run the application 
shinyApp(ui = ui, server = server)


