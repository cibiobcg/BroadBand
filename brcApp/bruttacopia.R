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
library(gridExtra)

# cercare modo per quelle cazzo di booleane

file <- read.table('sif_cbioportal_brca.tsv', header = TRUE)
file2<- read.delim('snv_freq_brca.tsv', header = TRUE, stringsAsFactors = FALSE) %>% group_by(class,type) 

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
    conditionalPanel( # questo per creare un slider aggiuntivo quando si passa alla secondo pannello 
              condition = 'input.tabs== 2',
                sliderInput(inputId = 'Gene_filter',label = 'Gene_filter', min = 5, max = 25, value = 25)
              )
        ),
        # Show a plot of the generated distribution
        mainPanel(
              tabsetPanel(type= 'tabs', id = 'tabs',                           #aggisutare nomi dei pannelli
                tabPanel(id='main', 'First screen',
                  fluidRow( # fluidrow serve come comnado per mettere dove vooglio i vari plot all'interno degli output
                    column(6,plotOutput(outputId = 'myplot')),
                    column(6,plotOutput(outputId = 'classplot' ))
                    ),
                  fluidRow(
                    column(3,dataTableOutput(outputId = 'classtable'))
                   )),
                tabPanel(id ='sec', value = 2,
                         'Second screen',
                   fluidRow(
                     column(6,plotOutput(outputId = 'barplot')),
                     column(6,plotOutput(outputId = 'heatmap'))
                   ),
                   fluidRow(
                     column(3,dataTableOutput(outputId = 'table2'))
                   )
                  ))
                           # downloadButton(
                            #  outputId = 'downloadPlot',
                             # label = 'Download Plot'),
                           
                         )
                      )
                    ))
                 
server <- function(input, output) {

#############################################################################################################
##################################### Per primo pannello ###################################################
  
  mw <- file2 %>%
    arrange(desc(w.mean))
  
  #in questo modo sono risucito a rendere reattivo il file rawdata per tutte e tre i chechboxinput, in questo modo posson modificare i plot in base agli input dei checkbox della ui
  
  newData <- reactive({
    rawdata <- file %>%       
      group_by(data,class,type,cna.data,snv.data) %>%
      summarise(n.samples=n_distinct(sample.id),n.patients=n_distinct(patient.id))
    data1 <- rawdata
    data1 <- subset(data1, data %in% input$Resources)
    data1 <- subset(data1, type %in% input$Types)
    data1 <- subset(data1, class %in% input$Class)
    if(input$ALL == TRUE){
      data1 <- data1 %>%
              filter(snv.data == TRUE | cna.data == TRUE)
    }else if (input$SNV == TRUE & input$CNA == FALSE){
      data1 <- data1 %>% 
              filter(snv.data == TRUE)
    }else if (input$CNA == TRUE & input$SNV == FALSE) {
      data1 <- data1 %>% 
              filter(cna.data == TRUE)
    }else if (input$SNV & input$CNA){
      data1 <- data1 %>% 
                filter(cna.data == TRUE & snv.data ==TRUE)
    }
    # data1 <- subset(data1, cna.data %in% input$CNA)
    # data1 <- subset(data1, snv.data %in% input$SNV)
    })
  #rendo interattivo solo per tutti e tre i gli input ed ? resa interattiva!!
  
  newData2 <- reactive({
    data2 <- file
    data2 <- subset(data2, type %in% input$Types)
    data2 <- subset(data2, data %in% input$Resources)
    data2 <- subset(data2, class %in% input$Class)
    if(input$ALL == TRUE){
      data2 <- data2 %>%
        filter(snv.data == TRUE | cna.data == TRUE)
    }else if (input$SNV == TRUE & input$CNA == FALSE){
      data2 <- data2 %>% 
        filter(snv.data == TRUE)
    }else if (input$CNA== TRUE & input$SNV == FALSE) {
      data2 <- data2 %>% 
        filter(cna.data == TRUE)
    }else if (input$SNV & input$CNA){
      data2 <- data2 %>% 
        filter(cna.data == TRUE & snv.data ==TRUE)
    }
    # data2 <- subset(data2, cna.data %in% input$CNA)
    # data2 <- subset(data2, snv.data %in% input$SNV)
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
  
  # !!! eliminare le due colonne di cnv e snv come output? nope
  newData_table <- reactive({
    rawdata <- file %>%       
      group_by(data,class,type,cna.data,snv.data) %>%
      summarise(n.samples=n_distinct(sample.id),n.patients=n_distinct(patient.id))
    data3 <- rawdata 
    data3 <- subset(data3, data %in% input$Resources) # problema con all perche' non e' nell'input ( quindi o non lo rendo reattivo, o inserisco nuovo input oopure boh?)
    data3 <- subset(data3, type %in% input$Types)
    data3 <- subset(data3, class %in% input$Class)
    if(input$ALL == TRUE){
      data3 <- data3 %>%
        filter(snv.data == TRUE | cna.data == TRUE)
    }else if (input$SNV == TRUE & input$CNA == FALSE){
      data3 <- data3 %>% 
        filter(snv.data == TRUE)
    }else if (input$CNA == TRUE & input$SNV == FALSE) {
      data3 <- data3 %>% 
        filter(cna.data == TRUE)
    }else if (input$SNV & input$CNA){
      data3 <- data3 %>% 
        filter(cna.data == TRUE & snv.data ==TRUE)
    }
    # data3 <- subset(data3, cna.data %in% input$CNA)
    # data3 <- subset(data3, snv.data %in% input$SNV)
    all2 <- data3 %>%
          group_by(class,type) %>%
          summarise(n.samples_sum=sum(n.samples),n.patients_sum=sum(n.patients))%>% # sum somma, che e' diverso da n_dsitincr il quale e' un equivalnte della funzionelenght
          add_column(data='all',.before = 'class')
    all2 <- all2 %>%
          rename(n.samples =n.samples_sum, n.patients = n.patients_sum)
    data3 <- data3 %>%
           rbind(all2)
    })
  
# quindi in questo caso riarriangiamo con arrange dati del file in modo decrescente con desc in base a w.mean (colonna)
  
 
  
  #data per barplot
  data_second_pannel <- reactive({
    data_pan_1 <- mw
    data_pan_1 <- subset(data_pan_1, class %in% input$Class)
    data_pan_2 <- data_pan_1 %>% 
                  filter(type == 'all')
    data_pan_1 <- data_pan_1 %>% 
                  filter(type != 'all')
    data_pan_1 <- subset(data_pan_1, type %in% input$Types)
    data_pan_3 <- data_pan_1 %>% 
                rbind(data_pan_2)
    data_pan_3 <- data_pan_3 %>%
                  slice_head(n = input$Gene_filter)  %>%
                  group_split()  
  })
  
  data_second_pannel_heatmap <- reactive({
    mw2 <- mw %>%
          slice_head(n = input$Gene_filter) %>%
          group_split()
    mat <- do.call(rbind,mw2)
    mat <- subset(mat, type %in% input$Types)
    mat <- subset(mat, class %in% input$Class)
    mat <- mat %>%
      filter(type != 'all')
  })
  
  data_second_pannel_table <- reactive({
    mw3<- mw %>%
      slice_head(n = input$Gene_filter) %>%
      group_split()
    mat1 <- do.call(rbind,mw3)  
    mat1 <- mat1[order(-mat1$w.mean),]  # ordinato per w.mean decrescente 
    mat1 <- subset(mat1, type %in% input$Types)
    mat1 <- subset(mat1, class %in% input$Class)
    mat1 <- mat1 %>%
      filter(type != 'all')
  })
############### gestione combinazione tasti ALL - SNV - CNA #####################
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
    output$classtable <- renderDataTable(newData_table()) 
    
#############################################################################################################
##################################### Per secondo pannello ###################################################
  
       output$barplot <- renderPlot({        
         
         plist <- list()
       for(i in 1:length(data_second_pannel())){
         dn <- data_second_pannel()[[i]]
         dn$Hugo_Symbol <- factor(dn$Hugo_Symbol,levels = rev(dn$Hugo_Symbol))
         # data_pan_1 <- data_pan_1[order(-data_pan_1$Hugo_Symbol),]
         # dn <- subset(dn, type %in% input$Types)
         # dn <- subset(dn, class %in% input$Class)
         p<- ggplot(data=dn, aes(x=Hugo_Symbol, y=w.mean)) +
           geom_bar(stat="identity") + coord_flip() +
           facet_wrap(type~class,scales = 'free')
         plist[[i]] <- ggplotGrob(p)
       }
         grid.arrange(grobs=plist,ncol=4) 
    })
    
    output$heatmap <- renderPlot({ 
      ggplot(data_second_pannel_heatmap(), aes(type, Hugo_Symbol)) +  
        geom_tile(aes(fill = w.mean)) + 
        geom_text(aes(label = round(w.mean, 3)),size=2) +
        scale_fill_gradient(low = "white", high = "red") +
        facet_wrap(~class)
    })
    output$table2 <- renderDataTable(data_second_pannel_table())
}

# Run the application 
shinyApp(ui = ui, server = server)


