library(shiny)
library(shinythemes)
library(tidyverse)
library(shinyjs) 
library(gridExtra)
library(data.table)
library(parallel)
library(wesanderson)
library(shinyBS)
library(UpSetR)

#loading files
file <- read.delim('sif_cbioportal_brca.tsv', header = TRUE,stringsAsFactors = FALSE)
file2 <- read.delim('snvs_raw_data.tsv', header = TRUE, stringsAsFactors = FALSE)
ensembl <- read.delim('mart_export_GRCh38p13.tsv',check.names = F,stringsAsFactors = F)
goi <- readLines('genes_of_interest.txt')
load('scna_data.RData')

# file <- read.delim('./Files_app/sif_cbioportal_brca.tsv', header = TRUE,stringsAsFactors = FALSE)
# file2 <- read.delim('./Files_app/snvs_raw_data.tsv', header = TRUE, stringsAsFactors = FALSE)
# ensembl <- read.delim('./Files_app/mart_export_GRCh38p13.tsv',check.names = F,stringsAsFactors = F)
# goi <- readLines('./Files_app/genes_of_interest.txt')
# load('./Files_app/scna_data.RData')


ui <- shinyUI(fluidPage(#shinythemes::themeSelector(),
  theme = shinytheme('superhero'),
  shinyjs::useShinyjs(), # questo mi serve per 'controllare' i vari click per gli input riguradanti snv e cna   
  # Application title
  tags$head(HTML("<title>BroadBand</title>")),
  titlePanel(
    div(img(height = 150, width = 250, src='logo_app.jpg',style = 'border-radius: 20%'),"BroadBand")),
  # Sidebar with a slider input for number of bins 
  sidebarLayout( position = 'left' ,# posso indicare al posizone dove mettere la sidebar
                 sidebarPanel('Options', width = 2, # posso anche mettere un sottotitolo nella sidebar
                              conditionalPanel( # questo per creare un slider aggiuntivo quando si passa alla secondo pannello 
                                condition = 'input.tabs== 1',
                                fileInput(
                                  inputId = 'otherfile1',
                                  label = 'choose tsv file',
                                  multiple = TRUE,
                                  accept = '.tsv'),
                                div(style = "margin-top:-550px;display:inline-block ",
                                  actionButton("cancel1", "Cancel",icon("recycle"),#paper-plane
                                             style="color: #fff; background-color: #c64c04; padding: 14px; border-radius: 20%"),
                                  actionButton('info1','',icon('question'))),
                                div(style = 'margin-top: 20px',
                                checkboxGroupInput(
                                  inputId = 'Resources1',
                                  label = 'Data resources',
                                  choices = c('brca_metabric','brca_igr_2015','brca_mbcproject_wagle_2017','breast_msk_2018','brca_tcga_pan_can_atlas_2018'),
                                  selected = c('brca_metabric','brca_igr_2015','brca_mbcproject_wagle_2017','breast_msk_2018','brca_tcga_pan_can_atlas_2018'))),
                                h6('Available somatic data'),
                                div(style = "margin-top: -5px ",
                                checkboxInput(
                                  inputId = 'ALL',
                                  label = 'all',
                                  value = TRUE)),
                                div(style = "margin-top: -10px ",
                                checkboxInput(
                                  inputId = 'CNA',
                                  label = 'Somatic copy number alterations',
                                  value = FALSE)),
                                div(style = "margin-top: -10px ",
                                checkboxInput(
                                  inputId = 'SNV',
                                  label = 'Somatic single nucleotide variants',
                                  value = FALSE)),
                                checkboxGroupInput(
                                  inputId = 'Types1',
                                  label = 'Breast cancer subtypes',
                                  choices = c('HR+','HER2+','TNBC'),
                                  selected = c('HR+','HER2+','TNBC')),
                                checkboxGroupInput(
                                  inputId = 'Class1',
                                  label = 'Tumor classification',
                                  choices = c('Primary','Metastasis'),
                                  selected = c('Primary','Metastasis')),
                                actionButton("updatebutton1", "Apply",icon("dna"),#paper-plane
                                             style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
                                div(style="display:inline-block;width:32%;text-align: center;",
                                actionButton('Reset1',"Reset",icon("trash-alt"),
                                             style = 'color: #fff; background-color: #BF0000; border-color:#BF0000'))),
                              conditionalPanel( # questo per creare un slider aggiuntivo quando si passa alla secondo pannello 
                                condition = 'input.tabs== 2',
                                fileInput(
                                  inputId = 'otherfile2',
                                  label = 'choose tsv file',
                                  multiple = TRUE,
                                  accept = '.tsv'),
                                div(style = "margin-top:-550px;display:inline-block ",
                                    actionButton("cancel2", "Cancel",icon("recycle"),#paper-plane
                                                 style="color: #fff; background-color: #c64c04; padding: 14px; border-radius: 20%"),
                                actionButton('info2','',icon('question'))),
                                div(style = 'margin-top: 20px',
                                checkboxGroupInput(
                                  inputId = 'Resources2',  
                                  label = 'Data resources',
                                  choices = c('brca_metabric','brca_igr_2015','brca_mbcproject_wagle_2017','breast_msk_2018','brca_tcga_pan_can_atlas_2018'),
                                  selected = c('brca_metabric','brca_igr_2015','brca_mbcproject_wagle_2017','breast_msk_2018','brca_tcga_pan_can_atlas_2018'))),
                                checkboxGroupInput(
                                  inputId = 'Types2',
                                  label = 'Breast cancer subtypes',
                                  choices = c('HR+','HER2+','TNBC'),
                                  selected = c('HR+','HER2+','TNBC')),
                                checkboxGroupInput(
                                  inputId = 'Class2',
                                  label = 'Tumor classification',
                                  choices = c('Primary','Metastasis'),
                                  selected = c('Primary','Metastasis')),
                                sliderInput(inputId = 'Gene_filter',label = 'Gene_filter', min = 5, max = 25, value = 25),
                                actionButton("updatebutton2", "Apply",icon("dna"),#paper-plane
                                             style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
                                div(style="display:inline-block;width:32%;text-align: center;",
                                    actionButton('Reset2',"Reset",icon("trash-alt"),
                                                 style = 'color: #fff; background-color: #BF0000; border-color:#BF0000'))),
                              conditionalPanel(
                                condition = 'input.tabs== 3',
                                fileInput(
                                  inputId = 'otherfile3',
                                  label = 'choose tsv file',
                                  multiple = TRUE,
                                  accept = '.tsv'),
                                div(style = "margin-top:-550px;display:inline-block ",
                                    actionButton("cancel3", "Cancel",icon("recycle"),#paper-plane3
                                                 style="color: #fff; background-color: #c64c04; padding: 14px; border-radius: 20%"),
                                    actionButton('info3','',icon('question'))),
                                div(style = 'margin-top: 20px',
                                checkboxGroupInput(
                                  inputId = 'Resources3',
                                  label = 'Data resources',
                                  choices = c('brca_metabric','brca_igr_2015','brca_mbcproject_wagle_2017','brca_tcga_pan_can_atlas_2018'),
                                  selected = c('brca_metabric','brca_igr_2015','brca_mbcproject_wagle_2017','brca_tcga_pan_can_atlas_2018'))),
                                checkboxGroupInput(
                                  inputId = 'Types3',
                                  label = 'Breast cancer subtypes',
                                  choices = c('HR+','HER2+','TNBC'),
                                  selected = c('HR+','HER2+','TNBC')),
                                checkboxGroupInput(
                                  inputId = 'Class3',
                                  label = 'Tumor classification',
                                  choices = c('Primary','Metastasis'),
                                  selected = c('Primary','Metastasis')),
                                radioButtons(inputId = 'copynumber_granularity',
                                             label = 'Copy number granularity',
                                             choices = c('2 classes ' = TRUE, '4 classes' = FALSE),  # vedere se va bene senza true e false , inoltre aggiunger e
                                             selected = TRUE),
                                radioButtons(inputId = 'Groups3',
                                             label = 'Choose',choices =c('deletion' = 'homodel','amplification' = 'ampl'),selected = 'ampl'),
                                sliderInput(inputId = 'filter_median_freq',
                                            label= 'Filter Median frequencing', min = 0, max= 1, value=0.02,step = 0.01),
                                selectInput(inputId ='Chromosomes',
                                            label = 'Chromosomes',
                                            choices = c('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X'),
                                            selected = '1'),
                                selectInput(inputId = 'Cytoband',
                                            label = 'Cytoband',
                                            choices = c(''),
                                            multiple = FALSE),   
                                actionButton("updatebutton3", "Apply",icon("dna"),#paper-plane
                                             style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
                                div(style="display:inline-block;width:32%;text-align: center;",
                                    actionButton('Reset3',"Reset",icon("trash-alt"),
                                                 style = 'color: #fff; background-color: #BF0000; border-color:#BF0000')))
                 ),
                 # Show a plot of the generated ditribution
                 mainPanel(
                   tabsetPanel(type= 'tabs', id = 'tabs',                           
                               tabPanel(id='main', value= 1, 
                                        'Overview',
                                        downloadButton(outputId = 'download_myplot',label = 'Download count plot'),
                                        downloadButton(outputId = 'download_classplot',label = 'Download primary-metastasis plot'),
                                        downloadButton(outputId = 'download_classtable',label = 'Download table'),
                                        tableOutput(outputId = 'out_data_1'),
                                        bsModal(id= 'information1', title = 'INFORMATION ABOUT FILES UPLOAD', trigger = 'info1', size = 'medium',textOutput('textinfo1')),
                                        fluidRow( style='height:60vh',# fluidrow serve come comnado per mettere dove vooglio i vari plot all'interno degli output
                                                  column(6, plotOutput(outputId = 'myplot', height = '700px')#, style = "height:400px; width:400 px"
                                                  ),
                                                  column(6, plotOutput(outputId = 'classplot', height = '700px')
                                                  )),
                                        fluidRow(
                                          column(12,dataTableOutput(outputId = 'classtable'), style = 'width:100.5%')
                                        )),
                               tabPanel(id ='sec', value = 2,
                                        'Somatic Single Nucleotide Variants (SNVs)',
                                        downloadButton(outputId = 'download_barplot',label = 'Download barplot'),
                                        downloadButton(outputId = 'download_heatmap',label = 'Download heatmap'),
                                        downloadButton(outputId = 'download_table2',label = 'Download table'),
                                        tableOutput(outputId = 'out_data_2'),
                                        bsModal(id= 'information2', title = 'INFORMATION ABOUT FILES UPLOAD', trigger = 'info2', size = 'medium',textOutput('textinfo2')),
                                        fluidRow(style='height:80vh',
                                                 column(12,plotOutput(outputId = 'barplot', width = '125%', height = '820px'))),
                                        fluidRow(style='height:80vh',
                                                 column(12,plotOutput(outputId = 'heatmap', width = '125%', height = '950px'))),
                                        fluidRow(style='height:80vh',
                                                 column(12,plotOutput(outputId = 'barplot_heat', width = '125%', height = '1000px'))),
                                        fluidRow(
                                          column(12,dataTableOutput(outputId = 'table2'), style = 'width:125%')
                                        )),
                               tabPanel(id = 'thrd', value = 3,
                                        'Somatic Copy Number Aberrations (SCNAs)',
                                        downloadButton(outputId = 'download_plot_All', label = 'download all chromosomes plot'),
                                        downloadButton(outputId = 'download_plots', label = 'downlaod chromosome plot'),
                                        downloadButton(outputId = 'download_cytoband',label = 'download cytoband plot'),
                                        downloadButton(outputId = 'download_table_chromosome', label = 'Download table chromosome'),
                                        downloadButton(outputId = 'download_table_cytoband', label = 'Download table cytoband'),
                                        tableOutput(outputId = 'out_data_3'),
                                        bsModal(id= 'information3', title = 'INFORMATION ABOUT FILES UPLOAD', trigger = 'info3', size = 'medium',textOutput('textinfo3')),
                                        fluidRow(style='height:80vh',
                                                 column(12,plotOutput(outputId = 'All_plot',width = '125%', height = '820px'))
                                        ),
                                        fluidRow(style='height:80vh',
                                                 column(12,plotOutput(outputId = 'plot',width = '125%', height = '820px'))
                                        ),
                                        fluidRow(
                                          column(12,plotOutput(outputId = 'cytoband', width = '125%', height = '820px'))
                                        ),
                                        fluidRow( 
                                          column(12,
                                                 tabsetPanel(
                                                   tabPanel('Chromosome',
                                                            fluidRow(
                                                              column(12,dataTableOutput(outputId = 'table_chromosome'), style = "width:125%"))),
                                                   tabPanel('Cytoband',
                                                            fluidRow(
                                                              column(12,dataTableOutput(outputId = 'table_cytoband'), style = "width:125%")))
                                                 )
                                          )
                                        )
                               )
                   )
                 )
  )
)
)

server <- function(input, output, session) {
  options(shiny.maxRequestSize=30*1024^2) 
  
  #############################################################################################################
  ##################################### Per primo pannello ###################################################
 observe({
   inFile <- input$otherfile1
   if(is.null(inFile)){
     file <- file  
     
   } else{
     numfiles = nrow(inFile)                
     kata = list()
     
     for(i in 1:numfiles){
       files = read.delim(input$otherfile1[[i,'datapath']],header = TRUE,stringsAsFactors = FALSE)
       kata[[i]] = files
     }
     finals <- do.call(rbind,kata)
     file <- file %>% 
       rbind(finals)
     
   }
    
    rawdata <- file %>%       
      group_by(data,class,type,cna.data,snv.data) %>%
      summarise(n.samples=n_distinct(sample.id),n.patients=n_distinct(patient.id))
    risorse <- rawdata$data[!duplicated(rawdata$data)]
    updateCheckboxGroupInput(session, 'Resources1', choices = risorse, selected = risorse)})
  
  observeEvent(input$Reset1,{
    shinyjs::reset("Resources1")
    shinyjs::reset("ALL")
    shinyjs::reset("CNA")
    shinyjs::reset("SNV")
    shinyjs::reset("Types1")
    shinyjs::reset("Class1")
    })
  
  observeEvent(input$cancel1,{
    updateCheckboxGroupInput(session, 'Resources1', choices = c('brca_metabric','brca_igr_2015','brca_mbcproject_wagle_2017','breast_msk_2018','brca_tcga_pan_can_atlas_2018'), selected = c('brca_metabric','brca_igr_2015','brca_mbcproject_wagle_2017','breast_msk_2018','brca_tcga_pan_can_atlas_2018'))
    shinyjs::reset("otherfile1") # non sembrano esserci rimasugli di file updated da test 
  })
  
  manipulation1 <-reactive({
    inFile1 <- input$otherfile1
    if(is.null(inFile1)){
      x1 <- file  
      
    } else{
      numfiles = nrow(inFile1)                
      kata2 = list()
      
      for(i in 1:numfiles){
        files = read.delim(input$otherfile1[[i,'datapath']],header = TRUE,stringsAsFactors = FALSE)
        kata2[[i]] = files
      }
      finals <- do.call(rbind,kata2)
      x1 <- file %>% 
        rbind(finals)
    }
    
    if(input$ALL == TRUE){
      x1 <- x1 %>%
        filter(snv.data == TRUE | cna.data == TRUE)
    }else if (input$SNV == TRUE & input$CNA == FALSE){
      x1 <- x1 %>% 
        filter(snv.data == TRUE)
    }else if (input$CNA == TRUE & input$SNV == FALSE) {
      x1 <- x1 %>% 
        filter(cna.data == TRUE)
    }else if (input$SNV & input$CNA){
      x1 <- x1 %>% 
        filter(cna.data == TRUE & snv.data ==TRUE)
    }
   
    x1 <- x1 %>%
      subset(data %in% input$Resources1) %>%
      subset(type %in% input$Types1) %>%
      subset(class %in% input$Class1)
    
    return(x1)
    
  })
  
  newData <-eventReactive(input$updatebutton1,ignoreNULL = F,ignoreInit = F,{ #eventReactive(input$updatebutton1,ignoreNULL = F,ignoreInit = F,
    data1 <- manipulation1()
    
    data1 <- data1%>%       
      group_by(data,class,type,cna.data,snv.data) %>%
      summarise(n.samples=n_distinct(sample.id),n.patients=n_distinct(patient.id))
    
    if(nrow(data1) == 0){
      df <- data.frame()
      ggplot(df) + geom_point() + xlim(0, 10) + ylim(0, 10) + 
        annotate("text", x=3.9, y=5.0, size=40, col="red", label="(" ) +
        annotate("text", x=5, y=5.6, size=12, col="red", label="o  o" ) +
        annotate("text", x=6.1, y=5.0, size=40, col="red", label=")" ) +
        annotate("text", x=5, y=5.1, size=12, col="red", label="|" ) +
        geom_segment(aes(x = 4.7, xend = 5.3, y = 4.4, yend = 4.4), size=2, color="red") +
        annotate("text", x=5, y=3, size=10, col="red", label="No Data") 
    }else{
      ggplot(data1, aes(x=factor(type, levels = c('HR+','HER2+','TNBC')), y=n.samples, fill=factor(class, levels = c('Primary','Metastasis')))) +
        geom_bar(stat="identity") + 
        scale_fill_manual(values=c('#999999','#E69F00'))+
        theme(text = element_text(size=18))+ 
        xlab('Type')+
        labs(fill = 'Classes')
      }
  })
  
  newData2 <-eventReactive(input$updatebutton1,ignoreNULL = F,ignoreInit = F,{
    data2 <- manipulation1()
    
    data2 <- data2 %>%
      group_by(class,type) %>%
      summarise(n.samples= n_distinct(sample.id),n.patients=n_distinct(patient.id)) 
    if(nrow(data2) == 0){
      df <- data.frame()
      ggplot(df) + geom_point() + xlim(0, 10) + ylim(0, 10) + 
        annotate("text", x=3.9, y=5.0, size=40, col="red", label="(" ) +
        annotate("text", x=5, y=5.6, size=12, col="red", label="o  o" ) +
        annotate("text", x=6.1, y=5.0, size=40, col="red", label=")" ) +
        annotate("text", x=5, y=5.1, size=12, col="red", label="|" ) +
        geom_segment(aes(x = 4.7, xend = 5.3, y = 4.4, yend = 4.4), size=2, color="red") +
        annotate("text", x=5, y=3, size=10, col="red", label="No Data") 
    }else{
      ggplot(data2, aes(x=factor(type, levels = c('HR+','HER2+','TNBC')), y=n.samples,fill=factor(class, levels = c('Primary','Metastasis')))) +
        geom_bar(stat="identity") + #theme(aspect.ratio = 1,legend.position = "none") +
        scale_fill_manual(values=c('#999999','#E69F00')) +
        geom_text(aes(label=n.samples),size =5) + 
        facet_wrap(~factor(class,levels = c('Primary','Metastasis')))+
      theme(text = element_text(size=18))+
      xlab('Type')+
      labs(fill = 'Classes')}
    
  })
  
  newData_table <-eventReactive(input$updatebutton1,ignoreNULL = F,ignoreInit = F,{
    data3 <- manipulation1()
    
    data3 <- data3%>%       
      group_by(data,class,type,cna.data,snv.data) %>%
      summarise(n.samples=n_distinct(sample.id),n.patients=n_distinct(patient.id))
    
    all2 <- data3 %>%
      group_by(class,type) %>%
      summarise(n.samples_sum=sum(n.samples),n.patients_sum=sum(n.patients))%>% 
      add_column(data='all',.before = 'class')
    all2 <- all2 %>%
      rename(n.samples =n.samples_sum, n.patients = n.patients_sum)
    data3 <- data3 %>%
      rbind(all2)
  })
  
  information_1 <- eventReactive(input$info1,{
    print('It is posssible to upload more files up to a maximum of 30MB, for uploading more files at the same time they need to be selected together')
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
  
  #############################################################################################################
  ##################################### Per secondo pannello ###################################################  
  
  observe({
  inFile2 <- input$otherfile2
  if(is.null(inFile2)){
    x <- file2
    
  } else{
    numfiles = nrow(inFile2)                
    kata = list()
    
    for(i in 1:numfiles){
      files = read.delim(input$otherfile2[[i,'datapath']],header = TRUE,stringsAsFactors = FALSE)
      kata[[i]] = files
    }
    finals <- do.call(rbind,kata)
    x <- file2 %>% 
      rbind(finals)
  }
    rawdatasecpan <- x %>%       
      group_by(data)
    risorse <- rawdatasecpan$data[!duplicated(rawdatasecpan$data)]
    updateCheckboxGroupInput(session, 'Resources2', choices = risorse, selected = risorse)})
  
  observeEvent(input$Reset2,{
    shinyjs::reset("Resources2")
    shinyjs::reset("Gene_filter")
    shinyjs::reset("Types2")
    shinyjs::reset("Class2")
  })
  
  observeEvent(input$cancel2,{
    updateCheckboxGroupInput(session, 'Resources2', choices = c('brca_metabric','brca_igr_2015','brca_mbcproject_wagle_2017','breast_msk_2018','brca_tcga_pan_can_atlas_2018'), selected = c('brca_metabric','brca_igr_2015','brca_mbcproject_wagle_2017','breast_msk_2018','brca_tcga_pan_can_atlas_2018'))
    shinyjs::reset("otherfile2") 
    })
  
  manipulation2 <- reactive({
    inFile2 <- input$otherfile2
    if(is.null(inFile2)){
      x <- file2
      
    } else{
      numfiles = nrow(inFile2)                
      kata = list()
      
      for(i in 1:numfiles){
        files = read.delim(input$otherfile2[[i,'datapath']],header = TRUE,stringsAsFactors = FALSE)
        kata[[i]] = files
      }
      finals <- do.call(rbind,kata)
      x <- file2 %>% 
        rbind(finals)
      
    }
    
    x <- filter(x, data %in% input$Resources2)  #input$Resources2
    m <- x %>% 
        group_by(Hugo_Symbol,class,type) %>% 
        summarise(max.freq=max(freq),w.mean=weighted.mean(x=freq, w = n.samples),median.freq = median(freq)) %>% 
        ungroup() %>% 
        group_by(class,type)

    mw <- m %>%
        arrange(desc(w.mean)) 
  
    return(mw)
  })
  #data per barplot
  data_second_pannel <- eventReactive(input$updatebutton2,ignoreNULL = F,ignoreInit = F,{
    
    data_pan_1 <- manipulation2()
    data_pan_1 <- filter(data_pan_1, class %in% input$Class2)
    data_pan_2 <- data_pan_1 %>% 
      filter(type == 'all')
    data_pan_1 <- data_pan_1 %>% 
      filter(type != 'all')
    data_pan_1 <- filter(data_pan_1, type %in% input$Types2)
    
    data_pan_3 <- data_pan_1 %>% 
      rbind(data_pan_2)
    data_pan_3 <- data_pan_3 %>%
      slice_head(n = input$Gene_filter)  %>%
      group_split()  
    
    plist <- list()
    if(length(data_pan_3) == 0){
      df <- data.frame()
      ggplot(df) + geom_point() + xlim(0, 10) + ylim(0, 10) +
        annotate("text", x=3.9, y=5.0, size=40, col="red", label="(" ) +
        annotate("text", x=5, y=5.6, size=12, col="red", label="o  o" ) +
        annotate("text", x=6.1, y=5.0, size=40, col="red", label=")" ) +
        annotate("text", x=5, y=5.1, size=12, col="red", label="|" ) +
        geom_segment(aes(x = 4.7, xend = 5.3, y = 4.4, yend = 4.4), size=2, color="red") +
        annotate("text", x=5, y=3, size=10, col="red", label="No Data")
    }else{
      for(i in 1:length(data_pan_3)){
        dn <- data_pan_3[[i]]
        dn$Hugo_Symbol <- factor(dn$Hugo_Symbol,levels = rev(dn$Hugo_Symbol))
        p<- ggplot(data=dn, aes(x=Hugo_Symbol, y=w.mean)) +
          geom_bar(stat="identity") + coord_flip() +
          facet_wrap(type~class,scales = 'free') +
          theme(text = element_text(size=18))
        
        plist[[i]] <- ggplotGrob(p)
      }
      grid.arrange(grobs=plist,ncol=4)
      
    }
  })
  
  data_second_pannel_heatmap <- eventReactive(input$updatebutton2,ignoreNULL = F,ignoreInit = F,{
    mw <- manipulation2()
    if(nrow(mw) == 0){
      df <- data.frame()
      ggplot(df) + geom_point() + xlim(0, 10) + ylim(0, 10) +
        annotate("text", x=3.9, y=5.0, size=40, col="red", label="(" ) +
        annotate("text", x=5, y=5.6, size=12, col="red", label="o  o" ) +
        annotate("text", x=6.1, y=5.0, size=40, col="red", label=")" ) +
        annotate("text", x=5, y=5.1, size=12, col="red", label="|" ) +
        geom_segment(aes(x = 4.7, xend = 5.3, y = 4.4, yend = 4.4), size=2, color="red") +
        annotate("text", x=5, y=3, size=10, col="red", label="No Data")
    }else{
      mw2 <- mw %>%
        slice_head(n = input$Gene_filter) %>%
        group_split()
      mat <- do.call(rbind,mw2)    
      mat <- filter(mat, type %in% input$Types2)
      mat <- filter(mat, class %in% input$Class2)
      
      ggplot(mat, aes(factor(type, levels = c('HR+','HER2+','TNBC')), Hugo_Symbol)) +
        geom_tile(aes(fill = w.mean)) + 
        geom_text(aes(label = round(w.mean, 3)),size=5) +
        scale_fill_gradient(low = "white", high = "red") +
        facet_wrap(~factor(class,levels = c('Primary','Metastasis'))) +
        theme(text = element_text(size=18))+
        xlab('Type')+
        theme(panel.background = element_rect(fill = '#599ad3')) #'#003f5c'
    }
  })
  
  
  data_barplot_heat <- eventReactive(input$updatebutton2, ignoreNULL= F, ignoreInit = F, {
    
    mw <- manipulation2()
    
    if(nrow(mw) == 0 ){
      df <- data.frame()
      ggplot(df) + geom_point() + xlim(0, 10) + ylim(0, 10) +
        annotate("text", x=3.9, y=5.0, size=40, col="red", label="(" ) +
        annotate("text", x=5, y=5.6, size=12, col="red", label="o  o" ) +
        annotate("text", x=6.1, y=5.0, size=40, col="red", label=")" ) +
        annotate("text", x=5, y=5.1, size=12, col="red", label="|" ) +
        geom_segment(aes(x = 4.7, xend = 5.3, y = 4.4, yend = 4.4), size=2, color="red") +
        annotate("text", x=5, y=3, size=10, col="red", label="No Data")
    }else{
      
      mw2 <- mw %>%
        slice_head(n = input$Gene_filter ) %>% #input$Gene_filter
        group_split()
      mat <- do.call(rbind,mw2)    
      mat <- filter(mat, type %in% input$Types2)
      mat <- filter(mat, class %in% input$Class2)
      mat <- mat %>% 
            arrange(desc(Hugo_Symbol))
      
      ggplot(data=mat, aes(x=Hugo_Symbol, y=w.mean, fill= factor(class, levels = c('Primary','Metastasis')))) +
        geom_bar(stat="identity", position = 'dodge') + coord_flip() +
        scale_fill_manual(values=c('#599ad3','#f9a65a')) +   #999999 #E69F00
        facet_wrap(~factor(type,levels = c('HR+','HER2+','TNBC')))+
        theme(text = element_text(size=18))+
        labs(fill = 'Classes')
    }
  })
  
  data_second_pannel_table <- eventReactive(input$updatebutton2,ignoreNULL = F,ignoreInit = F,{
    mw <- manipulation2()
    if(nrow(mw)==0){
      mw
    }else{
      mw3<- mw %>%
        slice_head(n = input$Gene_filter) %>%
        group_split()
      mat1 <- do.call(rbind,mw3)  
      mat1 <- mat1[order(-mat1$w.mean),]  
      mat1 <- subset(mat1, type %in% input$Types2)
      mat1 <- subset(mat1, class %in% input$Class2)
      mat1 <- mat1 %>%
        filter(type != 'all')}
  })
  
  information_2 <- eventReactive(input$info2,{
    print('It is posssible to upload more files up to a maximum of 30MB, for uploading more files at the same time they need to be selected together')
  })
  
  #############################################################################################################
  ##################################### Per terzo pannello ###################################################  
  
  observe({
    if(input$copynumber_granularity== FALSE){
      updateRadioButtons(session,'Groups3',choices =c('homozygous deletion' = 'homodel','hemizygous deletion' = 'hemidel', 'gain' = 'gain','amplification' = 'ampl'),selected = 'ampl')
    }else{
      updateRadioButtons(session,'Groups3',choices =c('deletion' = 'homodel','amplification' = 'ampl'),selected = 'ampl')
    }
  })
  
  colnames(ensembl) <- c('ensg','start','end','chr','band','Hugo_Symbol')
  ensembl$band <- paste0(ensembl$chr,ensembl$band)
  
  check.ensembl <- ensembl %>%
    group_by(band) %>% 
    summarise(n=n()) %>% 
    arrange(n)
  
  observe({
    inFile3 <- input$otherfile3
    if(is.null(inFile3)){
      
      scna_data <- scna_data
      
    }else{
      numfiles = nrow(inFile3)                
      kata = list()
      
      for(i in 1:numfiles){
        files = read.delim(input$otherfile3[[i,'datapath']],header = TRUE,stringsAsFactors = FALSE)
        kata[[i]] = files
      }
      finals <- do.call(rbind,kata)
      scna_data <- scna_data %>% 
                   rbind(finals)
    }
    
    rawdataterzpan <- scna_data %>%       
      group_by(Data)
    resources <- rawdataterzpan$Data[!duplicated(rawdataterzpan$Data)]
    updateCheckboxGroupInput(session, 'Resources3', choices = resources, selected = resources)})
  
  
  observeEvent(input$cancel3,{
    updateCheckboxGroupInput(session, 'Resources3', choices = c('brca_metabric','brca_igr_2015','brca_mbcproject_wagle_2017','brca_tcga_pan_can_atlas_2018'), selected = c('brca_metabric','brca_igr_2015','brca_mbcproject_wagle_2017','brca_tcga_pan_can_atlas_2018'))
    shinyjs::reset("otherfile3") 
    })
  
  manipulation3 <- eventReactive(input$updatebutton3,ignoreNULL = F,ignoreInit = F,{ #eventReactive(input$updatebutton3,ignoreNULL = F,ignoreInit = F,
    inFile3 <- input$otherfile3
    if(is.null(inFile3)){
      
      scna_data <- scna_data
      
    }else{
      numfiles = nrow(inFile2)                
      kata = list()
      
      for(i in 1:numfiles){
        files = read.delim(input$otherfile2[[i,'datapath']],header = TRUE,stringsAsFactors = FALSE)
        kata[[i]] = files
      }
      finals <- do.call(rbind,kata)
      scna_data <- scna_data %>% 
        rbind(finals)
    }
    
    datasets <- scna_data %>%
      na.omit() %>%
      filter(Data %in% input$Resources3) %>% 
      filter(Type %in% input$Types3) %>% #input$Types3
      filter(Class %in% input$Class3) #input$Class3
    
    if(input$copynumber_granularity == TRUE){  #input$copynumber_granularity
      datasets$scna[datasets$scna %in% c(-1,-2)] <- 'homodel'
      datasets$scna[datasets$scna %in% c(1,2)] <- 'ampl'
      datasets$scna[datasets$scna == 0] <- 'neutral'
    } else{ #input$copynumber_granularity
      datasets$scna[datasets$scna == -2] <- 'homodel'
      datasets$scna[datasets$scna == 2] <- 'ampl'
      datasets$scna[datasets$scna == 0] <- 'neutral'
      datasets$scna[datasets$scna == 1] <- 'gain'
      datasets$scna[datasets$scna == -1] <- 'hemidel'
    }
    
    all_brca <- datasets %>% 
      group_by(Hugo_symbol,scna) %>% 
      summarise(n= sum(n_samples)) %>% 
      mutate(Data = 'all_brca') 
    
    all_data <- datasets %>% 
      group_by(Hugo_symbol,scna,Data) %>% 
      summarise(n= sum(n_samples)) %>% 
      select(Hugo_symbol,scna,n,Data) 
    
    total <- all_data %>%      
      group_by(Hugo_symbol) %>% 
      summarise(tot = sum(n))   
    
    pre_final <- rbind(all_brca,all_data) 
    
    final <- full_join(pre_final,total,by= "Hugo_symbol") 
    
    out <- data.frame(Hugo_Symbol = final$Hugo_symbol,
                      scna = final$scna,
                      freq = final$n/final$tot,
                      data = final$Data
    ) 
    
    frq <- left_join(x = out, y=ensembl, by = 'Hugo_Symbol') %>%
      filter(!is.na(ensg)) %>%
      filter(chr != 'Y')
    
    if(input$Groups3 == 'ampl'){
      frq<- frq %>% 
        filter(scna == 'ampl')
    }else if(input$Groups3 == 'homodel'){
      frq <- frq %>% 
        filter(scna == 'homodel')
    }else if(input$Groups3 == 'gain'){
      frq <- frq %>% 
        filter(scna == 'gain')
    }else if(input$Groups3 == 'hemidel'){
      frq <- frq %>% 
        filter(scna == 'hemidel')
    }
    
    return(frq)
  })
  
  manipulation4 <-reactive({#eventReactive(input$updatebutton3,ignoreNULL = F,ignoreInit = F,
    frq <- manipulation3()
    
    d <- frq %>%
      group_by(data,scna,chr,band) %>% 
      summarise(n=n(),
                median.freq=median(freq,na.rm = TRUE),
                max=max(freq,na.rm = TRUE), 
                max.name=Hugo_Symbol[which.max(freq)]) %>% 
      mutate(is.goi = max.name %in% goi)
    
    d$arm <- rep('p',nrow(d))
    d$arm[grep(d$band,pattern =  'q')] <- 'q'
    
    br <-d %>% 
      filter(median.freq >= input$filter_median_freq) %>%       #input$filter_median_freq
      add_column(max.name.goi = NA) 
    
    br$max.name.goi[which(br$is.goi)] <- br$max.name[which(br$is.goi)]
    
    return(br)
    
  })
  
  plotting2 <- eventReactive(input$updatebutton3,ignoreNULL = F,ignoreInit = F,{#eventReactive(input$updatebutton3,ignoreNULL = F,ignoreInit = F,{
    
    br<- manipulation4() 
    if(nrow(br)==0){
      df <- data.frame()
      ggplot(df) + geom_point() + xlim(0, 10) + ylim(0, 10) +
        annotate("text", x=3.9, y=5.0, size=40, col="red", label="(" ) +
        annotate("text", x=5, y=5.6, size=12, col="red", label="o  o" ) +
        annotate("text", x=6.1, y=5.0, size=40, col="red", label=")" ) +
        annotate("text", x=5, y=5.1, size=12, col="red", label="|" ) +
        geom_segment(aes(x = 4.7, xend = 5.3, y = 4.4, yend = 4.4), size=2, color="red") +
        annotate("text", x=5, y=3, size=10, col="red", label="No Data")
    }else{
      ggplot(br%>% filter(data == 'all_brca')
             ,aes(x=band,y=median.freq,fill=arm)) +
        ylab(paste(input$Groups3 ,paste('median.freq by cytoband',collapse = ' '))) +   # come cambaire didascali con aggiornatmento
        geom_bar(stat = 'identity') +
        facet_wrap(~factor(chr,levels = c('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X')),scales = 'free_x')+
        scale_fill_manual('arm',values = wes_palette("Chevalier1",n = 2)) +
        theme(axis.text.x = element_blank()) +
        scale_color_manual('data',values = wes_palette("GrandBudapest1",n = 4)) +
        ggtitle(paste('class:',paste(input$Class3,collapse = ','),'\ntype: ',paste(input$Types3,collapse = ','))) +
        theme(text = element_text(size=25))
    } 
  })
  
  plotting <- eventReactive(input$updatebutton3,ignoreNULL = F,ignoreInit = F,{#eventReactive(input$updatebutton3,ignoreNULL = F,ignoreInit = F,{
    
    br<- manipulation4() 
  
    br <- br %>%
        filter(chr == input$Chromosomes)

    if(nrow(br)==0){
      df <- data.frame()
      ggplot(df) + geom_point() + xlim(0, 10) + ylim(0, 10) +
        annotate("text", x=3.9, y=5.0, size=40, col="red", label="(" ) +
        annotate("text", x=5, y=5.6, size=12, col="red", label="o  o" ) +
        annotate("text", x=6.1, y=5.0, size=40, col="red", label=")" ) +
        annotate("text", x=5, y=5.1, size=12, col="red", label="|" ) +
        geom_segment(aes(x = 4.7, xend = 5.3, y = 4.4, yend = 4.4), size=2, color="red") +
        annotate("text", x=5, y=3, size=10, col="red", label="No Data")
    }else{
      ggplot(br%>% filter(data == 'all_brca'),aes(x=band,y=median.freq,fill=arm)) +
        ylab(paste(input$Groups3 ,paste('median.freq by cytoband',collapse = ' '))) +   # come cambaire didascali con aggiornatmento
        geom_bar(stat = 'identity') +
        facet_wrap(~factor(chr,levels = c('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X')),scales = 'free_x')+
        scale_fill_manual('arm',values = wes_palette("Chevalier1",n = 2)) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size = 16)) +
        geom_point(data = br %>% filter(data != 'all_brca'),mapping = aes(x=band,y=median.freq,color=data),size=6) +  # fare controllo data con deleteion, che in quel caso al atogliere di daTA DA WIDJET NE SPUNTAVANO ALTRI 
        scale_color_manual('data',values = wes_palette("GrandBudapest1",n = 4)) +
        ggtitle(paste('class:',paste(input$Class3,collapse = ','),'\ntype: ',paste(input$Types3,collapse = ','))) +
        geom_point(data =br %>% filter( data == 'all_brca'),mapping = aes(x=band,y=max),shape=4,size=6) +
        geom_text(data = br %>% filter( data == 'all_brca'),mapping = aes(x=band,y=max,label=max.name.goi),size=10,angle=90,hjust=0,nudge_y=0.01)+
        theme(text = element_text(size=25))
    } 
  })
  
  manipulation_cytoband2 <- reactive({
    
    br <- manipulation4()
    
    if(input$Chromosomes != 'All'){
      br <- br %>%
        filter(chr == input$Chromosomes)
    }else{
      br <- br
    }
    
    br$max.name.goi[which(br$is.goi)] <- br$max.name[which(br$is.goi)]
    
    bsel <- br %>%
      pull(band) %>%
      str_sort(numeric = TRUE) %>% 
      unique() 
    
    return(bsel)
  })
  
  
  observe({
    updateSelectInput(session, 'Cytoband', choices = manipulation_cytoband2())
  })
  
  manipulation_cytoband <- reactive({
    
    frq <- manipulation3()
    
    gfrq <- frq %>%
      filter(band == input$Cytoband) %>%
      filter(scna == input$Groups3) %>%
      mutate(is.goi = Hugo_Symbol %in% goi) %>%
      arrange(chr,start,end)
    
    return(gfrq)
    
  })
  
  plotting3 <-eventReactive(input$updatebutton3,ignoreNULL = F,ignoreInit = F,{ #eventReactive(input$updatebutton3,ignoreNULL = F,ignoreInit = F,
    
    gfrq <- manipulation_cytoband() 
    if(nrow(gfrq)==0){
      df <- data.frame()
      ggplot(df) + geom_point() + xlim(0, 10) + ylim(0, 10) +
        annotate("text", x=3.9, y=5.0, size=40, col="red", label="(" ) +
        annotate("text", x=5, y=5.6, size=12, col="red", label="o  o" ) +
        annotate("text", x=6.1, y=5.0, size=40, col="red", label=")" ) +
        annotate("text", x=5, y=5.1, size=12, col="red", label="|" ) +
        geom_segment(aes(x = 4.7, xend = 5.3, y = 4.4, yend = 4.4), size=2, color="red") +
        annotate("text", x=5, y=3, size=10, col="red", label="No Data")
    }else{
    ggplot(gfrq %>%
               filter(data == 'all_brca') %>%
               arrange(start,end) %>% 
               distinct(Hugo_Symbol, .keep_all = TRUE) %>%
               mutate(Hugo_Symbol=factor(Hugo_Symbol, levels = Hugo_Symbol)),
             aes(x=Hugo_Symbol,y=freq,fill=is.goi)) +
        ylab(paste(input$Groups3 ,paste('freq', collapse=' '))) + 
        geom_bar(stat = 'identity') +
        facet_wrap(~band,scales = 'free_x') +
        scale_x_discrete(guide = guide_axis(n.dodge=2)) +
        scale_fill_manual('genes of interest',values = wes_palette("Royal1",n = 2)) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size = 18)) +
        ggtitle(paste('class:',paste(input$Class3,collapse = ','),'\ntype: ',paste(input$Types3,collapse = ','))) +
        geom_point(data = gfrq %>% filter(data != 'all_brca'),mapping = aes(x=Hugo_Symbol,y=freq,color=data),size=6) +
        scale_color_manual('data',values = wes_palette("GrandBudapest1",n = 4))+
        theme(text = element_text(size=25))
    }
  })
  
  
  chromosome_table <- eventReactive(input$updatebutton3,ignoreNULL = F,ignoreInit = F,{
    br <- manipulation4()
    
    if(input$Chromosomes != 'All'){
      br <- br %>%
        filter(chr == input$Chromosomes)
    }else{
      br <- br
    }
    
    br$max.name.goi[which(br$is.goi)] <- br$max.name[which(br$is.goi)]
    
    br <- br %>%
      filter(data == 'all_brca')
    
  })
  
  cytoband_table <- eventReactive(input$updatebutton3,ignoreNULL = F,ignoreInit = F,{
    
    gfrq <- manipulation_cytoband()
    
    gfrq <-gfrq %>% 
      filter(data == 'all_brca') %>%
      arrange(start,end) %>%
      distinct(Hugo_Symbol, .keep_all = TRUE) %>%
      mutate(Hugo_Symbol=factor(Hugo_Symbol, levels = Hugo_Symbol))
    
  })
  
  information_3 <- eventReactive(input$info3,{
    print('It is posssible to upload more files up to a maximum of 30MB, for uploading more files at the same time they need to be selected together')
  })
  
  ############################################################################################################
  ##################################### Per primo pannello ###################################################
  output$textinfo1 <- renderText({
    information_1()
  })
  
  #plotting
  output$myplot <- renderPlot({
    newData()
      })
  
  output$classplot <- renderPlot({
    newData2()
  })
  
  #DataTable
  output$classtable <- renderDataTable({newData_table()})
  
  #PER DOWNLOAD
  output$download_myplot <- downloadHandler(
    filename = function(){
      paste('Count_plot','.pdf',sep = '')},
    content = function(file){
      ggsave(file,newData(), width = 20, height = 20)
    }
  )
  
  output$download_classplot <- downloadHandler(
    filename = function(){
      paste('Metastis_primary_plot','.pdf',sep = '')
    },
    content = function(file){
      ggsave(file,newData2(), width = 20,height = 20)
    }
  )
  
  output$download_classtable <-  downloadHandler(
    filename = function(){
      paste('Count_table','.csv',sep = '')
    },
    content = function(file){
      write.csv(newData_table(),file, row.names = FALSE, col.names = T, sep = ',')
    }
  )
  
  #############################################################################################################
  ##################################### Per secondo pannello ###################################################
  output$textinfo2 <- renderText({
    information_2()
  })
  
  output$barplot <- renderPlot({
    data_second_pannel()
  }, height = 800)
  
  output$heatmap <- renderPlot({ data_second_pannel_heatmap()
  }, height = 900)
  
  output$barplot_heat <- renderPlot({
    data_barplot_heat()
  }, height = 1000)
  
  output$table2 <- renderDataTable({data_second_pannel_table()})
  
  output$download_barplot<-  downloadHandler(
    filename = function(){
      paste('SNV_barplot','.pdf',sep = '')
    },
    content = function(file){
      ggsave(file,data_second_pannel(), width = 30, height = 30)
      
    }
  )
  
  output$download_heatmap <-  downloadHandler(
    filename = function(){
      paste('SNV_heatmap','.pdf',sep = '')
    },
    content = function(file){
      ggsave(file,data_second_pannel_heatmap(), width = 30, height = 30)
    }
  )
  
  output$download_table2 <-  downloadHandler(
    filename = function(){
      paste('SNV_table','.csv',sep = '')
    },
    content = function(file){
      write.csv(data_second_pannel_table(),file, row.names = FALSE, col.names = T, sep = ',')
    }
  )
  #############################################################################################################
  ##################################### Per terzo pannello ###################################################
  output$textinfo3 <- renderText({
    information_3()
  })
  
  output$All_plot <- renderPlot({
    plotting2()},
    height = 800)
  
  output$plot <-renderPlot({
    plotting()
  },height = 800
  )
  output$cytoband <- renderPlot({
    plotting3()
  }, height = 800)
  
  output$table_chromosome <- renderDataTable({chromosome_table()})
  
  output$table_cytoband <- renderDataTable({cytoband_table()})
  
  output$download_plot_All <- downloadHandler(
    filename = function(){
      paste('All_chromosomes_plot','.pdf',sep='')
    },
    content = function(file){
      ggsave(file,plotting2(), width = 30, height = 30)
    }
  )
  
  output$download_plots<-  downloadHandler(
    filename = function(){
      paste('Chromosome_plot','.pdf',sep = '')
    },
    content = function(file){
      ggsave(file,plotting(), width = 30, height = 30)
    }
  )
  
  output$download_cytoband<-  downloadHandler(
    filename = function(){
      paste('Cytoband_plot','.pdf',sep = '')
    },
    content = function(file){
      ggsave(file,plotting3(),width = 30,height = 30)
    }
  )
  
  output$download_table_chromosome <-  downloadHandler(
    filename = function(){
      paste('CNA_table_chromosome','.csv',sep = '')
    },
    content = function(file){
      write.csv(chromosome_table(),file, row.names = FALSE, col.names = T, sep = ',')
    }
  )
  
  output$download_table_cytoband <-  downloadHandler(
    filename = function(){
      paste('CNA_table_cytoband','.csv',sep = '')
    },
    content = function(file){
      write.csv(cytoband_table(),file, row.names = FALSE, col.names = T, sep = ',')
    }
  )
  
}


# Run the application
shinyApp(ui = ui, server = server)



