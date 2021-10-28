library(shiny)
library(shinythemes)
library(tidyverse)
library(shinyjs) 
library(gridExtra)
library(data.table)
library(parallel)
library(wesanderson)
library(UpSetR)

#loading files
file <- read.delim('sif_cbioportal_brca.tsv', header = TRUE,stringsAsFactors = FALSE)
file2 <- read.delim('snvs_raw_data.tsv', header = TRUE, stringsAsFactors = FALSE)
ensembl <- read.delim('mart_export_GRCh38p13.tsv',check.names = F,stringsAsFactors = F)
goi <- readLines('genes_of_interest.txt')
# load('scna.RData') #aka dd
load('scna_data.RData')
# load('outaggregati.RData')
# load('outnoaggregat.RData')
# frq <- read.delim('data-esempio-implementazione.tsv',header = TRUE, stringsAsFactors = FALSE)
# d <- read.delim('data-esempio-implementazione2.tsv',header = TRUE, stringsAsFactors = FALSE)
# Define UI for application that draws a histogram

ui <- shinyUI(fluidPage(#shinythemes::themeSelector(),
                theme = shinytheme('superhero'),
                shinyjs::useShinyjs(), # questo mi serve per 'controllare' i vari click per gli input riguradanti snv e cna  
    # Application title
    titlePanel("App Breast Cancer"),
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
                inputId = 'Resources1',
                label = 'Data resources',
                choices = c('brca_metabric','brca_igr_2015','brca_mbcproject_wagle_2017','breast_msk_2018','brca_tcga_pan_can_atlas_2018'),
                selected = c('brca_metabric','brca_igr_2015','brca_mbcproject_wagle_2017','breast_msk_2018','brca_tcga_pan_can_atlas_2018')),
            checkboxGroupInput(
                inputId = 'Types1',
                label = 'Breast cancer subtypes',
                choices = c('HER2+','HR+','TNBC'),
                selected = c('HER2+','HR+','TNBC')),
            checkboxGroupInput(
                inputId = 'Class1',
                label = 'Tumor classification',
                choices = c('Primary','Metastasis'),
                selected = c('Primary','Metastasis')),
            actionButton("updatebutton1", "Proceed",icon("dna"),#paper-plane
                         style="color: #fff; background-color: #337ab7; border-color: #2e6da4")),
    conditionalPanel( # questo per creare un slider aggiuntivo quando si passa alla secondo pannello 
              condition = 'input.tabs== 2',
              fileInput(
                inputId = 'otherfile2',
                label = 'choose tsv file',
                multiple = TRUE,
                accept = '.tsv'),
              checkboxGroupInput(
                inputId = 'Resources2',  
                label = 'Data resources',
                choices = c('brca_metabric','brca_igr_2015','brca_mbcproject_wagle_2017','breast_msk_2018','brca_tcga_pan_can_atlas_2018'),
                selected = c('brca_metabric','brca_igr_2015','brca_mbcproject_wagle_2017','breast_msk_2018','brca_tcga_pan_can_atlas_2018')),
              checkboxGroupInput(
                inputId = 'Types2',
                label = 'Breast cancer subtypes',
                choices = c('HER2+','HR+','TNBC'),
                selected = c('HER2+','HR+','TNBC')),
              checkboxGroupInput(
                inputId = 'Class2',
                label = 'Tumor classification',
                choices = c('Primary','Metastasis'),
                selected = c('Primary','Metastasis')),
              sliderInput(inputId = 'Gene_filter',label = 'Gene_filter', min = 5, max = 25, value = 25),
              actionButton("updatebutton2", "Proceed",icon("dna"),#paper-plane
                           style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
              ),
    conditionalPanel(
      condition = 'input.tabs== 3',
      fileInput(
        inputId = 'otherfile3',
        label = 'choose tsv file',
        multiple = TRUE,
        accept = '.tsv'),
      checkboxGroupInput(
        inputId = 'Resources3',
        label = 'Data resources',
        choices = c('brca_metabric','brca_igr_2015','brca_mbcproject_wagle_2017','brca_tcga_pan_can_atlas_2018'),
        selected = c('brca_metabric','brca_igr_2015','brca_mbcproject_wagle_2017','brca_tcga_pan_can_atlas_2018')
        ),
      checkboxGroupInput(
        inputId = 'Types3',
        label = 'Breast cancer subtypes',
        choices = c('HER2+','HR+','TNBC'),
        selected = c('HER2+','HR+','TNBC')),
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
                         label = 'Choose',
                         choices = c('')),
      sliderInput(inputId = 'filter_median_freq',
                  label= 'Filter Median frequencing', min = 0, max= 1, value=0.02,step = 0.01),
      selectInput(inputId ='Chromosomes',
                         label = 'Chromosomes',
                         choices = c('All','1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X'),
                  selected = 'All'),
      selectInput(inputId = 'Cytoband',
                  label = 'Cytoband',
                  choices = c(''),
                  multiple = FALSE),
      actionButton("updatebutton3", "Proceed",icon("dna"),#paper-plane
                   style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
      )
        ),
        # Show a plot of the generated ditribution
        mainPanel(
              tabsetPanel(type= 'tabs', id = 'tabs',                           
                tabPanel(id='main', value= 1, 
                         'Count pannel',
                      downloadButton(outputId = 'download_myplot',label = 'Download count plot'),
                      downloadButton(outputId = 'download_classplot',label = 'Download primary-metastasis plot'),
                      downloadButton(outputId = 'download_classtable',label = 'Download table'),
                  fluidRow( # fluidrow serve come comnado per mettere dove vooglio i vari plot all'interno degli output
                    column(6,plotOutput(outputId = 'myplot')),
                    column(6,plotOutput(outputId = 'classplot' ))
                    ),
                  fluidRow(
                    column(6,dataTableOutput(outputId = 'classtable'))
                   )),
                tabPanel(id ='sec', value = 2,
                         'SNV pannel',
                         downloadButton(outputId = 'download_barplot',label = 'Download barplot'),
                         downloadButton(outputId = 'download_heatmap',label = 'Download heatmap'),
                         downloadButton(outputId = 'download_table2',label = 'Download table'),
                   fluidRow(
                     column(6,plotOutput(outputId = 'barplot')),
                     column(6,plotOutput(outputId = 'heatmap'))
                   ),
                   fluidRow(
                     column(6,dataTableOutput(outputId = 'table2'))
                   )),
                tabPanel(id = 'thrd', value = 3,
                         'CNA pannel',
                         downloadButton(outputId = 'download_plots', label = 'downlaod chromosome plot'),
                         downloadButton(outputId = 'download_cytoband',label = 'download cytoband plot'),
                         downloadButton(outputId = 'download_table_chromosome', label = 'Download table chromosome'),
                         downloadButton(outputId = 'download_table_cytoband', label = 'Download table cytoband'),
                  fluidRow(style='height:80vh',
                    column(12,plotOutput(outputId = 'plot',width = '125%'))
                  ),
                  fluidRow(
                    column(12,plotOutput(outputId = 'cytoband', width = '125%'))
                  ),
                  fluidRow( 
                    column(12,
                           tabsetPanel(
                             tabPanel('Chromosome',
                                      fluidRow(
                                        column(12,dataTableOutput(outputId = 'table_chromosome')))),
                             tabPanel('Cytoband',
                                      fluidRow(
                                        column(12,dataTableOutput(outputId = 'table_cytoband'))))
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

#############################################################################################################
##################################### Per primo pannello ###################################################
observe({
  if(is.null(input$otherfile1))
      {
        file <- file
      }
      else
      {
        uploaded <- input$otherfile1
        uploaded2 <- read.delim(file = uploaded$datapath, header = TRUE, stringsAsFactors = FALSE)
        file <- file %>%
          rbind(uploaded2)    }
  
  rawdata <- file %>%       
    group_by(data,class,type,cna.data,snv.data) %>%
    summarise(n.samples=n_distinct(sample.id),n.patients=n_distinct(patient.id))
  risorse <- rawdata$data[!duplicated(rawdata$data)]
  updateCheckboxGroupInput(session, 'Resources1', choices = risorse, selected = risorse)})

manipulation1 <-reactive({
  if(is.null(input$otherfile1))
  {
    x1 <- file
  }
  else
  {
    x1<- file
    uploaded <- input$otherfile1
    uploaded2 <- read.delim(file = uploaded$datapath, header = TRUE, stringsAsFactors = FALSE) 
    x1 <- x1 %>%
      rbind(uploaded2)  }
  return(x1)
})
  
manipualtion1 <- eventReactive(input$updatebutton1,{
  
  if(is.null(input$otherfile1))
  {
    x1 <- file
  }
  else
  {
    x1<- file
    uploaded <- input$otherfile1
    uploaded2 <- read.delim(file = uploaded$datapath, header = TRUE, stringsAsFactors = FALSE) 
    x1 <- x1 %>%
      rbind(uploaded2)  }
  
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
  x1 <- subset(file, data %in% input$Resources1)
  x1  <- subset(file, type %in% input$Types1)
  x1 <- subset(file, class %in% input$Class1)
 
  return(x1)
  
   })
  
newData <- reactive({
  data1 <- manipulation1()
  
  data1 <- data1%>%       
    group_by(data,class,type,cna.data,snv.data) %>%
    summarise(n.samples=n_distinct(sample.id),n.patients=n_distinct(patient.id))
    
    ggplot(data1, aes(x=type, y=n.samples, fill=class)) +
      geom_bar(stat="identity") + 
      scale_fill_manual(values=c('#999999','#E69F00'))
    })

newData2 <- reactive({
    data2 <- manipulation1()
    
    data2 <- data2 %>%
        group_by(class,type) %>%
        summarise(n.samples= n_distinct(sample.id),n.patients=n_distinct(patient.id)) 
   
     ggplot(data2, aes(x=type,y=n.samples,fill=class)) +
      geom_bar(stat="identity") + theme(aspect.ratio = 1,legend.position = "none") +
      scale_fill_manual(values=c('#999999','#E69F00')) +
      geom_text(aes(label=n.samples)) + 
      facet_wrap(~class)
    })
  
  newData_table <- reactive({
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
    if(is.null(input$otherfile2))
    {
      x <- file2
    }
    else
    {
      uploaded <- input$otherfile2
      uploaded2 <- read.delim(file = uploaded$datapath, header = TRUE, stringsAsFactors = FALSE)
      x <- file2 %>%
        rbind(uploaded2)
    }
    rawdatasecpan <- x %>%       
      group_by(data)
    risorse <- rawdatasecpan$data[!duplicated(rawdatasecpan$data)]
    updateCheckboxGroupInput(session, 'Resources2', choices = risorse, selected = risorse)})
  
  manipulation2 <- reactive({
    if(is.null(input$otherfile2))
    {
      x<- file2
      x <- subset(x, data %in% input$Resources2 )
      m <- x %>% 
        group_by(Hugo_Symbol,class,type) %>% 
        summarise(max.freq=max(freq),w.mean=weighted.mean(x = freq,w = n.samples),median.freq = median(freq)) %>% 
        ungroup() %>% 
        group_by(class,type) 
      mw <- m %>% 
        arrange(desc(w.mean))
    }
    else
    {
      x<- file2
      x <- subset(x, data %in% input$Resources2 )
      m <- x %>% 
        group_by(Hugo_Symbol,class,type) %>% 
        summarise(max.freq=max(freq),w.mean=weighted.mean(x = freq,w = n.samples),median.freq = median(freq)) %>% 
        ungroup() %>% 
        group_by(class,type) 
      
      uploaded3<- input$otherfile2
      uploaded4 <- read.delim(file = uploaded3$datapath, header = TRUE, stringsAsFactors = FALSE)
      uploaded5 <- uploaded4 %>% 
        group_by(Hugo_Symbol,class,type) %>% 
        summarise(max.freq=max(freq),w.mean=weighted.mean(x=freq, w = n.samples),median.freq = median(freq)) %>% 
        ungroup() %>% 
        group_by(class,type)
      
      mw <- m %>%
        rbind(uploaded5)
      mw <- mw %>%
        arrange(desc(w.mean)) 
      
      mw <- mw[-which(duplicated(mw)), ] 
    }
  })
  #data per barplot
  data_second_pannel <- reactive({
  
    data_pan_1 <- manipulation2()
    data_pan_1 <- subset(data_pan_1, class %in% input$Class2)
    data_pan_2 <- data_pan_1 %>% 
                  filter(type == 'all')
    data_pan_1 <- data_pan_1 %>% 
                  filter(type != 'all')
    data_pan_1 <- subset(data_pan_1, type %in% input$Types2)
    data_pan_3 <- data_pan_1 %>% 
                rbind(data_pan_2)
    data_pan_3 <- data_pan_3 %>%
                  slice_head(n = input$Gene_filter)  %>%
                  group_split()  
    
    plist <- list()
    for(i in 1:length(data_pan_3)){
      dn <- data_pan_3[[i]]
      dn$Hugo_Symbol <- factor(dn$Hugo_Symbol,levels = rev(dn$Hugo_Symbol))
      p<- ggplot(data=dn, aes(x=Hugo_Symbol, y=w.mean)) +
        geom_bar(stat="identity") + coord_flip() +
        facet_wrap(type~class,scales = 'free')
      plist[[i]] <- ggplotGrob(p)
    }
    grid.arrange(grobs=plist,ncol=4)
    
  })
  
data_second_pannel_heatmap <- reactive({
    mw <- manipulation2()
    mw2 <- mw %>%
          slice_head(n = input$Gene_filter) %>%
          group_split()
    mat <- do.call(rbind,mw2)
    mat <- subset(mat, type %in% input$Types2)
    mat <- subset(mat, class %in% input$Class2)
    mat <- mat %>%
      filter(type != 'all')
    
    ggplot(mat, aes(type, Hugo_Symbol)) +  
      geom_tile(aes(fill = w.mean)) + 
      geom_text(aes(label = round(w.mean, 3)),size=2) +
      scale_fill_gradient(low = "white", high = "red") +
      facet_wrap(~class)
  })
  
data_second_pannel_table <- reactive({
    mw <- manipulation2()
    mw3<- mw %>%
      slice_head(n = input$Gene_filter) %>%
      group_split()
    mat1 <- do.call(rbind,mw3)  
    mat1 <- mat1[order(-mat1$w.mean),]  
    mat1 <- subset(mat1, type %in% input$Types2)
    mat1 <- subset(mat1, class %in% input$Class2)
    mat1 <- mat1 %>%
      filter(type != 'all')
  })

#############################################################################################################
##################################### Per terzo pannello ###################################################  

observe({
  if(input$copynumber_granularity== TRUE){
    updateRadioButtons(session,'Groups3',choices =c('deletion' = 'homodel','amplification' = 'ampl'),selected = 'ampl')
  } else if(input$copynumber_granularity== FALSE){
    updateRadioButtons(session,'Groups3',choices =c('homozygous deletion' = 'homodel','hemizygous deletion' = 'hemidel', 'gain' = 'gain','amplification' = 'ampl'),selected = 'ampl')
  }
})
  
colnames(ensembl) <- c('ensg','start','end','chr','band','Hugo_Symbol')
ensembl$band <- paste0(ensembl$chr,ensembl$band)
  
check.ensembl <- ensembl %>%
    group_by(band) %>% 
    summarise(n=n()) %>% 
    arrange(n)

observe({
  if(is.null(input$otherfile3))
  {
    scna_data <- scna_data
  }
  else
  {
    scna_data<- scna_data
    uploaded <- input$otherfile3
    uploaded2 <- read.delim(file = uploaded$datapath, header = TRUE, stringsAsFactors = FALSE) 
    scna_data <- scna_data %>%
      rbind(uploaded2)  
  }
  rawdataterzpan <- scna_data %>%       
    group_by(Data)
  resources <- rawdataterzpan$Data[!duplicated(rawdataterzpan$Data)]
  updateCheckboxGroupInput(session, 'Resources3', choices = resources, selected = resources)})
  

# # selected_class <- c('Metastasis','Primary') # possible options: 'Metastasis' 'Primary'
# # 
# # selected_type <- c('HR+','HER2+','TNBC') # possible options: 'HR+','HER2+','TNBC'
# # 
# # selected_data <- c("brca_igr_2015","brca_mbcproject_wagle_2017", "brca_metabric","brca_tcga_pan_can_atlas_2018") #all without metabric test
# 

manipulation3 <- reactive({
  if(is.null(input$otherfile3))
  {
    scna_data <- scna_data
  }
  else
  {
    scna_data<- scna_data
    uploaded <- input$otherfile3
    uploaded2 <- read.delim(file = uploaded$datapath, header = TRUE, stringsAsFactors = FALSE) 
    scna_data <- scna_data %>%
      rbind(uploaded2)  
  }
  datasets <- scna_data %>%
              na.omit() %>%
              subset(Data %in% input$Resources3) %>% 
              subset(Type %in% input$Types3) %>% #input$Types3
              subset(Class %in% input$Class3) #input$Class3

   if(input$copynumber_granularity == TRUE){  #input$copynumber_granularity
     datasets$scna[datasets$scna %in% c(-1,-2)] <- 'homodel'
     datasets$scna[datasets$scna %in% c(1,2)] <- 'ampl'
     datasets$scna[datasets$scna == 0] <- 'neutral'
    } else if (input$copynumber_granularity  == FALSE){ #input$copynumber_granularity
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
 
manipulation4 <- reactive({
  frq <- manipulation3()
  
  d <- frq %>%
    group_by(data,scna,chr,band) %>% 
    summarise(n=n(),
              median.freq=median(freq,na.rm = TRUE),
              max=max(freq,na.rm = TRUE), 
              max.name=Hugo_Symbol[which.max(freq)]) %>% 
    mutate(is.goi = max.name %in% goi)
  
  d$arm <- rep('p',nrow(d))
  d$arm[grep(d$band,pattern = 'q')] <- 'q'
  
  br <-d %>% 
    filter(median.freq > input$filter_median_freq) %>%       #input$filter_median_freq
    add_column(max.name.goi = NA) 
    

  if(input$Chromosomes != 'All'){
    br <- br %>%
      filter(chr == input$Chromosomes)
  }else{
    br <- br
  }
  
  br$max.name.goi[which(br$is.goi)] <- br$max.name[which(br$is.goi)]

  return(br)
})
# observe({
#   br <- manipulation4()
#   selected<- br %>%
#     filter( data != 'all_brca') 
#   data_names <- selected$data%>%
#     unique() %>% 
#     droplevels()
#   
#   updateCheckboxGroupInput(session,'Resources3',selected = data_names )
# })

plotting <- reactive({
   
  br<- manipulation4() 
  
  ggplot(br%>% filter(data == 'all_brca'),aes(x=band,y=median.freq,fill=arm)) +
    ylab(paste(input$Groups3 ,paste('median.freq by cytoband',collapse = ' '))) +   # come cambaire didascali con aggiornatmento
    geom_bar(stat = 'identity') +
    facet_wrap(~factor(chr,levels = c('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X')),scales = 'free_x')+
    scale_fill_manual('arm',values = wes_palette("Chevalier1",n = 2)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size = 4)) +
    geom_point(data = br %>% filter(data != 'all_brca'),mapping = aes(x=band,y=median.freq,color=data),size=0.5) +  # fare controllo data con deleteion, che in quel caso al atogliere di daTA DA WIDJET NE SPUNTAVANO ALTRI 
    scale_color_manual('data',values = wes_palette("GrandBudapest1",n = 4)) +
    ggtitle(paste('class:',paste(input$Class3,collapse = ','),'\ntype: ',paste(input$Types3,collapse = ','))) +
    geom_point(data =br %>% filter( data == 'all_brca'),mapping = aes(x=band,y=max),shape=4,size=0.5) +
    geom_text(data = br %>% filter( data == 'all_brca'),mapping = aes(x=band,y=max,label=max.name.goi),size=1,angle=90,hjust=0,nudge_y=0.01)
  
})

manipulation_cytoband <- reactive({
  
  br <- manipulation4()
  
  bsel <- br %>%
          pull(band) %>%
          unique()
  return(bsel)
})

observe({
  updateSelectInput(session, 'Cytoband', choices = manipulation_cytoband() )
})

plotting2 <-reactive({
  
  frq <- manipulation3() 
  
  gfrq <- frq %>%
    filter(band == input$Cytoband) %>%
    filter(scna == input$Groups3) %>%
    mutate(is.goi = Hugo_Symbol %in% goi) %>%
    arrange(chr,start,end)

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
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    ggtitle(paste('class:',paste(input$Class3,collapse = ','),'\ntype: ',paste(input$Types3,collapse = ','))) +
    geom_point(data = gfrq %>% filter(data != 'all_brca'),mapping = aes(x=Hugo_Symbol,y=freq,color=data)) +
    scale_color_manual('data',values = wes_palette("GrandBudapest1",n = 4))
})


chromosome_table <- reactive({
  br <- manipulation4()
  
  br <- br %>%
    filter(data == 'all_brca')
})

cytoband_table <- reactive({
    
    frq <- manipulation3()
    gfrq <- frq %>%
      filter(band == input$Cytoband) %>%
      mutate(is.goi = Hugo_Symbol %in% goi) %>%
      arrange(chr,start,end) %>%''
      filter(data == 'all_brca') %>%
      arrange(start,end) %>%
      distinct(Hugo_Symbol, .keep_all = TRUE) %>%
      mutate(Hugo_Symbol=factor(Hugo_Symbol, levels = Hugo_Symbol))
  })
############################################################################################################
##################################### Per primo pannello ###################################################

      #plotting
    output$myplot <- renderPlot({
            newData()
        })
    output$classplot <- renderPlot({
      newData2()
       })

    #DataTable
    output$classtable <- renderDataTable(newData_table())

    #PER DOWNLOAD
    output$download_myplot <- downloadHandler(
                              filename = function(){'Count_plot.png'},
                              content = function(file){
                                ggsave(file,newData())
                              }
    )

    output$download_classplot <- downloadHandler(
                                  filename = function(){
                                    paste('Metastis_primary_plot.png')
                                  },
                                  content = function(file){
                                    ggsave(file,newData2())
      }
    )

    output$download_classtable <-  downloadHandler(
                                   filename = function(){
                                     paste('Count_table','.csv',sep = '')
                                   },
                                   content = function(file){
                                     write.csv(newData_table(),file)
      }
    )

#############################################################################################################
##################################### Per secondo pannello ###################################################

       output$barplot <- renderPlot({
          data_second_pannel()
    })

    output$heatmap <- renderPlot({ data_second_pannel_heatmap()
    })
    output$table2 <- renderDataTable(data_second_pannel_table())

    output$download_barplot<-  downloadHandler(
      filename = function(){
        paste('SNV_barplot','.pdf',sep = '')
      },
      content = function(file){
        ggsave(file,data_second_pannel())

      }
    )

    output$download_heatmap <-  downloadHandler(
      filename = function(){
        paste('SNV_heatmap','.pdf',sep = '')
      },
      content = function(file){
        ggsave(file,data_second_pannel_heatmap())
        }
    )

    output$download_table2 <-  downloadHandler(
      filename = function(){
        paste('SNV_table','.csv',sep = '')
      },
      content = function(file){
        write.csv(data_second_pannel_table(),file)
       }
    )
#############################################################################################################
##################################### Per terzo pannello ###################################################

    output$plot <-renderPlot({
      plotting()
      },height = 800
    )
    output$cytoband <- renderPlot({
      plotting2()
    })

    # problema codice riguardo tabelle e blocca anche pollting2 in qualche modo

   output$table_chromosome <- renderDataTable({chromosome_table()})

   output$table_cytoband <- renderDataTable({cytoband_table()})

    output$download_plots<-  downloadHandler(
      filename = function(){
        paste('Chromosome_plot','.pdf',sep = '')
      },
      content = function(file){
        ggsave(file,plotting())
      }
    )

  output$download_cytoband<-  downloadHandler(
    filename = function(){
      paste('Chromosome_plot','.pdf',sep = '')
    },
    content = function(file){
      ggsave(file,plotting2())
    }
    )

    output$download_table_chromosome <-  downloadHandler(
      filename = function(){
        paste('CNA_table_chromosome','.csv',sep = '')
      },
      content = function(file){
        write.csv(chromosome_table(),file)
      }
    )

    output$download_table_cytoband <-  downloadHandler(
      filename = function(){
        paste('CNA_table_cytoband','.csv',sep = '')
      },
      content = function(file){
        write.csv(cytoband_table(),file)
      }
    )

}


# Run the application
shinyApp(ui = ui, server = server)

