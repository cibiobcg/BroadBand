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
file2<- read.delim('snv_freq_brca.tsv', header = TRUE, stringsAsFactors = FALSE) %>% group_by(class,type) 
ensembl <- read.delim('mart_export_GRCh38p13.tsv',check.names = F,stringsAsFactors = F)
goi <- readLines('genes_of_interest.txt')
load('scna.RData')
frq <- read.delim('data-esempio-implementazione.tsv',header = TRUE, stringsAsFactors = FALSE)
d <- read.delim('data-esempio-implementazione2.tsv',header = TRUE, stringsAsFactors = FALSE)
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
                selected = c('Primary','Metastasis'))),
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
              sliderInput(inputId = 'Gene_filter',label = 'Gene_filter', min = 5, max = 25, value = 25)
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
        choices = c('brca_metabric','brca_igr_2015','brca_mbcproject_wagle_2017','breast_msk_2018','brca_tcga_pan_can_atlas_2018'),
        selected = c('brca_metabric','brca_igr_2015','brca_mbcproject_wagle_2017','breast_msk_2018','brca_tcga_pan_can_atlas_2018')),
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
      radioButtons(inputId = 'Groups',
                         label = 'Choose',
                         choices = c('deletion' = 'homodel','amplification' = 'ampl'),
                         selected = 'ampl'),
      sliderInput(inputId = 'filter_median_freq',
                  label= 'Filter Median frequencing', min = 0, max= 1, value=0.02,step = 0.01),
      selectInput(inputId ='Chromosomes',
                         label = 'Chromosomes',
                         choices = c('All','1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X'),
                  selected = 'All'),
      selectInput(inputId = 'Cytoband',
                  label = 'Cytoband',
                  choices = c(''),
                  multiple = FALSE)
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
        file<- file
        uploaded <- input$otherfile1
        uploaded2 <- read.delim(file = uploaded$datapath, header = TRUE, stringsAsFactors = FALSE)
        file <- file %>%
          rbind(uploaded2)
    }
  rawdata <- file %>%       
  group_by(data,class,type,cna.data,snv.data) %>%
  summarise(n.samples=n_distinct(sample.id),n.patients=n_distinct(patient.id))
risorse <- rawdata$data[!duplicated(rawdata$data)]
updateCheckboxGroupInput(session, 'Resources1', choices = risorse, selected = risorse)})
  
newData <- reactive({
  if(is.null(input$otherfile1))
  {
    file <- file
  }
  else
  {
    file<- file
    uploaded <- input$otherfile1
    uploaded2 <- read.delim(file = uploaded$datapath, header = TRUE, stringsAsFactors = FALSE) #Warning: Error in read.table: more columns than column names
    file <- file %>%
      rbind(uploaded2)# Warning: Error in read.table: more columns than column names
  }
  rawdata <- file%>%       
      group_by(data,class,type,cna.data,snv.data) %>%
      summarise(n.samples=n_distinct(sample.id),n.patients=n_distinct(patient.id))
    data1 <- rawdata
    data1 <- subset(data1, data %in% input$Resources1)
    data1 <- subset(data1, type %in% input$Types1)
    data1 <- subset(data1, class %in% input$Class1)
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
    ggplot(data1, aes(x=type, y=n.samples, fill=class)) +
      geom_bar(stat="identity") + 
      scale_fill_manual(values=c('#999999','#E69F00'))
    })

  newData2 <- reactive({
    if(is.null(input$otherfile1))
    {
      file <- file
    }
    else
    {
      file<- file
      uploaded <- input$otherfile1
      uploaded2 <- read.delim(file = uploaded$datapath, header = TRUE, stringsAsFactors = FALSE) #Warning: Error in read.table: more columns than column names
      file <- file %>%
        rbind(uploaded2)# Warning: Error in read.table: more columns than column names
    }
    data2 <- file
    data2 <- subset(data2, type %in% input$Types1)
    data2 <- subset(data2, data %in% input$Resources1)
    data2 <- subset(data2, class %in% input$Class1)
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
    data2 <- data2 %>%
        group_by(class,type) %>%
        summarise(n.samples= n_distinct(sample.id),n.patients=n_distinct(patient.id)) 
   
     ggplot(data2, aes(x=type,y=n.samples,fill=class)) +
      geom_bar(stat="identity") + theme(aspect.ratio = 1,legend.position = "none") +
      scale_fill_manual(values=c('#999999','#E69F00')) +
      geom_text(aes(label=n.samples)) + # questo comando serve per aggiungere la numerazione delle quantita', in qunato non si riesce a capire bene il numero reale senza indicazione
      facet_wrap(~class)# serve per fare la divisione per classi --> in questo modo passo da un istogramma a 2 gfafici ad istogramma divisi tra metastatici e primari 
    })
  
  # devo rendere reattivo all in modo che anche quest cambi con il cambiare della soruces, quindi cambino i valori man mano
  # all <- file %>% # in questo caso %>% indica al file sif di applicare l aseguente funzione che segue la 'pipe'
  #  group_by(class,type) %>% # in questo caso il comando group_by singifica che il file vien raggruppato per classe e tipo, poi una volta applicato usa la seconda funzione che segue 
  #  summarise(n.samples=n_distinct(sample.id),n.patients=n_distinct(patient.id)) %>% # summarise sommato creand onumer ocolonne pazienti e sample facendo la conta per i subtypes
  #  add_column(data='all',.before = 'class')
  
  newData_table <- reactive({
    if(is.null(input$otherfile1))
    {
      file <- file
    }
    else
    {
      file<- file
      uploaded <- input$otherfile1
      uploaded2 <- read.delim(file = uploaded$datapath, header = TRUE, stringsAsFactors = FALSE) #Warning: Error in read.table: more columns than column names
      file <- file %>%
        rbind(uploaded2)# Warning: Error in read.table: more columns than column names
    }
    rawdata <- file%>%       
      group_by(data,class,type,cna.data,snv.data) %>%
      summarise(n.samples=n_distinct(sample.id),n.patients=n_distinct(patient.id))
    data3 <- rawdata 
    data3 <- subset(data3, data %in% input$Resources1) 
    data3 <- subset(data3, type %in% input$Types1)
    data3 <- subset(data3, class %in% input$Class1)
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
    all2 <- data3 %>%
          group_by(class,type) %>%
          summarise(n.samples_sum=sum(n.samples),n.patients_sum=sum(n.patients))%>% # sum somma, che e' diverso da n_dsitincr il quale e' un equivalnte della funzionelenght
          add_column(data='all',.before = 'class')
    all2 <- all2 %>%
          rename(n.samples =n.samples_sum, n.patients = n.patients_sum)
    data3 <- data3 %>%
           rbind(all2)
    })
  
#############################################################################################################
##################################### Per secondo pannello ###################################################  
# quindi in questo caso riarriangiamo con arrange dati del file in modo decrescente con desc in base a w.mean (colonna)
  mw <- file2 %>%
    arrange(desc(w.mean))
  
  #data per barplot
  data_second_pannel <- reactive({
    data_pan_1 <- mw
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
    mw3<- mw %>%
      slice_head(n = input$Gene_filter) %>%
      group_split()
    mat1 <- do.call(rbind,mw3)  
    mat1 <- mat1[order(-mat1$w.mean),]  # ordinato per w.mean decrescente 
    mat1 <- subset(mat1, type %in% input$Types2)
    mat1 <- subset(mat1, class %in% input$Class2)
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
  # observeEvent(input$)
  
#############################################################################################################
##################################### Per terzo pannello ###################################################  
  colnames(ensembl) <- c('ensg','start','end','chr','band','Hugo_Symbol')
  ensembl$band <- paste0(ensembl$chr,ensembl$band)
  
  check.ensembl <- ensembl %>%
    group_by(band) %>% 
    summarise(n=n()) %>% 
    arrange(n)
  
plotting <-  reactive({
      selected_class <- input$Class3
      selected_type <- input$Types3
      selected_data <- input$Resources3
      # datalist <- dd[which(names(dd) %in% selected_data)]
      # 
      # FilterSCNA <- function(data_source_id, datalist, file, selected_type, selected_class){
      #   
      #   sel <- file %>%  
      #     filter(data == data_source_id, class %in% selected_class, type %in% selected_type) %>% 
      #     pull(sample.id)
      #   
      #   if(length(sel) == 0){
      #     return( NA )
      #   } else{
      #     m <- datalist[[data_source_id]][,c('Hugo_Symbol',sel)]
      #     return(m)
      #   }
      # }
      # 
      # filtered_datalist <- lapply(X = names(datalist), FUN = FilterSCNA, datalist, file, selected_type, selected_class)
      # names(filtered_datalist) <- names(datalist)
      # 
      # ds <- lapply(names(filtered_datalist), FUN = function(x) setdiff(colnames(filtered_datalist[[x]]),'Hugo_Symbol'))
      # names(ds) <- names(filtered_datalist)
      # 
      # # get freq by gene, by class and subtype
      # 
      # qq <- Reduce(function(...) full_join(...,by='Hugo_Symbol'), filtered_datalist)
      # 
      # getfreq <- function(id,ds,qq,agg=FALSE,flag=NA){
      #   
      #   ss <- which(colnames(qq) %in% as.character(unlist(ds[id])))
      #   
      #   if(!is.na(flag)){
      #     id <- flag
      #   }
      #   
      #   df <- qq %>% select(c(1,all_of(ss))) 
      #   nas <- which(rowSums(is.na(df %>% select(-Hugo_Symbol))) == ncol(df %>% select(-Hugo_Symbol)))
      #   if(length(nas) > 0){
      #     df %>% slice(-nas)
      #   } else{
      #     return(NA)
      #   }
      #   
      #   cna <- c(-2,-1,0,1,2)
      #   names(cna) <- c('homodel','hemidel','neutral','gain','ampl')
      #   
      #   if(agg){
      #     df[which(df == 1,arr.ind = T)] <- 2
      #     df[which(df == -1,arr.ind = T)] <- -2
      #   }
      #   
      #   compfreq <- function(n,cna,df,id){
      #     w <- data.frame(Hugo_Symbol=df$Hugo_Symbol,
      #                     scna=n,
      #                     freq=rowSums(df %>% select(-Hugo_Symbol) == as.numeric(cna[n]),na.rm = TRUE) / rowSums(!is.na(df %>% select(-Hugo_Symbol))),
      #                     data=paste(id,collapse = ';'),
      #                     stringsAsFactors = FALSE)
      #     return(w)
      #   }
      #   
      #   return(do.call(rbind,lapply(names(cna),compfreq,cna,df,id)))
      #   
      # }
      # 
      # out.not <-  do.call(rbind,lapply(names(ds),getfreq,ds,qq,agg=FALSE)) %>% 
      #   rbind(getfreq(id = setdiff(names(ds),"breast_msk_2018"),flag='all_brca',ds,qq,agg=FALSE)) %>% 
      #   add_column(agg=FALSE)
      # 
      # out.agg <- do.call(rbind,lapply(names(ds),getfreq,ds,qq,agg=TRUE)) %>% 
      #   rbind(getfreq(id = setdiff(names(ds),"breast_msk_2018"),flag='all_brca',ds,qq,agg=TRUE)) %>% 
      #   add_column(agg=TRUE)
      # 
      # out <- rbind(out.not, out.agg)
      # 
      # out
      # by gene
      # frq <- left_join(x = out, y=ensembl, by = 'Hugo_Symbol') %>%
      #   filter(!is.na(ensg)) %>%
      #   filter(chr != 'Y') %>%
      #   filter(scna != 'neutral')

      
      
      # d <- frq %>%
      #   group_by(data,scna,chr,band,agg) %>%
      #   summarise(n=n(),
      #             median.freq=median(freq,na.rm = T),
      #             max=max(freq,na.rm = TRUE),
      #             max.name=Hugo_Symbol[which.max(freq)]) %>%
      #   mutate(is.goi = max.name %in% goi)
      # 
      # d$arm <- rep('p',nrow(d))
      # d$arm[grep(d$band,pattern = 'q')] <- 'q'
      
      # plots 
      
      # partendo caricando frq invece che direttamente d da problemi di calolco, perché? Errore : Warning in max(freq, na.rm = TRUE) :
      # no non-missing arguments to max; returning -Inf
      br <- d %>% 
        filter(agg == TRUE) %>% 
        filter(median.freq > input$filter_median_freq) %>%
        filter(data != 'breast_msk_2018') %>% 
        add_column(max.name.goi = NA)
      
      if(input$Chromosomes != 'All'){
        br <- br %>% 
          filter(chr == input$Chromosomes)
      }else{
        br <- br
      }
     
       # non posso collegare qua le resources, devo farlo da qualche altra parte perchè plot usa all brca, quindi avrei meno campione
      br$max.name.goi[which(br$is.goi)] <- br$max.name[which(br$is.goi)]
      
      # per filtro cromosomi


      ggplot(br %>% filter(data == 'all_brca',scna == input$Groups), aes(x=band,y=median.freq,fill=arm)) +
        ylab(paste(input$Groups ,paste('median.freq by cytoband',collapse = ' '))) +   # come cambaire didascali con aggiornatmento
        geom_bar(stat = 'identity') +
        facet_wrap(~chr,scales = 'free_x') +
        scale_fill_manual('arm',values = wes_palette("Chevalier1",n = 2)) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size = 4)) +
        geom_point(data = br %>% filter(scna == input$Groups, data != 'all_brca'),mapping = aes(x=band,y=median.freq,color=data),size=0.5) +
        scale_color_manual('data',values = wes_palette("GrandBudapest1",n = 4)) +
        ggtitle(paste('class:',paste(selected_class,collapse = ','),'\ntype: ',paste(selected_type,collapse = ','))) +
        geom_point(data = br %>% filter(scna == input$Groups, data == 'all_brca'),mapping = aes(x=band,y=max),shape=4,size=0.5) +
        geom_text(data = br %>% filter(scna == input$Groups, data == 'all_brca'),mapping = aes(x=band,y=max,label=max.name.goi),size=1,angle=90,hjust=0,nudge_y=0.01)
    
   
      })
 

#con reactive e observer problema che selectinput non rimande 'fermo' quando seleziono

plotting2 <- reactive({
  br <- d %>%
    filter(agg == TRUE) %>%
    filter(median.freq > input$filter_median_freq) %>%
    filter(data != 'breast_msk_2018') %>%
    add_column(max.name.goi = NA)
  
  if(input$Chromosomes != 'All'){
    br <- br %>%
      filter(chr == input$Chromosomes)
  }else{
    br <- br
  }

  # non posso collegare qua le resources, devo farlo da qualche altra parte perchè plot usa all brca, quindi avrei meno campione
  br$max.name.goi[which(br$is.goi)] <- br$max.name[which(br$is.goi)]
  
  bsel <- br %>% 
          filter(agg == TRUE, scna == input$Groups) %>% 
          pull(band) %>% 
          unique()
})
  
observe({
  updateSelectInput(session, 'Cytoband', choices = plotting2() )
})

                        
plotting3 <-reactive({
  selected_class <- input$Class3
  selected_type <- input$Types3
  selected_data <- input$Resources3
  
  gfrq <- frq %>%
    filter(agg == TRUE, data != "breast_msk_2018") %>%
    filter(band == input$Cytoband) %>%
    mutate(is.goi = Hugo_Symbol %in% goi) %>%
    arrange(chr,start,end)

  ggplot(gfrq %>%                         
           filter(data == 'all_brca',scna == input$Groups) %>%
           arrange(start,end) %>%
           distinct(Hugo_Symbol, .keep_all = TRUE) %>%
           mutate(Hugo_Symbol=factor(Hugo_Symbol, levels = Hugo_Symbol)),
         aes(x=Hugo_Symbol,y=freq,fill=is.goi)) +
    ylab(paste(input$Groups ,paste('freq', collapse=' '))) +
    geom_bar(stat = 'identity') +
    facet_wrap(~band,scales = 'free_x') +
    scale_x_discrete(guide = guide_axis(n.dodge=2)) +
    scale_fill_manual('is.goi',values = wes_palette("Royal1",n = 2)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    ggtitle(paste('class:',paste(selected_class,collapse = ','),'\ntype: ',paste(selected_type,collapse = ','))) +
    geom_point(data = gfrq %>% filter(scna == input$Groups, data != 'all_brca'),mapping = aes(x=Hugo_Symbol,y=freq,color=data)) +
    scale_color_manual('data',values = wes_palette("GrandBudapest1",n = 4))

})


chromosome_table <- reactive({
  br <- d %>%
    filter(agg == TRUE) %>%
    filter(median.freq > input$filter_median_freq) %>%
    filter(data != 'breast_msk_2018') %>%
    add_column(max.name.goi = NA)
  
  if(input$Chromosomes != 'All'){
    br <- br %>%
      filter(chr == input$Chromosomes)
  }else{
    br <- br
  }
  br$max.name.goi[which(br$is.goi)] <- br$max.name[which(br$is.goi)]
  
  br <- br %>% 
    filter(data == 'all_brca', scna == input$Groups)
})
  
  
cytoband_table <- reactive({
    br <- d %>%
      filter(agg == TRUE) %>%
      filter(median.freq > input$filter_median_freq) %>%
      filter(data != 'breast_msk_2018') %>%
      add_column(max.name.goi = NA)
    
    if(input$Chromosomes != 'All'){
      br <- br %>%
        filter(chr == input$Chromosomes)
    }else{
      br <- br
    }
    
    # non posso collegare qua le resources, devo farlo da qualche altra parte perchè plot usa all brca, quindi avrei meno campione
    br$max.name.goi[which(br$is.goi)] <- br$max.name[which(br$is.goi)]
    
    gfrq <- frq %>%
      filter(agg == TRUE, data != "breast_msk_2018") %>%
      filter(band == input$Cytoband) %>%
      mutate(is.goi = Hugo_Symbol %in% goi) %>%
      arrange(chr,start,end) %>% 
      filter(data == 'all_brca',scna == input$Groups) %>%
      arrange(start,end) %>%
      distinct(Hugo_Symbol, .keep_all = TRUE) %>%
      mutate(Hugo_Symbol=factor(Hugo_Symbol, levels = Hugo_Symbol))
  })
#############################################################################################################
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
      plotting3()
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
      ggsave(file,plotting3())
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


