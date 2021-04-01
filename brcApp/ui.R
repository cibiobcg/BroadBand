#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

# Define UI for application that draws a histogram
ui <- fluidPage(
    
    # Application title
    titlePanel("App Breast Cance"),
    
    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            fileInput(
                inputId = 'otherfile',
                label = 'choose file'),
            checkboxGroupInput(
                inputId = 'Resources',
                label = 'Data resources',
                choices = c('brca_metabric','brca_igr_2015','brca_mbcproject_wagle_2017','breast_msk_2018','brca_tcga_pan_can_atlas_2018'),
                selected = c('brca_metabric','brca_igr_2015','brca_mbcproject_wagle_2017','breast_msk_2018','brca_tcga_pan_can_atlas_2018')
            ),
            checkboxGroupInput(
                inputId = 'Types',
                label = 'Breast cancer subtypes',
                choices = c('HER2+','HR+','TNBC'),
                selected = c('HER2+','HR+','TNBC')
            ),
            checkboxGroupInput(
                inputId = 'Class',
                label = 'Tumor classification',
                choices = c('Primary','Metastasis'),
                selected = c('Primary','Metastasis')
                
            ),
        ),
        
        # Show a plot of the generated distribution
        mainPanel(
            downloadButton(
                outputId = 'downloadData',
                label = 'Download Data'),
            plotOutput(
                outputId = 'myplot',
                width = '100%',
                height = '400px'),
            
        )
    )
)
