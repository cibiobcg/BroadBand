# Web based-resource in Shiny for breast cancer genomic data exploration 

This app leads to exploring data and producing a list of gene and/or genomic alterations frequently altered in breast cancer, both genes and genomic regions.
The default data are from cBioPortal, in particular: 
- Somatic single nucleotide variants (non-synonymous)
- Somatic copy number alterations

The first panel of the tool is a count panel for an overview of the samples and their frequency. The second panel is only about SNV, the plots inside this panel want 
to show the genes most frequently altered and if there are specific genes for class. The third and last panel is about CNA, the plots show the aberration frequency of chromosomal cytobands and also the aberration frequency of the genes inside the cytobands based on the choice of copy number granularity, the user can choose between 4 and 2 classes.
4 classes which are: 
- Homozygous deletion
- Hemizygous deletion
- Gain
- Amplification

2 classes which are  :
- Deletion
- Amplification 

All panels contain widgets about the resources of the data, the tumor classification, and the tumor subtypes but some panels can be present with other widgets but really intuitive.
In the panels by pressing the button 'Apply' the user is able to filter the data and see the results of the selected filters on plots and tables.

Every panel has the possibility to upload a single file that will bind directly to the file already present in the application. This widget shows the resources will be automatically updated with all the new resources of the new file generated, for seeing the updated plots and tables with the new file the user need to press even in this case the button 'Apply' and then the new file can undergo under the same analysis procedures. Any plots and tables in the application can be downloaded thanks to the buttons in the application.

For complete use of this app, there is the need to download the file 'app.R' and folders inside 'brcApp', which contain the files, the application code, and the temporary logo for this app. The user needs to set the working directory to the source file location.  
 
More information on the [GitHub] (https://cibiobcg.github.io/BroadBand/) 
