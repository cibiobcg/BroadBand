# Web based-resource in Shiny for breast cancer genomic data exploration 

This app lead to explore data and to produce a list of gene and/or genomic alterations frequentaly alteretaed in breast cancer, both geens and genomic regions.
The dafault data are from cbioportal, in particular: 
- Somatic single nucleotide varaints (non sinonimus)
- Somatic copy number alterations

Fist pannel of the tool is a count pannel for an overview of the samples and their frequency. Second pannel is only about SNV, the plots inside this pannel want 
to show the genes most freqeuntly alterated and if there are specific genes for class. The third and last pannel is about CNA, the plots shows the aberration freqeuncy of chromosomial cytoband and also the abberation freqeuncy of the genes inside the cytobands based to the choose of copy number granularity, the user can choose between 4 and 2 classes, where 4 correspond to homozygous deletion-hemizygous deletion-gain-amplification and for 2 classes correspond deletion-amplification. All pannels contains widgets about the resources of the data, the tumor cassification and the tumor subtypes but in some pannels can be present others widgets but realy intuitive. In all pannels to be able to see on the plots the wanted filters applied there is the need to press the button 'Apply'. 

Every pannel have the possibility to upload a single file that will bind diretcly the file already present in the application, the widget that shows the resources will be automaticaly be updated with all the new resources of the new file generated, for seeing the updated plots and tabels wit hthe new file the user need to press even in this case the button 'Apply' and then the new file can undergoes under the same analysis procedurs. Any plots and tables presents in the application can be downloaded thanks to the buttons presents in the application.

For a comeplete use of this app there is the need to download the file 'app.R' and folders inside 'brcApp', which contain the files, the application code and the temporary logo for this app. Also the user need to change the path of the files inside the app code.
 
