# Web based-resource in Shiny for breast cancer genomic data exploration 

This app lead to explore data and to produce a list of gene and/or genomic alterations frequentaly alteretaed in breast cancer, both gnees and genomic regions.
The dafault data are from cbioportal, in particular: 
- Somatic single nucleotide varaints (non sinonimus)
- Somatic copy number alterations

Fist pannel of the tool is a count pannel for an overview of the samples and their frequency. Second pannel is only about SNV, the plots inside this pannel want 
to show the genes most freqeuntly alterated and if there are specific genes for class. The third and last pannel is about CNA, the plots shows the aberration freqeuncy of amplification/deletion of chromosomial cytoband and also the abberation freqeuncy of the genes inside the cytobands. In all pannels for apply the wanted filter there is the need to press the button 'Apply'. 

Every pannel have the possibility to upload a single file that will bind tto the file alrady present in the application and thanks to that the widget that shows the resources will be automaticaly be updated with all the new resources of the new file generated, for seeing the updated plots and tabels the user need to press even in this case the button 'Apply' and then at the new file generated can undergoes the same analysis procedur. Any plots and tables presents in the application can be downloaded thanks to the special buttons presents in the application.

For the 
 
