# web based-resource in Shiny for breast cancer genomic data exploration 

This app lead to explore data and to produce a list of gene and/or  genomic alterations frequentaly alteretaed in breast cancer, both gnees and genomic regions.
The dafault data are from cbioportal, in particular : 
- somatic single nucleotide varaints (non sinonimus)
-  somatic copy number alterations

Fist pannel of the tool is a count pannel for an overview of the samples and their frequency. Second pannel is only about SNV, the plots inside this pannel want 
to show the genes most freqeuntly alterated and if there are specific genes for class. The third and last pannel is about CNA, the plots shows the aberration freqeuncy of amplification/deletion of chromosomial cytoband and also the abberation freqeuncy of the genes inside the cytobands. 

Every pannel have the possibility to upload a single file that will bind the deafult file and undergoes the same analysis procedur. Any plot and table can be downloaded.
