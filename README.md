# COVID_MOFA
This is the code for the Dissecting the Blood Ecosystem in SARS-CoV-2 Omicron Patients.
## File directory description

Script files are stored in the  **bio/** directory. 

Script file introduction.

Script files are stored in the **bio/** directory. The functions of each script are as follows:

+ COVID_MOFA.R contains all of the code from data import cleaning, model generation to basic analysis, and is the basis for the rest of the scripts.

+ MOFA_limma.R The main function is to analyze the difference between platelet transcriptome and platelet proteome.

+ MOFA_GSEA.R makes it convenient to use the data in the new Ensembl ID for the GSEA analysis.

+ MOFA_dotplot.R Visualizes GSEA results.

+ MOFA_RCircos.R uses a circle graph to visualize the top gene for each factor and their weight.

+ MOFA_factor_metadata_correlation.R Scripts used to analyze and visualize the correlation between factors and metadata.

+ MOFA_factor_Omicron_phases_scatter script is used to analyze and visualize the correlation between factors in different Omicron phases.

+ MOFA_top_gene.R in the correlation analysis and visualization of features of interest in the various Omicron phases.
 
