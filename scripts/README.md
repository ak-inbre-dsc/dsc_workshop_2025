# Scripts  

This folder contains R and shell scripts we used during the workshop. 

- [**install_bioinfo_tools.sh**](install_bioinfo_tools.sh) - Shell script to install bioinformatics tools on a Google Cloud Vertex AI instance.  

- [**as_processing_cloud_v2.sh**](as_processing_cloud_v2.sh) - Shell script for processing and analyzing nanopore adaptive sampling data. It splits reads by treatment, maps them to a reference, and generates summaries.  

- [**get_all_stats.sh**](get_all_stats.sh) - Shell script to collect SeqKit stats for various FASTQ outputs from the adaptive sampling analysis pipeline, including: the number of sequences, total yield, minimum, average, and maximum sequence lengths, N50, and the average quality score of all bases.  

- [**violinplot_treatments_isolates.R**](violinplot_treatments_isolates.R) - R script that generates a series of violin plots to visualize sequencing data. It compares the distributions of read length and quality between two different treatments, both for the overall dataset and for reads mapped to specific isolates.  

- [**plot_treatments_isolates.R**](plot_treatments_isolates.R) - R script that creates bar plots to compare overall sequencing metrics like read count and quality between the sequences derived from the control and adaptive sampling portion of the flow cells. It also visualizes how the proportion of reads mapped to different isolates varies between these two treatments.  


