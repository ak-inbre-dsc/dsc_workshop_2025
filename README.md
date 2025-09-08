![header](assets/images/git_header_v2.png)

# Materials for the 2025 Alaska INBRE Data Science Core Genomics Workshop 

Welcome! This repository contains the materials that we used during our 2025 workshop. The easiest way to access all the materials is to [download the whole repository](https://github.com/KGrond/test_workshop). A full description of the materials can be found [here](https://docs.google.com/document/d/14ThjTwR0dHlLXQQcQyTGRRvncRI6Mnz8ogHLNMDPjTs/edit?usp=sharing)


## ![](assets/images/ppt_icon.png) [Lectures](lectures)  

1. [Introduction to Nanopore Sequencing](lectures/DSC-Workshop2025_Mod-1-2_IntroNanopore.pdf)
2. [Cloud Computing and HPCs](lectures/DSC-Workshop2025_Mod-2_Cloud-Computing-and-HPCs.pdf)
3. [Pipelines, Containers, and Workflow Management Systems](lectures/DSC-Workshop2025_Mod-3_pipelines-workflows-containers.pdf)
4. [Principles of Reproducible Research](lectures/DSC-Workshop2025_Mod-3_Reproducible-Research-Principles.pdf)
5. [Version Control with Git](lectures/DSC-Workshop2025_Mod-3_VersionControl.pdf)

##  :notebook: [Jupyter Notebooks](notebooks/)

### 1. Nanopore Analysis
* [Epi2me setup](notebooks/Nanopore_Analysis/DSC_Module_2_Epi2meSetup.ipynb)  
* [Epi2me bacterial-genomes workflow](notebooks/Nanopore_Analysis/DSC_Module_2_wf-bacterial-genomes.ipynb)  
* [Epi2me metagenomics workflow](notebooks/Nanopore_Analysis/DSC_Module_2_wf-metagenomics.ipynb)  
* [Walkthrough of the get_all_stats.sh script](notebooks/Nanopore_Analysis/DSC_Module_4_get_stats.ipynb) 
* [Walkthrough of the install_bioinfo_tools.sh script](notebooks/Nanopore_Analysis/DSC_Module_4_install_bioinfo_tools.ipynb)
* [Adaptive sampling analysis pipeline](notebooks/Nanopore_Analysis/DSC_Module_4_process_AS.ipynb)


### 2. Version Control with Git
* Lesson 1: [Configure Git](notebooks/Git/Lesson1_Setting_up_git.ipynb)  
* Lesson 2: [Getting started](notebooks/Git/Lesson2_Create_Project_and_View_History_Ignoring_Files.ipynb)  
* Lesson 3: [Connecting to Github](notebooks/Git/Lesson3_Connecting_to_Github.ipynb)  
* Lesson 4: [Branching and merging](notebooks/Git/Lesson4_Branching_merging.ipynb)  
* Lesson 5: [Collaborating with Git](notebooks/Git/Lesson5_Collaborating.ipynb)  
* Lesson 6: [Solve a merge conflicts](notebooks/Git/Lesson6_Managing_Conflict.ipynb)  
* Lesson 7: [RStudio & Git](notebooks/Git/Lesson7_RStudioServer_Git.ipynb)  


## ![](assets/images/code_icon.png) [Scripts](scripts/)

* Visualize sequencing data in R through [violin plots](scripts/violinplot_treatments_isolates.R) 
* Create [bar plots](scripts/plot_treatments_isolates.R) in R to compare different sequencing metrics.
* [Install bioinformatics tools](scripts/install_bioinfo_tools.sh) on a Google Cloud Vertex AI instance.
* [Collect SeqKit stats](scripts/get_all_stats.sh) for various FASTQ outputs from the adaptive sampling analysis pipeline.
* [Process and analyze nanopore adaptive sampling data](scripts/as_processing_cloud_v2.sh)


## :bookmark_tabs: [Nanopore Handouts](handouts)  
Descriptions of the mock communities and different library prepapration methods 

* [Mock Communities](handouts/Mock_Communities)
* [Sequencing protocols & Information](handouts/Protocols)
\

*Alaska INBRE is an Institutional Development Award (IDeA) from the National Institute of General Medical Sciences of the National Institute of Health (NIH) under grant number P20GM103395.*



## License
This project is licensed under the [Creative Commons Attribution (CC BY)](LICENSE.md).
