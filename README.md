# Transcript_plotter

Transcript plotter is a tool for visualising predicted protein domains and variants in the context of genes and transcripts.

The tool is based around a python script that generates inputs that can be used in a shiny app. 
The shiny app is found here: https://maxecoulter.shinyapps.io/Transcript_plotter/.

## transcript_plotter.py

* The python script is written in python 3.8. It does not require any packages to be installed

* The script takes as a main transcriptome input either a .gtf file, or a .bed file 

  -gtf **/path/to/GTF file** 
  
  -bed **/path/to/Bed file**

  -c **Coordinates**, The region of interest. This is optional, but we suggest you use this otherwise you will have a huge list of genes to visualise. Format is <chromosome>,<start>,<end>
  
  -interpro **Interpro scan tsv file**
  
  -p **rtPCRInput**, the RT-PCR proportions file
  
  -s **PathToSamplesFolder**, the salmon input folder
  
  -w **Match window size**, this is the wondow size that the RT-PCR/transcript product matching algorithm uses. For BaRTv2.0 it was set to 6.
  
  -o $output, The output path and prefix
