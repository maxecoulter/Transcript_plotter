# Transcript plotter

Transcript plotter is a tool for visualising predicted protein domains and variants in the context of genes and transcripts. You can produce publication quality figures for understanding both how transcript variation is affecting downstream protein function and also how variants are affecting your transcripts. 

### Overview

Genes can produce multiple transcripts with multiple isoforms. Advances in Iso-seq and nanopore technology mean we can now sequence whole transcripts. This is great, but the resulting data can be very complex. It is helpful to visualise RNA-seq data in the context of genes and transcripts. Transcript plotter is a tool for visualising genes, transcripts, predicted protein domains and variant calling. This can help you better understand expression changes at the transcript level (AS, DTU etc.), and so the tool goes hand in hand with transcriptome based RNA-seq quantification using Salmon/Kallisto. 




![Figure 1](https://github.com/maxecoulter/Transcript_plotter/blob/master/Figures/BaRT2v18chr3HG123000_with_domains.png)

All transcripts in a gene can be visualised with predicted protein domains, to assess how transcript diversity is affecting protein diversity.



![Figure 2](https://github.com/maxecoulter/Transcript_plotter/blob/master/Figures/Cellulose%20synthase.png)

Single transcripts can be visualised as above. You can visualise SnpEff data on individual transcripts.





The tool is currently based around a python script that generates inputs that can be used in a shiny app. 
The shiny app is found here: https://maxecoulter.shinyapps.io/Transcript_plotter_v2/.



## transcript_plotter.py

* The python script is written in python 3.8. It does not require any packages to be installed


* The script takes as a main transcriptome input either a .gtf file, or a .bed file 


  -gtf **/path/to/mytranscriptome.gtf** 
  
  -bed **/path/to/mytranscriptome.bed**

  -c **Coordinates**, The region of interest. This is optional, but I suggest you use this otherwise you will have a huge list of genes to visualise. Format is <chromosome>,<start>,<end>
  
  -interpro **InterProScan tsv file**
  
  -snpeff **SnpEff vcf file**
  
  -o $output, **outfile (must specify as .R). Default is "gene_models.R". All other files will be written to the directory the R script is in.**
  
  
### Detailed information on inputs:
  
### **Gtf**
  
      chr1H PBRI	exon	72060	72399	.	+	.	transcript_id "BaRT2v18chr1HG000020.1"; gene_id "BaRT2v18chr1HG000020";
      chr1H	PBRI	five_prime_utr	72060	72209	.	+	.	transcript_id "BaRT2v18chr1HG000020.1"; gene_id "BaRT2v18chr1HG000020";
      chr1H	PBRI	CDS	72210	72399	.	+	.	transcript_id "BaRT2v18chr1HG000020.1"; gene_id "BaRT2v18chr1HG000020";
      chr1H	PBRI	start_codon	72210	72212	.	+	.	transcript_id "BaRT2v18chr1HG000020.1"; gene_id "BaRT2v18chr1HG000020";
      chr1H	PBRI	exon	72524	73355	.	+	.	transcript_id "BaRT2v18chr1HG000020.1"; gene_id "BaRT2v18chr1HG000020";
      chr1H	PBRI	CDS	72524	72570	.	+	.	transcript_id "BaRT2v18chr1HG000020.1"; gene_id "BaRT2v18chr1HG000020";
      chr1H	PBRI	stop_codon	72568	72570	.	+	.	transcript_id "BaRT2v18chr1HG000020.1"; gene_id "BaRT2v18chr1HG000020";

* Note transcript_id and then gene_id
  
### **Bed**
  
      chr1H	366145087	366145296	BaRT2v18chr1HG024650;BaRT2v18chr1HG024650.1	.	+	366145098	366145296	255,0,0	1	209	0
      chr1H	491055869	491064771	BaRT2v18chr1HG043310;BaRT2v18chr1HG043310.1	.	+	491055915	491064496	255,0,0	2	969,321	0,8581
      chr1H	436979920	436980794	BaRT2v18chr1HG032290;BaRT2v18chr1HG032290.1	.	+	436980051	436980567	255,0,0	1	874	0
      chr1H	101565607	101572273	BaRT2v18chr1HG011810;BaRT2v18chr1HG011810.3	.	+	101566769	101572061	255,0,0	7	1344,85,65,100,42,42,455	0,4355,5396,5571,5908,6080,6211
      chr1H	101565607	101572273	BaRT2v18chr1HG011810;BaRT2v18chr1HG011810.2	.	+	101566769	101572061	255,0,0	8	49,206,85,65,100,42,42,455	0,1138,4355,5396,5571,5908,6080,6211
  
* File must be in bed12 format. Note the difference between coordinates in columns 2/3 and 7/8. Columns 7/8 are coordinates for CDS, so it is important that if you want to display CDS information (and downstream Interproscan results) you will need to make sure CDS information is included in your .bed file
 
 * Inputs are both optional, you only need one or the other
  
 ### **InterProScan (Optional)**

        BaRT2v18chr3HG122880.8	e0c507447a8cc27375bf2ecb1307867f	613	Gene3D	G3DSA:1.10.510.10	Transferase(Phosphotransferase) domain 1	365	580	7.60E-59	T	01/11/2021	-	-
        BaRT2v18chr3HG122880.8	e0c507447a8cc27375bf2ecb1307867f	613	Pfam	PF13855	Leucine rich repeat	102	161	1.30E-07	T	01/11/2021	IPR001611	Leucine-rich repeat
        BaRT2v18chr3HG122880.8	e0c507447a8cc27375bf2ecb1307867f	613	Pfam	PF08263	Leucine rich repeat N-terminal domain	35	74	7.60E-11	T	01/11/2021	IPR013210	Leucine-rich repeat-containing N-terminal, plant-type
        BaRT2v18chr3HG122880.8	e0c507447a8cc27375bf2ecb1307867f	613	ProSitePatterns	PS00107	Protein kinases ATP-binding region signature.	294	316	-	T	01/11/2021	IPR017441	Protein kinase, ATP binding site
        BaRT2v18chr3HG122880.8	e0c507447a8cc27375bf2ecb1307867f	613	PANTHER	PTHR47988	SOMATIC EMBRYOGENESIS RECEPTOR KINASE 1	6	613	0	T	01/11/2021	-	-
        BaRT2v18chr3HG122880.8	e0c507447a8cc27375bf2ecb1307867f	613	ProSitePatterns	PS00108	Serine/Threonine protein kinases active-site signature.	411	423	-	T	01/11/2021	IPR008271	Serine/threonine-protein kinase, active site
        BaRT2v18chr3HG122880.8	e0c507447a8cc27375bf2ecb1307867f	613	PANTHER	PTHR47988:SF2	PROTEIN NSP-INTERACTING KINASE 3	6	613	0	T	01/11/2021	-	-
        BaRT2v18chr3HG122880.8	e0c507447a8cc27375bf2ecb1307867f	613	Pfam	PF00069	Protein kinase domain	289	558	1.30E-42	T	01/11/2021	IPR000719	Protein kinase domain
        BaRT2v18chr3HG122880.8	e0c507447a8cc27375bf2ecb1307867f	613	SUPERFAMILY	SSF52058	L domain-like	35	202	5.33E-35	T	01/11/2021	-	-
        BaRT2v18chr3HG122880.8	e0c507447a8cc27375bf2ecb1307867f	613	Gene3D	G3DSA:3.30.200.20	Phosphorylase Kinase; domain 1	259	364	6.50E-34	T	01/11/2021	-	-
        BaRT2v18chr3HG122880.8	e0c507447a8cc27375bf2ecb1307867f	613	ProSiteProfiles	PS50011	Protein kinase domain profile.	288	571	35.255	T	01/11/2021	IPR000719	Protein kinase domain


  * This file is required for visualising predicted protein domains on your transcripts
  
  * To create this file you will need to have a fasta file of your amino acid sequences for each transcript, and you can run using instructions found here: https://interproscan-docs.readthedocs.io/en/latest/HowToRun.html 
  
  * Create the tsv file needed with option **-f tsv**
  
  * There maybe some domains you do not wish to display. In which case, you can easily adjust the output csv file **genemodel_input_domain.csv**, which feeds into the shiny app
  
 ### **SnpEff vcf file (Optional)**
  
 * You will need a .vcf file as produced by snpeff (see http://pcingola.github.io/SnpEff/) 
  
 * Transcript_plotter uses output from version snpeff-5.0-1, it is not tested on later versions
  
 * You will need to create a custom SnpEff database if you do not already have one for your transcriptome, in order for your snp annotations to make sense
  
 * You will need to filter your vcf file based on quality, otherwise all variants will be displayed (transcript_plotter doesn't do that for you)
  
 * If you would like to include (or remove) text labels for your variants (high impact red variants are automatically labelled) you can adjust the **snpeff_out.csv** file accordingly
  
 * The snpeff analysis is currently made with RNA-seq variant calling in mind, so all intergenenic. upstream/ downstream variants are removed. I will include an option to retain these, not implemented currently
 
  
 ### Outputs:
  
 * All outputs will be written to the directory you specified to put your output R file (ie /path/to/Rout.R)
  
 * **The .R outfile**
  
 * This is an ordinary R script that can be used to produce plots. In most cases you can ignore this and use the shiny app. However for advanced R users who might like more control you can use this rather long and cumbersome R script. It may help with debugging as well
  
 * The R script is written by python. It should produce the same outputs as the shiny app
 
  * **Shiny inputs**
  
 * The follwing files are produced by transcript_plotter.py that are required for the shiny app:
  
  * **genemodel_input.csv**
  
  * **shiny_out.csv**
  
  * **genemodel_input_domain.csv**
  
  * **snpeff_out.csv**
  
  
  * Should be clear in shiny app what goes where. If you get an error in shiny, check you have downloaded the correct input files!
  
  **SNPeff_output.tsv**
  
  * This is a tsv file with a gene based overview of snpeff results
  
  * Use this if you would like to see snpeff results in tabulated form 
 
  
  
  ## Shiny app
  
  You can use the shiny app to generate a variety of gene/transcript plots. Simply download the input files. The interProScan and snpeff files are optional. You can adjust legend size, plot size, axis text size etc. You can download the plots to eps format. Have a play! Some known glitches:
  
  * If you change the plot size, it will spread out the domain key. Simply turn domains on and off again, should fix
  
  * If you move on to the next gene the domain key from the previous will still be there. Just switch domains on and off again...
  
  * Not possible to display snpeff data on transcripts without interpro domains, if you have **genemodel_input_domain.csv** loaded. If you want to do this, restart the app and don't load **genemodel_input_domain.csv**
  
  * Everytime you press the domain button, it will change the colours of the domains. This is a quick way of visualising the data. However, if you download the plot it will also change the colours, so they won't be the same. You will need to put the colours in manually in order to make sure they stay the same in your final plot!
  
  * You need to input the correct number of manual colours for the number of domains. This is the number of domains you see in the key excluding CDS and UTR (these always have the same colour). If you input the wrong number the plot will not display

  * The shiny.io server does not work currently with large files. So you can't currently run transcript_plotter.py on a whole transcriptome and upload the csv file, it won't work. This can be done locally though
  
  
  * Please let me know if you find any bugs or have any problems (mecoulter@dundee.ac.uk). Enjoy!
  
  

  
 
  
  
  
  

 
  

  
