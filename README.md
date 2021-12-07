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
  
  -snpeff **SnpEff vcf file**
  
  -o $output, **outfile (must specify as .R). Default is "gene_models.R". All other files will be written to the directory the R script is in.**
  
  
### Detailed information on inputs:
  
* **gtf example**
  
      chr1H PBRI	exon	72060	72399	.	+	.	transcript_id "BaRT2v18chr1HG000020.1"; gene_id "BaRT2v18chr1HG000020";
      chr1H	PBRI	five_prime_utr	72060	72209	.	+	.	transcript_id "BaRT2v18chr1HG000020.1"; gene_id "BaRT2v18chr1HG000020";
      chr1H	PBRI	CDS	72210	72399	.	+	.	transcript_id "BaRT2v18chr1HG000020.1"; gene_id "BaRT2v18chr1HG000020";
      chr1H	PBRI	start_codon	72210	72212	.	+	.	transcript_id "BaRT2v18chr1HG000020.1"; gene_id "BaRT2v18chr1HG000020";
      chr1H	PBRI	exon	72524	73355	.	+	.	transcript_id "BaRT2v18chr1HG000020.1"; gene_id "BaRT2v18chr1HG000020";
      chr1H	PBRI	CDS	72524	72570	.	+	.	transcript_id "BaRT2v18chr1HG000020.1"; gene_id "BaRT2v18chr1HG000020";
      chr1H	PBRI	stop_codon	72568	72570	.	+	.	transcript_id "BaRT2v18chr1HG000020.1"; gene_id "BaRT2v18chr1HG000020";

* Note transcript_id and then gene_id
  
* **Bed example**
  
      chr1H	366145087	366145296	BaRT2v18chr1HG024650;BaRT2v18chr1HG024650.1	.	+	366145098	366145296	255,0,0	1	209	0
      chr1H	491055869	491064771	BaRT2v18chr1HG043310;BaRT2v18chr1HG043310.1	.	+	491055915	491064496	255,0,0	2	969,321	0,8581
      chr1H	436979920	436980794	BaRT2v18chr1HG032290;BaRT2v18chr1HG032290.1	.	+	436980051	436980567	255,0,0	1	874	0
      chr1H	101565607	101572273	BaRT2v18chr1HG011810;BaRT2v18chr1HG011810.3	.	+	101566769	101572061	255,0,0	7	1344,85,65,100,42,42,455	0,4355,5396,5571,5908,6080,6211
      chr1H	101565607	101572273	BaRT2v18chr1HG011810;BaRT2v18chr1HG011810.2	.	+	101566769	101572061	255,0,0	8	49,206,85,65,100,42,42,455	0,1138,4355,5396,5571,5908,6080,6211
