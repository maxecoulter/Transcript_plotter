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
  
* -gtf
  
      chr1H PBRI	exon	72060	72399	.	+	.	transcript_id "BaRT2v18chr1HG000020.1"; gene_id "BaRT2v18chr1HG000020";
      chr1H	PBRI	five_prime_utr	72060	72209	.	+	.	transcript_id "BaRT2v18chr1HG000020.1"; gene_id "BaRT2v18chr1HG000020";
      chr1H	PBRI	CDS	72210	72399	.	+	.	transcript_id "BaRT2v18chr1HG000020.1"; gene_id "BaRT2v18chr1HG000020";
      chr1H	PBRI	start_codon	72210	72212	.	+	.	transcript_id "BaRT2v18chr1HG000020.1"; gene_id "BaRT2v18chr1HG000020";
      chr1H	PBRI	exon	72524	73355	.	+	.	transcript_id "BaRT2v18chr1HG000020.1"; gene_id "BaRT2v18chr1HG000020";
      chr1H	PBRI	CDS	72524	72570	.	+	.	transcript_id "BaRT2v18chr1HG000020.1"; gene_id "BaRT2v18chr1HG000020";
      chr1H	PBRI	stop_codon	72568	72570	.	+	.	transcript_id "BaRT2v18chr1HG000020.1"; gene_id "BaRT2v18chr1HG000020";

