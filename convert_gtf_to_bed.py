

import argparse
import sys
#chr1H_part_1    scallop transcript      5085369 5088493 1000    +       .       gene_id "chr1H_part_1G000800"; transcript_id "chr1H_part_1G000800.RTD.4"; RPKM "0.6720"; cov "18.2501";
#chr1H_part_1    scallop exon    5085369 5085836 1000    +       .       gene_id "chr1H_part_1G000800"; transcript_id "chr1H_part_1G000800.RTD.4"; exon "1";

#bed12 format: https://genome.ucsc.edu/FAQ/FAQformat.html

#First import gtf file. what is needed is a dictionary of transcripts with features. Suggested data structure:
#{transcript_id:[transcript_object,[exon1_object,exon2_object]]
#Next create bed12 file
#Go through each transcript in dictionary (as each line in bed file)
#Get features from exons: block count (number of exons), block sizes , block starts


class Gtfinput:
	def __init__(self,line):
		line1 = line.rstrip("\n").split("\t")
		self.chromosome = line1[0]
		self.assembler = line1[1]
		self.type = line1[2]
		self.start = line1[3]
		self.end = line1[4]
		self.score = line1[5]
		self.strand = line1[6]
		self.information = line1[8]
		if self.type == "transcript":
			transcript_id, gene_id = line1[8].split(";")[0:2]
			self.gene_id = gene_id.split(" ")[2].replace('"','')
			self.transcript_id = transcript_id.split(" ")[1].replace('"','')
		elif self.type == "exon":
			transcript_id, gene_id, exon_number = line1[8].split(";")[0:3]
			self.gene_id = gene_id.split(" ")[2].replace('"','')
			self.transcript_id = transcript_id.split(" ")[1].replace('"','')
			#self.exon_number = exon_number.split(" ")[2].replace('"','')
	def process(self,gtf_dict):
		if self.type == "transcript":
			gtf_dict[self.transcript_id] = [self,[]]
		elif self.type == "exon":
			gtf_dict[self.transcript_id][1].append(self)
		else:
			print("Line error! Should only be transcript or exon types")
		return gtf_dict

def bedmaker(gtf_dict,output):
	outputfile = open(output,"w")
	print(f"Number of transcripts found: {len(gtf_dict.keys())}")
	for transcript in gtf_dict.keys():
		transcript_object, exons = gtf_dict[transcript]
		bedstart = int(transcript_object.start) - 1#500, 499
		exon_sizes = []
		exon_starts = []
		for exon in exons:
			exon_sizes.append(str(int(exon.end) + 1 - int(exon.start))) # + 1
			exon_starts.append(str(int(exon.start) - int(transcript_object.start)))
		#print(transcript)
		outputfile.write(transcript_object.chromosome + "\t" + str(bedstart) + "\t" + transcript_object.end + "\t" + transcript_object.gene_id + ";" + transcript_object.transcript_id + "\t" + transcript_object.score + "\t" + transcript_object.strand + "\t" + str(bedstart) + "\t" + transcript_object.end + "\t" + "255,0,0" + "\t" + str(len(exons)) + "\t" + ",".join(exon_sizes) + "\t" + ",".join(exon_starts) + "\n")
	outputfile.close()

def main():
	parser = argparse.ArgumentParser(description='Convert gtf to .bed')
	parser.add_argument('-i', dest = 'gtf', type = str, help = 'The gtf input file')
	parser.add_argument('-o', dest = 'output',type = str, help = 'Name of output file')
	args = parser.parse_args()
	gtf_dict = {}
	for line in open(args.gtf):
		gtfline = Gtfinput(line)
		try:
			gtf_dict = gtfline.process(gtf_dict)
		except:
			print(f"Error in line: {line}")
			sys.exit()
	bedmaker(gtf_dict,args.output)

if __name__ == "__main__":
	main()


