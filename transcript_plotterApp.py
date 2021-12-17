

import random
import sys
import pdb
import argparse



class Bed:
	b_Instances = {}
	gene_transcripts = {}
	genestrand = {}
	genestarts = {}
	geneends = {}
	genechromosome = {}
	def __init__(self, line, locus=[]):
		self.line = line
		self.fields = line.strip("\n").split("\t")
		self.chromosome = self.fields[0]
		self.start = int(self.fields[1])
		self.end = int(self.fields[2])
		#
		if locus:#assumes start and end != None
			chromosome, start, end = locus
			if self.chromosome != chromosome:
				return
			if self.start < int(start) or self.end > int(end):
				return
		#
		transcriptinfo=self.fields[3].split(";")#transcript id is transcriptinfo[1]
		self.gene = transcriptinfo[0]
		self.transcript = transcriptinfo[1]
		self.strand = self.fields[5]
		self.CDS_start = int(self.fields[6])
		self.CDS_end = int(self.fields[7])
		self.exon_no = int(self.fields[9])
		self.list_exon_lengths = self.fields[10].split(",")
		self.list_exon_starts = self.fields[11].split(",")
		self.b_Instances[self.transcript] = self
		self.gene_transcripts.setdefault(self.gene, []).append(self.transcript)
		self.genestrand[self.gene] = self.strand
		self.genechromosome[self.gene] = self.chromosome
		try:
			gene_start = self.genestarts[self.gene]
			if self.start < gene_start:
				self.genestarts[self.gene] = self.start 
		except KeyError:
			self.genestarts[self.gene] = self.start
		#
		try:
			gene_end = self.geneends[self.gene]
			if self.end > gene_end:
				self.geneends[self.gene] = self.end
		except KeyError:
			self.geneends[self.gene] = self.end
	#	 
	#
	def __getitem__(self, transcript):
		return self.b_Instances[transcript]
	#
	def __repr__(self):
		return self.line
	#
	def summarise(self):
		genes = set()
		transcripts = set()
		for bed in self.b_Instances.values():
			genes.add(bed.gene)
			transcripts.add(bed.transcript)
		print(f"Bed file has {len(genes)} genes and {len(transcripts)} transcripts in selected region")
	@staticmethod
	def get_f_coordinates(start, end, exon_coordinates):
		"""Get feature coordinates"""
		coordinates = []
		for exon in exon_coordinates:
			e_start, e_end = exon.split("-")
			if int(e_end) < start:#if exon completely before cds
				continue
			if start >= int(e_start) and int(e_end) < end:#if start within exon and end not
				coordinates.append(f"{start}-{e_end}")
			elif end > int(e_start) and end <= int(e_end) and start < int(e_start):#if end within the exon
				coordinates.append(f"{e_start}-{end}")
				break
			elif start >= int(e_start) and end <= int(e_end): #if cds start and end in one exon
				coordinates.append(f"{start}-{end}")
				break
			elif start < int(e_start) and end > int(e_end): #If exon within cds completely
				coordinates.append(f"{e_start}-{e_end}")
			else:
				pass
		return coordinates
	#
	@staticmethod
	def convert_coordinates(coordinates, transcript_start, transcript_end):
		"""Convert coordinates to relative to strand start if on '-' strand"""
		conversion = transcript_end - transcript_start
		new_coordinates = []
		for coordinate in coordinates:
			end, start = [conversion - int(co) + 1 for co in coordinate.split("-")]
			new_coordinates.append(f"{start}-{end}")
		#
		new_coordinates.sort()
		#
		return new_coordinates
	#
	def write_genemodel(self, gene, transcript, output):
		#
		with open(output, "a") as out:
			features_dict = {}
			bed = self.b_Instances[transcript]
			if bed.strand == "+":
				orientation = "forward"
			elif bed.strand == "-":
				orientation = "reverse"
			else:
				print("fatal strand error")
				sys.exit()
			gene_start = self.genestarts[gene]
			gene_end = self.geneends[gene]
			
			
			#
			#First generate exon coordinates (relative start)
			features_dict["exon"] = [f"{int(bed.list_exon_starts[i]) + 1}-{int(bed.list_exon_starts[i]) + int(bed.list_exon_lengths[i])}" for i in range(0, len(bed.list_exon_lengths))]
			#
			#intron coordinates
			features_dict["intron"] = [f"{int(bed.list_exon_lengths[i]) + int(bed.list_exon_starts[i]) + 1}-{bed.list_exon_starts[i + 1]}" for i in range(0, len(bed.list_exon_lengths) - 1)]
			#
			#if CDS_start != start (has CDS)
			
			if bed.CDS_start == bed.start and bed.CDS_end == bed.end:#Assumes no CDS
				for feature, coordinates in features_dict.items():
					[out.write(f"{gene},{transcript},{bed.start + 1}-{bed.end},{feature},{coordinate},{gene_start},{gene_end},{orientation},{bed.chromosome}\n") for coordinate in coordinates]
				return bed.start, bed.end
			else:
				#5' UTR coordinates (with gaps for introns)
				if bed.strand == "+":
					f_utr_start = 1 
					f_utr_end = bed.CDS_start - bed.start#if +
				else:
					f_utr_start = bed.CDS_end - bed.start + 1
					f_utr_end = bed.end
				#
				#
				features_dict["5' utr"] = bed.get_f_coordinates(f_utr_start, f_utr_end, features_dict["exon"])
				#
				#CDS coordinates (with gaps for introns)
				#
				cds_start = bed.CDS_start - bed.start + 1
				cds_end = bed.CDS_end - bed.start
				features_dict["coding_region"] = bed.get_f_coordinates(cds_start, cds_end, features_dict["exon"])
				#
				#3' utr coordinates
				if bed.strand == "+":
					thr_utr_start = bed.CDS_end - bed.start + 1
					thr_utr_end = bed.end - bed.start
				else:
					thr_utr_start = 1
					thr_utr_end = bed.CDS_start - bed.start
				features_dict["3' utr"] = bed.get_f_coordinates(thr_utr_start, thr_utr_end, features_dict["exon"])
				#
				#- strand
				if bed.strand == "-":#if '-' strand flip everything around
					features_dict2 = {}

					for feature, coordinates in features_dict.items():
						features_dict2[feature] = bed.convert_coordinates(coordinates, bed.start, bed.end)
					features_dict = features_dict2
				#
				#Write outputs
				
				for feature, coordinates in features_dict.items():
					[out.write(f"{gene},{transcript},{bed.start + 1}-{bed.end},{feature},{coordinate},{gene_start},{gene_end},{orientation},{bed.chromosome}\n") for coordinate in coordinates]
			#gene,transcript,transcript_coordinates,type,coordinates,start, end, orientation, chromosome
			
			#
		#return bed.start, bed.end

			


			


class SnpEffAnnotation:
	def __init__(self, annotation):
		try:
			self.allele, self.p_annotation, self.putative_impact, self.position, self.gene, self.feature_type, self.transcript, self.transcript_type  = annotation.split("|")[:8]
			self.all = annotation.split("|")[:8]
		except ValueError:
			print(f"Error at position {self.position}! Annotation is {annotation}")
			sys.exit()
	def __repr__(self):
		return self.p_annotation

class SnpEff:
	gene_impact = {}
	gene_SNP_annotation = {}
	Instances = {}#Gene: {position1: [annotations], position2: [annotations]}
	Instances_clean = {}
	position_info = {}
	dn_ds = {}
	#
	def __init__(self, line):
		self.fields = line.rstrip().split("\t")
		self.INFO = self.fields[7]
		self.SNP_position = int(self.fields[1])
		self.reference_allele = self.fields[3]
		self.alternative_allele = self.fields[4]
		annotations = self.INFO.split("ANN=")[1].split("=")[0].split(",")#List of all possible annotations
		self.annotations = []
		for a in annotations:
			annotation = SnpEffAnnotation(a)
			self.annotations.append(annotation)
				
			try:
				position_annotations = self.Instances[annotation.gene]
				position_annotations.setdefault(self.SNP_position, []).append(annotation)
				self.Instances[annotation.gene] = position_annotations
			except KeyError:
				positions_annotations = {}
				positions_annotations[self.SNP_position] = [annotation]
				self.Instances[annotation.gene] = positions_annotations
				
		self.position_info[self.SNP_position] = self
				
	def __getitem__(self, item):
		if self.Instances_clean:
			return self.Instances_clean[item]
		else:
			return self.Instances[item]
	#
	def items(self):
		if self.Instances_clean:
			return self.Instances_clean.items()
		else:
			return self.Instances.items()
	#
	def clean_up_annotations_transcript(self):
		"""As below but keeps transcript level information"""
		not_interesting = {"downstream_gene_variant", "upstream_gene_variant", \
		"intergenic_region", "intron_variant"}
		for gene, position_annotation in self.Instances.items():
			position_annotation_clean = {}
			for position, annotations in position_annotation.items():
				annotations_cleaned = [annotation for annotation in annotations if annotation.p_annotation not in not_interesting] 
				if annotations_cleaned:
					position_annotation_clean[position] = annotations_cleaned
			self.Instances_clean[gene] = position_annotation_clean
	#
	def clean_up_annotations(self):
		"""Many positions have multiple annotations.
		This is because of multiple transcripts. Most transcripts will be similar, 
		so will effectively have the same annotation. So remove duplicates. Will
		also remove annotations not of interest e.g upstream/ downstream"""
		
		not_interesting = {"downstream_gene_variant", "upstream_gene_variant", \
		"intergenic_region", "intron_variant"}
		
		for gene, position_annotation in self.Instances.items():
			position_annotation_clean = {}
			for position, annotations in position_annotation.items():
				annotations_cleaned = []
				all_p_annotations = set()
				for annotation in annotations:
					if annotation.p_annotation not in all_p_annotations and annotation.p_annotation not in not_interesting:
						#if annotation.putative_impact == "HIGH" or annotation.putative_impact == "MODERATE":
						annotations_cleaned.append(annotation)
						all_p_annotations.add(annotation.p_annotation)
				if annotations_cleaned:#Only keep high/moderate impact snps
					position_annotation_clean[position] = annotations_cleaned
			self.Instances_clean[gene] = position_annotation_clean
	def __repr__(self):
		string = ""
		print_dict = self.Instances_clean if self.Instances_clean else self.Instances
		
		for gene, position_annotation in print_dict.items():
			string += f"{gene}: \n"
			for position, annotations in position_annotation.items():
				string += f"{position} : {','.join([str(annotation) for annotation in annotations])}\n"
			string += "\n"  
		return string
	#
	def __len__(self):
		return(len(self.Instances.keys()))
	#
	def write_tsv(self, outfile):
		"""Write a .csv file with results of cleaned up annotation. 
		Gene: Position: Annotation"""
		with open(outfile, "w") as out:
			out.write("Gene\tPosition\tAnnotation\tref allele\talt allele\tdn/ds ratio\n")
			for gene, position_annotation in self.Instances_clean.items():
				
				index = 0
				
					
				try:
					dn_ds = self.dn_ds[gene]
				except KeyError:
					dn_ds = ""
					
				for position, annotations in position_annotation.items():
					
					ref = self.position_info[position].reference_allele
					alt = self.position_info[position].alternative_allele
					
					out.write(f"{gene}\t") if not index else out.write(f"\t")
					
					out.write(f"{position}\t{','.join([str(annotation) for annotation in annotations])}\t\
					{ref}\t{alt}\t")#
					
					#Only write annotation once
					out.write(f"{dn_ds}\n") if not index else out.write("\t\n")
					
					index += 1
	#		
	def get_dn_ds(self):
		"""Calculates dn/ds ratio for each gene where possible
		missense_variant = n, synonymous_variant"""
		#
		for gene, position_annotation in self.Instances_clean.items():
			n = []
			s = []
			for position, annotations in position_annotation.items():
				SNP_anotations = set([str(annotation) for annotation in annotations])
				if "missense_variant" in SNP_anotations and "synonymous_variant" in SNP_anotations:
					continue#A case where SNP is annotated as both, these cannot be used
				n += [annotation for annotation in SNP_anotations if annotation == "missense_variant"]
				s += [annotation for annotation in SNP_anotations if annotation == "synonymous_variant"]
			if n and s:
				self.dn_ds[gene] = len(n)/len(s)
			else:
				self.dn_ds[gene] = "NA"
	#
	def write_snp_toRscript(self, gene, transcript, snpeffout):
		"""Appends line of code to Rscript to enable visualisation of snpeff annotation"""
		if gene not in self.Instances.keys():
			return
		for position, snp_annotations in self[gene].items():
			transcript_snp_annotations = [snp_annotation for snp_annotation in snp_annotations if snp_annotation.transcript == transcript]
			if not len(transcript_snp_annotations):#No snp_annotations for this transcript
				continue
			elif len(transcript_snp_annotations) == 1:#In most cases this is true
				snp_annotation = transcript_snp_annotations[0]
				colour, label, depth = labels(snp_annotation.putative_impact, snp_annotation.p_annotation)
				label = label.replace("_variant", "")
				#
			else:#This is a case where one SNP has multiple annotations due SNP annotation being ambiguous
				previous_score = -1#take the highest impacting snp_annotation
				scores = {"HIGH": 3, "MODERATE" : 2, "LOW": 1, "MODIFIER" : 0}#Order of priority
				for snp_annotation in transcript_snp_annotations:
					if scores[snp_annotation.putative_impact] > previous_score:
						colour, label, depth = labels(snp_annotation.putative_impact, snp_annotation.p_annotation)
			#		
			#
			#depth = -0.10
			#
			ref = self.position_info[position].reference_allele
			alt = self.position_info[position].alternative_allele
			snp_length = max([len(ref), len(alt)])
			#Rout.write(f"mutation.plot({position}, {position + snp_length - 1}, text='{label}', col='black', drop={depth}, haplotypes=c('{colour}'))\n\n")
			snpeffout.write(f"{transcript},{position},{position + snp_length - 1}, {label},{depth}, {colour}\n")




class Gtf(Bed):
	gene_transcripts = {}
	gtfentries = []
	Instances = {}
	Transcript_Instances = {}
	genestarts = {}
	geneends = {}
	genestrand = {}
	genechromosome = {}
	transcript_starts = {}
	transcript_ends = {}
	transcript_CDS_starts = {}
	transcript_CDS_ends = {}
	def __init__(self, line, locus=[]):
		self.line = line
		fields = line.rstrip().split("\t")
		self.chromosome = fields[0]
		self.start = int(fields[3])
		self.end = int(fields[4])
		#
		if locus:
			chromosome, start, end = locus
			if self.chromosome != chromosome:
				return
			if self.start < int(start) or self.end > int(end):
				return
		#
		self.assembler = fields[1]
		self.type = fields[2]
		self.strand = fields[6]
		self.info = fields[8]
		self.info_split = self.info.split(";")
		for i in range(0,len(self.info_split)):
			if "transcript" in self.info_split[i]:
				self.transcript = self.info_split[i].split('"')[-2]
				break
		self.gene = self.transcript.split(".")[0]
		self.genechromosome[self.gene] = self.chromosome
		self.line = line
		self.gtfentries.append(self)
		self.Instances.setdefault(self.gene,[]).append(self)
		self.Transcript_Instances.setdefault(self.transcript,[]).append(self)
		self.genestrand[self.gene] = self.strand
		self.gene_transcripts.setdefault(self.gene,set()).add(self.transcript)
		#
		try:
			transcript_start = self.transcript_starts[self.transcript]
			#
			if self.start < transcript_start:
				transcript_start = self.start #Not a stranded start
				self.transcript_starts[self.transcript] = transcript_start
		#
		except KeyError:
			self.transcript_starts[self.transcript] = self.start
		#
		#
		try:
			transcript_end = self.transcript_ends[self.transcript]
			#
			if self.end > transcript_end:
				transcript_end = self.end #Not a stranded start
				self.transcript_ends[self.transcript] = transcript_end
		#
		except KeyError:
			self.transcript_ends[self.transcript] = self.end
		#
		try:
			gene_start = self.genestarts[self.gene]
			#
			if self.start < gene_start:
				gene_start = self.start #Not a stranded start
				self.genestarts[self.gene] = gene_start
		#
		except KeyError:
			self.genestarts[self.gene] = self.start
		#
		#
		try:
			gene_end = self.geneends[self.gene]
			#
			if self.end > gene_end:
				gene_end = self.end #Not a stranded start
				self.geneends[self.gene] = gene_end
		#
		except KeyError:
			self.geneends[self.gene] = self.end
	#
	#
	def __getitem__(self,gene):
		return self.Instances[gene]
	#
	def __repr__(self):
		return self.line
	#
	def get_CDS_start(self):
		for transcript, gtfs in self.Transcript_Instances.items():
			for gtf in gtfs:
				if gtf.type == "CDS":
					try:
						CDS_start = self.transcript_CDS_starts[transcript]
						if gtf.start < CDS_start:
							CDS_start = gtf.start
							self.transcript_CDS_starts[transcript] = CDS_start
					except KeyError:
						self.transcript_CDS_starts[transcript] = gtf.start
					try:
						CDS_end = self.transcript_CDS_ends[transcript]
						if gtf.end > CDS_end:
							CDS_end = gtf.end
							self.transcript_CDS_ends[transcript] = CDS_end
					except KeyError:
						self.transcript_CDS_ends[transcript] = gtf.end
	#
	def summarise(self):
		genes = set()
		transcripts = set()
		for gtf in self.gtfentries:
			genes.add(gtf.gene)
			transcripts.add(gtf.transcript)
		print(f"GTF file has {len(genes)} genes and {len(transcripts)} transcripts")
	#
	def get_gtf_region(self,chromosome,start,end,outfile,zerobased=False):
		with open(f"{outfile}_{chromosome}_{start}_{end}.gtf","w") as out:
			genes = set()
			transcripts = set()
			for gtf in self.gtfentries:
				if gtf.chromosome == chromosome and gtf.start >= int(start) and gtf.end <= int(end):
					genes.add(gtf.gene)
					transcripts.add(gtf.transcript)
					if zerobased:
						newstart = gtf.start - int(start)
						newend = gtf.end - int(start)
						out.write("\t".join(gtf.line.split("\t")[:3]) + "\t" + str(newstart) + "\t" + str(newend) + "\t" + "\t".join(gtf.line.split("\t")[5:]))
					else:
						out.write(gtf.line)
			print(f"Filtered GTF file has {len(genes)} genes and {len(transcripts)} transcripts")
	#
	def write_genemodel(self, gene, transcript, output):
		genemodel_conversion = {"five_prime_utr" : "5' utr", "CDS": "coding_region", "exon" : "exon", "three_prime_utr" : "3' utr"}
		with open(output, "a") as out:
			transcript_start = self.transcript_starts[transcript] - 1
			genestrand = self.genestrand[gene]
			#
			if genestrand == "+":
				orientation = "forward"
			elif genestrand == "-":
				orientation = "reverse"
			else:
				print("Fatal strand error")
				sys.exit()
			chromosome = self.genechromosome[gene]
			gene_start = self.genestarts[gene]
			gene_end = self.geneends[gene]
			transcript_end = self.transcript_ends[transcript] + 1
			strand_start = transcript_start if genestrand == "+" else transcript_end
			exons = []
			for gtf_line in self[gene]:#Convert coordinates to 1 based coordinates with start of gene as reference
				if gtf_line.transcript != transcript:
					continue
				#
				start = gtf_line.start - strand_start if genestrand == "+" else strand_start - gtf_line.end
				end = gtf_line.end - strand_start if genestrand == "+" else strand_start - gtf_line.start
				#
				if gtf_line.type == "exon":
					exons.append((start, end))
				try:
					feature = genemodel_conversion[gtf_line.type]
					out.write(f"{gene},{transcript},{transcript_start}-{transcript_end},{feature},{start}-{end},{gene_start},{gene_end},{orientation},{chromosome}\n")
				except KeyError:
					print(f"No availability for category {gtf_line.type}")
					continue
			exons.sort()
			previous_exon = 0
			for exon in exons:
				if not previous_exon:
					previous_exon = exon[1]
				else:
					out.write(f"{gene},{transcript},{transcript_start}-{transcript_end},intron,{previous_exon + 1}-{exon[0] - 1},{gene_start},{gene_end},{orientation},{chromosome}\n")
					previous_exon = exon[1]
		return transcript_start, transcript_end
		##gene,transcript,transcript_coordinates,type,coordinates,start, end, orientation, chromosome



class Interproscan(Gtf):#Need gtf.Instances attribute
	I_Instances = {}
	def __init__(self,line):#See https://interproscan-docs.readthedocs.io/en/latest/UserDocs.html#output-formats
		self.line = line
		self.fields = line.rstrip().split("\t")
		self.transcript = self.fields[0]
		self.sequence_length = int(self.fields[2])
		self.analysis = self.fields[3]
		self.sig_description = self.fields[5]
		self.start = int(self.fields[6])
		self.end = int(self.fields[7])
		self.score = self.fields[8]
		self.Status = self.fields[9]
		self.Interpro_annotation_description = self.fields[12]
		self.I_Instances.setdefault(self.transcript, []).append(self)
	#
	def __getitem__(self, transcript):
		return self.I_Instances[transcript]
	#
	@staticmethod
	def fix_interpro_co(CDS_start, interpro_co, exon_coordinates):
		"""Interpro coordinates do not take into account introns. So fix"""
		previous_end = 1
		CDS_past = False
		for exon in exon_coordinates:
			e_start, e_end = [int(co) for co in exon.split("-")]
			intron = e_start - previous_end
			if CDS_past:#CDS may start after many introns. Only add introns to interpro after cds start
				interpro_co += intron
			if interpro_co <= e_end and interpro_co >= e_start:
				return interpro_co
			previous_end = e_end
			if CDS_start <= e_end and CDS_start >= e_start:
				CDS_past = True
		print("Coordinate adjustment error!")
		pdb.set_trace()
		#sys.exit()
	#
	@staticmethod
	def sort_exons(exon_coordinates):
		"""sort list of exon coordinates in ascending order"""
		exons = []
		for exon in exon_coordinates:
			s, e = exon.split("-")
			exons.append((int(s), int(e)))
		exons.sort()
		return ["-".join((str(exon[0]), str(exon[1]))) for exon in exons]
	#
	def write_csv(self, filetype, transcript, output):
		"""Writes a csv file in the format required for putting domains in R gene model"""
		gene = transcript.split(".")[0]
		with open(output, "a") as out:
			for interpro in self[transcript]:
				if interpro.Interpro_annotation_description == "-":
					continue
				if filetype == "Bed":#If the input annotation is bed
					bed = self.b_Instances[transcript]
					transcript_start = bed.start + 1
					transcript_end = bed.end
					
					if bed.strand == "+":#Start is the start!
						CDS_start1 = bed.CDS_start - bed.start
						exon_coordinates = [f"{int(bed.list_exon_starts[i]) + 1}-{int(bed.list_exon_starts[i]) + int(bed.list_exon_lengths[i])}" for i in range(0, len(bed.list_exon_lengths))]
						
					elif bed.strand == "-":
						CDS_start1 = bed.end - bed.CDS_end
						exon_coordinates_reverse = [f"{int(bed.list_exon_starts[i]) + 1}-{int(bed.list_exon_starts[i]) + int(bed.list_exon_lengths[i])}" for i in range(0, len(bed.list_exon_lengths))]
						exon_coordinates = self.convert_coordinates(exon_coordinates_reverse,bed.start, bed.end)
						#exon_coordinates = [e.split("-") for e in exon_coordinates]
					else:
						print("Major strand error")
				elif filetype == "Gtf":
					transcript_start = self.transcript_starts[transcript]
					transcript_end = self.transcript_ends[transcript]
					if self.genestrand[gene] == "+":#Start is the start!
						transcript_start = self.transcript_starts[transcript]
						CDS_start1 = self.transcript_CDS_starts[transcript] - transcript_start #CDS start relative to transcript start
						exon_coordinates = [f"{gtf.start - (transcript_start - 1)}-{gtf.end - (transcript_start - 1)}" for gtf in self.Transcript_Instances[transcript] if gtf.type == "exon"]
						
					else:
						transcript_start = self.transcript_ends[transcript]#
						CDS_start1 = transcript_start - self.transcript_CDS_ends[transcript]  #CDS
						exon_coordinates = [f"{(transcript_start + 1) - gtf.end}-{(transcript_start + 1) - gtf.start}" for gtf in self.Transcript_Instances[transcript]if gtf.type == "exon"] #
				else:
					print("Filetype error!")
				#
				interpro_start_raw = (interpro.start * 3) + CDS_start1
				interpro_end_raw = (interpro.end * 3) + CDS_start1
				exon_coordinates = self.sort_exons(exon_coordinates)
				interpro_start = self.fix_interpro_co(CDS_start1, interpro_start_raw, exon_coordinates)
				interpro_end = self.fix_interpro_co(CDS_start1, interpro_end_raw, exon_coordinates)
				#
				interpro_exons = self.get_f_coordinates(interpro_start, interpro_end, exon_coordinates)
				domain = interpro.Interpro_annotation_description.split(",")[0]
				#pdb.set_trace()
				for n, exon in enumerate(interpro_exons):
					out.write(f"{gene},{transcript},{transcript_start}-{transcript_end},{domain},{exon}\n")
					
					















def labels(impact,label):
	if impact == "MODERATE":
		colour = "orange"
		label = ""
		depth = random.uniform(-.3, -.4)
	elif impact == "HIGH":
		colour = "red"
		depth = random.uniform(-.5, -.6)
	else:
		colour = "blue"
		label = ""
		depth = random.uniform(-.1, -.2)
	return colour, label, depth




def write_R_script(outdir, gene_output, domain_output, snpeff, annotation, interproscan, filetype):
	"""Put entries in format for R genemodel"""
	with open(outdir + "snpeff_out.csv","w") as snpeffout:
		#
		
		snpeffout.write("Transcript,position_start, position_end, label, depth, colour\n")
		

		with open(gene_output, "w") as out:
			out.write("gene,transcript,transcript_coordinates,type,coordinates,start, end, orientation, chromosome\n")
		#
		if interproscan:
			with open(domain_output, "w") as out:
				out.write("gene,transcript,transcript_coordinates,type,coordinates\n")
		#
		#
		for gene, transcripts in annotation.gene_transcripts.items():
			#pdb.set_trace()
			transcripts = list(transcripts) #with gtf transcripts are a set
			transcripts.sort()#So we go through .1 -> .n
			for transcript in transcripts:
				#
				annotation.write_genemodel(gene, transcript, gene_output)
				#
				
				if interproscan:
					try:
						interproscan.write_csv(filetype, transcript, domain_output)
						if transcript == "BaRT2v18chr3HG122970.1":
							interproscan.write_csv(filetype, transcript, "BaRT2v18chr3HG122970.1.bed.csv")
						
					except KeyError:
						pass
				#
				if snpeff:
					snpeff.write_snp_toRscript(gene, transcript, snpeffout)
				


def main():

	parser = argparse.ArgumentParser(description='Plot transcripts and annotate with Interproscan domains and snpEff. Output is and R script for visualisation. All arguments are optional, but you will need a gtf or a bed to start.')

	parser.add_argument('-gtf', dest = 'gtf_input', help = 'gtf annotation file', type = str, default = "")

	parser.add_argument('-bed', dest = 'bed_input', help = 'bed annotation file', type = str, default = "")

	parser.add_argument('-c', dest = 'coordinates', help = 'Physical location of interest in format <chromosome>,<start>,<end>', type = str, default = "")

	parser.add_argument('-snpeff', dest = 'snpeff_input_file', help = 'SnpEff vcf file', type = str, default = "")

	parser.add_argument('-interpro', dest = 'interproscan_input', help = 'Interproscan input. Needs to be in tsv format', type = str, default = "")

	parser.add_argument('-o', dest = 'Routfile', help = 'outfile (must specify as .R). Default is "gene_models.R". All other files will be written to the directory the R script is in.', type = str, default = "gene_models.R")

	args = parser.parse_args()

	if args.gtf_input and args.bed_input:
		print("You only need one input annotation (bed or gtf)! Try again")
		sys.exit()


	if args.coordinates:
		locus = args.coordinates.split(",")#3H locus boundaries on 3H
	else:
		locus = []
	#
	#Outfiles:
	if len(args.Routfile.split("/")) == 1:#No /
		outdir = ""
		gene_output = "genemodel_input.csv"
		domain_output = "genemodel_input_domain.csv"

	else:
		outdir = "/".join(args.Routfile.split("/")[:-1]) + "/"
		gene_output = outdir + "genemodel_input.csv"
		domain_output = outdir + "genemodel_input_domain.csv"
	
	#
	#
	#
	if args.snpeff_input_file:
		print("Parsing SNPeff input ...")
		for line in open(args.snpeff_input_file):
			#print(line)
			if line.startswith("#"):
				continue
			#print(line)
			
			if locus:
				chromosome, position = line.rstrip().split("\t")[0:2]
				if chromosome != locus[0]:
					continue
				if int(position) < int(locus[1]) or int(position) > int(locus[2]):#ONly include snps in 3H locus
					continue
			snpeff = SnpEff(line)
		#
		snpeff.clean_up_annotations_transcript()
		snpeff.write_tsv(outdir + "SNPeff_output.tsv")
	else:
		snpeff = {}
	#
	if args.gtf_input:
		print("Parsing gtf transript annotation...")
		filetype = "Gtf"
		for line in open(args.gtf_input):
			annotation = Gtf(line, locus)
		#
		annotation.get_CDS_start()
	elif args.bed_input:
		print("Parsing bed...")
		filetype = "Bed"
		for line in open(args.bed_input):
			annotation = Bed(line, locus)#Will just take in coordinates within the locus
	else:
		print("Input error! Needs bed or gtf input")

	#
	annotation.summarise()
	#
	
	#
	if args.interproscan_input:
		print("Parsing interproscan tsv...")
		for line in open(args.interproscan_input):
			interproscan = Interproscan(line)
	else:
		interproscan = {}
	#
	print("Writing R script and csvs files")
	#
	write_R_script(outdir, gene_output, domain_output, snpeff, annotation, interproscan, filetype)
	#	
	print("Script finished")


if __name__ == "__main__":
	main()


		








