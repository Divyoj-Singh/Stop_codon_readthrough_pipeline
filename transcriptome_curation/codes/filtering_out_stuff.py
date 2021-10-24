# initializing variables

from Bio import SeqIO

file_type = "ISR"

def return_good_genes():
	gene_all = []
	gene_dict = {}
	degenerate_isr = []
	matching_isr = []
	with open("./../output_files/annotation/annotation_file_"+file_type+".txt") as file:
		next(file)
		for line in file:
			a = line[:-1].split("\t")
			gene_all.append(a[0])
			gene_id = a[0]
			start = int(a[1])
			stop1 = int(a[2])
			stop2 = int(a[3])
			length = int(a[4])
			cds_len = stop1 - start + 1
			isr_len = stop2 - stop1
			utr3_len = length - stop2
			gene_dict[gene_id] = [start,stop1,stop2,length,cds_len,isr_len,utr3_len]
	gene_all = set(gene_all)
	with open("./../output_files/blast_results/"+file_type+"_CDS.txt") as file:
		for line in file:
			a = line[:-1].split("\t")
			degenerate_isr.append(a[0])
			#print(a[0])
	degenerate_isr = set(degenerate_isr)
	with open("./../output_files/blast_results/"+file_type+"_"+file_type+".txt") as file:
		for line in file:
			a = line[:-1].split("\t")
			if a[0].split(".")[0] != a[1].split(".")[0]:
				matching_isr.append(a[0])
				matching_isr.append(a[1])
				#print(a[0],"\t",a[1])
	matching_isr = set(matching_isr)
	good_genes = list(gene_all - degenerate_isr - matching_isr)
	print("# genes before blasting",len(gene_all))
	print("# genes in which ISR matched with some part of CDS", len(degenerate_isr))
	print("# genes in which ISR matched with some part of other ISR", len(matching_isr))
	print("# genes remaining after removing degenrate ISRs", len(good_genes))
	#for i in good_genes:
	#	print(i)
	return gene_dict, good_genes

gene_dict, good_genes = return_good_genes()
with open ("./../output_files/annotation/annotation_file_"+file_type+"_after_removing_degenerate_isr.txt","w") as file:
	for gene in good_genes:
		file.write(str(gene)+"\t"+str(gene_dict[gene][0])+"\t"+str(gene_dict[gene][1])+"\t"+str(gene_dict[gene][2])+"\t"+str(gene_dict[gene][3])+"\t"+str(gene_dict[gene][4])+"\t"+str(gene_dict[gene][5])+"\t"+str(gene_dict[gene][6])+"\n")

for record in SeqIO.parse("./../cdna/Arabidopsis_thaliana.TAIR10.cdna.all.fa", "fasta"):
	if str(record.id) in good_genes:
		with open ("./../output_files/nucleotide/Arabidopsis_thaliana_fasta_file_after_removing_degenerate_isr.fasta","a+") as file:
			file.write(">"+str(record.description)+"\n")
			file.write(str(record.seq)+"\n")
