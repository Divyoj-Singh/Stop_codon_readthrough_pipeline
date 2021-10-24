from Bio import SeqIO

dict_detail={}

path_to_gene_details_file="./../gene_details/cds.txt"
path_isr_fasta_file="./../output_files/protein/plant_transcriptome_ISR_protein.fasta"
path_to_mass_spec_file="./../output_files/protein/ISR_protein_mass_spec.fasta"

with open(path_to_gene_details_file) as file1:
	for line in file1:
		a=line.split("\t")
		if a[0] not in dict_detail:
			dict_detail[a[0]]=a[1]

for record in SeqIO.parse(path_isr_fasta_file, "fasta"):
	isr=record.seq	
	new_isr=isr[1:-1]
	with open (path_to_mass_spec_file,"a+") as file:
			file.write(">"+"sp|"+str(record.id)+"|"+str(dict_detail[record.id])+" OS=Arabidopsis thaliana OX=3702 GN=GENE PE=2 SV=1"+"\n")
			file.write(str(new_isr)+"\n")

'''
>sp|Q9MAJ7|BGAL5_ARATH Beta-galactosidase 5 OS=Arabidopsis thaliana OX=3702 GN=BGAL5 PE=2 SV=1
'''