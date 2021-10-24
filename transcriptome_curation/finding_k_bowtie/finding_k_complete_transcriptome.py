from Bio import SeqIO

dic={}
for seq_record in SeqIO.parse("combined_transcriptome_ISR_FSP_FSM.fasta", "fasta"):
	gene_id = seq_record.id
	a = str(gene_id.split(".")[0])
	if a not in dic:
		dic[a]={"freq":0}
	dic[a]["freq"]+=1

for i in dic:
	print(i,"\t",dic[i]["freq"])
