from Bio import SeqIO

gene_list = []

with open("combined_transcriptome_ISR_FSP_FSM.fasta","a+") as f:
    for i in ["ISR","FSP","FSM"]:
        for seq_record in SeqIO.parse("./Arabidopsis_thaliana_fasta_file_"+i+".fasta", "fasta"):
            gene_id = seq_record.id
            gene_desc = seq_record.description
            gene_seq = seq_record.seq
            if gene_id not in gene_list:
                gene_list.append(gene_id)
                f.write(">"+gene_desc+"\n"+str(gene_seq)+"\n")
print("Lenght of transcriptome = ",len(gene_list))

