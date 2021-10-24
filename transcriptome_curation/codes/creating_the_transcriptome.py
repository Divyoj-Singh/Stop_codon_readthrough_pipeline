import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

def function_to_filter_out_gene_details(gene_description):
 #### this function reads the description of the gene and separates all the relevant details from it.    
    gene_description = gene_description.split(" ")
    
    gene = "NA"
    gene_symbol = "NA"
    gene_biotype = "NA"
    transcript_biotype = "NA"
    transcript_broad_type = gene_description[1]

    for i in gene_description:
        if "gene:" in i:
            gene = i[len("gene:"):]
        elif "gene_symbol:" in i:
            gene = i[len("gene_symbol:"):]
        elif "gene_biotype:" in i:
            gene_biotype = i[len("gene_biotype:"):]
        elif "transcript_biotype:" in i:
            transcript_biotype = i[len("transcript_biotype:"):]

    return gene,transcript_broad_type,gene_symbol,gene_biotype,transcript_biotype


def reading_fasta_file(path_to_the_fasta_file,path_to_the_gene_details_folder,broad_type):
# filtering out and creating a specific file for the rRNA and tRNA sequences from the ncRNA file
### this function reads creats the gene details folder and writes the 4 files in to it which have details  from description.
    if not os.path.exists(path_to_the_gene_details_folder):
        os.mkdir(path_to_the_gene_details_folder)
    
    with open (path_to_the_gene_details_folder+broad_type+".txt","w") as file_details:
        
        file_details.write("gene_id\tgene_name\tgene_symbol\tgene_biotype\ttranscript_biotype\tgene_description\n")
        
        for seq_record in SeqIO.parse(path_to_the_fasta_file, "fasta"):
            
            gene_id = seq_record.id
            gene_description = seq_record.description
            
            gene,transcript_broad_type,gene_symbol,gene_biotype,transcript_biotype = function_to_filter_out_gene_details(gene_description)
            if transcript_broad_type != broad_type:
                print("Warning! Mismatch in transcipt broad type: ",broad_type)

            if "description:" in gene_description:
                gene_description = gene_description.split("description:")[-1]
                if "[Source:" in gene_description:
                    gene_description = gene_description.split("[Source:")[0]
            else:
                gene_description = "NA"

            file_details.write(gene_id+"\t"+gene+"\t"+gene_symbol+"\t"+gene_biotype+"\t"+transcript_biotype+"\t"+gene_description+"\n")

def creating_rRNA_tRNA_snRNA_snoRNA_miRNA_file_and_ncRNA_file(path_to_the_ncRNA_fasta_file,path_to_rRNA_tRNA_snRNA_snoRNA_miRNA_file,path_to_ncRNA_only_file):
    
    counter_removal = 0
    counter_ncRNA = 0

    with open (path_to_ncRNA_only_file,"w") as file_for_ncRNA:
        with open (path_to_rRNA_tRNA_snRNA_snoRNA_miRNA_file,"w") as file_for_removal:
            for seq_record in SeqIO.parse(path_to_the_ncRNA_fasta_file, "fasta"):

                    gene_id = str(seq_record.id)
                    gene_description = seq_record.description
                    gene_sequence = str(seq_record.seq)
            
                    gene,transcript_broad_type,gene_symbol,gene_biotype,transcript_biotype = function_to_filter_out_gene_details(gene_description)

                    if gene_biotype in ["snRNA","rRNA","tRNA","snoRNA","miRNA"]:
                        file_for_removal.write(">"+gene_id+"\n"+str(gene_sequence)+"\n")
                        counter_removal += 1
                    elif gene_biotype in ["lncRNA","ncRNA"]:
                        file_for_ncRNA.write(">"+gene_description+"\n"+gene_sequence+"\n")
                        counter_ncRNA += 1

    print("Created the file for preprocessing of Ribosomal profiling data: ")
    print("\tthat contains ",counter_removal," gene sequences which includes rRNA tRNA snRNA snoRNA and miRNA. ")

    print("Created the file for non-coding RNA: ")
    print("\tthat contains ",counter_ncRNA," gene sequences which includes ncRNA and lncRNA. ")

path_to_the_gene_details_folder = "./../gene_details/"
broad_types = ["cdna","cds","ncrna","pep"]
####################### take care of path here:##############################
for transcript_broad_type in broad_types:
    print("Extracting the gene details for: ",transcript_broad_type)
    path_to_the_fasta_file = "./../"+transcript_broad_type+"/Arabidopsis_thaliana.TAIR10."+transcript_broad_type+".all.fa"
    reading_fasta_file(path_to_the_fasta_file,path_to_the_gene_details_folder,transcript_broad_type)

path_to_the_ncRNA_fasta_file = "./../ncrna/Arabidopsis_thaliana.TAIR10.ncrna.all.fa"

if not os.path.exists("./../gene_files/"):
    os.mkdir("./../gene_files/")

path_to_rRNA_tRNA_snRNA_snoRNA_miRNA_file = "./../gene_files/rRNA_tRNA_snRNA_snoRNA_and_miRNA.fasta"
path_to_ncRNA_only_file = "./../gene_files/ncRNA_and_lncRNA.fasta"
creating_rRNA_tRNA_snRNA_snoRNA_miRNA_file_and_ncRNA_file(path_to_the_ncRNA_fasta_file,path_to_rRNA_tRNA_snRNA_snoRNA_miRNA_file,path_to_ncRNA_only_file)
############################################################################## 
def reading_the_cdna_and_cds_file(path_to_cdna_file,path_to_cds_file):
    dict_cdna_sequences = {}
    dict_cds_sequences = {}

    for seq_record in SeqIO.parse(path_to_cdna_file, "fasta"):
            gene_id = seq_record.id
            gene_sequence = seq_record.seq
            if gene_id not in dict_cdna_sequences:
                dict_cdna_sequences[gene_id] = gene_sequence
            else: 
                print("Warning! Duplication found in cdna sequence: ",gene_id)
    
    for seq_record in SeqIO.parse(path_to_cds_file, "fasta"):
            gene_id = seq_record.id
            gene_sequence = seq_record.seq
            if gene_id not in dict_cds_sequences:
                dict_cds_sequences[gene_id] = gene_sequence
            else: 
                print("Warning! Duplication found in cds sequence: ",gene_id)
    
    print("Number of cdna sequences: ",len(dict_cdna_sequences))
    print("Number of cds sequences: ",len(dict_cds_sequences))
    
    cds_that_have_a_cdna = list(set(dict_cdna_sequences.keys()).intersection(set(dict_cds_sequences.keys())))
    print("Number of CDS sequences that have a corresponding cDNA: ",len(cds_that_have_a_cdna))

    return dict_cdna_sequences,dict_cds_sequences,cds_that_have_a_cdna

def computing_the_stop_codon_positions(utr3_sequence):
    inframe = -1
    for i,j in enumerate(utr3_sequence):
        if str(utr3_sequence[i:i+3]) in ["TAG","TGA","TAA"]:
            if i%3 == 0:    
                if inframe == -1:
                    inframe = i+2
        
    return inframe

def initialising_output_files(path_to_output_directory_CDS_ISR_nucleotide_and_protein_sequences):
    
    if not os.path.exists(path_to_output_directory_CDS_ISR_nucleotide_and_protein_sequences+"nucleotide/"):
        os.mkdir(path_to_output_directory_CDS_ISR_nucleotide_and_protein_sequences+"nucleotide/")
    if not os.path.exists(path_to_output_directory_CDS_ISR_nucleotide_and_protein_sequences+"protein/"):
        os.mkdir(path_to_output_directory_CDS_ISR_nucleotide_and_protein_sequences+"protein/")
    if not os.path.exists(path_to_output_directory_CDS_ISR_nucleotide_and_protein_sequences+"annotation/"):
        os.mkdir(path_to_output_directory_CDS_ISR_nucleotide_and_protein_sequences+"annotation/")
    
    with open(path_to_output_directory_CDS_ISR_nucleotide_and_protein_sequences+"annotation/annotation_file_ISR.txt","w") as f_isr:
        f_isr.write("gene_id\tstart\tstop1\tstop2\tlen_cdna\n")
    with open(path_to_output_directory_CDS_ISR_nucleotide_and_protein_sequences+"nucleotide/plant_transcriptome_CDS_nucleotide.fasta","w"):
        pass
    with open(path_to_output_directory_CDS_ISR_nucleotide_and_protein_sequences+"nucleotide/plant_transcriptome_ISR_nucleotide.fasta","w"):
        pass
    with open(path_to_output_directory_CDS_ISR_nucleotide_and_protein_sequences+"protein/plant_transcriptome_CDS_protein.fasta","w"):
        pass
    with open(path_to_output_directory_CDS_ISR_nucleotide_and_protein_sequences+"protein/plant_transcriptome_ISR_protein.fasta","w"):
        pass

def creating_the_transciptome_annotation_file_for_ISR(path_to_output_directory_CDS_ISR_nucleotide_and_protein_sequences,gene_id,start_0_based,stop_0_based,length_of_cdna,new_stop,type_of_file):
    with open(path_to_output_directory_CDS_ISR_nucleotide_and_protein_sequences+"/annotation/annotation_file_"+type_of_file+".txt","a+") as file_type:
        file_type.write(gene_id+"\t"+str(start_0_based)+"\t"+str(stop_0_based)+"\t"+str(new_stop)+"\t"+str(length_of_cdna)+"\n")

def writing_to_the_fasta_file(path,type_file,seq_type,gene_id,sequence):
    with open(path+seq_type+"/plant_transcriptome_"+type_file+"_"+seq_type+".fasta","a+") as f:
        f.write(">"+gene_id+"\n"+str(sequence)+"\n")

def finding_the_inframe_stop_codons(dict_cdna_sequences,dict_cds_sequences,cds_that_have_a_cdna,path_to_output_directory_CDS_ISR_nucleotide_and_protein_sequences,minimum_length_of_3utr=45):
    
    initialising_output_files(path_to_output_directory_CDS_ISR_nucleotide_and_protein_sequences)

    num_genes_not_a_multiple_of_3 = 0
    num_genes_no_valid_stop_codon = 0
    num_genes_that_have_small_utr = 0
    num_genes_that_have_small_isr = 0
    num_of_valid_genes = 0

    for gene_id in cds_that_have_a_cdna:
        if dict_cds_sequences[gene_id] not in dict_cdna_sequences[gene_id]:
            print("Warning! The following cds doesn't lie in the corresponding cdna: ",gene_id)
        else:
            CDS_nucleotide_sequence = dict_cds_sequences[gene_id]
            if len(CDS_nucleotide_sequence) % 3 != 0:
                num_genes_not_a_multiple_of_3 += 1
                continue
            if len(CDS_nucleotide_sequence) < 6 or CDS_nucleotide_sequence[-3:] not in ["TGA","TAG","TAA"]:
                num_genes_no_valid_stop_codon += 1
                continue
            CDS_protein_sequence = str(dict_cds_sequences[gene_id].translate())
            writing_to_the_fasta_file(path_to_output_directory_CDS_ISR_nucleotide_and_protein_sequences,"CDS","nucleotide",gene_id,CDS_nucleotide_sequence)
            writing_to_the_fasta_file(path_to_output_directory_CDS_ISR_nucleotide_and_protein_sequences,"CDS","protein",gene_id,CDS_protein_sequence)
            
            ###### Making the start, stop annotation here.
            length_of_cdna = len(dict_cdna_sequences[gene_id])
            start_0_based = dict_cdna_sequences[gene_id].find(dict_cds_sequences[gene_id]) + 1 ## i based
            
            stop_0_based = start_0_based + len(dict_cds_sequences[gene_id])-1# 1 based
            
            utr3_sequence = dict_cdna_sequences[gene_id][stop_0_based:]
            if len(utr3_sequence) <= minimum_length_of_3utr:
                num_genes_that_have_small_utr += 1
                continue
            else:
                inframe = computing_the_stop_codon_positions(utr3_sequence)
                #print(gene_id,"\t",inframe,"\t",outframe_minus,"\t",outframe_plus,"\n")
                
                if inframe < 45:
                    num_genes_that_have_small_isr += 1
                else:
                    inframe_stop_codon = stop_0_based + inframe +1
                    #print(dict_cdna_sequences[gene_id][inframe_stop_codon-3:inframe_stop_codon+5])
                    ISR_nucleotide_sequence = dict_cdna_sequences[gene_id][stop_0_based-3:inframe_stop_codon]
                    ISR_protein_sequence = str(ISR_nucleotide_sequence.translate())
                    creating_the_transciptome_annotation_file_for_ISR(path_to_output_directory_CDS_ISR_nucleotide_and_protein_sequences,gene_id,start_0_based,stop_0_based,length_of_cdna,inframe_stop_codon,"ISR")
                    num_of_valid_genes += 1
                    writing_to_the_fasta_file(path_to_output_directory_CDS_ISR_nucleotide_and_protein_sequences,"ISR","nucleotide",gene_id,ISR_nucleotide_sequence)
                    writing_to_the_fasta_file(path_to_output_directory_CDS_ISR_nucleotide_and_protein_sequences,"ISR","protein",gene_id,ISR_protein_sequence)

    print("Number of genes that do not have a valid stop codon: ",num_genes_no_valid_stop_codon)
    print("Number of genes that have a length of CDS that is not a multiple of 3: ",num_genes_not_a_multiple_of_3)
    print("Number of genes that have a small utr(<45): ",num_genes_that_have_small_utr)
    print("Number of genes that have a small isr(<45): ",num_genes_that_have_small_isr)
    print("Number of valid ISRs: ",num_of_valid_genes)

def creating_the_transciptome_for_RT_alignment(path_to_cdna_file,path_to_cds_file,path_to_output_directory_CDS_ISR_nucleotide_and_protein_sequences):
    dict_cdna_sequences,dict_cds_sequences,cds_that_have_a_cdna = reading_the_cdna_and_cds_file(path_to_cdna_file,path_to_cds_file)
    finding_the_inframe_stop_codons(dict_cdna_sequences,dict_cds_sequences,cds_that_have_a_cdna,path_to_output_directory_CDS_ISR_nucleotide_and_protein_sequences)

path_to_output_directory_CDS_ISR_nucleotide_and_protein_sequences = "./../output_files/"
if not os.path.exists("./../output_files/"):
    os.mkdir("./../output_files/")
########################################3 take care of path here:#####################    
creating_the_transciptome_for_RT_alignment("./../cdna/Arabidopsis_thaliana.TAIR10.cdna.all.fa","./../cds/Arabidopsis_thaliana.TAIR10.cds.all.fa",path_to_output_directory_CDS_ISR_nucleotide_and_protein_sequences)
######################################### ############################################
