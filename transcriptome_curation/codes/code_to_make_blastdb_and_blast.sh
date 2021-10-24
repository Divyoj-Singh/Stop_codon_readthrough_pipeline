path_to_makeblastdb="/home/csb/ncbi-blast-2.9.0+/bin/makeblastdb"
path_to_blastn="/home/csb/ncbi-blast-2.9.0+/bin/blastn"
$path_to_makeblastdb -dbtype nucl -in ./../output_files/nucleotide/plant_transcriptome_CDS_nucleotide.fasta -hash_index -out ./../output_files/blast_results/blast_db/plant_transcriptome_CDS_nucleotide
$path_to_makeblastdb -dbtype nucl -in ./../output_files/nucleotide/plant_transcriptome_ISR_nucleotide.fasta -hash_index -out ./../output_files/blast_results/blast_db/plant_transcriptome_ISR_nucleotide

$path_to_blastn -query ./../output_files/nucleotide/plant_transcriptome_ISR_nucleotide.fasta -db ./../output_files/blast_results/blast_db/plant_transcriptome_CDS_nucleotide -out ./../output_files/blast_results/ISR_CDS.txt -strand plus -word_size 24 -outfmt 6 -num_threads 10 -perc_identity 100
$path_to_blastn -query ./../output_files/nucleotide/plant_transcriptome_ISR_nucleotide.fasta -db ./../output_files/blast_results/blast_db/plant_transcriptome_ISR_nucleotide -out ./../output_files/blast_results/ISR_ISR.txt -strand plus -word_size 24 -outfmt 6 -num_threads 10 -perc_identity 100
