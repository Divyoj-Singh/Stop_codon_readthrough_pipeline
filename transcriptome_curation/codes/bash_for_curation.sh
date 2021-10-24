### After making 4 folders of ncrna, cds, cdna and pep and their fasta files with proper name
# eg: ./../cds/Arabidopsis_thaliana.TAIR10.cds.all.fa

### run the python code to make the first version of transcriptome, rrna and annotation:
python3 ./creating_the_transcriptome.py

## now install blast beforehand and add the path to blast in the following lines:

path_to_makeblastdb="/home/csb/ncbi-blast-2.9.0+/bin/makeblastdb"
path_to_blastn="/home/csb/ncbi-blast-2.9.0+/bin/blastn"

### Running blast to find degenerate ISRs  
$path_to_makeblastdb -dbtype nucl -in ./../output_files/nucleotide/plant_transcriptome_CDS_nucleotide.fasta -hash_index -out ./../output_files/blast_results/blast_db/plant_transcriptome_CDS_nucleotide
$path_to_makeblastdb -dbtype nucl -in ./../output_files/nucleotide/plant_transcriptome_ISR_nucleotide.fasta -hash_index -out ./../output_files/blast_results/blast_db/plant_transcriptome_ISR_nucleotide

$path_to_blastn -query ./../output_files/nucleotide/plant_transcriptome_ISR_nucleotide.fasta -db ./../output_files/blast_results/blast_db/plant_transcriptome_CDS_nucleotide -out ./../output_files/blast_results/ISR_CDS.txt -strand plus -word_size 24 -outfmt 6 -num_threads 10 -perc_identity 100
$path_to_blastn -query ./../output_files/nucleotide/plant_transcriptome_ISR_nucleotide.fasta -db ./../output_files/blast_results/blast_db/plant_transcriptome_ISR_nucleotide -out ./../output_files/blast_results/ISR_ISR.txt -strand plus -word_size 24 -outfmt 6 -num_threads 10 -perc_identity 100
### After this run the following python code to filter out degenerate ISRs and writing new transcriptome.
python3 ./filtering_out_stuff.py

### finallly making another folder in which we will keep the rRANA fasta, transcriptome fasta and annotation.
mkdir ./../final_output_files
mv ./../gene_files/rRNA_tRNA_snRNA_snoRNA_and_miRNA.fasta ./../final_output_files/rRNA_tRNA_snRNA_snoRNA_and_miRNA.fasta
cp ./../output_files/annotation/annotation_file_ISR_after_removing_degenerate_isr.txt ./../final_output_files/annotation_file.txt
cp ./../output_files/nucleotide/Arabidopsis_thaliana_fasta_file_after_removing_degenerate_isr.fasta ./../final_output_files/transcriptome.fasta
