###############################################################################
# Counter for no of files done
counter=0
# declaring count to be an integer
declare -i counter
###############################################################################

### Declaring paths to hard disks
###### paths to ssd drive 
path_to_ssd='/media/csb/D/Plant_Readthrough'
###### paths to external drive
path_to_hdd='/media/csb/D/Plant_Readthrough'

# defining various paths

path_to_compressed_adapter_rrna_removed_fasta_files=$path_to_hdd'/compressed_adapter_rrna_removed_files/'

path_to_aligned_bam_file=$path_to_hdd'/bam_files/'

path_to_transcriptome=$path_to_hdd'/transcriptome/transcriptome'

path_to_stats_folder=$path_to_hdd'/stats'

#this will remain same for all species
path_to_transc_aligned_sam_file=$path_to_ssd'/transc_align_sam_files/'

path_to_decompressed_adapter_rrna_removed_files=$path_to_ssd'/decompressed_adapter_rrna_removed_files/'

# main loop
while IFS=" " read -r word1 remainder; 
do 
#word1=${word1:0:-1}

echo "running on file: "$word1

cd $path_to_compressed_adapter_rrna_removed_fasta_files


################### ###########################
echo "decompress rrna_adapter_rem_file"
time 7z e $word1.7z -o$path_to_decompressed_adapter_rrna_removed_files

cd ~

#echo "Aligning with transciptome:"
#time bowtie2 -p 11 --very-sensitive-local -k 15 --no-unal -U $path_to_decompressed_adapter_rrna_removed_files$word1.fastq -x $path_to_transcriptome -S $path_to_transc_aligned_sam_file$word1.sam |& tee $path_to_stats_folder/transc_align/fout5_$word1.txt

echo " Step 2:Aligning with ISR transcriptome- end to end "
time bowtie2 -p 11 -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 --end-to-end -k 15 --no-unal -U $path_to_decompressed_adapter_rrna_removed_files$word1.fastq -x $path_to_transcriptome -S $path_to_transc_aligned_sam_file$word1.sam |& tee $path_to_stats_folder/transc_align/$word1.txt


echo "deleting fastq file: "
rm -rf $path_to_decompressed_adapter_rrna_removed_files$word1.fastq

echo "converting sam to bam"
time samtools view -S -@ 10 -b $path_to_transc_aligned_sam_file$word1.sam > $path_to_aligned_bam_file$word1.bam

echo "deleting sam file: "
rm -f $path_to_transc_aligned_sam_file$word1.sam

# updating the counter to keep track of files that are done
counter=$counter+1

echo "NUMBER OF FILES DONE                                                       "$counter


done <./../../filenames/all_files.txt