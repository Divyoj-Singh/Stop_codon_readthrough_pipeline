# counter for number of files after which compression starts
count=0
# declaring count to be an integer
declare -i count
# declaring pids to be an array storing the process ids of compression
declare -a pids
###############################################################################
# Counter for no of files done
counter=0
# declaring count to be an integer
declare -i counter
###############################################################################
declare -a list_of_failed_runs
###############################################################################
###############################################################################
# main loop
while IFS=" " read -r word1 word2 remainder; 
do 

# removing the mysterious NULL
#word2=${word2:0:-1}## this needs to be checked

echo "Adapter seq is: "$word2
###########################################3
# defining various paths

### Declaring paths to hard disks
###### paths to ssd drive 
path_to_ssd='/media/csb/D/Plant_Readthrough'
###### paths to external drive
path_to_hdd='/media/csb/D/Plant_Readthrough'

### paths to programs  fastp and bowtie2 in the code later
path_to_prefetch='prefetch'
#path_to_fastqdump='/home/smelab/sratoolkit/bin/fastq-dump'
path_to_parallel_fastq_dump="/home/csb/parallel-fastq-dump-master/parallel-fastq-dump"
### path to folders
path_to_sra_folder='/home/csb/ncbi/public/sra/'

path_to_adapter_removed_files_folder=$path_to_ssd'/adapter_removed_files/'
path_to_folder_for_sam_files_to_be_deleted=$path_to_ssd

path_to_uncompressed_adapter_rrna_removed_fasta_files=$path_to_hdd'/uncompressed_adapter_rrna_removed_files/'
path_to_stats_folder=$path_to_hdd'/stats'
path_to_rrna_index=$path_to_hdd'/rRNA_index/rRNA_index'
path_to_compressed_adapter_rrna_removed_fasta_files=$path_to_hdd'/compressed_adapter_rrna_removed_files/'
##################################################################################################################
echo "Running on file: "$word1

echo "Downloading sra file using prefetch:"
time $path_to_prefetch $word1 

while [ ! -f $path_to_sra_folder$word1.sra ]; do
    echo "file not downloaded properly, doing it again!"
    rm $path_to_sra_folder*
    prefetch $word1
done


echo "Coverting to fastq format using fastqdump"
#time $path_to_fastqdump $path_to_sra_folder$word1.sra --outdir $path_to_adapter_removed_files_folder$word1.fastq |& tee $path_to_stats_folder/sra_fastq/fout1_$word1.txt
time $path_to_parallel_fastq_dump --threads 10 --sra-id $word1 --tmpdir $path_to_sra_folder --outdir $path_to_adapter_removed_files_folder$word1.fastq |& tee $path_to_stats_folder/sra_fastq/fout1_$word1.txt

echo "Deleting sra file:"
rm -f $path_to_sra_folder$word1.sra

echo "removing rrna: "
time bowtie2 -p 11 --very-sensitive-local -U $path_to_adapter_removed_files_folder$word1.fastq/$word1.fastq -x $path_to_rrna_index  --un=$path_to_uncompressed_adapter_rrna_removed_fasta_files$word1.fastq -S  $path_to_folder_for_sam_files_to_be_deleted/delete.sam |& tee $path_to_stats_folder/adapter_rrna_removed/fout3_$word1.txt

echo "deleting sam file:"
rm -f $path_to_folder_for_sam_files_to_be_deleted/delete.sam

echo "deleting adapterremoved file: "
rm -rf $path_to_adapter_removed_files_folder$word1.fastq

# updating the counter to keep track of files that will be compressed
count=$count+1

# updating the counter to keep track of files that are done
counter=$counter+1

echo "NUMBER OF FILES DONE                                                       "$counter
########################################################################################################################################
# enter if equal to 24 files
if [ $count -eq 12 ]
then

# defining pids array to be NULL
pids=()

echo "compressing these many output files at a time             "$count

# setting counter back to 0
count=0

cd $path_to_uncompressed_adapter_rrna_removed_fasta_files

# for each file in folder compressing
for i in *.fastq 
do 
7z a -t7z  $path_to_compressed_adapter_rrna_removed_fasta_files${i:0:-6}.7z $i &
pids+=($!)
echo "started compressing ${pids[*]} "
done

# wait till all done
for pid in ${pids[*]};
do
echo "waiting for $pid"
wait $pid
done

# delete done files
for i in *.fastq
do 
rm $i
echo "deleted $i"
done

cd ~
fi

########################################################################################################
done <./../../filenames/no_adapter_filenames.txt


