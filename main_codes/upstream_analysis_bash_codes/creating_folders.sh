#####Change the species name for other species####

#declaring paths to hard disks
## path to ssd drive 
path_to_ssd='/media/csb/D/Plant_Readthrough/'
## path to external drive
path_to_hdd='/media/csb/D/Plant_Readthrough/'

mkdir $path_to_hdd
mkdir $path_to_ssd

# Now creating separate folders
cd $path_to_hdd
mkdir -v bam_files compressed_adapter_rrna_removed_files stats 
ls 

cd $path_to_ssd
mkdir -v adapter_removed_files uncompressed_adapter_rrna_removed_files fastq_files transc_align_sam_files
ls

cd $path_to_hdd
cd ./stats
mkdir -v adapter_rrna_removed fastq_adaptrem sra_fastq transc_align