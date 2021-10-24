# importing pyhton packages

import os
import sys
import math
import pysam
import collections
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
from collections import defaultdict



# global variables
readthrough_abs_reads = 0 # no of reads in isr combining all lengths 
min_reads_in_1_length_of_isr_of_rt_genes=0
ratio_of_den_isr_to_3utr=0

min_no_of_reads_in_cds_for_calculating_cds_data =0
coverage_threshold_isr = 0
coverage_threshold_3utr = 1
condition = "ISR"
fold_std = 50

## paths:
path_to_density = "./../../downstream_analysis_files_relaxed_copy/density_filtered_genes/"+condition+"/"
if not os.path.exists(path_to_density):
    os.makedirs(path_to_density)

path_to_coverage = "./../../downstream_analysis_files_relaxed_copy/coverage_filtered_genes/"+condition+"/"
if not os.path.exists(path_to_coverage):
    os.makedirs(path_to_coverage)

path_to_frame_info = "./../../downstream_analysis_files_relaxed_copy/periodicity_filtered_genes/"+condition+"/"
if not os.path.exists(path_to_frame_info):
    os.makedirs(path_to_frame_info)

path_to_histogram_data = "./../../downstream_analysis_files_relaxed_copy/length_dist_genes/"+condition+"/"
if not os.path.exists(path_to_histogram_data):
    os.makedirs(path_to_histogram_data)

path_to_final_RT_list = "./../../downstream_analysis_files_relaxed_copy/consolidated_filtered_genes/"+condition+"/"
if not os.path.exists(path_to_final_RT_list):
    os.makedirs(path_to_final_RT_list)

path_to_RT_gene_list = "./../../downstream_analysis_files_relaxed_copy/filtered_genes_frequencies/"+condition+"/"
if not os.path.exists(path_to_RT_gene_list):
    os.makedirs(path_to_RT_gene_list)

path_to_CDS_stats = "./../../downstream_analysis_files_relaxed_copy/cds_stats/"
if not os.path.exists(path_to_CDS_stats):
    os.makedirs(path_to_CDS_stats)

path_to_heatmaps = "./../../downstream_analysis_files_relaxed_copy/heatmaps/"+condition+"/"
if not os.path.exists(path_to_heatmaps):
	os.makedirs(path_to_heatmaps)
	os.makedirs(path_to_heatmaps+"/start/")
	os.makedirs(path_to_heatmaps+"/stop/")

path_to_bam_files = "./../../bam_files/"

path_of_RT_graphs = "./../../downstream_analysis_files_relaxed_copy/gene_profiles/"+condition+"/"
if not os.path.exists(path_of_RT_graphs):
	os.makedirs(path_of_RT_graphs)
	os.makedirs(path_of_RT_graphs+"/whole/")
	os.makedirs(path_of_RT_graphs+"/zoomed/")
	os.makedirs(path_of_RT_graphs+"/data/")

path_to_files_to_run = "./../../filenames/files_to_run_downstream/files_to_run_main.txt"
path_to_transcriptome_annotation = "./../../transcriptome/annotation_file.txt"


## functions:
def div(dividend,divisor):

    ## function to get rid od divided by zero error: 
    ## if divided by zero it will give an arbitrary value. 
	if divisor!=0:	
		q=float((dividend/divisor))
	if divisor==0:
		q=1000000000000000000000
		print('divsion by zero!!!')
	return q

def normalize(array):
    ## it just normalises the array u give it. 
	# ###  ## ###########################3######  adding a pseudo-count 1
	array[0] += 1
	array[1] += 1
	array[2] += 1
	summy = sum(array)
	if summy != 0:
		return np.divide(array,summy)

def function_to_initialise_dictionaries(length_range_start,length_range_stop):

	length_for_start = len(range(-24,63))
	length_for_stop = len(range(-62,1))

	start_region = {}
	stop_region = {}
    
	for i in range(length_range_start,length_range_stop+1):
		start_region[i] = [0] * length_for_start       # creating an array for storing heat map values for each length
		stop_region[i] = [0] * length_for_stop         # creating an array for storing heat map values for each length

	reads_in_frames = defaultdict(lambda:defaultdict(lambda:{"5utr":[0,0,0],"cds":[0,0,0],"isr":[0,0,0],"3utr":[0,0,0]}))
	len_hist = defaultdict(lambda:defaultdict(lambda:{"5utr":[0,0,0],"cds":[0,0,0],"isr":[0,0,0],"3utr":[0,0,0]}))
    # dictionary to store the reads in each gene and each length
	return start_region, stop_region, reads_in_frames,len_hist

def function_to_make_an_annotation_dictionary():

	transcriptome_annotation = {}

	with open (path_to_transcriptome_annotation) as file1: 
		for line in file1:
			array = line[:-1].split("\t")
			if len(array)==1:
				continue
			gene_id = array[0]
			start = int(array[1]) + 1
			stop1 = int(array[2]) + 1
			stop2 = int(array[3]) + 3
			length = int(array[4])
			len_utr5 = length - start
			len_cds = stop1 - start + 1
			len_isr = stop2 - stop1
			len_utr3 = length - stop2
			if len_utr3 < 45 or len_isr < 45:
				continue

			if gene_id not in transcriptome_annotation:
				transcriptome_annotation[gene_id] = [start, stop1, stop2, length, len_utr5, len_cds, len_isr, len_utr3]
			else:
				print("Something is wrong!")
				continue

	return transcriptome_annotation

def function_to_collect_frame_data(iid, pos, length, transcriptome_annotation, reads_in_frames):
	### function to collect frame data gene-length wise for a file
    ## 13 and 22 are window size for ribosome size  

	# CDS
	if pos-transcriptome_annotation[iid][0] >= -13 and pos-transcriptome_annotation[iid][1] <= -22:
		reads_in_frames[iid][length]["cds"][(transcriptome_annotation[iid][0]-pos+2) % 3] += 1

	# ISR
	elif pos-transcriptome_annotation[iid][1]-1 >= -13 and pos-transcriptome_annotation[iid][2] <= -22:
		reads_in_frames[iid][length]["isr"][(transcriptome_annotation[iid][0]-pos+2) % 3] += 1

	# 3UTR
	elif pos-transcriptome_annotation[iid][2]-1 >= -13:
		reads_in_frames[iid][length]["3utr"][(transcriptome_annotation[iid][0]-pos+2) % 3] += 1

	# 5UTR
	if pos-transcriptome_annotation[iid][0]+1 <= -22:
		reads_in_frames[iid][length]["5utr"][(transcriptome_annotation[iid][0]-pos+2) % 3] += 1

	return reads_in_frames

def function_to_print_length_distribution(len_hist,path_to_histogram_data,srr):
	
	with open(path_to_histogram_data+srr+".txt","w") as file:

		for gene in len_hist:
			for length in len_hist[gene]:
				file.write(gene+"\t"+str(length)+"\t")
				file.write(str(len_hist[gene][length]["5utr"][0])+"\t"+str(len_hist[gene][length]["5utr"][1])+"\t"+str(len_hist[gene][length]["5utr"][2])+"\t")
				file.write(str(len_hist[gene][length]["cds"][0])+"\t"+str(len_hist[gene][length]["cds"][1])+"\t"+str(len_hist[gene][length]["cds"][2])+"\t")
				file.write(str(len_hist[gene][length]["isr"][0])+"\t"+str(len_hist[gene][length]["isr"][1])+"\t"+str(len_hist[gene][length]["isr"][2])+"\t")
				file.write(str(len_hist[gene][length]["3utr"][0])+"\t"+str(len_hist[gene][length]["3utr"][1])+"\t"+str(len_hist[gene][length]["3utr"][2])+"\n")

		

def function_to_read_a_BAM_file(BAM_filename, transcriptome_annotation, start_region, stop_region, reads_in_frames, len_hist, length_range_start,length_range_stop):
    
	bam_file = pysam.AlignmentFile(BAM_filename,"rb") # using pysam for reading the bam file.

	for read in bam_file.fetch(until_eof = True): ## reading the bam file here

		iid = read.reference_name
		if iid not in transcriptome_annotation:
			continue
		pos = read.pos
		length = read.qlen

		len_hist = function_to_collect_frame_data(iid, pos, length, transcriptome_annotation, len_hist)

		cigar_criteria = len(read.cigarstring) == 3

		if length >= length_range_start and length <= length_range_stop and cigar_criteria: ### taking account of no mismatch rn 
			
			reads_in_frames = function_to_collect_frame_data(iid, pos, length, transcriptome_annotation, reads_in_frames)
            ## read the bam file and store the frame wise data in dictionary

            ## storing data in heatmap array,   62 is the window length
			
			if pos - transcriptome_annotation[iid][0] >= -24 and pos - transcriptome_annotation[iid][0] <= 62:
				start_region[length][pos - transcriptome_annotation[iid][0] + 24] += 1
			
			if pos - transcriptome_annotation[iid][1] >= -62 and pos - transcriptome_annotation[iid][1] <= 0:
				stop_region[length][pos - transcriptome_annotation[iid][1] + 62] += 1
    
    ### sorting the heat map array so according to length:
	start_region = collections.OrderedDict(sorted(start_region.items()))
	stop_region = collections.OrderedDict(sorted(stop_region.items()))

	return start_region, stop_region, reads_in_frames, len_hist

def function_to_plot_heat_map(dictionary,startorstop,range_provided,srr,length_range_start,length_range_stop):
    
	on_x = range_provided
	# setting the range on the x axis
    	
	on_y = range(length_range_start,length_range_stop+1)
	# setting the range on the y axis
    
	data_for_map = []
	# all data for the heat map will be stored here

	dictionary = collections.OrderedDict(sorted(dictionary.items()))

	for i in dictionary:
		# compiling all the data into one single array of arrays (2D array)
		#print(i)
		data_for_map.append(dictionary[i])
    
	data_for_map = np.array(data_for_map)
	# making the 2D array a numpy array so that it can be used for plotting
    
	fig, ax = plt.subplots()
	# initializing a subplot
    
	im = ax.imshow(data_for_map,cmap="magma")
	# putting in data into the heat map
    
	ax.set_xticks(np.arange(len(on_x)))
	# marking the X ticks
    
	ax.set_yticks(np.arange(len(on_y)))
	# marking the Y ticks
    
	ax.set_xticklabels(on_x,fontsize=5)
	# labelling the X ticks
    
	ax.set_yticklabels(on_y,fontsize=5)
	# labelling the Y ticks

	#for i in range(len(on_y)):# for loop for writing the corresponding values into the heat map
		#for j in range(len(on_x)):
			#text = ax.text(j, i, data_for_map[i,j],ha="center", va="center", color="w")
    	
	plt.setp(ax.get_xticklabels(), rotation=90, ha="right",rotation_mode="anchor")
	# setting the format for labelling the x ticks
    
	ax.set_title("Heat map for periodicity "+startorstop+": "+srr)
	# setting the title of the plot
    
	fig.tight_layout()
	# setting a tight layout for the figure

	#plt.show()
	plt.savefig(path_to_heatmaps+startorstop+"/"+srr+".png")
	plt.close()
	pass

def function_to_check_robustness_of_a_ribo_seq_file(BAM_filename,iid,srr,million,transcriptome_annotation,heatmap_start_region, heatmap_stop_region, reads_in_frames,len_hist,length_range_start,length_range_stop):
	## i think the name should be changed, but i dont know to what.
    ## anyways the function just take inputs from previous step and passes it to other function, and yah also prints theat the bam file is read.
	start_region, stop_region, reads_in_frames, len_hist = function_to_read_a_BAM_file(BAM_filename, transcriptome_annotation, heatmap_start_region, heatmap_stop_region, reads_in_frames, len_hist, length_range_start,length_range_stop)
	print("Done reading the BAM file!")

	return start_region, stop_region, reads_in_frames,len_hist

def function_to_get_RT_candidates_first_list(reads_in_frames,iid,transcriptome_annotation,million,length_range_start,length_range_stop):

	with open (path_to_density+iid+".txt","w") as file2:
		pass
	for gene in reads_in_frames:
		#if transcriptome_annotation[gene][7] == 0: ## if 3_UTR is 0 in length, then leave that gene.
		#	continue
		reads_in_cds = 0
		density_in_cds = 0
		reads_in_isr = 0
		density_in_isr = 0
		reads_in_utr3 = 0
		density_in_utr3 = 0
		# calculating RPKM density: reads per kilo base per million reads.
		for i in range(length_range_start,length_range_stop+1):
			reads_in_cds += sum(reads_in_frames[gene][i]["cds"])            
			reads_in_isr += sum(reads_in_frames[gene][i]["isr"])
			reads_in_utr3 += sum(reads_in_frames[gene][i]["3utr"]) 


		density_in_cds = div((reads_in_cds*1000),(transcriptome_annotation[gene][5]*million))
		density_in_isr = div((reads_in_isr*1000),(transcriptome_annotation[gene][6]*million))
		density_in_utr3 = div((reads_in_utr3*1000),(transcriptome_annotation[gene][7]*million))

        ## putting criteria for density comparison and no of reads in isr.
		##  REMOVING THE IF CONDITION TO PRINT FOR ALL GENES:
		#if (density_in_cds > density_in_isr) and (density_in_isr >= ratio_of_den_isr_to_3utr*density_in_utr3) and (reads_in_isr > readthrough_abs_reads):
			#print(filename+"\t"+gene+"\t"+str(reads_in_cds)+"\t"+str(reads_in_isr)+"\t"+str(reads_in_utr3)+"\t"+str(round(density_in_cds,4))+"\t"+str(round(density_in_isr,4))+"\t"+str(round(density_in_utr3,4))+"\n")
			# column names: gene_name   transcript_legth    cds_length  isr_length  3utr_length reads_in_cds    reads_in_isr    reads_in_3utr   density_in_cds  density_in_isr  density_in_3utr
		
        ## for gene which pass through this criteria are written into a text file along with their annotation.
		with open (path_to_density+iid+".txt","a+") as file2:
			file2.write(gene+"\t"+str(transcriptome_annotation[gene][3])+"\t"+str(transcriptome_annotation[gene][5])+"\t"+str(transcriptome_annotation[gene][6])+"\t"+str(transcriptome_annotation[gene][7])+"\t"+str(reads_in_cds)+"\t"+str(reads_in_isr)+"\t"+str(reads_in_utr3)+"\t"+str(round(density_in_cds,4))+"\t"+str(round(density_in_isr,4))+"\t"+str(round(density_in_utr3,4))+"\n")

def function_to_plot_bar_plot(reads_in_frames,iid,length_range_start,length_range_stop):

	with open (path_to_CDS_stats+iid+".txt","w") as file1:
		pass

	data_bar_plot = defaultdict(lambda:[[],[],[]])
	
	for i in range(length_range_start,length_range_stop+1):
		for gene in reads_in_frames:
			if sum(reads_in_frames[gene][i]["cds"])>=min_no_of_reads_in_cds_for_calculating_cds_data:
				# if no of reads in cds is more than a threshold, accept it for cds stats.

				data_bar_plot[i][0].append(div(reads_in_frames[gene][i]["cds"][0],sum(reads_in_frames[gene][i]["cds"])))
				data_bar_plot[i][1].append(div(reads_in_frames[gene][i]["cds"][1],sum(reads_in_frames[gene][i]["cds"])))
				data_bar_plot[i][2].append(div(reads_in_frames[gene][i]["cds"][2],sum(reads_in_frames[gene][i]["cds"])))
                #appending the list for one length.
	## writing the mean and standard deviation for each length for each frame	
		with open (path_to_CDS_stats+iid+".txt","a+") as file1:
			file1.write(str(i)+"\t"+str(len(data_bar_plot[i][0]))+"\t"+str(round(np.mean(data_bar_plot[i][0]),5))+"\t"+str(round(np.mean(data_bar_plot[i][1]),5))+"\t"+str(round(np.mean(data_bar_plot[i][2]),5))+"\t"+str(round(np.std(data_bar_plot[i][0]),5))+"\t"+str(round(np.std(data_bar_plot[i][1]),5))+"\t"+str(round(np.std(data_bar_plot[i][2]),5))+"\n")

def best_among_the_first_genes(srr):
    # for a set calculating the no of times  gene has occured in readthrough candidadates.   
    d = {}
    with open(path_to_density+srr+".txt") as file:
        for line in file:
            a = line[:-1].split("\t")
            gene = a[0]
            if gene not in d:
                d[gene] = 1
            else:
                d[gene] += 1
    return d

def function_to_give_potential_readthrough_list(d):
	return d.keys() 
    ## just returns the list of candidate genes on the basis of density.

def initialise_RT_dictionary(potn_readthrough_list,transcriptome_annotation):
	RT_profiles = {} ## dictionary to store gene profiles.
	RT_frame_info = {} ## frame info dictionary
	
	for gene in potn_readthrough_list:
		RT_profiles[gene] = [0]*transcriptome_annotation[gene][3]
		RT_frame_info[gene] = defaultdict(lambda:{"5utr":[0,0,0],"cds":[0,0,0],"isr":[0,0,0],"3utr":[0,0,0]})
	
	return RT_profiles,RT_frame_info

def function_to_read_BAM_file_again(BAM_filename, transcriptome_annotation, RT_profiles,RT_frame_info,length_range_start,length_range_stop):
	bam_file = pysam.AlignmentFile(BAM_filename,"rb")
	for read in bam_file.fetch(until_eof = True):

		iid = read.reference_name
		if iid not in RT_profiles or iid not in transcriptome_annotation:
			continue
		pos = read.pos
		length = read.qlen
		#sequence = read.seq
		cigar_criteria = len(read.cigarstring) == 3

		if length >= length_range_start and length <= length_range_stop and cigar_criteria: ### taking account of no mismatch rn 
			## collecting frame data for the rt genes
			RT_frame_info = function_to_collect_frame_data(iid, pos, length, transcriptome_annotation, RT_frame_info)
			## making profile for rt genes 
			for i in range(pos,pos+length):
				RT_profiles[iid][i] += 1

	return RT_profiles,RT_frame_info

def function_to_calculate_coverage(RT_profiles,transcriptome_annotation,iid):
	coverage = {}
	with open(path_to_coverage+iid+".txt","w") as file:
		pass
    ## calculating coverage from the profile 
	for gene in RT_profiles:
		cov_utr5 = np.count_nonzero(RT_profiles[gene][:transcriptome_annotation[gene][0]-1])/transcriptome_annotation[gene][4]
		cov_cds = np.count_nonzero(RT_profiles[gene][transcriptome_annotation[gene][0]-1:transcriptome_annotation[gene][1]])/transcriptome_annotation[gene][5]
		cov_isr = np.count_nonzero(RT_profiles[gene][transcriptome_annotation[gene][1]:transcriptome_annotation[gene][2]])/transcriptome_annotation[gene][6]
		cov_utr3 = np.count_nonzero(RT_profiles[gene][transcriptome_annotation[gene][2]:])/transcriptome_annotation[gene][7]
		coverage[gene] = [cov_utr5,cov_cds,cov_isr,cov_utr3]
		if cov_isr <= coverage_threshold_isr:
			continue
		if cov_utr3 >= coverage_threshold_3utr:
			continue
		# do not calculate coverage for those genes that do not have at least one read spanning the entire canonical stop codon
		## REMOVING THIS CRITERA AS WELL FOR PRINTING COV AND DENSITY OF ALL GENES.
		#if (RT_profiles[gene][int(transcriptome_annotation[gene][1])-3] == 0) or (RT_profiles[gene][transcriptome_annotation[gene][1]-2] == 0) or (RT_profiles[gene][transcriptome_annotation[gene][1]-1] == 0):
		#	continue
		with open(path_to_coverage+iid+".txt","a+") as file:
			file.write(gene+"\t"+str(cov_utr5)+"\t"+str(cov_cds)+"\t"+str(cov_isr)+"\t"+str(cov_utr3)+"\n")
	
	return coverage

def function_to_report_readthrough_genes(RT_frame_info,CDS_stats_dict,iid,good_coverage_genes):
    ## this function prints the genes
	with open(path_to_frame_info+iid+".txt","w") as file:
		pass
	for gene in RT_frame_info:
		if gene not in good_coverage_genes:
			continue
		for i in CDS_stats_dict:
			if "nan" in CDS_stats_dict[i][1:4]:
				continue
			## putting threshold on min no of reads in one length of rt genes
			if sum(RT_frame_info[gene][i]["isr"]) < min_reads_in_1_length_of_isr_of_rt_genes:
				continue
            ## checking if cds and then isr are within the std dev for 3 nucleotide periodicity..
            ## first checking for cds
			if np.array(CDS_stats_dict[i][1:4]).astype(float)[0]-fold_std*np.array(CDS_stats_dict[i][4:]).astype(float)[0]<=normalize(RT_frame_info[gene][i]["cds"])[0]<=np.array(CDS_stats_dict[i][1:4]).astype(float)[0]+fold_std*np.array(CDS_stats_dict[i][4:]).astype(float)[0]:
				if np.array(CDS_stats_dict[i][1:4]).astype(float)[1]-fold_std*np.array(CDS_stats_dict[i][4:]).astype(float)[1]<=normalize(RT_frame_info[gene][i]["cds"])[1]<=np.array(CDS_stats_dict[i][1:4]).astype(float)[1]+fold_std*np.array(CDS_stats_dict[i][4:]).astype(float)[1]:
					if np.array(CDS_stats_dict[i][1:4]).astype(float)[2]-fold_std*np.array(CDS_stats_dict[i][4:]).astype(float)[2]<=normalize(RT_frame_info[gene][i]["cds"])[2]<=np.array(CDS_stats_dict[i][1:4]).astype(float)[2]+fold_std*np.array(CDS_stats_dict[i][4:]).astype(float)[2]:
                        ## now checking for isr    
						if np.array(CDS_stats_dict[i][1:4]).astype(float)[0]-fold_std*np.array(CDS_stats_dict[i][4:]).astype(float)[0]<=normalize(RT_frame_info[gene][i]["isr"])[0]<=np.array(CDS_stats_dict[i][1:4]).astype(float)[0]+fold_std*np.array(CDS_stats_dict[i][4:]).astype(float)[0]:
							if np.array(CDS_stats_dict[i][1:4]).astype(float)[1]-fold_std*np.array(CDS_stats_dict[i][4:]).astype(float)[1]<=normalize(RT_frame_info[gene][i]["isr"])[1]<=np.array(CDS_stats_dict[i][1:4]).astype(float)[1]+fold_std*np.array(CDS_stats_dict[i][4:]).astype(float)[1]:
								if np.array(CDS_stats_dict[i][1:4]).astype(float)[2]-fold_std*np.array(CDS_stats_dict[i][4:]).astype(float)[2]<=normalize(RT_frame_info[gene][i]["isr"])[2]<=np.array(CDS_stats_dict[i][1:4]).astype(float)[2]+fold_std*np.array(CDS_stats_dict[i][4:]).astype(float)[2]:
									with open(path_to_frame_info+iid+".txt","a+") as file:
										file.write(gene+"\t"+str(i)+"\t"+str(np.array(CDS_stats_dict[i][1:4]).astype(float)[0])+"\t"+str(np.array(CDS_stats_dict[i][1:4]).astype(float)[1])+"\t"+str(np.array(CDS_stats_dict[i][1:4]).astype(float)[2])+"\t"+str(normalize(RT_frame_info[gene][i]["cds"])[0])+"\t"+str(normalize(RT_frame_info[gene][i]["cds"])[1])+"\t"+str(normalize(RT_frame_info[gene][i]["cds"])[2])+"\t"+str(normalize(RT_frame_info[gene][i]["isr"])[0])+"\t"+str(normalize(RT_frame_info[gene][i]["isr"])[1])+"\t"+str(normalize(RT_frame_info[gene][i]["isr"])[2])+"\t"+str(normalize(RT_frame_info[gene][i]["3utr"])[0])+"\t"+str(normalize(RT_frame_info[gene][i]["3utr"])[1])+"\t"+str(normalize(RT_frame_info[gene][i]["3utr"])[2])+"\n")

def function_to_read_CDS_stats(iid):
	CDS_stats_dict = {}
	with open (path_to_CDS_stats+iid+".txt") as file:
			for line in file:
				a = line[:-1].split("\t")
				CDS_stats_dict[int(a[0])] = a[1:]
	return CDS_stats_dict

def reading_bam_files_again(BAM_filename,iid,srr,million,RT_profiles, RT_frame_info,length_range_start,length_range_stop):

	RT_profiles,RT_frame_info = function_to_read_BAM_file_again(BAM_filename, transcriptome_annotation, RT_profiles,RT_frame_info,length_range_start,length_range_stop)
	print("Done reading the bam file again!")
	
	return RT_profiles,RT_frame_info


with open (path_to_files_to_run) as file3:
	for line in file3:

		
		a = line[:-1].split("\t")
		iid = a[0]
		srr_list = a[1].split(",")
		million = int(a[2])/1000000
		length_range_start= int(a[3])
		length_range_stop= int(a[4])
		print("Running on the set: ",iid)

		
		transcriptome_annotation = function_to_make_an_annotation_dictionary()
		print("Done reading the transcriptome!")
	
		heatmap_start_region, heatmap_stop_region, reads_in_frames,len_hist = function_to_initialise_dictionaries(length_range_start,length_range_stop)
		print("Done initialising dictionaries!")
		
		for i in srr_list:
			print(i)
			start_region, stop_region, reads_in_frames,len_hist = function_to_check_robustness_of_a_ribo_seq_file(path_to_bam_files+i+".bam",iid,i,million,transcriptome_annotation,heatmap_start_region, heatmap_stop_region, reads_in_frames,len_hist,length_range_start,length_range_stop)

		function_to_print_length_distribution(len_hist,path_to_histogram_data,iid)

		function_to_plot_heat_map(start_region,"start",range(-24,63),iid,length_range_start,length_range_stop)
		function_to_plot_heat_map(stop_region,"stop",range(-62,1),iid,length_range_start,length_range_stop)
		print("Done plotting both the heat maps!")

		function_to_get_RT_candidates_first_list(reads_in_frames,iid,transcriptome_annotation,million,length_range_start,length_range_stop)
		print("Done printing the initial list of RT candidates!")

		function_to_plot_bar_plot(reads_in_frames,iid,length_range_start,length_range_stop)
		print("Done creating the CDS stats for the files!")
		
		d = best_among_the_first_genes(iid)

		potn_readthrough_list = function_to_give_potential_readthrough_list(d)
		RT_profiles, RT_frame_info = initialise_RT_dictionary(potn_readthrough_list,transcriptome_annotation)
		CDS_stats_dict = function_to_read_CDS_stats(iid)
		for i in srr_list:
			print(i)
			RT_profiles,RT_frame_info = reading_bam_files_again(path_to_bam_files+i+".bam",iid,i,million,RT_profiles, RT_frame_info,length_range_start,length_range_stop)
		
		coverage = function_to_calculate_coverage(RT_profiles,transcriptome_annotation,iid)
		print("Done calculating Coverage values!")


