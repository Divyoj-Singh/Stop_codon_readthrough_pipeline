import matplotlib.pyplot as plt
import numpy as np
import pysam

path_to_files_to_run = "./../../filenames/files_to_run_downstream_13_sept/files_to_run_main.txt"
gene_id=["AT4G15510.3","AT1G26600.1"]
length_g=[1229,900]
for a , b in enumerate (gene_id):
        
	RT_profiles = [0]*length_g[a] 
	RT_profiles_start = [0]*length_g[a]
	with open (path_to_files_to_run) as file3:
		for line in file3:
			a = line[:-1].split("\t")
			srr_id = a[0]
			srr_list = a[1].split(",")
			million = int(a[2])/1000000
			length_range_start= int(a[3])
			length_range_stop= int(a[4])
			print("Running on the set: ",srr_id)

			for i in srr_list:
				print(i)
			
				BAM_filename = "./../../bam_files/"+i+".bam"
				bam_file = pysam.AlignmentFile(BAM_filename,"rb")
				for read in bam_file.fetch(until_eof = True):
					iid = read.reference_name
					if iid != b:
						continue
					pos = read.pos
					length = read.qlen
            		## making profile for rt genes
					RT_profiles_start[pos] += 1
					for i in range(pos,pos+length):
						RT_profiles[i] += 1
            
			plt.plot(np.array(RT_profiles))
			plt.show()
			plt.plot(np.array(RT_profiles_start))
			plt.show()
			with open(b+"_"+srr_id+"_profiles.txt","w") as f:
				for i in RT_profiles:
					f.write(str(i)+"\n")
            
			with open(b+"_"+srr_id+"_start_position.txt","w") as f:
				for i in RT_profiles_start:
					f.write(str(i)+"\n")