import collections
from collections import defaultdict
import os
import sys
import matplotlib.pyplot as plt


condition="ISR"

path_to_length_dist = "./../../downstream_analysis_files_relaxed/length_dist_genes/"+condition+"/"



path_to_length_histogram_data = "./../../downstream_analysis_files_relaxed/length_histogram/"+condition+"/"
if not os.path.exists(path_to_length_histogram_data):
    os.makedirs(path_to_length_histogram_data)


for filename in os.listdir(path_to_length_dist):
    if filename.endswith(".txt"): 
         srr_id=filename[:-4]

    else:
        continue    

    orf_utr5=0
    orf_utr3=0

    if not os.path.exists(path_to_length_histogram_data+"data/"):
        os.makedirs(path_to_length_histogram_data+"data/")
        os.makedirs(path_to_length_histogram_data+"combined/")
        os.makedirs(path_to_length_histogram_data+"isr/")
        os.makedirs(path_to_length_histogram_data+"3utr/")
        os.makedirs(path_to_length_histogram_data+"cds/")
        os.makedirs(path_to_length_histogram_data+"5utr/")
        

    with open (path_to_length_histogram_data+"data/"+srr_id+".txt","w") as file:
        pass

    with open (path_to_length_dist+srr_id+".txt") as file:
        length_dic= defaultdict(lambda:{"combined":0,"5utr":0,"cds":0,"isr":0,"3utr":0})

        for line in file:
            a=line[:-1].split("\t")
            length=int(a[1])
            gene=str(a[0])
            utr5 = int(a[2])+int(a[3])+ int(a[4])  
            cds  = int(a[5])+int(a[6])+ int(a[7])
            isr  = int(a[8])+int(a[9])+ int(a[10])
            utr3 = int(a[11])+int(a[12])+ int(a[13])

            if utr5 > cds:
                orf_utr5+=1
                #print(gene)
                continue     
            if utr3 > cds:
                orf_utr3+=1
                #print(gene,"3utr")
                continue
        
            length_dic[length]["combined"]= length_dic[length]["combined"]+utr5+cds+isr+utr3
            length_dic[length]["5utr"]+=utr5
            length_dic[length]["cds"]+=cds
            length_dic[length]["isr"]+=isr
            length_dic[length]["3utr"]+=utr3

        
        length_dic = collections.OrderedDict(sorted(length_dic.items()))
        with open(path_to_length_histogram_data+"data/"+srr_id+".txt","a+") as file1:
            for i in length_dic:
                file1.write(str(i)+"\t"+str(length_dic[i]["combined"])+"\t"+str(length_dic[i]["5utr"])+"\t"+str(length_dic[i]["cds"])+"\t"+str(length_dic[i]["isr"])+"\t"+str(length_dic[i]["3utr"])+"\n")        

        x_axis=[]
        combined=[]
        utr5=[]
        cds=[]
        isr=[]
        utr3=[]
        
        for i in length_dic:
            x_axis.append(i)
            combined.append(length_dic[i]['combined'])
            utr5.append(length_dic[i]['5utr'])
            utr3.append(length_dic[i]['3utr'])
            cds.append(length_dic[i]['cds'])
            isr.append(length_dic[i]['isr'])
        
        plt.bar(x_axis,combined)
        plt.title("Combined reads for "+srr_id)
        plt.xlabel("Length of fragment")
        plt.ylabel("Number of reads")
        plt.savefig(path_to_length_histogram_data+"combined/"+"combined_"+srr_id+".png",dpi = 1000)
        plt.close()
        
        plt.bar(x_axis,utr5)
        plt.title("5'UTR reads for "+srr_id)
        plt.xlabel("Length of fragment")
        plt.ylabel("Number of reads")
        plt.savefig(path_to_length_histogram_data+"5utr/"+"5utr_"+srr_id+".png",dpi = 1000)
        plt.close()

        plt.bar(x_axis,cds)
        plt.title("CDS reads for "+srr_id)
        plt.xlabel("Length of fragment")
        plt.ylabel("Number of reads")
        plt.savefig(path_to_length_histogram_data+"cds/"+"cds_"+srr_id+".png",dpi = 1000)
        plt.close()

        plt.bar(x_axis,isr)
        plt.title("ISR reads for "+srr_id)
        plt.xlabel("Length of fragment")
        plt.ylabel("Number of reads")
        plt.savefig(path_to_length_histogram_data+"isr/"+"isr_"+srr_id+".png",dpi = 1000)
        plt.close()

        plt.bar(x_axis,utr3)
        plt.title("3'UTR reads for "+srr_id)
        plt.xlabel("Length of fragment")
        plt.ylabel("Number of reads")
        plt.savefig(path_to_length_histogram_data+"3utr/"+"3utr"+srr_id+".png",dpi = 1000)
        plt.close()

		
        #print("Number of orf in 5utr in file ", srr_id, "is", orf_utr5)
        #print("Number of orf in 3utr in file ", srr_id, "is", orf_utr3)
   