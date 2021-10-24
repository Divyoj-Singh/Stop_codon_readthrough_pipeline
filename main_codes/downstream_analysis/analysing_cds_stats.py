import collections
from collections import defaultdict
import os
import sys
import matplotlib.pyplot as plt
import array as arr


std_fold=1.5

path_to_cds_stats = "./../../downstream_analysis_files/cds_stats/"



path_to_cds_overlap_data = "./../../downstream_analysis_files/cds_overlap/"
if not os.path.exists(path_to_cds_overlap_data):
    os.makedirs(path_to_cds_overlap_data)

def function_to_find_overlap(start_1,end_1,start_2,end_2):
    overlap=0
    if start_1 >= start_2:
        if start_2 < end_1:
            overlap = 0
        elif end_2 > end_1:
            overlap = 1  
        elif end_2 <= end_1:
            overlap = (start_2-end_1)

    if start_1 < start_2:
        if start_1 < end_2:
            overlap =0
        elif end_1 > end_2:
            overlap = 1  
        elif end_1 <= end_2:
            overlap = (start_1-end_2)

    return overlap

for filename in os.listdir(path_to_cds_stats):
    if filename.endswith(".txt"): 
         srr_id=filename[:-4]
    else:
        continue    


    #with open (path_to_cds_overlap_data+srr_id+".txt","w") as file:
    #    pass

    with open (path_to_cds_stats+srr_id+".txt") as file:
        for line in file:
            start=arr.array('d',[0.0,0.0,0.0])
            end=arr.array('d',[0.0,0.0,0.0])

            a=line[:-1].split("\t")
            length=int(a[0])
            
            mean_0 = float(a[2])
            mean_1 = float(a[3])
            mean_2 = float(a[4])

            std_0 = float(a[5])
            std_1 = float(a[6])
            std_2 = float(a[7])

            start[0] =  mean_0 + std_fold*std_0   
            start[1] =  mean_1 + std_fold*std_1 
            start[2] =  mean_2 + std_fold*std_2 

            end[0] =  mean_0 - std_fold*std_0 
            end[1] =  mean_1 - std_fold*std_1 
            end[2] =  mean_2 - std_fold*std_2 
            
            
            overlap_0_1 = function_to_find_overlap(start[0],end[0],start[1],end[1])
            overlap_0_2 = function_to_find_overlap(start[0],end[0],start[2],end[2])
            overlap_1_2 = function_to_find_overlap(start[1],end[1],start[2],end[2])
            print (srr_id, length,overlap_0_1,overlap_0_2,overlap_1_2)