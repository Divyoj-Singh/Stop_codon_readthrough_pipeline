import os
files = []
#following line prints the first line for the file i.e. the names for the stats
print('filename'+'\t'+'total#reads'+'\t'+'#aligned0times'+'\t'+'%aligned0times'+'\t'+'#aligned1times'+'\t'+'%aligned1times'+'\t'+'# aligned>1times'+'\t'+'%aligned>1times'+'\t'+'overallalignment')
#	
for i in os.listdir('./'):
	if i.endswith('.txt'):
		with open (i) as file:
			linenumber=0
			filename=i[6:-4]
			for line in file:
				a=line.split(' ')
				linenumber +=1

				if linenumber ==1:
					total_read=a[0]  
				if linenumber ==3:
					zero_times_read=a[4]
					zero_times_readptg=a[5]	  
				if linenumber ==4:
					one_times_read=a[4]
					one_times_readptg=a[5] 
				if linenumber ==5:
					more1times_read=a[4]
					more1times_readptg=a[5] 
				if linenumber ==6:
					overall=a[0]  
			
			print (filename +'\t'+total_read+'\t'+zero_times_read+'\t'+zero_times_readptg+'\t'+one_times_read+'\t'+one_times_readptg+'\t'+more1times_read+'\t'+more1times_readptg+'\t'+overall)

# print the output using >./../align_stats_ and give a name in command line.
