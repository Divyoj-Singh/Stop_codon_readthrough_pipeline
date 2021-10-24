file_type = "FSM"
dic={}
with open ("annotation_file_"+file_type+"_latest.txt") as file:
	for line in file:
		a=line[:-1].split("\t")[0].split(".")[0]
		if a not in dic:
			dic[a]={"freq":0}
		dic[a]["freq"]+=1

for i in dic:
	print(i,"\t",dic[i]["freq"])
