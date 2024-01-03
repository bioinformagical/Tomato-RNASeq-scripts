#!/usr/bin/python

import re
import sys
#import die

#########################       add2MasterTable.py     ----Adding to a big table, another column of EBSEQ-HMM results        #########################################
##########################################################################################################
######
######		This will read in a master file  and an EBSEQ result file.  
	# 	We are going to:
	#		1.  Read the EBSEQ file into a HASH, key = geneid, value = HMM pattern value
	#		1.5 Read data file in to get the description file and raw data too!!!
	#		2.  Iterate through the master file. Check if geneID exists, add a new column to the line with the hash value in af hash.
	#		3.  Otherwise add a "0" to denote no result.
	#		4.  
	#		3. Write results to table.
	#	USAGE:		 python ./add2MasterTable.py ./masterTable-ebseqhmm.txt ./a34gc_ebseqGeneCalls-SL4.txt ./muday-144-SL4_counts-salmon.txt
Prefix = ""
print(sys.argv[1])
print("The number of sys.argv given were; ",len(sys.argv),"\n")
#if (len(sys.argv) > 4 or len(sys.argv) < 4 ):
#	sys.exit("Exiting due to wrong # of arguments !!! ")
	#die.Die("Must specify 1 ARG. LaST ARG WAS ",str(sys.argv)," ")
	#Prefix = sys.argv[2]

###  making some OUTPUT FILES.
print(sys.argv[2])
tmp=sys.argv[2]
tmparray=tmp.split("gc_")
front=tmparray[0].split("/")[1]
print(front)
### OUtfile
outname = "masterTable-ebseqhmm.new"
#sys.exit()
outfile = open(outname, 'w')


###  A dict of A data / desc file 
deschash={}
with open(sys.argv[3]) as afile:
	for line in afile:
		line =  line.strip()
		if not line:
			continue
		pattern=re.compile("^gene")
		if pattern.match(line):
			continue
		#elif "exon" not in line:
		#	continue
		else:
			linearray = line.split()
			#print(line)
			geneid=linearray[0]
			desc=linearray[73]
			lengtharray=len(linearray)
			count=73
			while (count < lengtharray):
				desc = desc+" "+linearray[count]
				count = count + 1
			print(geneid)
			deschash[geneid]=desc
			
#sys.exit()

###  A dict of the EBSEQ result file
ebseqhash={}
with open(sys.argv[2]) as afile:
        for line in afile:
                line =  line.strip()
                if not line:
                        continue
                pattern=re.compile("\"Most")
                if pattern.match(line):
                        continue
                #elif "exon" not in line:
                #       continue
                else:
                        linearray = line.split()
                        #print(line)
                        geneid=linearray[0].replace('"', '')
                        pattern=linearray[1]
                        print(geneid)
                        ebseqhash[geneid]=pattern

#sys.exit()

### Read in the big table and check for each coord in gene range exists in the vcfhash. COunt # of snps and report that. Add as a column.
with open(sys.argv[1]) as tablefile:
	for line in tablefile:
		line =  line.strip()
		if not line:
			continue
		pattern=re.compile("^Gene_ID")
		if pattern.match(line):
			header=line+"\t"+front
			print(header)
			outfile.write(header)
			outfile.write("\n")
			continue
		#elif "Gene_id" in line:
		#	header=line+"\t"+front
		#	print(header)
		#	outfile.write(header)
		#	outfile.write("\n")
		else:
			outfile.write(line+"\t")
			linearray = line.split()
			geneid=linearray[0].replace('"', '')
			
			###  This section only gets used on the first time though to add a description to the master file. Not needed after.
			#if geneid in deschash:
			#	outfile.write(deschash[geneid])
			#	outfile.write("\t")
			#else:
			#	sys.exit("We failed to find a description entry in hash.")
			
			if geneid in ebseqhash:
				outfile.write(ebseqhash[geneid])
			else:
				outfile.write("-")

			outfile.write("\n")
				
#sys.exit()

print("The Final is done ")
