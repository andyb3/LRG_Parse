import xml.etree.ElementTree as etree
import os
import csv
import sys

mydict = {'foo': 1, 'bar': 2}

myvar = 'foo'

if myvar in mydict:
	print("YAS")




# usergene = sys.argv[1]

# print('Usergene is %s' % usergene)

# LRG_dir = '../LRG_Data/'

# code to take user's gene and find the corresponding LRG
# def findLRG(usergene, LRG_dir):
	
# 	for files in os.listdir(LRG_dir):
# 		root = etree.parse(LRG_dir + files)
# 		for annotation_set in root.findall('updatable_annotation/'):
# 			if annotation_set.attrib['type'] == 'lrg':
# 				genename = annotation_set.find('lrg_locus').text
# 				if genename == usergene:
# 					return(files)

# print(findLRG(usergene, LRG_dir))




# if os.path.isfile(filepath) != True:
# 	print("Please specify a valid LRG number")
# else:
	# run the function

# put the stuff below this line into a function

# root = etree.parse(filepath)

# gDNAseq = root.find('./fixed_annotation/sequence').text

# mylist = []


# for annotation_set in root.findall('updatable_annotation/'):
# 	if annotation_set.attrib['type'] == 'lrg':
# 		for mapping in annotation_set.findall('mapping'):
# 			if mapping.attrib['type'] == 'main_assembly':
# 				build = mapping.attrib['coord_system']
# 				chrom = mapping.attrib['other_name']
# 				other_start = int(mapping.find('mapping_span').attrib['other_start'])
# 				conv_offset = other_start - 1

# for transcript in root.findall('./fixed_annotation/transcript'):
# 	tx = transcript.attrib['name']
# 	for exon in transcript.findall('./exon'):
# 		exonnumber = exon.attrib['label']
# 		coordinates = exon.find('coordinates')
# 		coord_sys = coordinates.attrib['coord_system']
# 		start = coordinates.attrib['start']
# 		end = coordinates.attrib['end']
# 		seq = gDNAseq[int(start)-1:int(end)]
# 		tmp = [build, tx, exonnumber, coord_sys, chrom, str(int(start)+conv_offset), str(int(end)+conv_offset), seq]
# 		mylist.append(tmp)

# for row in mylist:
# 	print(row)

# Create the headers
# headers = ['Transcript_ID', 'Exon_no', 'Coordinate_system', 'Start', 'End', 'Sequence']

# Add headers to the output
# mylist.insert(0, headers)

# Define the output file
# myoutput = os.getcwd() + "/output.csv"

# Write the output file
# with open(myoutput, 'wb') as f:
# 	writer = csv.writer(f)
# 	writer.writerows(mylist)

# print("Done")