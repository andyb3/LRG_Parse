import xml.etree.ElementTree as etree
import os
import csv
import sys

lrg = sys.argv[1]

filepath = '../LRG_Data/' + lrg + '.xml'

if os.path.isfile(filepath) != True:
	print("Please specify a valid LRG number")
else:
	# run the function


# root = etree.parse(filepath)

# gDNAseq = root.find('./fixed_annotation/sequence').text

# mylist = []

# for transcript in root.findall('./fixed_annotation/transcript'):
# 	tx = transcript.attrib['name']
# 	for exon in transcript.findall('./exon'):
# 		exonnumber = exon.attrib['label']
# 		coordinates = exon.find('coordinates')
# 		coord_sys = coordinates.attrib['coord_system']
# 		start = coordinates.attrib['start']
# 		end = coordinates.attrib['end']
# 		seq = gDNAseq[int(start)-1:int(end)]
# 		tmp = [tx, exonnumber, coord_sys, start, end, seq]
# 		mylist.append(tmp)

# for row in mylist:
# 	print(row)

# headers = ['Transcript_ID', 'Exon_no', 'Coordinate_system', 'Start', 'End', 'Sequence']

# mylist.insert(0, headers)

# myoutput = os.getcwd() + "/output.csv"

# with open(myoutput, 'wb') as f:
# 	writer = csv.writer(f)
# 	writer.writerows(mylist)

# print("Done")