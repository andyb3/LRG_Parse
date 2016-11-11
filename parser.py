import xml.etree.ElementTree as etree
import os
import csv
import sys


def readLRG(lrg):
    filepath = '../LRG_Data/' + lrg + '.xml'
    if os.path.isfile(filepath) == False:
         print("Please specify a valid LRG number")
         return
    else:
        tree = etree.parse(filepath)
        if tree.getroot().tag == 'lrg' and tree.getroot().attrib['schema_version'] == "1.9":
            return tree
        else:
            print("File does not conform to lrg version 1.9 specification.")
            return

def getBuildChrom(tree):
    for annotation_set in tree.findall('updatable_annotation/'):
        if annotation_set.attrib['type'] == 'lrg':
            for mapping in annotation_set.findall('mapping'):
                if mapping.attrib['type'] == 'main_assembly':
                    build = mapping.attrib['coord_system']
                    chrom = mapping.attrib['other_name']
                    other_start = int(mapping.find('mapping_span').attrib['other_start'])
                    other_end = int(mapping.find('mapping_span').attrib['other_end'])
                    strand = mapping.find('mapping_span').attrib['strand']
                    return build, chrom, other_start, other_end, strand

def getExons(tree, dnaSeq, build, chrom, other_start, other_end, strand):
    exonList = []
    for transcript in tree.findall('./fixed_annotation/transcript'):
        tx = transcript.attrib['name']
        for exon in transcript.findall('./exon'):
            exonnumber = exon.attrib['label']
            coordinates = exon.find('coordinates')
            coord_sys = coordinates.attrib['coord_system']
            start = coordinates.attrib['start']
            end = coordinates.attrib['end']
            if strand == "1":
                gDNA_start = (int(start) + int(other_start)) - 1
                gDNA_end = (int(end) + int(other_start)) - 1
            elif strand == "-1":
                gDNA_start = (int(other_end) - int(start)) + 1
                gDNA_end = (int(other_end) - int(end)) + 1
            seq = dnaSeq[int(start)-1:int(end)]
            tmp = [coord_sys, build, tx, exonnumber, chrom, gDNA_start, gDNA_end, seq]
            exonList.append(tmp)
    return exonList

def writeCSV(headers, exonList):
    exonList.insert(0, headers)
    outputFile = os.getcwd() + "/output.csv"
    # Write the output file
    with open(outputFile, 'wb') as f:
        writer = csv.writer(f)
        writer.writerows(exonList)

# for i in range(250):
#     lrg = 'LRG_' + str(i)
#     print(lrg)
#     tree = readLRG(lrg)
#     if tree:
#         dnaSeq = tree.find('./fixed_annotation/sequence').text
#         build, chrom, conv_offset = getBuildChrom(tree)
#         exonList = getExons(tree, dnaSeq, build, chrom, conv_offset)
#         headers = ['LRG_Number', 'Build', 'Transcript_ID', 'Exon_no', 'Chrom', 'Start', 'End', 'Sequence']
#         # writeCSV(headers, exonList)
#         for row in exonList:
#             print(row)


tree = readLRG(sys.argv[1])
if tree:
    dnaSeq = tree.find('./fixed_annotation/sequence').text
    build, chrom, other_start, other_end, strand = getBuildChrom(tree)
    exonList = getExons(tree, dnaSeq, build, chrom, other_start, other_end, strand)
    headers = ['LRG_Number', 'Build', 'Transcript_ID', 'Exon_no', 'Chrom', 'Start', 'End', 'Sequence']
    # writeCSV(headers, exonList)
    for row in exonList:
        print(row)




# dnaSeq = tree.find('./fixed_annotation/sequence').text

# mylist = []

# build, chrom, conv_offset = getBuildChrom()




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