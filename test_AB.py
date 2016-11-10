'''
Created on 9 Nov 2016

@author: andy
'''
import xml.etree.ElementTree as etree

root = etree.parse('../LRG_Data/LRG_1.xml')

for exon in root.findall('./fixed_annotation/transcript/exon'):
    print exon.attrib['label']
    for coordinates in exon.findall('coordinates'):
        print coordinates.attrib 

#===============================================================================
# for child in root.iter('cdna'):
#     for sequence in child.iter('sequence'):
#         print sequence.text
#       
#   
# for child in root.iter('fixed_annotation'):
#     for exon in child.iter('exon'):
#         print exon.attrib['label']
#         for coordinate in exon.iter('coordinates'):
#             print coordinate.attrib 
#===============================================================================
