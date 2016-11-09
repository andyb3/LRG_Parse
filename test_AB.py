'''
Created on 9 Nov 2016

@author: andy
'''
from lxml import etree

#help(etree.XMLParser)

context = etree.iterparse('LRG_1.xml')
#print help(context)
for action, elem in context:
    if elem.tag == 'exon':
        for child in elem:
            print elem.tag
            print elem.attrib
            print child.tag
            print child.attrib
        #print elem.tag
        #print elem.attrib
    #print elem.text
    #print elem.attrib


    
