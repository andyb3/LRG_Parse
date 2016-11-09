import xml.etree.ElementTree as ET
import os

mydir = os.getcwd() + "/LRG_public_xml_files/"

tree = ET.parse(mydir + "LRG_1.xml")

root = tree.getroot()	

print(root[0][0].text)
print(root[0][1].text)
print(root[0][2].text)
print(root[0][3].text)
print(root[0][4][0].text)
print(root[0][4][1].text)
print(root[0][4][2][0].text)