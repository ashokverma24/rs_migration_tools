import sys, os
import numpy as np

dataPath = "./context-pds4/target/"

for path, subdirs, files in os.walk(dataPath):
    lid = {}
    for f in files:
        p = path+f
        titles =[]
        for p in open(p).readlines():
            if "<type>" in p: type = p.strip().split("type")[1][1:-2]
            if "logical_identifier" in p: ld = p.strip().split('logical_identifier')[1][1:-2]
            if "<title>" in p: titles.append(p.strip().split("title")[1][1:-2])
            if "<alternate_title>" in p: titles.append(p.strip().split("alternate_title")[1][1:-2])
        lid.update({ld: {'type': type, 'title' : titles}})


fmt = "{:<50} =  {:}\n"
lidFile = open("targLid.txt", "w")
typFile = open("targType.txt", "w")

for k in lid.keys():
    vals = lid[k]
    for t in vals['title']:
        if 'DEPRECATED' not in t.upper() and len(t.strip()):
           lidFile.write(fmt.format(t.title(), k))
           typFile.write(fmt.format(t.title(), vals['type']))
     
lidFile.close()
typFile.close()
