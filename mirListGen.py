#!/usr/local/bin/python
#sys.argv list of command line arguments

#lists all human miRNAs
import sys
sys.path.append('/home/ogan/-PModules')
import ogbox as o
import rpy2.robjects as r


mirBase = o.mirbaseGFF('/home/ogan/Desktop/PORTAL/hsa.gff3')


mirName = [None] * len(mirBase.name)
mirPre = [None] * len(mirBase.name)
s=0
step=1
i=0
while i < len(mirBase.name):
    while mirBase.mirType[i+step] == 'miRNA':
        mirName[s] = mirBase.name[i+step]
        mirPre[s] = mirBase.name[i]
        step = step + 1
        if i + step == len(mirBase.mirType):
            break
        s=s+1
    i = i + step
    step = 1   
    if i + step == len(mirBase.mirType):
        break
mirName = o.trimNone(mirName)
mirPre = o.trimNone(mirPre)     
with open('/home/ogan/Desktop/PORTAL/mirName&Origin','w') as daFile:
    for n,p  in zip(mirName, mirPre):
        daFile.write(n+' '+p+'\n')
    
    
    
