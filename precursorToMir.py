#to get microRNAs from a list of precursors
import sys
sys.path.append('/home/ogan/-PModules')
import ogbox as o

lines=o.readtxt('/home/ogan/Desktop/PORTAL/FANTOM5/newCAPPED')

mirBase=o.mirbaseGFF('/home/ogan/Desktop/PORTAL/hsa.gff3')

mir=[None]*len(mirBase.name)
mirStep=0
for i in lines:
    index=o.find(mirBase.name,i)[0]
    step=1
    while mirBase.mirType[index+step]=='miRNA':
        mir[mirStep]=mirBase.name[index+step]
        mirStep=mirStep+1
        step=step+1
    
mir=o.trimNone(mir)

o.writeLines(mir,'/home/ogan/Desktop/PORTAL/FANTOM5/precursors')
