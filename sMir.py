#!/raid2/local/python2.7.3/bin/python2.7
#$ -N
from pybedtools import BedTool
import os
import time
import sys

fileno=sys.argv[1]

def cuffToBed(address,output):
    f=open(address)
    #f=open('CufflinksOut/transcripts.gtf')
    cuffToBed=open(output,'w')
    while True:
        line=f.readline()
        if line=='':
            break
        lineSplit=line.split('\t')
        features=lineSplit[8].split(';')
        if lineSplit[2]=='transcript':
            print >>cuffToBed, lineSplit[0]+'\t'+str(int(lineSplit[3])-1)+'\t'+\
                  lineSplit[4]+'\t'+features[1][16:len(features[1])-1]\
                  +'\t'+lineSplit[5]+'\t'+lineSplit[6]
    f.close()

#def cuffToBed2(address,output):
#    f=open(address)
#    #f=open('CufflinksOut/transcripts.gtf')
#    cuffToBed=open(output,'w')
#    while True:
#        line=f.readline()
#        if line=='':
#            break
#        lineSplit=line.split('\t')
#        features=lineSplit[8].split(';')
#        if lineSplit[2]=='exon':
#            print >>cuffToBed, lineSplit[0]+'\t'+str(int(lineSplit[3])-1)+'\t'+\
#                  lineSplit[4]+'\t'+features[1][16:len(features[1])-1]\
#                  +'\t'+lineSplit[5]+'\t'+lineSplit[6]
#    f.close()


miRNA=BedTool('hsa.bed')
exome=BedTool('exome.bed')
print 'bed files loaded'
intergenic=miRNA.subtract(exome)

windowed=intergenic.slop(b=50000,g='human.hg19.genome')


cuffToBed('/raid6/ogan/cuffLinksOutput/'+str(fileno)+'/transcripts.gtf','/raid6/ogan/cuffBedT/'+str(fileno))
#cuffToBed2('/raid6/ogan/cuffLinksOutput/'+str(fileno)+'/transcripts.gtf','/raid6/ogan/cuffBedE/'+str(fileno))

cuff=BedTool('/raid6/ogan/cuffBedT/'+str(fileno))
#cuff2=BedTool('/raid6/ogan/cuffBedT/'+str(fileno))
cuffmiRNA=cuff.intersect(intergenic,wa=True,s=True,wb=True)
#cuff2miRNA=cuff2.intersect(intergenic,wa=True,s=True,wb=True)



cuffmiRNA.saveas('/raid6/ogan/cuffInt/cuffmiRNA'+str(fileno))
#cuff2miRNA.saveas('/raid6/ogan/cuffInt/cuffmiRNA2'+str(fileno))



reads=BedTool("/raid6/ogan/RNASeq_Tuxedo3_Cufflink2_Mapping.RDhi"+str(fileno)+".tophat2.accepted_hits.bed")
inter=reads.intersect(windowed,wa=True,s=True)
print 'intersection complete'

contigs=inter.merge()
select=contigs.intersect(intergenic,wa=True,wb=True)
select.saveas('/raid6/ogan/select/select'+str(fileno))
print 'contigs are selected'

f=open('/raid6/ogan/maxLength/maxLength'+str(fileno),'w')
selectFile=open('/raid6/ogan/select/select'+str(fileno))
while True:
    line=selectFile.readline()
    if line=='':
        break
    lineSplit=line.split('\t')
    name=lineSplit[6]
    start=lineSplit[1]
    end=lineSplit[2]
    chrom=lineSplit[0]
    strand=lineSplit[8][0]
    f.write(chrom+'\t'+start+'\t'+end+'\t'+name+'\t10\t'+strand+'\n')
f.close()
selectFile.close()

f=open('/raid6/ogan/cuffInt/cuffmiRNA'+str(fileno))
cuffmiRNA=open('/raid6/ogan/cuffmiRNAfix/cuffmiRNAfix'+str(fileno),'w')
while True:
    line=f.readline()
    if line=='':
        break
    lineSplit=line.split('\t')
    name=lineSplit[9]
    start=lineSplit[1]
    end=lineSplit[2]
    chrom=lineSplit[0]
    strand=lineSplit[5][0]
    score=lineSplit[4]
    cuffID=lineSplit[3]
    cuffmiRNA.write(chrom+'\t'+start+'\t'+end+'\t'+name+";"+cuffID+'\t'+score+'\t'+strand+'\n')
f.close()
cuffmiRNA.close()


f=open('/raid6/ogan/cuffInt/cuffmiRNA'+str(fileno))
cuffmiRNA=open('/raid6/ogan/cuffmiRNAfix/cuff2miRNAfix'+str(fileno),'w')
while True:
    line=f.readline()
    if line=='':
        break
    lineSplit=line.split('\t')
    name=lineSplit[9]
    start=lineSplit[1]
    end=lineSplit[2]
    chrom=lineSplit[0]
    strand=lineSplit[5][0]
    score=lineSplit[4]
    cuffID=lineSplit[3]
    cuffmiRNA.write(chrom+'\t'+start+'\t'+end+'\t'+name+";"+cuffID+'\t'+score+'\t'+strand+'\n')
f.close()
cuffmiRNA.close()



cuffmiRNA=BedTool('/raid6/ogan/cuffmiRNAfix/cuffmiRNAfix'+str(fileno))
#cuff2miRNA=Bedtool('/raid6/ogan/cuffmiRNAfix/cuff2miRNAfix'+str(fileno))

maxlength=BedTool('/raid6/ogan/maxLength/maxLength'+str(fileno))
coverageMax=newInter.coverage(maxlength,d=True)
coverageMax.saveas('/raid6/ogan/coverageMax/coverageMax'+str(fileno))

coverageCuff=newInter.coverage(cuffmiRNA,d=True)
coverageCuff2=newInter.coverage(cuff2miRNA,d=True)

coverageCuff.saveas('/raid6/ogan/coverageCuff/coverageCuff'+str(fileno))
coverageCuff2.saveas('/raid6/ogan/coverageCuff/coverageCuff2'+str(fileno))

print 'library complete'


