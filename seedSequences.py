import sys
sys.path.append('/home/ogan/-PModules')
import ogbox
import seqH19 as seqH
#import seqMM10 as seqM
import re
from Bio import SeqIO
from Bio.Seq import Seq
import numpy as np
import os

outdir          = '/home/ogan/Desktop/PORTAL/Ubuntu Python 2.7/results3/ALLMICRORNA/'
hMirbaseAd      = '/home/ogan/Desktop/PORTAL/hsa.gff3'
#nameAd          = outdir+'selection'
mouseCappedAd   = '/home/ogan/Desktop/PORTAL/Ubuntu Python 2.7/namesMIR'
mMirbaseAd      = '/home/ogan/Desktop/PORTAL/Ubuntu Python 2.7/mouseCluster/mmu.gff3.txt'

hMir=ogbox.mirbaseGFF(hMirbaseAd)
#alternating between this part and the next is to target all mouse microRNA or 
'''
#opening list of RNAs to get the seeds of
with open(nameAd) as cappedHFile:
    cappedLines = cappedHFile.readlines()
 
def replaceUnderscore(inp):
    return inp.replace('_', '-').replace('\r\n', '').replace('\n','')

cappedLines = map(replaceUnderscore, cappedLines)

seedNames = [None] * 2 * len(cappedLines)
seedSeqs  = [None] * 2 * len(cappedLines)
seedStep = 0


for i in xrange(0,len(cappedLines)):
    index=ogbox.find(hMir.name, cappedLines[i])[0]
    if type(index)!=type(1):
        index=index[0]
    step=1
    while hMir.mirType[index+step] == 'miRNA':
        if hMir.strand[index] == '-':
            seedEnd = hMir.end[index + step]
            seedStart = hMir.end[index + step] - 8

        if hMir.strand[index] == '+':
            seedStart = hMir.start[index + step]
            seedEnd = seedStart + 8
            
        seedSeq = seqH.giveSeq(hMir.chro[index], seedStart - 1, seedEnd-1, \
                               strand = hMir.strand[index], string = True)

        
        seedSeqs[seedStep] = seedSeq
        seedNames[seedStep] = hMir.name[index+step]
        seedStep = seedStep + 1
        step = step + 1
        
        if index+step+1>len(hMir.mirType):
            break

'''
seedNames = [None] * 2 * len(hMir.name)
seedSeqs  = [None] * 2 * len(hMir.name)
seedStep = 0

def replaceUnderscore(inp):
    return inp.replace('_', '-').replace('\r\n', '').replace('\n','')

for i in xrange(0,len(hMir.name)):
    if hMir.mirType[i] == 'miRNA':
        if hMir.strand[i] == '-':
            seedEnd = hMir.end[i]
            seedStart = hMir.end[i] - 8

        if hMir.strand[i] == '+':
            seedStart = hMir.start[i]
            seedEnd = seedStart + 8
            
        seedSeq = seqH.giveSeq(hMir.chro[i], seedStart - 1, seedEnd-1, \
                               strand = hMir.strand[i], string = True)
        
        seedStep = seedStep + 1
        seedSeqs[seedStep] = seedSeq
        seedNames[seedStep] = hMir.name[i]                       
                               
#'''

seedSeqs=ogbox.trimNone(seedSeqs)
seedNames=ogbox.trimNone(seedNames)

seqH.clearModule()

###################################################################3
#Mouse part
import seqMM10 as seqM
with open(mouseCappedAd) as mCappedHFile:
    mCappedLines = mCappedHFile.readlines()
 
mCappedLines = map(replaceUnderscore, mCappedLines)

mMir=ogbox.mirbaseGFF(mMirbaseAd)

mSeedNames = [None]*2*len(mCappedLines)
mSeedSeqs = [None]*2*len(mCappedLines)
seedStep=0

for i in xrange(0,len(mCappedLines)):
    index=ogbox.find(mMir.name, mCappedLines[i])[0]
    step=1
    while mMir.mirType[index+step] == 'miRNA':
        if mMir.strand[index] == '-':
            seedEnd = mMir.end[index + step]
            seedStart = mMir.end[index + step] - 8

        if mMir.strand[index] == '+':
            seedStart = mMir.start[index + step]
            seedEnd = seedStart+8
            
        seedSeq = seqM.giveSeq(mMir.chro[index], seedStart - 1, seedEnd-1, \
                               strand = mMir.strand[index], string = True)
        
        mSeedSeqs[seedStep] = seedSeq
        mSeedNames[seedStep] = mMir.name[index+step]
        seedStep = seedStep + 1
        step = step + 1
        if index+step+1>len(hMir.mirType):
            break

mSeedSeqs=ogbox.trimNone(mSeedSeqs)
mSeedNames=ogbox.trimNone(mSeedNames)
seqM.clearModule()
#seed is taken from the very beginning. for correct comparison start from 2nd index
def similarity7(seq1,seq2):
    o=0
    for i in range(1,8):
       o=o+int(seq1[i]==seq2[i])
    return o
        
        
def similarity6(seq1,seq2):
    o=0
    for i in range(1,7):
        o=o+int(seq1[i]==seq2[i])
    return o

#def similarity(seq1,seq2):
#    for i in range(8,0,-1):
#        if seq1[0:i]==seq2[0:i]:
#            return i
#    return 0
            
#simMat=np.zeros((len(mSeedSeqs),len(seedSeqs)))
#simMat=simMat.astype(int)

oneDfile=open('1DSeedSimilarity7','w')
oneDfile2=open('1DSeedSimilarity6','w')
for i in range(0,len(mSeedSeqs)):
    for t in range(0,len(seedSeqs)):
        #simMat[i,t]=similarity7(mSeedSeqs[i],seedSeqs[t])
        if similarity7(mSeedSeqs[i],seedSeqs[t])>0:
            print >>oneDfile, mSeedNames[i] +' '+seedNames[t]+' '+str(similarity7(mSeedSeqs[i],seedSeqs[t])) 
            print >>oneDfile2,mSeedNames[i] +' '+seedNames[t]+' '+str(similarity6(mSeedSeqs[i],seedSeqs[t]))
#np.savetxt('seedSimilarity',simMat,'%int')        
oneDfile.close()
oneDfile2.close()

os.system("sort -rnk3 1DSeedSimilarity7 > '"+outdir+"1DSorted7'")
os.system("sort -rnk3 1DSeedSimilarity6 > '"+outdir+"1DSorted6'")

#for i in seedSeqs:
    




