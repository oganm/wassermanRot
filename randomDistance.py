import re
import random
import numpy as np
from Bio import SeqIO
from Bio import SeqUtils
import ogbox
import math
repeats=500
cappedFile='RelaxedQPCR'
#######################################################
#find the capped ones in all of them
with open('humanS/names') as name:
    names=name.readlines()
    numberList=range(1,len(names)+1)
    
def idremove(string):
    string=re.sub(r'id.*?_','',string)
    return string

propNames=map(idremove,names)

with open(cappedFile) as name:
    cappedNames=name.readlines()

#remove the \r character
def rRemove(string):
    string=re.sub(r'\r','',string)
    return string

cappedNames=map(rRemove,cappedNames)

def lookinPropNames(string):
    try:
        return propNames.index(string.lower())+1
    except ValueError:
        return propNames.index(string[0:len(string)-1].lower()+'_1\n')

cappedLocations=map(lookinPropNames,cappedNames)
allSeqs=ogbox.readFasta('hallMirShortM.fa')
mirSeqs=allSeqs[0:19]
allSeqs=allSeqs[19:len(allSeqs)]

    


cappedSeqs=[allSeqs[i] for i in cappedLocations]
#sequences of capped ones
sequences=[x.seq for x in cappedSeqs]
GCcontent=map(SeqUtils.GC,sequences)

sample=cappedLocations[:]
print 'gc contents of capped ones are acquired'
#################################################################
#the data is loaded below so operation will be continued there

#seqs=ogbox.readFasta('cappedMirBaseOut/seqs.fasta')
#sequences=[x.seq for x in seqs]
#GCcontent=map(SeqUtils.GC,sequences)



#######################################################
#read capped ones. capped ones only, the distance value depends on other distances too
'''with open('cappedMirBaseOut/dist-list') as distances:
    dc=[None]*171
    step=0
    while True:
        line=distances.readline()
        if line=='':
            break
        lineSplit=line.split(' ')
        dc[step]=lineSplit[2]
        step=step+1
         
dc=map(int,dc)
meanCapped=np.mean(dc)
stdCapped=np.std(dc)
#save the files
with open('distancesOCapped','w') as f:
    print >>f,'\n'.join(map(str,dc))

print 'capped distances saved'   '''  
##############################################
#get all distances
#preallocate
lineNo=ogbox.flines('humanS/dist-list')
f=[None]*lineNo
s=[None]*lineNo
d=[None]*lineNo
step=0
with open('humanS/dist-list') as distances:
    while True:
        line=distances.readline()
        if line=='':
            break
        lineSplit=line.split(' ')
        f[step]=lineSplit[0]
        s[step]=lineSplit[1]
        d[step]=lineSplit[2]
        step=step+1
f=map(int,f)
s=map(int,s)  
d=map(int,d)

print 'info of all mir is calculated'
#put tuples in a seperate list for fast search
def tupler(foo,bar):
    return (foo,bar)
fs=map(tupler,f,s)

np.std(d)
np.mean(d)
np.median(d)
#######################################################
##looking at the capped data from full data
'''with open('allDistances','w') as file:
    print >>file,'\n'.join(map(str,d))

def nCr(n,r):
    f = math.factorial
    return f(n) / f(r) / f(n-r)

sample=cappedLocations[:]
with open('RelaxedQPCRDistances','w') as file:
    sample.sort()
    sample.reverse()
    infSize=int(nCr(len(cappedLocations),2))
    inf=[None]*2016
    ins=[None]*2016
    step=0
    #structure data to look for
    while True:
        for i in range(1,len(sample)):
            inf[step]=sample[0]
            ins[step]=sample[i]
            step+=1
        del sample[0]
        if sample==[]:
            break
        
    #look for the data
    sampDistances=[None]*len(inf)
    print 'reached hard part'
    step=0
    for t in range(0,len(inf)):
        location=fs.index((inf[t],ins[t]))
        sampDistances[t]=d[location]
        step=step+1
        print step
        
    print 'finished hard part '
    print np.mean(sampDistances)
    print np.std(sampDistances)
    print np.median(sampDistances)
    print >>file, '\n'.join(map(str,sampDistances))

'''
    
#####################################################################
#match GC content to read

allMirDistLocations=[i for i,val in enumerate(s) if val in range(1,20)]
arrayD=np.array(d)

allMirDist=arrayD[allMirDistLocations]


###meanwhile lets take the distances of human capped to mir capped###
print 'meanwhile lets take the distances of human capped to mir capped'
cappedMirDistLocations=[i for i,val in enumerate(f) if val in cappedLocations]
mirToCappedDistLocations=list(set(cappedMirDistLocations)&set(allMirDistLocations))


with open('humanCappedToMouseCapped','w') as file:
    print >>file,'\n'.join(map(str,arrayD[mirToCappedDistLocations].tolist()))
print 'done'

allSequences=[x.seq for x in allSeqs]
AllGC=map(SeqUtils.GC,allSequences)
AllGC=[999]*19+AllGC

def nCr(n,r):
    f = math.factorial
    return f(n) / f(r) / f(n-r)

infSize=int(nCr(len(cappedLocations),2))

stdList=[None]*repeats
meanList=[None]*repeats
medianList=[None]*repeats
distanceList=[None]*repeats
print 'matched sample calculations'
for zaytung in xrange(0,repeats):
    sample=[None]*len(GCcontent)
    arrayGC=np.array(AllGC)
    for i in range(0,len(GCcontent)):
        matches=(arrayGC>GCcontent[i]-5)&(arrayGC<GCcontent[i]+5)
        possible=np.where(matches)
        possible=possible[0].tolist()
        random.shuffle(possible)
        sample[i]=possible[0]+1 #fix for index idiot
        arrayGC[possible[0]]=999


    #mir's distances to sample
    print 'mouse mirs distances to chosen sample'
    temp=sample[:]
    allSempDistLocations=[i for i,val in enumerate(f) if val in temp]
    mirToSempDistLocations=list(set(allSempDistLocations)&set(allMirDistLocations))
    print 'locations are found'
    distanceList[zaytung]=np.mean(arrayD[mirToSempDistLocations])
    with open('matchedDistanceToMirCap/'+str(zaytung),'w') as file:
        print >>file,'\n'.join(map(str,arrayD[mirToSempDistLocations].tolist()))

    print 'distances to mouse capped are written'

    sample.sort()
    sample.reverse()
    inf=[None]*infSize
    ins=[None]*infSize
    step=0

    

    #structure data to look for
    while True:
        for i in range(1,len(sample)):
            inf[step]=sample[0]
            ins[step]=sample[i]
            step+=1
        del sample[0]
        if sample==[]:
            break
        
    #look for the data
    sampDistances=[None]*len(inf)
    print 'reached hard part'


    for t in range(0,len(inf)):
        location=fs.index((inf[t],ins[t]))
        sampDistances[t]=d[location]
    
    
    print 'finished hard part '+str(zaytung)
    print np.mean(sampDistances)
    print np.std(sampDistances)
    print np.median(sampDistances)
    stdList[zaytung]=np.std(sampDistances)
    meanList[zaytung]=np.mean(sampDistances)
    medianList[zaytung]=np.median(sampDistances)
    with open('matchedDistr/repeat'+str(zaytung),'w') as file:
       print >>file,'\n'.join(map(str,sampDistances))      

with open('matchedDistancetoMouseCapped','w') as file:
    print >>file,'\n'.join(map(str,distanceList))

with open('matchedMeans','w') as file:
    print >>file,'\n'.join(map(str,meanList))

with open('mathecMedians','w') as file:
    print >>file,'\n'.join(map(str,medianList))

with open('mathecSTD','w') as file:
    print >>file,'\n'.join(map(str,stdList))





stdList=[None]*repeats
meanList=[None]*repeats
medianList=[None]*repeats
#################################################################
#random sampling

'''print 'random'
for zaytung in xrange(0,repeats):
    print zaytung
    random.shuffle(numberList)
    sample=numberList[1:len(GCcontent)+1]
    sample.sort()
    sample.reverse()
    inf=[None]*infSize
    ins=[None]*infSize
    step=0
    while True:
        for i in range(1,len(sample)):
            inf[step]=sample[0]
            ins[step]=sample[i]
            step+=1
        del sample[0]
        if sample==[]:
            break
        
        
    sampDistances=[None]*len(inf)
    print 'beginning hard part'
    for t in xrange(0,len(inf)):
        location=fs.index((inf[t],ins[t]))
        sampDistances[t]=d[location]
    print 'finished hard part '+str(zaytung)
    stdList[zaytung]=np.std(sampDistances)
    meanList[zaytung]=np.mean(sampDistances)
    medianList[zaytung]=np.median(sampDistances)
    with open('randomDistr/repeat'+str(zaytung),'w') as file:
       print >>file,'\n'.join(map(str,sampDistances)) 


with open('randomMeans','w') as file:
    print >>file,'\n'.join(map(str,meanList))

with open('randomMedians','w') as file:
    print >>file,'\n'.join(map(str,medianList))

with open('randomSTD','w') as file:
    print >>file,'\n'.join(map(str,stdList))'''


#################################################################
#Distances to mouse capped miRNA (of all)


cappedLocations
sArray=np.array(s)


allMirDistLocations=[i for i,val in enumerate(s) if val in range(1,20)]

arrayD=np.array(d)

allMirDist=arrayD[allMirDistLocations]

np.mean(allMirDist)



