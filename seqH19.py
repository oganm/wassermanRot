import sys
sys.path.append('/home/ogan/-PModules')
from Bio import SeqIO
from Bio.Seq import Seq
from ogbox import trimNone
#have whole genome in a single fa document labeled as chr1 chr2 etc.
#use cat folderAddress/* regular exp
#changed to take in string name of the chromosome
folder='/home/ogan/-PModules/'
filename=folder+'humangenome.fa'
print(filename)

chrID = [None] * 50
chrSeq = [None] * 50

with open(filename) as fastaFile:
    for record in SeqIO.parse(fastaFile, "fasta") :
        print(record.id)
        print(repr(record.seq))
        print(len(record))
#    try:
#        temp=int(record.id[3:len(record.id)])
#    except ValueError:
#        temp=(record.id[3:len(record.id)])
    
        
        chrID=chrID+[record.id]
        chrSeq=chrSeq+[record.seq]

chrID=trimNone(chrID)
chrSeq=trimNone(chrSeq)


def trial():
    print(chrID)

#enter chr as string or number
def giveSeq(chro,start,end,strand='+',string=False):
    try:
        sequence=chrSeq[chrID.index(chro)][start:end]
    except ValueError:
        sequence=chrSeq[chrID.index('chr'+str(chro))][start:end]
    if strand=='-':
        sequence=sequence.reverse_complement()
    sequence=sequence.transcribe()
    if string==True:
        sequence=sequence.tostring()
    return sequence
    
def clearModule():
    global chrSeq
    global chrID
    del chrSeq
    del chrID
    
    
    
    

