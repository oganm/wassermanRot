from Bio import SeqIO
from os import listdir
from os.path import isfile, join
import math
import numpy as np
import re

print 'official ogbox copy'


#return whole fasta sequences as a list
def readFasta(filename):
    with open(filename, 'rU') as f:
        i = 0
        l = flines(filename) / 2
        seqs = [None] * l
        for record in SeqIO.parse(f, 'fasta'):
            seqs[i] = record
            i = i + 1
        return seqs
     

#return number of lines in a file
def flines(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

#list files of a directory
def getFiles(directory):
    onlyfiles = [f for f in listdir(directory) if isfile(join(directory, f))]
    return onlyfiles

#combination
def nCr(n, r):
    f = math.factorial
    return f(n) / f(r) / f(n - r)

#removes None values from a list
def trimNone(L):
    return [x for x in L if x is not None]

#find object in list. can be slow
def find(L, obj):
    return [i for i, x in enumerate(L) if x == obj]


#read file output list of lines. removes endline characters. just because
def removeNewLine(string):
    return string.replace('\n','').replace('\r','')
    
#outputs a list of lines from a file
def readtxt(fileName):
    with open(fileName) as leFile:
        out=leFile.readlines()
        return map(removeNewLine,out)

def writeLines(lineList,fileName):
    with open(fileName,'w') as leFile:
        for i in lineList:
            print >>leFile,i

# a class to read GFF files from mirbase. does not ignore the comments. manually
#delete them.
#TURN IT INTO A PROPER LIST OF OBJECTS BEFORE RE-USE
#AND IGNORE THE BLOODY comments
class mirbaseGFF:
    def __init__(self, filename):
        with open(filename) as mirbase:
            lines = mirbase.readlines()
        self.chro    = []
        self.mirType = []
        self.start   = []
        self.end     = []
        self.strand  = []
        self.name    = []
        self.edge    = []
        self.derives = []
        self.ID      = []


        for i in range(0, len(lines)):
            mirna = lines[i].split('\t')
            try:
                temp = int(mirna[0][3:len(mirna[0])])
            except ValueError:
                temp = (mirna[0][3:len(mirna[0])])

            self.chro = self.chro + [temp]
            self.mirType = self.mirType + [mirna[2]]

            self.start = self.start + [int(mirna[3])]
            self.end = self.end + [int(mirna[4])]

            self.strand = self.strand + [mirna[6]]
            temp = re.findall(r'Name=.*?\n', mirna[8])[0]
            self.name = self.name + [temp[5:len(temp) - 1]]
            tempID = re.findall(r'ID=.*?;', mirna[8])[0]
            self.ID = self.ID + [tempID[3:len(tempID) - 1]]
            
            if self.mirType[i] == 'miRNA_primary_transcript':
                self.edge = self.edge + ['NA']
                self.derives = self.derives + ['NA']
            else:
                temp = re.findall(r'[5|3]p', self.name[i])
                tempD = re.findall(r'MI.......', self.name[i])
                self.derives = self.derives + tempD
        
                tempN = re.findall(r'.*?;', self.name[i])
                self.name[i] = tempN[0][0:len(tempN[0]) - 1]
        
                if temp == []:
                    temp = ['?']
            
                self.edge = self.edge + temp



#def readTable(filename,split='\t')
#while True:
#    line=selectFile.readline()
#    if line=='':
#        break
#    lineSplit=line.split('\t')
    
