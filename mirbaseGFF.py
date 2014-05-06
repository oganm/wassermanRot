import time
timer = time.time()
import seq #comment out for fast test
import re
import os

#instead of loading the whole human genome, index it to avoid excessive RAM usage

file='C:/Users/Ogan/Documents/Rotation 2/files/hsa.gff3.txt'

#delete the comments in the file before using

f=open(file)


lines=f.readlines()
f.close()

#parse gff
#use objects next time or at least comprehension
chr=[]
mirType=[]
str=[]
end=[]
strand=[]
name=[]
edge=[]
derives=[]
ID=[]
for i in range(0,len(lines)):
    mirna=lines[i].split('\t')
    try:
        temp=int(mirna[0][3:len(mirna[0])])
    except ValueError:
        temp=(mirna[0][3:len(mirna[0])])

    chr=chr+[temp]
    mirType=mirType+[mirna[2]]

    str=str+[int(mirna[3])]
    end=end+[int(mirna[4])]

    strand=strand+[mirna[6]]
    temp=re.findall(r'Name=.*?\n',mirna[8])[0]
    name=name+[temp[5:len(temp)-1]]
    tempID=re.findall(r'ID=.*?;',mirna[8])[0]
    ID=ID+[tempID[3:len(tempID)-1]]

    if mirType[i]=='miRNA_primary_transcript':
        edge=edge+['NA']
        derives=derives+['NA']
    else:
        temp=re.findall(r'[5|3]p',name[i])
        tempD=re.findall(r'MI.......',name[i])
        derives=derives+tempD
        
        tempN=re.findall(r'.*?;',name[i])
        name[i]=tempN[0][0:len(tempN[0])-1]
        
        if temp==[]:
            temp=['?']
            
        edge=edge+temp

    
    
#requires seq after this point
unkowns=[i for i, x in enumerate(edge) if x == "?"]

toFold=open('C:/Users/Ogan/Documents/Rotation 2/files/toFold.txt','w')
foldedStart=[]
foldedEnd=[]
for i in range(0,len(unkowns)):
    tempIndex=[t for t, x in enumerate(ID) if x == derives[unkowns[i]]][0]
    print('>'+name[tempIndex],file=toFold)
    foldedStart=foldedStart+[str[tempIndex]]
    foldedEnd=foldedEnd+[end[tempIndex]]
    sequence=seq.giveSeq(chr[tempIndex],str[tempIndex],end[tempIndex],strand[tempIndex],True)
    print(sequence,file=toFold)

toFold.close()

# note to self: \ is for line continuation
os.system("rnafold < 'C:\\Users\\Ogan\\Documents\\Rotation 2\\files\\toFold.txt' \
> 'C:\\Users\\Ogan\\Documents\\Rotation 2\\files\\output.txt'")


foldOut=open("C:\\Users\\Ogan\\Documents\\Rotation 2\\files\\output.txt")

foldLines=foldOut.readlines()
foldedName=[]

foldedName=[foldLines[i][1:len(foldLines[i])-1] for i in range(0,len(foldLines),3)]
foldedFold=[foldLines[i] for i in range(2,len(foldLines),3)]

#use of compiled regular expression
side=[]
results=open('C:/Users/Ogan/Documents/Rotation 2/files/sides.txt','w')
for i in range(0,len(foldedFold)):
    m=re.search(r'\([.]*?\)',foldedFold[i])
    hairpin=int((m.span()[0]+m.span()[1])/2+foldedStart[i])
    tI=unkowns[i] #temp index of RNA
    tIhp=[t for t, x in enumerate(ID) if x == derives[unkowns[i]]][0] #temp index of hairpin
    
    startMir=str[tI]
    endMir=end[tI]
    if ((startMir<hairpin) & (endMir<hairpin)):
        side=side+['5p']*(strand[tIhp]=='+')+['3p']*(strand[tIhp]=='-')
    elif ((startMir>hairpin) & (endMir>hairpin)):
        side=side+['3p']*(strand[tIhp]=='+')+['5p']*(strand[tIhp]=='-')
    elif (startMir-hairpin+endMir-hairpin)<0:
        side=side+['?5p']*(strand[tIhp]=='+')+['?3p']*(strand[tIhp]=='-')
    elif (startMir-hairpin+endMir-hairpin)>0:
        side=side+['?3p']*(strand[tIhp]=='+')+['?5p']*(strand[tIhp]=='-')
    else:
        side=side+['?']

        
    print(name[tI]+'\t'+name[tIhp]+'\t'+side[i],file=results)


results.close()

wtf=[i for i, x in enumerate(side) if (x == "?3p")|(x=='?5p')]

elapsed = time.time() - timer
minutes=elapsed/60

