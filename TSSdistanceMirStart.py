import urllib.request
import re
mirStart=urllib.request.urlopen('http://mirstart.mbc.nctu.edu.tw/browse.php')

toParse=mirStart.readall()
toParse=toParse.decode()
newlinks=re.findall(r'href="mirna.php.........',toParse)
print('website loaded')
theyAreClose=[]
theyAreReallyClose=[]
for i in range(0,len(newlinks)):
    print('miRNA no: '+str(i+1))
    newLink=newlinks[i][6:len(newlinks[i])]
    mirLink=urllib.request.urlopen('http://mirstart.mbc.nctu.edu.tw/'+newLink)
    toParse=mirLink.readall()
    toParse=toParse.decode()
    name=re.findall(r'hsa-.*?<',toParse)
    name=name[0][0:len(name[0])-1]
    distances=re.findall(r'[(][0-9]*?[)]',toParse)

    

    distances=[int(distances[x][1:len(distances[x])-1]) for x in range(1,len(distances))]
    #remove the first one as it also the second one

    for j in range(0,len(distances)):
        if (distances[j]<=10):
            theyAreReallyClose.append(name)
            print('really close detected')
        if (distances[j]>10)&(distances[j]<=20):
            theyAreClose.append(name)
            print('close detected')
        if distances[j]>20:
            print('meh')
            
            
    
    


    



#toParse=mirStart.readline()

#toParse=toParse.decode()


#toParse=toParse.strip() #removes whitespace

#regular expression search for link
#newlink=re.findall(r'href="mirna.php.........',toParse)
#newlink=re.findall(r' ',toParse)

#miRNA=urllib.request.urlopen('http://mirstart.mbc.nctu.edu.tw/'+
#                             newlink
