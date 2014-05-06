#master table

require('VennDiagram')


names = read.table('/home/ogan/Desktop/PORTAL/mirName&Origin')

colnames(names)=c('mirName','precName')
########################################################################### to finalframe 1
names$mirName=tolower(names$mirName)
############################################################################
listCaged=read.table('/home/ogan/Desktop/PORTAL/FANTOM5/newCapped')
listCaged=tolower(listCaged$V1)
#length(listCaged) #188
######################################################################### to final frame 2
isCaged=names$precName %in% listCaged
#######################################################################
shareSeed6=read.table('/home/ogan/Desktop/PORTAL/Ubuntu Python 2.7/results3/ALLMICRORNA/1DSorted6')
shareSeed7=read.table('/home/ogan/Desktop/PORTAL/Ubuntu Python 2.7/results3/ALLMICRORNA/1DSorted7')

tolerant6=unique(shareSeed6[shareSeed6$V3>=5,])
shareSeed6=unique(shareSeed6[shareSeed6$V3==6,])

tolerant7=unique(shareSeed7[shareSeed7$V3>=6,])
shareSeed7=unique(shareSeed7[shareSeed7$V3==7,])




#capital for final shared seeds
ShareSeed6=rep(NA,length(names$mirName))
Tolerant6=rep(NA,length(names$mirName))
ShareSeed7=rep(NA,length(names$mirName))
Tolerant7=rep(NA,length(names$mirName))
############################################################### to final frame 3
for (i in 1:length(shareSeed6$V1)){
  ShareSeed6[which(names$mirName %in% tolower(shareSeed6$V2[i]))]=as.character(shareSeed6$V1[i])
}

for (i in 1:length(shareSeed7$V1)){
  ShareSeed7[which(names$mirName %in% tolower(shareSeed7$V2[i]))]=as.character(shareSeed7$V1[i])
}
for (i in 1:length(tolerant6$V1)){
  Tolerant6[which(names$mirName %in% tolower(tolerant6$V2[i]))]=as.character(tolerant6$V1[i])
}
for (i in 1:length(tolerant7$V1)){
  Tolerant7[which(names$mirName %in% tolower(tolerant7$V2[i]))]=as.character(tolerant7$V1[i])
}
####
longRNAclust1p=read.table('/home/ogan/Desktop/PORTAL/Ubuntu Python 2.7/results3/longlocarna/selectionWithMouse')
longRNAclust1p$V1=gsub('_','-',as.character(longRNAclust1p$V1))

shortRNAclust1p=read.table('/home/ogan/Desktop/PORTAL/Ubuntu Python 2.7/results3/shortlocarna/selectionWithMouse')
shortRNAclust1p$V1=gsub('_','-',as.character(shortRNAclust1p$V1))

longRNAdist1p=read.table('/home/ogan/Desktop/PORTAL/Ubuntu Python 2.7/results3/longRNAdist/selectionWithMouse')
longRNAdist1p$V1=gsub('_','-',as.character(longRNAdist1p$V1))

shortRNAdist1p=read.table('/home/ogan/Desktop/PORTAL/Ubuntu Python 2.7/results3/shortRNAdist/selectionWithMouse')
shortRNAdist1p$V1=gsub('_','-',as.character(shortRNAdist1p$V1))


################################################################## to final frame 4
isLongRNAclust1p = rep(NA,length(names$mirName))
isShortRNAclust1p = rep(NA,length(names$mirName))
isLongRNAdist1p = rep(NA,length(names$mirName))
isShortRNAdist1p = rep(NA,length(names$mirName))

for (i in 1:length(longRNAclust1p$V1)){
  isLongRNAclust1p[which(names$precName %in% tolower(longRNAclust1p$V1[i]))]=as.character(longRNAclust1p$V2[i])
}
for (i in 1:length(shortRNAclust1p$V1)){
  isShortRNAclust1p[which(names$precName %in% tolower(shortRNAclust1p$V1[i]))]=as.character(shortRNAclust1p$V2[i])
}
for (i in 1:length(longRNAdist1p$V1)){
  isLongRNAdist1p[which(names$precName %in% tolower(longRNAdist1p$V1[i]))]=as.character(longRNAdist1p$V2[i])
}
for (i in 1:length(shortRNAdist1p$V1)){
  isShortRNAdist1p[which(names$precName %in% tolower(shortRNAdist1p$V1[i]))]=as.character(shortRNAdist1p$V2[i])
}

#variant that includes full cage coordinates
cagePeaks=read.table('/home/ogan/Desktop/PORTAL/FANTOM5/CAGEHits2')

#there are are multiple matches. reports the one closest to start site. dependent on interval of seach
#set during the creation of preMirStart.bed
cagePeakLoc=rep(NA,length(names$mirName))

for (i in 1:length(names$precName)){
  peaks=which(cagePeaks$V10 %in% names$precName[i])  
  if (length(peaks)>0){
    if (cagePeaks$V6[peaks[1]]=='+'){
      mirStart=cagePeaks$V8[peaks[1]]+30
      cagePeakLoc[i]=cagePeaks$V2[peaks[which.min(abs(cagePeaks$V2[peaks]-mirStart))]]      
    }
    if (cagePeaks$V6[peaks[1]]=='-'){
      mirStart=cagePeaks$V9[peaks[1]]-30
      cagePeakLoc[i]=cagePeaks$V3[peaks[which.min(abs(cagePeaks$V3[peaks]-mirStart))] ]     
    }
  }
  ############################################################# to final frame 5
}


unkownSides=read.table('/home/ogan/Desktop/PORTAL/sides.txt')
sides=rep(NA,length(names$mirName))

for (i in 1:length(names$mirName)){
  if ((tolower(names$mirName[i]) %in% tolower(unkownSides$V1))&(tolower(names$precName[i]) %in% tolower(unkownSides$V2))){
    sides[i]=as.character(unkownSides$V3[which((tolower(unkownSides$V1) %in% tolower(names$mirName[i]))
                                  &(tolower(unkownSides$V2) %in% tolower(names$precName[i])))])
  }
  
  if (length(grep('-5p',names$mirName[i]))>0){
    sides[i]='5p'
  }
  if (length(grep('3p',names$mirName[i]))>0){
    sides[i]='3p'
  }
  
  
}


finalFrame=data.frame(mirName=names$mirName,precName=names$precName,isCaged=isCaged,seed6=ShareSeed6,seed7=ShareSeed7,
                      longRNAclust=isLongRNAclust1p,shortRNAClust=isShortRNAclust1p, longRNAdist=isLongRNAdist1p, shortRNAdist=isShortRNAdist1p,cageLocation=cagePeakLoc,sides=sides )


write.table(finalFrame,row.name=FALSE,quote=FALSE,file='/home/ogan/Desktop/PORTAL/Ubuntu Python 2.7/results3/allInclusive')

plot.new()
venList=list(shortRNAclust=which(!is.na(finalFrame$shortRNAClust)),longRNAdist=which(!is.na(finalFrame$longRNAdist)),
             longRNAclust=which(!is.na(finalFrame$longRNAclust)),shortRNAdist=which(!is.na(finalFrame$shortRNAdist)))

venn.plot=venn.diagram(venList,filename=NULL,fill=c('red','blue','green','orange'))
grid.draw(venn.plot)



plot.new()
venList=list(shortRNAclust=which(!is.na(finalFrame$shortRNAClust)),longRNAdist=which(!is.na(finalFrame$longRNAdist)),
             longRNAclust=which(!is.na(finalFrame$longRNAclust)),shortRNAdist=which(!is.na(finalFrame$shortRNAdist)),
             CAGED=which(finalFrame$isCaged))

venn.plot=venn.diagram(venList,filename=NULL,fill=c('red','blue','green','orange','cyan'))
grid.draw(venn.plot)

