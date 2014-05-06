#to get %1 interval from the locarna results
#folder='shortlocarna'
folder='longlocarna'

#distances=read.table('/home/ogan/Downloads/RNAclust-1.3/examples/humanMirOut/dist-list')
distances=read.table('/home/ogan/Desktop/NewRNAclust/out/dist-list')

#names=read.table('/home/ogan/Downloads/RNAclust-1.3/examples/humanMirOut/names')
names=read.table('/home/ogan/Desktop/NewRNAclust/out/names')

select=distances[(distances$V2 %in% 1:19)&(distances$V1 %in% 20:max(distances$V1)),]
####select=distances[(distances$V2 %in% c(1:7,9:15))&(distances$V1 %in% 20:max(distances$V1)),]


names=as.character(names$V1)
names=gsub('id.*?_','',names)



q=quantile(select$V3,probs=0.01)

q01=select$V1[select$V3<=q]

fullq01=select[select$V3<=q,]

uniqueQ01=unique(q01)

findMouseMir = function(mirNo){
  index=which(fullq01$V1 %in% mirNo)
  out=''
  for (i in 1:length(index)){
    out=paste(out,names[fullq01$V2[index[i]]],sep=',')
  } 
  return(substr(out,2,nchar(out)))
}


write.table(names[uniqueQ01],row.name=FALSE,col.name=F,quote=F,file=paste('/home/ogan/Desktop/PORTAL/Ubuntu Python 2.7/results3/',folder,'/selection',sep=''))

write.table(data.frame(names=names[uniqueQ01],
                       cappedMouse = as.character(lapply(uniqueQ01,findMouseMir))),row.name=FALSE,col.name=F,quote=F,file=paste('/home/ogan/Desktop/PORTAL/Ubuntu Python 2.7/results3/',folder,'/selectionWithMouse',sep=''))


