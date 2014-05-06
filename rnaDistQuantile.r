require(ctc)
require(ape)
RNAnames=read.table('names')
RNAnames=as.character(RNAnames$V1)
RNAnames=gsub('>','',RNAnames)
RNAnames=make.names(RNAnames,unique=T)
txt=readChar('justDistances',file.info('justDistances')$size)

mouseDot=gsub('_','.',as.character(mouse$V1))


# hclust(dist(data.frame(a=c(1,2,3,4),b=c(1,2,3,4))))

distances=read.table('justDistances',fill=T,col.names=1:1908)
colnames(distances)=RNAnames[1:length(RNAnames)-1]
rownames(distances)=RNAnames[2:length(RNAnames)]

distanceMat=as.dist(distances)

#for tree stuff. didnt work with this data.
#cluster=hclust(distanceMat)
#newick=hc2Newick(cluster)
#tree=read.tree(text=newick)

#colors=rep('black',length(tree$tip.label))
#colors[which(tolower(tree$tip.label) %in% mouseDot)]='red'


#plot.phylo(tree,no.margin=T,direction='downwards',tip.color=colors)
#plot.phylo(mytree,no.margin=T,tip.color=colors,direction='downwards')


#human=read.table(file.choose())
#mouse=read.table(file.choose())

toMouse=distances[38:nrow(distances),1:38]
toShortMouse=(distances[38:nrow(distances),20:38])
toLongMouse=distances[38:nrow(distances),1:19] #distances human to 
#toLongMouse=distances[38:nrow(distances),c(1:7,9:15)] #distances human to 
shortMouseMouse=distances[20:38,20:38] #distances of short to short mouse
longMouseMouse=distances[1:19,1:19] #distances of long to long mouse
#longMouseMouse=distances[c(1:7,9:15),c(1:7,9:15)] #distances of long to long mouse

frame=data.frame(distances=c(unlist(t(shortMouseMouse)),unlist(t(toShortMouse))),
                 label=c(rep('mousetomouse',length(unlist(t(shortMouseMouse)))),rep('humantomouse',length(unlist(t(toShortMouse))))))


require(ggplot2)

ggplot(frame, aes(x=distances, fill=label)) + geom_density(alpha=.3)

frame=data.frame(distances=c(unlist(t(longMouseMouse)),unlist(t(toLongMouse))),
                 label=c(rep('mousetomouse',length(unlist(t(longMouseMouse)))),rep('humantomouse',length(unlist(t(toLongMouse))))))

ggplot(frame, aes(x=distances, fill=label)) + geom_density(alpha=.3)


longQ=quantile(unlist(t(toLongMouse)),probs=0.01)


shortQ=quantile(unlist(t(toShortMouse)),probs=0.01)

selectLong=toLongMouse<longQ
selectShort=toShortMouse<shortQ
length(which(apply(selectLong,1,any)))
length(which(apply(selectShort,1,any)))

length(intersect(which(apply(selectLong,1,any)),which(apply(selectShort,1,any))))

rownames(selectLong)[intersect(which(apply(selectLong,1,any)),which(apply(selectShort,1,any)))]

rownames(longQ)

giveMirNameLong=function(index){
  return(paste(gsub('[.]','_',gsub('[.]LONG[.]','',colnames(selectLong)[which(selectLong[index,])])),collapse=','))
}

giveMirNameShort=function(index){
  return(paste(gsub('[.]','_',gsub('[.]LONG[.]','',colnames(selectShort)[which(selectShort[index,])])),collapse=','))
  
}


write.table(gsub('[.]','_',rownames(selectLong)[intersect(which(apply(selectLong,1,any)),which(apply(selectShort,1,any)))]),file="onePercentIntersect",quote=F,col.names=F,row.names=F)

write.table(data.frame(V1=gsub('[.]','_',rownames(selectLong)[which(apply(selectLong,1,any))]),V2=as.character(lapply(as.numeric(which(apply(selectLong,1,any))),giveMirNameLong))),file='longSelectionWithMouse',row.names=F,col.names=F,quote=F)
write.table(data.frame(V1=gsub('[.]','_',rownames(selectLong)[which(apply(selectLong,1,any))])),file='longSelection',row.names=F,col.names=F,quote=F)

write.table(data.frame(V1=gsub('[.]','_',rownames(selectShort)[which(apply(selectShort,1,any))]),V2=as.character(lapply(as.numeric(which(apply(selectShort,1,any))),giveMirNameShort))),file='shortSelectionWithMouse',row.names=F,col.names=F,quote=F)
write.table(data.frame(V1=gsub('[.]','_',rownames(selectShort)[which(apply(selectShort,1,any))])),file='shortSelection',row.names=F,col.names=F,quote=F)

