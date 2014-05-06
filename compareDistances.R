require('ggplot2')
cnames=c('Distance','label')
setwd('C:/Users/Ogan/Desktop/PORTAL/For R/human')
cappedDist=read.table('RelaxedQPCRDistances')

allDist=read.table('allDistances')
allFrame=data.frame(allDist,label='all')
colnames(allFrame)=c(cnames)

matchedMean=read.table('matchedMeans')
matchedFrame=data.frame(matchedMean,label='GCmatchedMeans')
colnames(matchedFrame)=c(rnames)
#remove outlier 0s
fineCapped=cappedDist[-which(t(cappedDist) %in% 0),]

#wilcoxon test between all distances and 
wilcox.test(fineCapped,as.vector(allDist$V1),alternative='less')



frame=data.frame(fineCapped,label='capped')
colnames(frame)=c(cnames)


againstAll=rbind(frame,allFrame)
colnames(againstAll)=c(cnames)

(p=ggplot(againstAll, aes(label,Distance))
 +geom_violin(#position = position_jitter(height = 0.1)
   )
 )



n=0
allMatched=vector(length=500*2016)
for (i in 0:499){
  matchedExp=read.table(paste('matchedDistr/repeat',i,sep=''))
  matchedExp=as.vector(matchedExp)
  
  allMatched[((i+1)*2016-2015):((i+1)*2016)]=matchedExp$V1

    result=wilcox.test(fineCapped,matchedExp$V1,alternative='less')
    if (result$p.value<0.05){
      n=n+1
    }
  if (result$p.value>=0.05){
    print(i)
  }
    
  
}

matchedExp=read.table(paste('matchedDistr/repeat',i,sep=''))

exampleFrame=data.frame(matchedExp,label='GCmatched')
fullFrame=data.frame(allMatched,label='GCmatched')
colnames(exampleFrame)=c(cnames)
colnames(fullFrame)=c(cnames)
colnames(frame)=c(cnames)
(p=ggplot(rbind(frame,exampleFrame), aes(label,Distance))
 +geom_violin(#position = position_jitter(height = 0.1)
 )
)

(p=ggplot(rbind(frame,fullFrame), aes(label,Distance))
 +geom_violin(#position = position_jitter(height = 0.1)
 )
)

mean(frame$Distance)

mean(fullFrame$Distance)

mlv(fullFrame$Distance, method = "mfv")
#from matched
wilcox.test(fineCapped,allMatched,alternative='less')
#from all
wilcox.test(fineCapped,allDist$V1,alternative='less')

#######################################################################
#for distances from mouse mir
cappedDist2=read.table('humanCappedToMouseCapped')
fineCapped2=cappedDist2$V1
1216
n=0
allMatched2=vector(length=500*1216)
for (i in 0:499){
  matchedExp=read.table(paste('matchedDistanceToMirCap/',i,sep=''))
  matchedExp=as.vector(matchedExp)
  
  allMatched2[((i+1)*1216-1215):((i+1)*1216)]=matchedExp$V1
  
  result=wilcox.test(fineCapped2,matchedExp$V1,alternative='less')
  if (result$p.value<0.05){
    n=n+1
  }
  if (result$p.value>=0.05){
    print(i)
  }
  
  
}

result=wilcox.test(fineCapped2,allMatched2,alternative='less')

result$p.value<0.05


(p=ggplot(as.data.frame(Distance=c(fineCapped2,allMatched2),label=(rep(fineCapped))), aes(label,Distance))
 +geom_violin(#position = position_jitter(height = 0.1)
 )
)


plot(density(allMatched2))
plot(density(allMatched))

################################################
#random stuff

n=0
allRandom=vector(length=1000*171)
for (i in 0:999){
  randomExp=read.table(paste('randomDistr/repeat',i,sep=''))
  randomExp=t(as.vector(matchedExp))
  
  allRandom[((i+1)*171-170):((i+1)*171)]=randomExp
  if (mean(matchedExp)>mean(fineCapped)){
    result=t.test(fineCapped)
    if (result$p.value<0.05){
      n=n+1  
    }
    if (result$p.value>=0.05){
      
    }
    
  }
}
randomFrame=data.frame(allRandom,label='random')
colnames(randomFrame)=cnames
wilcox.test(fineCapped,allRandom,alternative='less')

(p=ggplot(rbind(frame,randomFrame), aes(label,Distance))
 +geom_violin(#position = position_jitter(height = 0.1)
 )
)
wilcox.test(fineCapped,allRandom,alternative='less')




mouseDist=read.table('mouseDistDistr')
plot(density(mouseDist$V3))
