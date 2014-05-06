setwd('C:/Users/Ogan/Desktop/PORTAL/For R/')

miRNA=read.table('datanorm.txt')
miRNA2=read.table('datanorm.txt')
sides=read.table('sides.txt')
rnalist=read.table('intersect_pre-mirANDdpi_maturemirna_list.txt')

names=vector(length=nrow(miRNA))
noName=0
for (i in 1:nrow(miRNA)){
  s=unlist(strsplit(rownames(miRNA)[i],'-'))
  s=unlist(strsplit(s[length(s)],'_'))
  s[1]
  
  url=paste ('https://www.lifetechnologies.com/order/genome-database/browse/mirna/keyword/',s[1],'?ICID=uc-mirna-',s[1], sep ="", collapse = NULL)
  lines=readLines(url)
  code=paste(lines,collapse='\n')
  m=regexpr(".*?Homo sapiens.*?</td>\n.*\n.*\n.*\n.*\n.*?/", code, perl=TRUE)
  miRName=regmatches(code, m)
  m=regexpr("hsa-.*?[;|<]",miRName,perl=TRUE)
  #m=regexpr(">.*?</$", miRName, perl=TRUE)
  name=regmatches(miRName, m)
  name=substring(name,1,nchar(name)-1)
  if (length(name)==0){
    name=paste("NN",noName)
    noName=noName+1
  }
  print(name)
  print(i)
  names[i]=name
}
m=regexpr("RNU.*", rownames(miRNA), perl=TRUE)
aa=regmatches(rownames(miRNA),m)
controls=which(rownames(miRNA) %in% aa)
controlV=apply(miRNA[controls,],2,FUN=mean)

rownames(miRNA)=names

miRNA=miRNA[order(rownames(miRNA)),]

miRNA=miRNA-controlV


miRNA=cbind(miRNA,LMBchange=miRNA$LMB-miRNA$untreated,PHAXchange=miRNA$PHAX_rnai-miRNA$untreated)
miRNA=cbind(miRNA,LMBfold=2^-miRNA$LMBchange,PHAXfold=2^-miRNA$PHAXchange)

miRNA3p=miRNA[which(grepl('3p',rownames(miRNA))),]

miRNA5p=miRNA[which(grepl('5p',rownames(miRNA))),]

p3=gsub("-3p","",rownames(miRNA3p))
p5=gsub("-5p","",rownames(miRNA5p))
rownames(miRNA3p)=p3
rownames(miRNA5p)=p5

#5p 
miRNA5pC=miRNA5p[-which(match(p5,p3) %in% NA),]
miRNA3pC=miRNA3p[-which(match(p3,p5) %in% NA),]

LMBratio=miRNA3pC$LMBfold/miRNA5pC$LMBfold
PHAXratio=miRNA3pC$PHAXfold/miRNA5pC$PHAXfold

versus=as.data.frame(cbind(LMBratio,PHAXratio))
rownames(versus)=rownames(miRNA3pC)

versus[with(versus,order(LMBratio)),]
versus[with(versus,order(PHAXratio)),]


rnalist=read.table('intersect_pre-mirANDdpi_maturemirna_list.txt')
temp=gsub("-3p","",as.character(t(rnalist)))
temp=gsub("-5p","",temp)

interesting=versus[(versus$LMBratio<1)|(versus$PHAXratio<1),]
evenMoreInterestin=interesting[rownames(interesting) %in% temp,]

write.table(gsub('-','_',temp),file=CAGEnames,quote=FALSE)


rownames(interesting) %in% temp


just3p=miRNA[which(rownames(miRNA) %in% sides[which(sides$V3 %in% '3p'),]$V1),]

interestingJ3P=just3p[just3p$LMBfold<1,]
evenMoreInterestinJ3p=interestingJ3P[rownames(interestingJ3P) %in% temp,]

just3p2=miRNA3p[which(match(p3,p5) %in% NA),]
interestingJ3P2=just3p2[just3p2$LMBfold<1,]
evenMoreInterestinJ3p2=interestingJ3P2[rownames(interestingJ3P2) %in% temp,]


#filtering ####################################################################
miRNA2=read.table('datanorm.txt')
rownames(miRNA2)=names
miRNA
miRNA2=miRNA2[order(rownames(miRNA2)),]

miRNA2$LMB[miRNA2$LMB>35]=40
miRNA2$untreated[miRNA2$untreated>35]=40
miRNA2$PHAX_rnai[miRNA2$PHAX_rnai>35]=40
miRNA2=miRNA2[order(rownames(miRNA2)),]

rnalist=read.table('hitpremir')
temp=data.frame(lapply(rnalist, as.character), stringsAsFactors=FALSE)
temp=temp$V1
rownames(miRNA2)=tolower(rownames(miRNA2))
#miRNA2=cbind(miRNA2,LMBchange=miRNA2$LMB-miRNA2$untreated,PHAXchange=miRNA2$PHAX_rnai-miRNA2$untreated)
#miRNA2=cbind(miRNA2,LMBfold=2^-miRNA2$LMBchange,PHAXfold=2^-miRNA2$PHAXchange)

miRNA2=miRNA2-controlV


miRNA2f=cbind(miRNA2,LMBchange=miRNA2$LMB-miRNA2$untreated,PHAXchange=miRNA2$PHAX_rnai-miRNA2$untreated)
miRNA2f=cbind(miRNA2f,LMBfold=2^-miRNA2f$LMBchange,PHAXfold=2^-miRNA2f$PHAXchange)

miRNA3pf=miRNA2f[which(grepl('3p',rownames(miRNA2f))),]
miRNA5pf=miRNA2f[which(grepl('5p',rownames(miRNA2f))),]


p3f=gsub("-3p","",rownames(miRNA3pf))
p5f=gsub("-5p","",rownames(miRNA5pf))
rownames(miRNA3pf)=p3f
rownames(miRNA5pf)=p5f

#5p 
miRNA5pCf=miRNA5pf[-which(match(p5f,p3f) %in% NA),]
miRNA3pCf=miRNA3pf[-which(match(p3f,p5f) %in% NA),]




LMBratiof=miRNA3pCf$LMBfold/miRNA5pCf$LMBfold
PHAXratiof=miRNA3pCf$PHAXfold/miRNA5pCf$PHAXfold

versusf=as.data.frame(cbind(LMBratiof,PHAXratiof))
rownames(versusf)=rownames(miRNA3pCf)

versusf[with(versusf,order(LMBratiof)),]
versusf[with(versusf,order(PHAXratiof)),]
versusf[(!is.nan(versusf$LMBratiof))|(versusf$LMBratiof!=Inf),]

interestingf=versusf[(versusf$LMBratio<0.8)|(versusf$PHAXratio<0.8),]
evenMoreInterestinf=interestingf[tolower(rownames(interestingf)) %in% temp,]
sum(tolower(rownames(miRNA3pCf)) %in% temp)


phyper(nrow(evenMoreInterestin),sum(tolower(rownames(miRNA3pC)) %in% temp),nrow(miRNA3pC)-sum(tolower(rownames(miRNA3pC)) %in% temp),nrow(interesting))

phyper(nrow(evenMoreInterestinf),sum(rownames(miRNA3pCf) %in% temp),nrow(miRNA3pCf)-sum(rownames(miRNA3pCf) %in% temp),nrow(interestingf))


