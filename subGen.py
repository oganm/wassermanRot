for i in range(0,100):
    if i<10:
        i='0'+str(i)
    elif i>=10:
        i=str(i)

    string=\
          "#!/bin/sh\n\
#$ -N cufflinks"+i+"\n#$ -cwd\n\
#$ -j y\n\
#$ -S /bin/bash\n\
cuffAdd=/raid2/local/exome/software/cufflinks-2.1.1.Linux_x86_64/cufflinks\n\
echo $cuffAdd\n\
echo "+i+"\n\
$cuffAdd -o ./cuffOut/100"+i+" -b  /home/ogan/humangenome.fa /raid6/ogan/RNASeq_Tuxedo3_Cufflink2_Mapping.RDhi100"+i+".tophat2.accepted_hits.bam"
    f=open('cuffSubs/cuff'+i,'w')
    f.write(string)
    f.close()
    

    
