Wasserman Lab Rotation
=====================
Analysis of human and mouse capped microRNA. For manipulation of files and non structured data python is used. R is used for advanced plotting functions and manipulation of structured data. Most of the scripts are directory dependent. Meaning that they will need modification of the code to work in different systems.

General Purpose Tools
--------------------
**SeqXXX.py:** They are used to give the sequence of a said interval from a fasta files. Used to get miRNA sequences of mouse and humans. H19 and mm10 assemblies are used for human and mice respectively. They require around 4gb of RAM per genome as it loads the whole sequence. Use `clearModule` to clean up the sequence files. It will be modified to use indexing. **Modification of the filename and address is necessary before further use (9th and 10th lines. Direct it to the location of the fasta file.).**. Syntax for `giveSeq` is `giveSeq(chro, start, end, strand = '+', string = False)`. `chro` is either a string with the name of the sequence in the fasta [eg. 'X', 'chr1'] or a number if sequences are named as chr1 chr2 etc. `start` is the **0 based** starting coordinates of the desired location while `end` is **1 based**. If string is False it outputs a piopython seq object. Otherwise the output is a plain string. Keep in mind that it defaults to False.


**ogbox.py:** Small re-occuring functions. `readFasta` outputs a list of biophyton seq objects from a fasta file. `flines` gives the number of lines in a file. `mirbaseGFF` reads a GFF file downloaded from mirbase and and returns an object that carries lists storing the relevant information about the file. The comments need to be deleted from the GFF before processing. The rest is self explanatory.

Analysis Tools
--------------
**TSSdistanceMirStart.py:** From mirStart database, takes miRNA precursors that are close to their transcription site. Written in python 3.3. Not reverse compatible.

**mirbaseGFF.py:** import seq imports an old version of the SeqXXX.py. Just change the name to the current version to reuse it. It also requires the gff file from mirbase. It invokes RNA fold to fold RNAs of unkown sides and output their correct side to sides.txt.

**subGen.py:** Used to generate variants of sMir.py executors in clustdell. It generates 100 files in its native form while there are only 70 files. Extra files do not cause problems in downstream analysis. A similar script was used to generate the commands for cufflinks but that script was also this script so it was overwritten.

**sMir.py:** Requires pybedtools. Responsible analysis of cufflinks results. Requires fileno of the RNAseq data from FANTOM as input. Takes in the transcript information sends it to cufflinks. It discludes exomic miRNAs from the analysis. The file exome.bed is downloaded from USCS table browser (USCSgenes-> knownGenes). After that it takes a 5 kb long window around miRNAs to look for transcripts. This requires human.hg19.genome file to be present the working directory. Alternatively, the file can be generated with the commands below.


`fout = open('human.hg19.genome','w')
chromdict = pybedtools.get_chromsizes_from_ucsc('hg19')
for chrom, size in chromdict.items():
fout.write("%s\t%s\n" % (chrom, size[1]))
fout.close()`


Required cufflinks output is at clustdell raid6/ogan/cuffLinksOutput. Inside the numbered folders you can find the transcript files. This output converted to a BED file by function `cuffToBed` defined inside the script.

Commented out sections are to include individual exons to the analysis but these are not used in the end as we are looking for miRNAs that are oustide exomic regions. 



**lifeCrawler.R:** For analysis of QPCR data. Gets the names of probes from Life Technologies. Looks for an increase in 3p/5p ratio.
 
 **randomDistance.py:** Using the output of RNAclust, it picks groups of pre-miRNA randomly and with matching GC content. Puts their distance information in a file. Repeats.
 
 **compareDistances.R:** Analyses the output of randomDistance.py. Plotting and wilcoxon rank sum test.
 
**locarn.R:** Using RNAclust distances and names, outputs 1% quantile of the ones closest to capped miRNA precursors. Alternating commented out sections are used to do it both to long and short versions. There is also a difference between selecting microRNA as few of the mouse miRNAs at the long set were actually of the normal size. This also decreased the number of unique results when selecting precursors.

**rnaDistQuantile.r:** Includes commented out code to create hierarchical clusters which did not provide valuable information. It takes %1 quantile of distances calculated by RNAdistance of ViennaRNA. RNAclust and ViennaRNA has different outputs hence there are two scripts who does the same thing. Also RNAdistance has a single output for both long and short sequences as inputs were given together. The input file justDistances was created manually.

**seedSequences.py:** Compares human miRNA seeds with capped mouse ones. Uses two seed definitions, 2-8, 2-7. Outputs 2 files listing the matching nucleotides for every seed. Commented out parts used to deal with comparing selected miRNAs from quantiles but this was unnecessary and removed. Its output is 1DSorted6 and 1DSorted7

**plotdistr.R:** Draws distribution plots of distances in between capped mouse miRNA precursors and and mouse capped to human miRNA precursors. Commented out sections are for alternating between long and short versions. It does not save the plot.

**startFinder.awk:**Finds start sites of microRNA precursors.

**mirListGen.py:** From a mirbase GFF file generates a list of all miRNAs along with their precursors. Generator of mirName&Origin

**precursorToMir.py:** Given a list of precursors, generates a list of all miRNAs.

**allInclusive.R:** Gathers all the data and forms the final table and Venn diagrams.

Files and Folders
----
These files and folders are placed inside files folder. They are generated by various scripts.

**Sides**: Sides designated to human microRNA based on their folding by RNAfold. Output of RNA fold is put into mirbaseGFF.py (python3.3)

**quantile.png:** venn diagram of microRNA that were found at the1 1% quantile with different methods. 

**quantileWithCage.png:** With caged ones added.

**RNAdistanceR:** Inputs and outputs of RNAdistance and RNAfold. hallMirVShortincludes all miRNA sequences.  


**mirName&Origin:** Generated by mirListGen.py

**preMirStart.bed:** generated by startFinder.awk. Start sites of precursor miRNAs. 

**selectionWithMouse(xxxx):** miRNA precursors that lie within the %1 range. Generated by rnaDistQuantile.r and locarn.r

**cageHits:** miRNAs that have intersecting CAGE reads. Generated by:

`intersectBed -b CAGE_dpi.bed -a preMirStart.bed -s -u`

**cageHits2:** CAGE sequences that coincide with miRNAs. Used to get their starting location. Generated by:

`intersectBed -b CAGE_dpi.bed -a preMirStart.bed -s -u >CAGEHits2`

**CAGE_dpi.bed:** Cage reads. Taken from Charles.

