awk 'BEGIN {}
$6=="+" {print $1"\t"$2-30"\t"$2+15"\t"$4"\t"$5"\t"$6}
$6=="-" {print $1"\t"$3-15"\t"$3+30"\t"$4"\t"$5"\t"$6}' pre-miRNAs.bed >preMirStart.bed
