#/usr/bin
for file in *_unmapped.bam;
	do filename=`echo $file | cut -d "." -f 1`;
	bedtools bamtofastq -i $file -fq ../bwamem/${filename}.fastq;	
	done;
exit
