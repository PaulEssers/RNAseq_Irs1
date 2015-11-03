java -jar ~/../../software/picard-1.130/picard.jar AddOrReplaceReadGroups I=wt20038_accepted.sorted.dedup.bam O=wt20038_accepted.sorted.dedup.gr.bam LB=lane6 PL=illumina PU=lane6 SM=lane6 >AddorReplaceReadGroups2clippeddata.log 2>&1
#java -jar ~/../../software/picard-1.130/picard.jar AddOrReplaceReadGroups I=wt20040_accepted.sorted.dedup.bam O=wt20040_accepted.sorted.dedup.gr.bam LB=lane6 PL=illumina PU=lane6 SM=lane6 >AddorReplaceReadGroups2clippeddata.log 2>&1
#java -jar ~/../../software/picard-1.130/picard.jar AddOrReplaceReadGroups I=wt20041_accepted.sorted.dedup.bam O=wt20041_accepted.sorted.dedup.gr.bam LB=lane6 PL=illumina PU=lane6 SM=lane6 >AddorReplaceReadGroups2clippeddata.log 2>&1
#java -jar ~/../../software/picard-1.130/picard.jar AddOrReplaceReadGroups I=wt20042_accepted.sorted.dedup.bam O=wt20042_accepted.sorted.dedup.gr.bam LB=lane6 PL=illumina PU=lane6 SM=lane6 >AddorReplaceReadGroups2clippeddata.log 2>&1
java -jar ~/../../software/picard-1.130/picard.jar AddOrReplaceReadGroups I=irs1ko619_accepted.sorted.dedup.bam O=irs1ko619_accepted.sorted.dedup.gr.bam LB=lane6 PL=illumina PU=lane6 SM=lane6 >AddorReplaceReadGroups2clippeddata.log 2>&1
java -jar ~/../../software/picard-1.130/picard.jar AddOrReplaceReadGroups I=irs1ko620_accepted.sorted.dedup.bam O=irs1ko620_accepted.sorted.dedup.gr.bam LB=lane6 PL=illumina PU=lane6 SM=lane6 >AddorReplaceReadGroups2clippeddata.log 2>&1
java -jar ~/../../software/picard-1.130/picard.jar AddOrReplaceReadGroups I=irs1ko621_accepted.sorted.dedup.bam O=irs1ko621_accepted.sorted.dedup.gr.bam LB=lane6 PL=illumina PU=lane6 SM=lane6 >AddorReplaceReadGroups2clippeddata.log 2>&1

#java -jar ~/../../software/picard-1.130/picard.jar ReorderSam I=wt20040_accepted.sorted.dedup.gr.bam O=wt20040_accepted.sorted.dedup.gr.reorder.bam R=~/references/GRCm38/Sequence/WholeGenomeFasta/genome.fa &
#java -jar ~/../../software/picard-1.130/picard.jar ReorderSam I=wt20041_accepted.sorted.dedup.gr.bam O=wt20041_accepted.sorted.dedup.gr.reorder.bam R=~/references/GRCm38/Sequence/WholeGenomeFasta/genome.fa &
#java -jar ~/../../software/picard-1.130/picard.jar ReorderSam I=wt20042_accepted.sorted.dedup.gr.bam O=wt20042_accepted.sorted.dedup.gr.reorder.bam R=~/references/GRCm38/Sequence/WholeGenomeFasta/genome.fa &
#java -jar ~/../../software/picard-1.130/picard.jar ReorderSam I=irs1ko619_accepted.sorted.dedup.gr.bam O=irs1ko619_accepted.sorted.dedup.gr.reorder.bam R=~/references/GRCm38/Sequence/WholeGenomeFasta/genome.fa &
#java -jar ~/../../software/picard-1.130/picard.jar ReorderSam I=irs1ko620_accepted.sorted.dedup.gr.bam O=irs1ko620_accepted.sorted.dedup.gr.reorder.bam R=~/references/GRCm38/Sequence/WholeGenomeFasta/genome.fa &
#java -jar ~/../../software/picard-1.130/picard.jar ReorderSam I=irs1ko621_accepted.sorted.dedup.gr.bam O=irs1ko621_accepted.sorted.dedup.gr.reorder.bam R=~/references/GRCm38/Sequence/WholeGenomeFasta/genome.fa &

#samtools index wt20038_accepted.sorted.dedup.gr.reorder.bam &
#samtools index wt20040_accepted.sorted.dedup.gr.reorder.bam &
#samtools index wt20041_accepted.sorted.dedup.gr.reorder.bam &
#samtools index wt20042_accepted.sorted.dedup.gr.reorder.bam &
#samtools index irs1ko619_accepted.sorted.dedup.gr.reorder.bam &
#samtools index irs1ko620_accepted.sorted.dedup.gr.reorder.bam &
#samtools index irs1ko621_accepted.sorted.dedup.gr.reorder.bam &

#echo -e "SampleID\tBamFile\tNotes\nwt20038\twt20038_accepted.sorted.dedup.gr.reorder.bam\tnone" > wt20038.txt
#echo -e "SampleID\tBamFile\tNotes\nwt20040\twt20040_accepted.sorted.dedup.gr.reorder.bam\tnone" > wt20040.txt
#echo -e "SampleID\tBamFile\tNotes\nwt20041\twt20041_accepted.sorted.dedup.gr.reorder.bam\tnone" > wt20041.txt
#echo -e "SampleID\tBamFile\tNotes\nwt20042\twt20042_accepted.sorted.dedup.gr.reorder.bam\tnone" > wt20042.txt
#echo -e "SampleID\tBamFile\tNotes\nirs1ko619\tirs1ko619_accepted.sorted.dedup.gr.reorder.bam\tnone" > irs1ko619.txt
#echo -e "SampleID\tBamFile\tNotes\nirs1ko620\tirs1ko620_accepted.sorted.dedup.gr.reorder.bam\tnone" > irs1ko620.txt
#echo -e "SampleID\tBamFile\tNotes\nirs1ko621\tirs1ko621_accepted.sorted.dedup.gr.reorder.bam\tnone" > irs1ko621.txt

#java -jar ~/../../software/RNA-SeQC_v1.1.8.jar -o ./wt20038_RNAseqQC -r ~/references/GRCm38/Sequence/WholeGenomeFasta/genome.fa -singleEnd -s wt20038.txt -t ~/references/GRCm38/Annotation/Genes/genes.gtf -ttype 2 &
#java -jar ~/../../software/RNA-SeQC_v1.1.8.jar -o ./wt20040_RNAseqQC -r ~/references/GRCm38/Sequence/WholeGenomeFasta/genome.fa -singleEnd -s wt20040.txt -t ~/references/GRCm38/Annotation/Genes/genes.gtf -ttype 2 &
#java -jar ~/../../software/RNA-SeQC_v1.1.8.jar -o ./wt20041_RNAseqQC -r ~/references/GRCm38/Sequence/WholeGenomeFasta/genome.fa -singleEnd -s wt20041.txt -t ~/references/GRCm38/Annotation/Genes/genes.gtf -ttype 2 &
#java -jar ~/../../software/RNA-SeQC_v1.1.8.jar -o ./wt20042_RNAseqQC -r ~/references/GRCm38/Sequence/WholeGenomeFasta/genome.fa -singleEnd -s wt20042.txt -t ~/references/GRCm38/Annotation/Genes/genes.gtf -ttype 2 &
#java -jar ~/../../software/RNA-SeQC_v1.1.8.jar -o ./irs1ko619_RNAseqQC -r ~/references/GRCm38/Sequence/WholeGenomeFasta/genome.fa -singleEnd -s irs1ko619.txt -t ~/references/GRCm38/Annotation/Genes/genes.gtf -ttype 2 &
#java -jar ~/../../software/RNA-SeQC_v1.1.8.jar -o ./irs1ko620_RNAseqQC -r ~/references/GRCm38/Sequence/WholeGenomeFasta/genome.fa -singleEnd -s irs1ko620.txt -t ~/references/GRCm38/Annotation/Genes/genes.gtf -ttype 2 &
#java -jar ~/../../software/RNA-SeQC_v1.1.8.jar -o ./irs1ko621_RNAseqQC -r ~/references/GRCm38/Sequence/WholeGenomeFasta/genome.fa -singleEnd -s irs1ko621.txt -t ~/references/GRCm38/Annotation/Genes/genes.gtf -ttype 2 &
