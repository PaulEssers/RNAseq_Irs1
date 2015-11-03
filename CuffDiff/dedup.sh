#samtools sort wt20040_accepted.bam wt20040_accepted.sorted.bam &
#samtools sort wt20041_accepted.bam wt20041_accepted.sorted.bam &
#samtools sort wt20042_accepted.bam wt20042_accepted.sorted.bam &
#samtools sort irs1ko619_accepted.bam irs1ko619_accepted.sorted.bam &
#samtools sort irs1ko620_accepted.bam irs1ko620_accepted.sorted.bam &
#samtools sort irs1ko621_accepted.bam irs1ko621_accepted.sorted.bam &

#java -jar  ~/../../software/picard-1.130/picard.jar MarkDuplicates INPUT=wt20040_accepted.sorted.bam OUTPUT=wt20040_accepted.sorted.dedup.bam METRICS_FILE=wt20040_dedup_metrics.txt REMOVE_DUPLICATES=true
java -jar  ~/../../software/picard-1.130/picard.jar MarkDuplicates INPUT=wt20041_accepted.sorted.bam OUTPUT=wt20041_accepted.sorted.dedup.bam METRICS_FILE=wt20041_dedup_metrics.txt REMOVE_DUPLICATES=true &
java -jar  ~/../../software/picard-1.130/picard.jar MarkDuplicates INPUT=wt20042_accepted.sorted.bam OUTPUT=wt20042_accepted.sorted.dedup.bam METRICS_FILE=wt20042_dedup_metrics.txt REMOVE_DUPLICATES=true &
java -jar  ~/../../software/picard-1.130/picard.jar MarkDuplicates INPUT=irs1ko619_accepted.sorted.bam OUTPUT=irs1ko619_accepted.sorted.dedup.bam METRICS_FILE=irs1ko619_dedup_metrics.txt REMOVE_DUPLICATES=true &
java -jar  ~/../../software/picard-1.130/picard.jar MarkDuplicates INPUT=irs1ko620_accepted.sorted.bam OUTPUT=irs1ko620_accepted.sorted.dedup.bam METRICS_FILE=irs1ko620_dedup_metrics.txt REMOVE_DUPLICATES=true &
java -jar  ~/../../software/picard-1.130/picard.jar MarkDuplicates INPUT=irs1ko621_accepted.sorted.bam OUTPUT=irs1ko621_accepted.sorted.dedup.bam METRICS_FILE=irs1ko621_dedup_metrics.txt REMOVE_DUPLICATES=true &
