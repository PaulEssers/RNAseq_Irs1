#nice tophat2 -o tophat_fusion/wt20038/ -p 2 --fusion-search --keep-fasta-order --bowtie1 --no-coverage-search  ~/references/mm10/Mus_musculus/UCSC/mm10/Sequence/BowtieIndex/genome fastq/wt20038_unmapped.fastq
#nice tophat2 -o tophat_fusion/wt20040/ -p 2 --fusion-search --keep-fasta-order --bowtie1 --no-coverage-search  ~/references/mm10/Mus_musculus/UCSC/mm10/Sequence/BowtieIndex/genome fastq/wt20040_unmapped.fastq
nice tophat2 -o tophat_fusion/wt20041/ -p 16 --fusion-search --keep-fasta-order --bowtie1 --no-coverage-search  ~/references/mm10/Mus_musculus/UCSC/mm10/Sequence/BowtieIndex/genome fastq/wt20041_unmapped.fastq
nice tophat2 -o tophat_fusion/wt20042/ -p 16 --fusion-search --keep-fasta-order --bowtie1 --no-coverage-search  ~/references/mm10/Mus_musculus/UCSC/mm10/Sequence/BowtieIndex/genome fastq/wt20042_unmapped.fastq
nice tophat2 -o tophat_fusion/irs1ko619/ -p 16 --fusion-search --keep-fasta-order --bowtie1 --no-coverage-search  ~/references/mm10/Mus_musculus/UCSC/mm10/Sequence/BowtieIndex/genome fastq/irs1ko619_unmapped.fastq
nice tophat2 -o tophat_fusion/irs1ko620/ -p 16 --fusion-search --keep-fasta-order --bowtie1 --no-coverage-search  ~/references/mm10/Mus_musculus/UCSC/mm10/Sequence/BowtieIndex/genome fastq/irs1ko620_unmapped.fastq
nice tophat2 -o tophat_fusion/irs1ko621/ -p 16 --fusion-search --keep-fasta-order --bowtie1 --no-coverage-search  ~/references/mm10/Mus_musculus/UCSC/mm10/Sequence/BowtieIndex/genome fastq/irs1ko621_unmapped.fastq

#python ~/../../software/CIRCexplorer-1.1.1/CIRCexplorer.py -f tophat_fusion/wt20038/accepted_hits.bam -g ~/references/mm10/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa -r ~/references/mm10/Mus_musculus/UCSC/mm10/Annotation/Genes/refFlat.txt.gz
python ~/../../software/CIRCexplorer-1.1.1/CIRCexplorer.py -f tophat_fusion/wt20040/accepted_hits.bam -g ~/references/mm10/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa -r ~/references/mm10/Mus_musculus/UCSC/mm10/Annotation/Genes/refFlat.txt
mv CIRCexplorer_circ.txt wt20040_circ.txt
python ~/../../software/CIRCexplorer-1.1.1/CIRCexplorer.py -f tophat_fusion/wt20041/accepted_hits.bam -g ~/references/mm10/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa -r ~/references/mm10/Mus_musculus/UCSC/mm10/Annotation/Genes/refFlat.txt
mv CIRCexplorer_circ.txt wt20041_circ.txt
python ~/../../software/CIRCexplorer-1.1.1/CIRCexplorer.py -f tophat_fusion/wt20042/accepted_hits.bam -g ~/references/mm10/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa -r ~/references/mm10/Mus_musculus/UCSC/mm10/Annotation/Genes/refFlat.txt
mv CIRCexplorer_circ.txt wt20042_circ.txt
python ~/../../software/CIRCexplorer-1.1.1/CIRCexplorer.py -f tophat_fusion/irs1ko619/accepted_hits.bam -g ~/references/mm10/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa -r ~/references/mm10/Mus_musculus/UCSC/mm10/Annotation/Genes/refFlat.txt
mv CIRCexplorer_circ.txt irs1ko619_circ.txt
python ~/../../software/CIRCexplorer-1.1.1/CIRCexplorer.py -f tophat_fusion/irs1ko620/accepted_hits.bam -g ~/references/mm10/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa -r ~/references/mm10/Mus_musculus/UCSC/mm10/Annotation/Genes/refFlat.txt
mv CIRCexplorer_circ.txt irs1ko620_circ.txt
python ~/../../software/CIRCexplorer-1.1.1/CIRCexplorer.py -f tophat_fusion/irs1ko621/accepted_hits.bam -g ~/references/mm10/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa -r ~/references/mm10/Mus_musculus/UCSC/mm10/Annotation/Genes/refFlat.txt
mv CIRCexplorer_circ.txt irs1ko621_circ.txt
