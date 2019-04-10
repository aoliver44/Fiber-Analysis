#!/bin/bash


BASEDIR=/dfs3/bio/aoliver2/fiber_metagenome/plus_new_data/
QUALITY=/dfs3/bio/aoliver2/fiber_metagenome/plus_new_data/quality_filtered/
CROSSASSEMBLY=/dfs3/bio/aoliver2/fiber_metagenome/plus_new_data/megahit/final.contigs.fa
ANVIO=/dfs3/bio/aoliver2/fiber_metagenome/plus_new_data/anvio/
MERGED=/dfs3/bio/aoliver2/fiber_metagenome/plus_new_data/quality_filtered/merged/


## if you need to build the bowtie database again make sure you do!
module load bowtie2
bowtie2-build ${CROSSASSEMBLY} ${ANVIO}contigs-fixed

## make bam files for each individual (high - low) for anvio profile db
## ap = anvio profile

## For merging the decon fasta files, do this:

cd ${QUALITY}
counter=0
for f in *_qc_decon_pe.1.fastq.gz; do
sample=$(basename $f _qc_decon_pe.1.fastq.gz)

echo "merging sample ${counter} out of $(ls -l *_qc_decon_pe.1.fastq.gz | wc -l) samples."
echo "hang tight"

gunzip -c ${sample}_qc_decon* | cat - >> ${MERGED}${sample}_qc_merged.fastq
counter=$((counter+1))
done

while read ind fiber readgroup;
do
echo "#!/bin/bash
#$ -N fiber_${fiber}.ap
#$ -q bio*,abio*,pub*
#$ -pe openmp 4
#$ -R y
#$ -j y
#$ -ckpt restart

module load bowtie2/2.2.7
module load samtools/1.3


cd ${QUALITY}

bowtie2 -x ${ANVIO}contigs-fixed --very-fast-local -k 1 -t -p 4 --reorder --mm \
-U ${MERGED}${readgroup} | \
samtools view -S -bh -T ${CROSSASSEMBLY} - > ${ANVIO}${fiber}.bowtie2.bam

" > ${ANVIO}/ap.${fiber}.sh

qsub ${ANVIO}/ap.${fiber}.sh

done < ${ANVIO}/sample_map.txt

