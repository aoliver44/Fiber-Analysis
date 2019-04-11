#!/bin/bash
#$ -N anvio
#$ -q bio,abio*,pub*
#$ -pe openmp 12
#$ -R y
#$ -ckpt restart

ANVIO=/dfs3/bio/aoliver2/fiber_metagenome/plus_new_data/anvio/

module load anaconda
source activate anvio5

module purge

cd ${ANVIO}

## The below is for generating the contigs database for Anvio5
anvi-gen-contigs-database -f contigs-fixed.fa -o contigs.db -n 'Contigs for fiber plus extra samples'

anvi-run-hmms -c contigs.db --num-threads 12

anvi-run-ncbi-cogs -c contigs.db --num-threads 12

anvi-get-sequences-for-gene-calls -c contigs.db -o gene_calls.fa


## Adding gene level taxonomy calls 
/dfs3/bio/aoliver2/database/kaiju/bin/kaiju \
-t /dfs3/bio/aoliver2/database/kaiju/kaijudb/nodes.dmp \
-f /dfs3/bio/aoliver2/database/kaiju/kaijudb/kaiju_db.fmi \
-i gene_calls.fa \
-o gene_calls_nr.out \
-z 12 \
-v

addTaxonNames -t /dfs3/bio/aoliver2/database/kaiju/kaijudb/nodes.dmp \
              -n /dfs3/bio/aoliver2/database/kaiju/kaijudb/names.dmp \
              -i gene_calls_nr.out \
              -o gene_calls_nr.names \
              -r superkingdom,phylum,order,class,family,genus,species
              
anvi-import-taxonomy-for-genes -i gene_calls_nr.names \
                               -c contigs.db \
                               -p kaiju --just-do-it


## making profiles for each sample

cd ${ANVIO}bams

for sample in *.bowtie2.bam; do
newname=$(basename $sample .bowtie2.bam)
anvi-init-bam ${sample} \
-o ${newname}_init.bam
done

mkdir anvio_profiles

for sample in *_init.bam; do
name=$(basename $sample _init.bam)
anvi-profile -i ${sample} -c contigs.db \
--skip-SNV-profiling \
--min-contig-length 2000 \
--sample-name ${name} \
--output-dir ${ANVIO}bams/anvio_profiles
done

## merge all the samples
cd ${ANVIO}bams/anvio_profiles

anvi-merge */PROFILE.db -o ${ANVIO}SAMPLES-MERGED -c ${ANVIO}contigs.db \
--skip-hierarchical-clustering

