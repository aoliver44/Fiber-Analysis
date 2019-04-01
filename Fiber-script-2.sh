#!/bin/bash
#$ -N script2-fiber
#$ -q bio,abio*
#$ -pe openmp 16
#$ -ckpt restart
## Fiber script 2, for use after Full_fiber_basic.sh

BASEDIR=/dfs3/bio/aoliver2/fiber_metagenome/plus_new_data/
QUALITY=/dfs3/bio/aoliver2/fiber_metagenome/plus_new_data/quality_filtered/
MIDAS=/dfs3/bio/aoliver2/fiber_metagenome/plus_new_data/midas/
METAPHLAN=/dfs3/bio/aoliver2/fiber_metagenome/plus_new_data/metaphlan/
MEGAHIT=/dfs3/bio/aoliver2/fiber_metagenome/plus_new_data/megahit/
ANVIO=/dfs3/bio/aoliver2/fiber_metagenome/plus_new_data/anvio/

## Clean and merge the metaphlan data
cd ${METAPHLAN}
# mkdir bowtie_outs
# for f in *.fastq.bowtie2out.txt; do mv $f bowite_outs/${f}; done

python /dfs3/bio/aoliver2/merge_metaphlan.py *.txt > merged_fiber_metaphlan.txt

## Clean and merge the midas data

module purge
module load enthought_python/7.3.2
export PYTHONPATH=$PYTHONPATH:/data/users/wengland/software/MIDAS
export PATH=$PATH:/data/users/wengland/software/MIDAS/scripts

cd ${BASEDIR}
## Merge the species for MIDAS
/data/users/wengland/software/python3/bin/python3 /data/users/wengland/software/MIDAS/scripts/merge_midas.py species \
${MIDAS}MERGED \
-i ${MIDAS} \
-t dir \
-d /dfs3/bio/wengland/MIDAS-dbs/midas_db_v1.2/

## Merge the genes for MIDAS
/data/users/wengland/software/python3/bin/python3 /data/users/wengland/software/MIDAS/scripts/merge_midas.py genes \
${MIDAS}MERGED \
-i ${MIDAS} \
-t dir \
-d /dfs3/bio/wengland/MIDAS-dbs/midas_db_v1.2

## Assemble all the reads into contigs, a la how Dr. Alex Chase does it

module purge
module load megahit/1.1.1

# --kmin-1pass == memory efficient for ultra low-depth datasets, such as soil metagenomics data


READ1=${QUALITY}WF_10_qc_decon_pe.1.fastq.gz,${QUALITY}WF_11_qc_decon_pe.1.fastq.gz,${QUALITY}WF_12_qc_decon_pe.1.fastq.gz,${QUALITY}WF_13_qc_decon_pe.1.fastq.gz,${QUALITY}WF_14_qc_decon_pe.1.fastq.gz,${QUALITY}WF_15_qc_decon_pe.1.fastq.gz,${QUALITY}WF_16_qc_decon_pe.1.fastq.gz,${QUALITY}WF_17_qc_decon_pe.1.fastq.gz,${QUALITY}WF_18_qc_decon_pe.1.fastq.gz,${QUALITY}WF_19_qc_decon_pe.1.fastq.gz,${QUALITY}WF_1_qc_decon_pe.1.fastq.gz,${QUALITY}WF_201_qc_decon_pe.1.fastq.gz,${QUALITY}WF_202_qc_decon_pe.1.fastq.gz,${QUALITY}WF_203_qc_decon_pe.1.fastq.gz,${QUALITY}WF_204_qc_decon_pe.1.fastq.gz,${QUALITY}WF_20_qc_decon_pe.1.fastq.gz,${QUALITY}WF_21_qc_decon_pe.1.fastq.gz,${QUALITY}WF_22_qc_decon_pe.1.fastq.gz,${QUALITY}WF_23_qc_decon_pe.1.fastq.gz,${QUALITY}WF_24_qc_decon_pe.1.fastq.gz,${QUALITY}WF_25_qc_decon_pe.1.fastq.gz,${QUALITY}WF_26_qc_decon_pe.1.fastq.gz,${QUALITY}WF_27_qc_decon_pe.1.fastq.gz,${QUALITY}WF_28_qc_decon_pe.1.fastq.gz,${QUALITY}WF_29_qc_decon_pe.1.fastq.gz,${QUALITY}WF_2_qc_decon_pe.1.fastq.gz,${QUALITY}WF_30_qc_decon_pe.1.fastq.gz,${QUALITY}WF_31_qc_decon_pe.1.fastq.gz,${QUALITY}WF_32_qc_decon_pe.1.fastq.gz,${QUALITY}WF_33_qc_decon_pe.1.fastq.gz,${QUALITY}WF_34_qc_decon_pe.1.fastq.gz,${QUALITY}WF_35_qc_decon_pe.1.fastq.gz,${QUALITY}WF_36_qc_decon_pe.1.fastq.gz,${QUALITY}WF_37_qc_decon_pe.1.fastq.gz,${QUALITY}WF_38_qc_decon_pe.1.fastq.gz,${QUALITY}WF_39_qc_decon_pe.1.fastq.gz,${QUALITY}WF_3_qc_decon_pe.1.fastq.gz,${QUALITY}WF_40_qc_decon_pe.1.fastq.gz,${QUALITY}WF_41_qc_decon_pe.1.fastq.gz,${QUALITY}WF_42_qc_decon_pe.1.fastq.gz,${QUALITY}WF_43_qc_decon_pe.1.fastq.gz,${QUALITY}WF_44_qc_decon_pe.1.fastq.gz,${QUALITY}WF_45_qc_decon_pe.1.fastq.gz,${QUALITY}WF_46_qc_decon_pe.1.fastq.gz,${QUALITY}WF_47_qc_decon_pe.1.fastq.gz,${QUALITY}WF_48_qc_decon_pe.1.fastq.gz,${QUALITY}WF_49_qc_decon_pe.1.fastq.gz,${QUALITY}WF_4_qc_decon_pe.1.fastq.gz,${QUALITY}WF_50_qc_decon_pe.1.fastq.gz,${QUALITY}WF_51_qc_decon_pe.1.fastq.gz,${QUALITY}WF_52_qc_decon_pe.1.fastq.gz,${QUALITY}WF_53_qc_decon_pe.1.fastq.gz,${QUALITY}WF_54_qc_decon_pe.1.fastq.gz,${QUALITY}WF_55_qc_decon_pe.1.fastq.gz,${QUALITY}WF_56_qc_decon_pe.1.fastq.gz,${QUALITY}WF_57_qc_decon_pe.1.fastq.gz,${QUALITY}WF_58_qc_decon_pe.1.fastq.gz,${QUALITY}WF_59_qc_decon_pe.1.fastq.gz,${QUALITY}WF_5_qc_decon_pe.1.fastq.gz,${QUALITY}WF_60_qc_decon_pe.1.fastq.gz,${QUALITY}WF_61_qc_decon_pe.1.fastq.gz,${QUALITY}WF_62_qc_decon_pe.1.fastq.gz,${QUALITY}WF_63_qc_decon_pe.1.fastq.gz,${QUALITY}WF_64_qc_decon_pe.1.fastq.gz,${QUALITY}WF_65_qc_decon_pe.1.fastq.gz,${QUALITY}WF_66_qc_decon_pe.1.fastq.gz,${QUALITY}WF_67_qc_decon_pe.1.fastq.gz,${QUALITY}WF_68_qc_decon_pe.1.fastq.gz,${QUALITY}WF_69_qc_decon_pe.1.fastq.gz,${QUALITY}WF_6_qc_decon_pe.1.fastq.gz,${QUALITY}WF_70_qc_decon_pe.1.fastq.gz,${QUALITY}WF_71_qc_decon_pe.1.fastq.gz,${QUALITY}WF_72_qc_decon_pe.1.fastq.gz,${QUALITY}WF_73_qc_decon_pe.1.fastq.gz,${QUALITY}WF_74_qc_decon_pe.1.fastq.gz,${QUALITY}WF_75_qc_decon_pe.1.fastq.gz,${QUALITY}WF_76_qc_decon_pe.1.fastq.gz,${QUALITY}WF_77_qc_decon_pe.1.fastq.gz,${QUALITY}WF_78_qc_decon_pe.1.fastq.gz,${QUALITY}WF_79_qc_decon_pe.1.fastq.gz,${QUALITY}WF_7_qc_decon_pe.1.fastq.gz,${QUALITY}WF_80_qc_decon_pe.1.fastq.gz,${QUALITY}WF_81_qc_decon_pe.1.fastq.gz,${QUALITY}WF_82_qc_decon_pe.1.fastq.gz,${QUALITY}WF_83_qc_decon_pe.1.fastq.gz,${QUALITY}WF_84_qc_decon_pe.1.fastq.gz,${QUALITY}WF_85_qc_decon_pe.1.fastq.gz,${QUALITY}WF_86_qc_decon_pe.1.fastq.gz,${QUALITY}WF_87_qc_decon_pe.1.fastq.gz,${QUALITY}WF_88_qc_decon_pe.1.fastq.gz,${QUALITY}WF_89_qc_decon_pe.1.fastq.gz,${QUALITY}WF_8_qc_decon_pe.1.fastq.gz,${QUALITY}WF_90_qc_decon_pe.1.fastq.gz,${QUALITY}WF_91_qc_decon_pe.1.fastq.gz,${QUALITY}WF_92_qc_decon_pe.1.fastq.gz,${QUALITY}WF_93_qc_decon_pe.1.fastq.gz,${QUALITY}WF_94_qc_decon_pe.1.fastq.gz,${QUALITY}WF_95_qc_decon_pe.1.fastq.gz,${QUALITY}WF_96_qc_decon_pe.1.fastq.gz,${QUALITY}WF_97_qc_decon_pe.1.fastq.gz,${QUALITY}WF_9_qc_decon_pe.1.fastq.gz
READ2=${QUALITY}WF_10_qc_decon_pe.2.fastq.gz,${QUALITY}WF_11_qc_decon_pe.2.fastq.gz,${QUALITY}WF_12_qc_decon_pe.2.fastq.gz,${QUALITY}WF_13_qc_decon_pe.2.fastq.gz,${QUALITY}WF_14_qc_decon_pe.2.fastq.gz,${QUALITY}WF_15_qc_decon_pe.2.fastq.gz,${QUALITY}WF_16_qc_decon_pe.2.fastq.gz,${QUALITY}WF_17_qc_decon_pe.2.fastq.gz,${QUALITY}WF_18_qc_decon_pe.2.fastq.gz,${QUALITY}WF_19_qc_decon_pe.2.fastq.gz,${QUALITY}WF_1_qc_decon_pe.2.fastq.gz,${QUALITY}WF_201_qc_decon_pe.2.fastq.gz,${QUALITY}WF_202_qc_decon_pe.2.fastq.gz,${QUALITY}WF_203_qc_decon_pe.2.fastq.gz,${QUALITY}WF_204_qc_decon_pe.2.fastq.gz,${QUALITY}WF_20_qc_decon_pe.2.fastq.gz,${QUALITY}WF_21_qc_decon_pe.2.fastq.gz,${QUALITY}WF_22_qc_decon_pe.2.fastq.gz,${QUALITY}WF_23_qc_decon_pe.2.fastq.gz,${QUALITY}WF_24_qc_decon_pe.2.fastq.gz,${QUALITY}WF_25_qc_decon_pe.2.fastq.gz,${QUALITY}WF_26_qc_decon_pe.2.fastq.gz,${QUALITY}WF_27_qc_decon_pe.2.fastq.gz,${QUALITY}WF_28_qc_decon_pe.2.fastq.gz,${QUALITY}WF_29_qc_decon_pe.2.fastq.gz,${QUALITY}WF_2_qc_decon_pe.2.fastq.gz,${QUALITY}WF_30_qc_decon_pe.2.fastq.gz,${QUALITY}WF_31_qc_decon_pe.2.fastq.gz,${QUALITY}WF_32_qc_decon_pe.2.fastq.gz,${QUALITY}WF_33_qc_decon_pe.2.fastq.gz,${QUALITY}WF_34_qc_decon_pe.2.fastq.gz,${QUALITY}WF_35_qc_decon_pe.2.fastq.gz,${QUALITY}WF_36_qc_decon_pe.2.fastq.gz,${QUALITY}WF_37_qc_decon_pe.2.fastq.gz,${QUALITY}WF_38_qc_decon_pe.2.fastq.gz,${QUALITY}WF_39_qc_decon_pe.2.fastq.gz,${QUALITY}WF_3_qc_decon_pe.2.fastq.gz,${QUALITY}WF_40_qc_decon_pe.2.fastq.gz,${QUALITY}WF_41_qc_decon_pe.2.fastq.gz,${QUALITY}WF_42_qc_decon_pe.2.fastq.gz,${QUALITY}WF_43_qc_decon_pe.2.fastq.gz,${QUALITY}WF_44_qc_decon_pe.2.fastq.gz,${QUALITY}WF_45_qc_decon_pe.2.fastq.gz,${QUALITY}WF_46_qc_decon_pe.2.fastq.gz,${QUALITY}WF_47_qc_decon_pe.2.fastq.gz,${QUALITY}WF_48_qc_decon_pe.2.fastq.gz,${QUALITY}WF_49_qc_decon_pe.2.fastq.gz,${QUALITY}WF_4_qc_decon_pe.2.fastq.gz,${QUALITY}WF_50_qc_decon_pe.2.fastq.gz,${QUALITY}WF_51_qc_decon_pe.2.fastq.gz,${QUALITY}WF_52_qc_decon_pe.2.fastq.gz,${QUALITY}WF_53_qc_decon_pe.2.fastq.gz,${QUALITY}WF_54_qc_decon_pe.2.fastq.gz,${QUALITY}WF_55_qc_decon_pe.2.fastq.gz,${QUALITY}WF_56_qc_decon_pe.2.fastq.gz,${QUALITY}WF_57_qc_decon_pe.2.fastq.gz,${QUALITY}WF_58_qc_decon_pe.2.fastq.gz,${QUALITY}WF_59_qc_decon_pe.2.fastq.gz,${QUALITY}WF_5_qc_decon_pe.2.fastq.gz,${QUALITY}WF_60_qc_decon_pe.2.fastq.gz,${QUALITY}WF_61_qc_decon_pe.2.fastq.gz,${QUALITY}WF_62_qc_decon_pe.2.fastq.gz,${QUALITY}WF_63_qc_decon_pe.2.fastq.gz,${QUALITY}WF_64_qc_decon_pe.2.fastq.gz,${QUALITY}WF_65_qc_decon_pe.2.fastq.gz,${QUALITY}WF_66_qc_decon_pe.2.fastq.gz,${QUALITY}WF_67_qc_decon_pe.2.fastq.gz,${QUALITY}WF_68_qc_decon_pe.2.fastq.gz,${QUALITY}WF_69_qc_decon_pe.2.fastq.gz,${QUALITY}WF_6_qc_decon_pe.2.fastq.gz,${QUALITY}WF_70_qc_decon_pe.2.fastq.gz,${QUALITY}WF_71_qc_decon_pe.2.fastq.gz,${QUALITY}WF_72_qc_decon_pe.2.fastq.gz,${QUALITY}WF_73_qc_decon_pe.2.fastq.gz,${QUALITY}WF_74_qc_decon_pe.2.fastq.gz,${QUALITY}WF_75_qc_decon_pe.2.fastq.gz,${QUALITY}WF_76_qc_decon_pe.2.fastq.gz,${QUALITY}WF_77_qc_decon_pe.2.fastq.gz,${QUALITY}WF_78_qc_decon_pe.2.fastq.gz,${QUALITY}WF_79_qc_decon_pe.2.fastq.gz,${QUALITY}WF_7_qc_decon_pe.2.fastq.gz,${QUALITY}WF_80_qc_decon_pe.2.fastq.gz,${QUALITY}WF_81_qc_decon_pe.2.fastq.gz,${QUALITY}WF_82_qc_decon_pe.2.fastq.gz,${QUALITY}WF_83_qc_decon_pe.2.fastq.gz,${QUALITY}WF_84_qc_decon_pe.2.fastq.gz,${QUALITY}WF_85_qc_decon_pe.2.fastq.gz,${QUALITY}WF_86_qc_decon_pe.2.fastq.gz,${QUALITY}WF_87_qc_decon_pe.2.fastq.gz,${QUALITY}WF_88_qc_decon_pe.2.fastq.gz,${QUALITY}WF_89_qc_decon_pe.2.fastq.gz,${QUALITY}WF_8_qc_decon_pe.2.fastq.gz,${QUALITY}WF_90_qc_decon_pe.2.fastq.gz,${QUALITY}WF_91_qc_decon_pe.2.fastq.gz,${QUALITY}WF_92_qc_decon_pe.2.fastq.gz,${QUALITY}WF_93_qc_decon_pe.2.fastq.gz,${QUALITY}WF_94_qc_decon_pe.2.fastq.gz,${QUALITY}WF_95_qc_decon_pe.2.fastq.gz,${QUALITY}WF_96_qc_decon_pe.2.fastq.gz,${QUALITY}WF_97_qc_decon_pe.2.fastq.gz,${QUALITY}WF_9_qc_decon_pe.2.fastq.gz

megahit \
-1 $READ1 \
-2 $READ2 \
-t 16 \
--min-count 3 \
--k-list 31,41,51,61,71,81,91,95,101,105,111 \
--kmin-1pass \
--min-contig-len 1000 \
--memory 0.95 \
--out-dir ${MEGAHIT} \
--continue

