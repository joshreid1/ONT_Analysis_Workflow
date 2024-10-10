#!/bin/bash
#SBATCH --partition gpuq     # submit to the gpuq partition
#SBATCH --cpus-per-task 6   # request 6 CPUs
#SBATCH --mem 64G         # request 64GB
#SBATCH --gres gpu:A30:4     # requesting 4 x A30 GPU
#SBATCH --job-name <sample_basecalling>

module load dorado/0.7.3

dorado basecaller /stornext/System/data/nvidia/dorado/models/dna_r10.4.1_e8.2_400bps_sup@v5.0.0 \
--modified-bases 5mCG_5hmCG --reference /vast/projects/bahlo_epilepsy/ref_genomes/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna \
./path_to_all_pod5 > <sample>_sup_v5.0.0_5mCG_5hmCG_aligned.bam

samtools sort -@ 6 -o <sample>_sup_v5.0.0_5mCG_5hmCG_sorted.bam <sample>_sup_v5.0.0_5mCG_5hmCG_aligned.bam
