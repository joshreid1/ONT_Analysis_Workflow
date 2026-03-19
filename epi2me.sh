#!/bin/bash
#SBATCH --cpus-per-task 6   # request 6 CPUs
#SBATCH --mem 64G         	# request 64GB
#SBATCH --job-name epi2me

samtools sort -@ 6 -o <sample>_sup_v5.0.0_5mCG_5hmCG_sorted.bam <sample>_sup_v5.0.0_5mCG_5hmCG_aligned.bam

nextflow run epi2me-labs/wf-human-variation -r v2.8.0 -c epi2me.config -w ./work \
 --snp --sv --str --cnv --mod --phased --bam <sorted.bam> \
 --ref /vast/projects/bahlo_epilepsy/ref_genomes/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna \
 --sample_name <sample_name> -with-report --bam_min_coverage 5 -resume