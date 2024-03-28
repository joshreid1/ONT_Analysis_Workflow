# running-wf-human-variation

**fast5 to pod5**
pod5 convert fast5 ./input/*.fast5 --output converted.pod5


**Basecalling (required GPU)**
srun --partition=gpuq -n1 -c6 --mem=100GB --gres=gpu:A30:1 --pty bash
module load dorado/0.5.2
dorado basecaller /stornext/System/data/nvidia/dorado/models/dna_r9.4.1_e8_sup@v3.6 --reference /vast/projects/bahlo_ukbiobank/GRCh38_full_analysis_set_plus_decoy_hla.fa <pod5> > T21583_dorado_v0.5.2_dna_r9.4.1_e8_sup@v3.6_aligned.bam


samtools sort aligned_bam > sorted_bam


/stornext/System/data/tools/nextflow/nextflow-23.04.2/nextflow-23.04.2-all run /home/users/allstaff/reid.j/bahlo_reidj/analysis/wehi-wf-human-variation/wehi-wf-human-variation -profile apptainer -w /vast/scratch/users/reid.j/wf-human-variation/workspace --snp --sv --str --cnv --bam /vast/scratch/users/reid.j/pod5/T21583_dorado_v0.5.2_dna_r9.4.1_e8_sup@v3.6_sorted.bam --ref /vast/projects/bahlo_epilepsy/ref_genomes/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna --sample_name T21583_dorado_dna_r9.4.1_e8_supv3.6 -with-report --basecaller_cfg dna_r9.4.1_e8_sup@v3.3 --bam_min_coverage 5 -resume



**Common Issues:**
apptainer pull

bam_ingress:minimap2_alignment
samtools sort: truncated file. Aborting

Increase process memory!




# running nf-core sarek
/stornext/System/data/tools/nextflow/nextflow-23.10.0/nextflow-23.10.0-all run nf-core/sarek -revision 3.4.0 -profile wehi --tools freebayes,mutect2,strelka,manta,cnvkit,mpileup,deepvariant --outdir ./BGI_240222 --input /stornext/Bioinf/data/lab_bahlo/users/reid.j/nf-sarek/20240222_bgi_wes.csv  --fasta /stornext/Bioinf/data/lab_bahlo/projects/epilepsy/hg38/reference/fasta/Homo_sapiens_assembly38.fasta --wes --email joshua.reid@unimelb.edu.au -c custom.config --save_mapped --save_output_as_bam --joint_mutect2 -c /stornext/Bioinf/data/lab_bahlo/users/reid.j/nf-sarek/amazon.config -resume

