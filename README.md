# ONT Analysis Workflow  
Start point = 
A) Raw pod5 files (basecalling required)
B) Basecalled bam files (alignment required)

_Please note: The P2i does not run the most up to date version of Dorado. Therefore, standalone Dorado basecalling on HPC (start point A) recommended over proceeding with live basecalling output (start point B)_  

**1) Request and download files via _Mediaflux Download Share_ (contact Josh Reid)**  
Email will provide a token string. Run following commands in a vast scratch directory  
```
module load mediaflux-data-mover
mediaflux-data-mover -download -exists skip <token> ./
```
_**Note: Request all pod5 files (pass, fail and recovered (if applicable) when basecalling)**_ 

_**Note: Complete steps below if fast5 files (legacy format) are received instead of pod5**_
```
pod5 convert fast5 ./input/*.fast5 --output converted.pod5
```

**2) Dorado (requires GPU)**
_See [Dorado](https://github.com/nanoporetech/dorado?tab=readme-ov-file#dna-models) link for available DNA models_  
> Current sup models (as at 03/06/2024):  
> LSK114 = dna_r10.4.1_e8.2_400bps_sup@v5.0.0
> LSK110 = dna_r9.4.1_e8_sup@v3.6  

_Note: See [link here](https://lh3.github.io/2017/11/13/which-human-reference-genome-to-use) for notes on recommended reference sequence_

_Option A: Run Basecalling_

Submit bash script below:
```
#!/bin/bash
#SBATCH --partition gpuq     # submit to the gpuq partition
#SBATCH --cpus-per-task 6   # request 6 CPUs
#SBATCH --mem 64G         # request 64GB
#SBATCH --gres gpu:A30:4     # requesting 4 x A30 GPU
#SBATCH --job-name <sample_basecalling>

module load dorado

dorado basecaller /stornext/System/data/nvidia/dorado/models/dna_r10.4.1_e8.2_400bps_sup@v5.0.0 --min-qscore 10 \
--modified-bases 5mCG_5hmCG --reference /vast/projects/bahlo_epilepsy/ref_genomes/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna \
./path_to_all_pod5 > <sample>_sup_v5.0.0_5mCG_5hmCG_aligned.bam
```
_Note: Approx real runtime for basecalling:_
- 5mC_5hmC (all methylated C's_ ~80 hours
- 5mCG_5hmCG (methylated 'CG' C's only) ~40 hours

_Note: 4 x A30 GPUs requested per job. 8 total available per user at once (i.e. 2 samples can be run simultaneuously)_

_Note: If basecalling is interrupted, command can be resumed by adding ```--resume-from <incomplete.bam>``` command_

_Option B: Run alignment_

Submit bash script below:
```
#!/bin/bash
#SBATCH --cpus-per-task 6   # request 6 CPUs
#SBATCH --mem 64G         # request 64GB
#SBATCH --job-name <sample_alignment>

module load dorado

samtools merge -@6 bam_pass/*bam -o <sample>_sup_v5.0.0_5mCG_5hmCG_unaligned.bam

dorado aligner /vast/projects/bahlo_epilepsy/ref_genomes/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna <sample>_sup_v5.0.0_5mCG_5hmCG_unaligned.bam > <sample>_sup_v5.0.0_5mCG_5hmCG_aligned.bam
```

Note: Alignment only does not require GPU

**3) Run epi2me wf-human-variation**  
_See details [here](https://github.com/epi2me-labs/wf-human-variation)_

https://github.com/epi2me-labs/wf-human-variation/releases
```
#!/bin/bash
#SBATCH --cpus-per-task 6   # request 6 CPUs
#SBATCH --mem 64G         # request 64GB
#SBATCH --job-name epi2me

samtools sort -@ 6 -o <sample>_sup_v5.0.0_5mCG_5hmCG_sorted.bam <sample>_sup_v5.0.0_5mCG_5hmCG_aligned.bam

nextflow run epi2me-labs/wf-human-variation -r v2.6.0 -c /vast/projects/reidj-project/nextflow-config/epi2me.config -w ./work --snp --sv --str --cnv --mod --phased --bam <sorted.bam> --ref /vast/projects/bahlo_epilepsy/ref_genomes/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna --sample_name <sample_name> -with-report --bam_min_coverage 5 -resume
```

**Common Issues:**  
> apptainer pull unsuccessful inside Nextflow
```
apptainer pull /vast/scratch/users/reid.j/nextflow/singularity_cache/ontresearch-wf-human-variation-snp-sha0d7e7e8e8207d9d23fdf50a34ceb577da364373e.img docker://ontresearch/wf-human-variation-snp:sha0d7e7e8e8207d9d23fdf50a34ceb577da364373e
apptainer pull /vast/scratch/users/reid.j/nextflow/singularity_cache/ontresearch-wf-cnv-sha428cb19e51370020ccf29ec2af4eead44c6a17c2.img docker://ontresearch/wf-cnv:sha428cb19e51370020ccf29ec2af4eead44c6a17c2
apptainer pull /vast/scratch/users/reid.j/nextflow/singularity_cache/ontresearch-wf-human-variation-str-sha28799bc3058fa256c01c1f07c87f04e4ade1fcc1.img docker://ontresearch/wf-human-variation-str:sha28799bc3058fa256c01c1f07c87f04e4ade1fcc1
apptainer pull /vast/scratch/users/reid.j/nextflow/singularity_cache/ontresearch-wf-human-variation-sha0800eade05e4cbb75d45421633c78c4f6320b2f6.img docker://ontresearch/wf-human-variation:sha0800eade05e4cbb75d45421633c78c4f6320b2f6
apptainer pull /vast/scratch/users/reid.j/nextflow/singularity_cache/ontresearch-wf-human-variation-sv-shabc3ac908a14705f248cdf49f218956ec33e93ef9.img docker://ontresearch/wf-human-variation-sv:shabc3ac908a14705f248cdf49f218956ec33e93ef9
apptainer pull /vast/scratch/users/reid.j/nextflow/singularity_cache/ontresearch-wf-common-sha0a6dc21fac17291f4acb2e0f67bcdec7bf63e6b7.img docker://ontresearch/wf-common:sha0a6dc21fac17291f4acb2e0f67bcdec7bf63e6b7
apptainer pull /vast/scratch/users/reid.j/nextflow/singularity_cache/ontresearch-snpeff-sha4f289afaf754c7a3e0b9ffb6c0b5be0f89a5cf04.img docker://ontresearch/snpeff:sha4f289afaf754c7a3e0b9ffb6c0b5be0f89a5cf04
```

>bam_ingress:minimap2_alignment  
Increase process memory!  

>samtools sort: truncated file. Aborting

**4) Additional structural variant annotation
```
conda activate svafotate-env
svafotate annotate --vcf <input.wf_sv.vcf.gz> --target /vast/projects/reidj-project/Genes4Epilepsy/EpilepsyGenes_NamesOnly_v2024-09.bed --out <output.vcf> \
-b /vast/projects/reidj-project/software/SVAFotate/svafotate_data/SVAFotate_core_SV_popAFs.GRCh38.v4.1.bed.gz
```




