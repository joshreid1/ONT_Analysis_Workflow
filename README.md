# running-wf-human-variation  

**1) Request and download pod5 files via _Mediaflux Download Share_ (contact Josh Reid)**  
Email will provide a token string. Run following commands in a vast scratch directory  
```
module load mediaflux-data-mover
mediaflux-data-mover -download <token> ./
jar xf pod5.zip
tar -xf pod5/pod5_pass.tar
```

_**Note: Complete steps below if fast5 files (legacy format) are received instead of pod5**_
```
pod5 convert fast5 ./input/*.fast5 --output converted.pod5
```

**2) Run Basecalling (requires GPU)**  
_Note: Run via interactive GPU node. Recommended to run in a screen session due to long run-time. See steps below_  
```
screen -S <session-name>
srun --partition=gpuq -n1 -c6 --mem=100GB --gres=gpu:A30:1 --pty bash
```
_See [Dorado](https://github.com/nanoporetech/dorado?tab=readme-ov-file#dna-models) link for available DNA models_  
> Current sup models (as at 04/04/2024):  
> LSK114 = dna_r10.4.1_e8.2_400bps_sup@v4.3.0  
> LSK110 = dna_r9.4.1_e8_sup@v3.6  
```
module load dorado/0.5.2
dorado basecaller /stornext/System/data/nvidia/dorado/models/<model version> --reference /vast/projects/bahlo_epilepsy/ref_genomes/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna <pod5_pass> > <sample id_model version_aligned>.bam
samtools sort aligned.bam > sorted.bam
```
  
_Note: If basecalling is interrupted, command can be resumed by adding ```--resume-from <incomplete.bam>``` command_

**3) Run wf-human-variation**  
```
/stornext/System/data/tools/nextflow/nextflow-23.04.2/nextflow-23.04.2-all run /home/users/allstaff/reid.j/bahlo_reidj/analysis/wehi-wf-human-variation/wehi-wf-human-variation -profile apptainer -w /vast/scratch/users/reid.j/wf-human-variation/workspace --snp --sv --str --cnv --bam <sorted.bam> --ref /vast/projects/bahlo_epilepsy/ref_genomes/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna --sample_name <sample id_model version> -with-report --basecaller_cfg <insert version> --bam_min_coverage 5 -resume
```


**Common Issues:**  
> apptainer pull unsuccessful inside Nextflow
```
apptainer pull /vast/scratch/users/reid.j/nextflow/singularity_cache/ontresearch-wf-human-variation-snp-sha0d7e7e8e8207d9d23fdf50a34ceb577da364373e.img docker://ontresearch/wf-human-variation-snp:sha0d7e7e8e8207d9d23fdf50a34ceb577da364373e
apptainer pull /vast/scratch/users/reid.j/nextflow/singularity_cache/ontresearch-wf-cnv-sha428cb19e51370020ccf29ec2af4eead44c6a17c2.img docker://ontresearch/wf-cnv:sha428cb19e51370020ccf29ec2af4eead44c6a17c2
apptainer pull /vast/scratch/users/reid.j/nextflow/singularity_cache/ontresearch-wf-human-variation-str-sha28799bc3058fa256c01c1f07c87f04e4ade1fcc1.img docker://ontresearch/wf-human-variation-str:sha28799bc3058fa256c01c1f07c87f04e4ade1fcc1
apptainer pull /vast/scratch/users/reid.j/nextflow/singularity_cache/ontresearch-wf-human-variation-sha0800eade05e4cbb75d45421633c78c4f6320b2f6.img docker://ontresearch/wf-human-variation:sha0800eade05e4cbb75d45421633c78c4f6320b2f6
apptainer pull /vast/scratch/users/reid.j/nextflow/singularity_cache/ontresearch-wf-human-variation-sv-shabc3ac908a14705f248cdf49f218956ec33e93ef9.img docker://ontresearch/wf-human-variation-sv:shabc3ac908a14705f248cdf49f218956ec33e93ef9
```


bam_ingress:minimap2_alignment  
samtools sort: truncated file. Aborting  

Increase process memory!  




# running nf-core sarek  
```
/stornext/System/data/tools/nextflow/nextflow-23.10.0/nextflow-23.10.0-all run nf-core/sarek -revision 3.4.0 -profile wehi --tools freebayes,mutect2,strelka,manta,cnvkit,mpileup,deepvariant --outdir ./BGI_240222 --input /stornext/Bioinf/data/lab_bahlo/users/reid.j/nf-sarek/20240222_bgi_wes.csv  --fasta /stornext/Bioinf/data/lab_bahlo/projects/epilepsy/hg38/reference/fasta/Homo_sapiens_assembly38.fasta --wes --email joshua.reid@unimelb.edu.au -c custom.config --save_mapped --save_output_as_bam --joint_mutect2 -c /stornext/Bioinf/data/lab_bahlo/users/reid.j/nf-sarek/amazon.config -resume
```
