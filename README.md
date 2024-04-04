# running-wf-human-variation  

**Request pod5 files via Mediaflux Download Share (contact Josh Reid)**  
Email will provide a token string. Run following commands in a vast scratch directory  
```
module load mediaflux-data-mover
mediaflux-data-mover -download <token> ./
tar -xf pod5_pass.tar
```

_**Note: Complete steps below if fast5 files (legacy format) are received instead of pod5**_
```
pod5 convert fast5 ./input/*.fast5 --output converted.pod5
```

**Basecalling (requires GPU)**  
```
screen -S <session-name>
srun --partition=gpuq -n1 -c6 --mem=100GB --gres=gpu:A30:1 --pty bash
module load dorado/0.5.2
dorado basecaller /stornext/System/data/nvidia/dorado/models/<model version>dna_r9.4.1_e8_sup@v3.6 --reference /vast/projects/bahlo_epilepsy/ref_genomes/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna <pod5_pass> > <sample id_model version_aligned>.bam
samtools sort aligned.bam > sorted.bam
```

**Running wf-human-variation**  
```
/stornext/System/data/tools/nextflow/nextflow-23.04.2/nextflow-23.04.2-all run /home/users/allstaff/reid.j/bahlo_reidj/analysis/wehi-wf-human-variation/wehi-wf-human-variation -profile apptainer -w /vast/scratch/users/reid.j/wf-human-variation/workspace --snp --sv --str --cnv --bam <sorted.bam> --ref /vast/projects/bahlo_epilepsy/ref_genomes/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna --sample_name <sample id_model version> -with-report --basecaller_cfg <insert version> --bam_min_coverage 5 -resume
```


**Common Issues:**  
apptainer pull  

bam_ingress:minimap2_alignment  
samtools sort: truncated file. Aborting  

Increase process memory!  




# running nf-core sarek  
```
/stornext/System/data/tools/nextflow/nextflow-23.10.0/nextflow-23.10.0-all run nf-core/sarek -revision 3.4.0 -profile wehi --tools freebayes,mutect2,strelka,manta,cnvkit,mpileup,deepvariant --outdir ./BGI_240222 --input /stornext/Bioinf/data/lab_bahlo/users/reid.j/nf-sarek/20240222_bgi_wes.csv  --fasta /stornext/Bioinf/data/lab_bahlo/projects/epilepsy/hg38/reference/fasta/Homo_sapiens_assembly38.fasta --wes --email joshua.reid@unimelb.edu.au -c custom.config --save_mapped --save_output_as_bam --joint_mutect2 -c /stornext/Bioinf/data/lab_bahlo/users/reid.j/nf-sarek/amazon.config -resume
```
