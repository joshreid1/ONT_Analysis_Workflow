# ONT Analysis Workflow

A workflow for running ONT basecalling (via Dorado) and downstream analysis (via epi2me wf-human-variation).

## 🚀 Overview
This repository provides guidance and example `sbatch` scripts for:

- **Option A:** Basecalling from raw **pod5** files using **Dorado** (recommended)
- **Option B:** Aligning already basecalled **BAM** files using **Dorado Aligner**
- Running **epi2me wf-human-variation** (Nextflow) on the resulting aligned BAM

> **Note:** The PromethION 2 Integrated uses the Dorado server (essentially a wrapper on standalone Dorado) that often lags behind the latest releases. Versioning is different and can not be directly compared. We recommend basecalling on HPC (Option A) instead of using on-device output (Option B).

---

## ✅ Prerequisites

- Access to WEHI HPC with **GPU** access to run Dorado basecalling (Option A)
- `dorado` (basecaller + aligner), `samtools`, and `nextflow`
- Access to the reference genome(s) used for alignment (e.g. GRCh38)

---

## 1) Download data (Mediaflux)

Download raw **pod5** files via Mediaflux. For most workflows you should request **all pod5 files** (pass, fail, recovered) to ensure complete coverage.

Move all pod5 files to a single directory for processing. Example:

```bashbash
mkdir ./pod5_all
mv pod5_pass/* pod5_all/
mv pod5_fail/* pod5_all/
mv pod5_recovered/* pod5_all/ 

**If you are given legacy FAST5 files** instead of pod5, convert them:

```bash
pod5 convert fast5 ./input/*.fast5 --output converted.pod5
```

📄 Pod5 docs: https://pod5-file-format.readthedocs.io/en/0.1.21/docs/tools.html

---

## 2) Run Dorado (basecalling + optional alignment)

Dorado is the current ONT basecaller. It can also align reads to a reference in the same run (recommended for epi2me wf-human-variation).

### 🧠 Model notes
- Available models: https://software-docs.nanoporetech.com/dorado/latest/models/list/
- (Example from March 2026):
  - `dna_r10.4.1_e8.2_400bps_sup@v5.2.0` (LSK114)
  - `dna_r9.4.1_e8_sup@v3.6` (LSK110)

### 📌 Reference genome notes
- We typically use **GRCh38 no-alt** for human analyses.
- See guidance: https://lh3.github.io/2017/11/13/which-human-reference-genome-to-use
- Optional alternative for certain analyses: T2T CHM13 v2.0 (`/vast/projects/bahlo_epilepsy/ref_genomes/chm13v2.0.fa.gz`)

---

### Option A — Basecalling (recommended)

Example `sbatch` script (GPU job):

```bash
#!/bin/bash
#SBATCH --partition gpuq
#SBATCH --cpus-per-task 6
#SBATCH --mem 64G
#SBATCH --gres gpu:A30:4
#SBATCH --job-name <sample_basecalling>

module load dorado

dorado basecaller /stornext/System/data/nvidia/dorado/models/dna_r10.4.1_e8.2_400bps_sup@v5.0.0 --min-qscore 10 \
  --modified-bases 5mCG_5hmCG \
  --reference /vast/projects/bahlo_epilepsy/ref_genomes/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna \
  ./<path_to_all_pod5> \
  > <sample>_sup_v5.0.0_5mCG_5hmCG_aligned.bam
```

> **Runtime estimates (approx):**
> - `5mCG_5hmCG` (methylated CG only): ~40 hours *Recommended*
> - `5mC_5hmC` (all methylated Cs): ~80 hours

> **Tip:** If your job is interrupted, you can resume with:
> ```bash
> --resume-from <incomplete.bam>
> ```

---

### Option B — Align already-basecalled BAMs (no GPU required)

Example `sbatch` script:

```bash
#!/bin/bash
#SBATCH --cpus-per-task 6
#SBATCH --mem 64G
#SBATCH --job-name <sample_alignment>

module load dorado

samtools merge -@6 bam_pass/*bam \
  -o <sample>_sup_v5.0.0_5mCG_5hmCG_unaligned.bam

dorado aligner \
  /vast/projects/bahlo_epilepsy/ref_genomes/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna \
  <sample>_sup_v5.0.0_5mCG_5hmCG_unaligned.bam \
  > <sample>_sup_v5.0.0_5mCG_5hmCG_aligned.bam
```

---

## 3) Run epi2me wf-human-variation (Nextflow)

See the workflow docs: https://github.com/epi2me-labs/wf-human-variation

Example `sbatch`:

```bash
#!/bin/bash
#SBATCH --cpus-per-task 6
#SBATCH --mem 64G
#SBATCH --job-name epi2me

samtools sort -@ 6 \
  -o <sample>_sup_v5.0.0_5mCG_5hmCG_sorted.bam \
  <sample>_sup_v5.0.0_5mCG_5hmCG_aligned.bam

nextflow run epi2me-labs/wf-human-variation -r v2.8.0 \
  -c epi2me.config \
  -w ./work \
  --snp --sv --str --cnv --mod --phased \
  --bam <sample>_sup_v5.0.0_5mCG_5hmCG_sorted.bam \
  --ref /vast/projects/bahlo_epilepsy/ref_genomes/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna \
  --sample_name <sample_name> \
  -with-report \
  --bam_min_coverage 5 \
  -resume
```
> NOTE: nextflow run can be resumed in an interactive session if the sbatch job is interrupted. Just make sure to use the same `-w` work directory and `-c` config file. Recommended to run in a screen session to avoid accidental disconnection.

> NOTE: `epi2me.config` in this repo contains some minor resource adjustments to the default configuration.

---

## 🛠️ Common issues & troubleshooting

### Image pull fails in Nextflow

**Fix 1:** Increase Singularity/Apptainer pull timeout in `nextflow.config`:

```groovy
singularity {
  pullTimeout = '2h'  // or longer
}
```

**Fix 2:** Manually pull the failing image(s) to the cache directory:

```bash
apptainer pull $APPTAINER_CACHEDIR/ontresearch-wf-human-variation-snp-sha0d7e7e8e8207d9d23fdf50a34ceb577da364373e.img \
  docker://ontresearch/wf-human-variation-snp:sha0d7e7e8e8207d9d23fdf50a34ceb577da364373e
```

### Processes fail with insufficient memory/time

Increase resources for the failing process in `nextflow.config`.
Example (for `sniffles2`):

```groovy
process {
  withName: sniffles2 {
    memory = '30.GB'
    cpus   = 8
    time   = '12h'
  }
}
```

---

## 📚 References

- epi2me wf-human-variation: https://github.com/epi2me-labs/wf-human-variation
- Dorado models: https://software-docs.nanoporetech.com/dorado/latest/models/list/
