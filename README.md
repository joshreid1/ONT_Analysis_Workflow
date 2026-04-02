# ONT Analysis Workflow

A workflow for running ONT basecalling (via Dorado) and downstream analysis (via epi2me wf-human-variation).

## 🚀 Overview

This repository provides guidance and example `sbatch` scripts for:

- **Option A:** Basecalling from raw **pod5** files using **Dorado** (recommended)
- **Option B:** Aligning already basecalled **BAM** files using **Dorado Aligner**
- Running **epi2me wf-human-variation** (Nextflow) on the resulting aligned BAM

> **Note:** The PromethION 2 Integrated uses the Dorado server (essentially a wrapper on standalone Dorado) that often lags behind the latest releases. Versioning is different and cannot be directly compared. We recommend basecalling on HPC (Option A) instead of using on-device output (Option B).

---

## 📁 Repository Structure

<!-- TODO: Update this tree to reflect the actual repo layout once scripts are finalised -->

```
.
├── README.md
├── epi2me.config              # Nextflow config with WEHI HPC resource adjustments
├── scripts/
│   ├── basecalling_optionA.sh # Dorado basecalling from pod5 (GPU)
│   ├── alignment_optionB.sh   # Dorado aligner for pre-basecalled BAMs
│   └── epi2me_run.sh          # epi2me wf-human-variation sbatch script
```

---

## ✅ Prerequisites

- Access to WEHI HPC with **GPU** access to run Dorado basecalling (Option A)
- `dorado` (basecaller + aligner), `samtools`, and `nextflow` (all available as modules on WEHI HPC)

---

## 1) Download data (Mediaflux)

Download raw **pod5** files via Mediaflux. For most workflows you should request **all pod5 files** (pass, fail, recovered) to ensure complete coverage.

Move all pod5 files to a single directory for processing. Example:

```bash
mkdir ./pod5_all
mv pod5_pass/* pod5_all/
mv pod5_fail/* pod5_all/
mv pod5_recovered/* pod5_all/
```

**If you are given legacy FAST5 files** instead of pod5, convert them:

```bash
pod5 convert fast5 ./input/*.fast5 --output converted.pod5
```

📄 Pod5 docs: https://pod5-file-format.readthedocs.io/en/0.1.21/docs/tools.html

---

## 2) Run Dorado (basecalling + optional alignment)

Dorado is the current ONT basecaller. It can also align reads to a reference in the same run (alignment via minimap2).

📄 Minimap2 docs: https://github.com/lh3/minimap2

### 🧠 Model notes

- Available models: https://software-docs.nanoporetech.com/dorado/latest/models/list/
- (Example from March 2026):
  - `dna_r10.4.1_e8.2_400bps_sup@v5.2.0` (LSK114)
  - `dna_r9.4.1_e8_sup@v3.6` (LSK110)

### 📌 Reference genome notes

- We typically use **GRCh38 no-alt** for human analyses.
- See guidance: https://lh3.github.io/2017/11/13/which-human-reference-genome-to-use
- WEHI HPC path: `/vast/projects/bahlo_epilepsy/ref_genomes/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna`
- Optional alternative for certain analyses: T2T CHM13 v2.0 (`/vast/projects/bahlo_epilepsy/ref_genomes/chm13v2.0.fa.gz`)

> **Note:** The `--str` subworkflow in epi2me wf-human-variation is **only compatible with hg38/GRCh38**, not CHM13.

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

module load samtools
module load nextflow

# Sort and index the aligned BAM — required before running epi2me
samtools sort -@ 6 \
  -o <sample>_sup_v5.0.0_5mCG_5hmCG_sorted.bam \
  <sample>_sup_v5.0.0_5mCG_5hmCG_aligned.bam

samtools index -@ 6 <sample>_sup_v5.0.0_5mCG_5hmCG_sorted.bam

# Specify --sex XX or --sex XY for accurate STR chrX genotyping.
# If omitted, the workflow will attempt to infer sex from allosome coverage.
nextflow run epi2me-labs/wf-human-variation -r v2.8.0 \
  -c epi2me.config \
  -w ./work \
  --snp --sv --str --cnv --mod --phased \
  --sex XY \
  --bam <sample>_sup_v5.0.0_5mCG_5hmCG_sorted.bam \
  --ref /vast/projects/bahlo_epilepsy/ref_genomes/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna \
  --sample_name <sample_name> \
  -with-report \
  --bam_min_coverage 5 \
  -resume
```

> **NOTE:** `nextflow run` can be resumed in an interactive session if the sbatch job is interrupted. Just make sure to use the same `-w` work directory and `-c` config file. Recommended to run in a `screen` session to avoid accidental disconnection.

> **NOTE:** `epi2me.config` in this repo contains some minor resource adjustments to the default configuration.

---

## 4) epi2me wf-human-variation Outputs

Output files are prefixed with `<sample_name>`. Each subworkflow produces distinct outputs:

| Flag | Tool | Key output files | Description |
|------|------|-----------------|-------------|
| `--snp` | [Clair3](https://github.com/HKU-BAL/Clair3) | `<sample>.wf_snp.vcf.gz`<br>`<sample>.wf_snp_clinvar.vcf.gz`<br>`<sample>.wf-human-snp-report.html` | SNPs and small indels, annotated with SnpEff and ClinVar |
| `--sv` | [Sniffles2](https://github.com/fritzsedlazeck/Sniffles) | `<sample>.wf_sv.vcf.gz`<br>`<sample>.wf_sv.snf`<br>`<sample>.wf-human-sv-report.html` | Structural variants (DEL, INS, DUP, INV, BND); `.snf` enables multi-sample SV calling |
| `--str` | [Straglr](https://github.com/bcgsc/straglr) | `<sample>.wf_str.vcf.gz`<br>`<sample>.wf_str.tsv`<br>`<sample>.haplotagged.cram`<br>`<sample>.wf-human-str-report.html` | STR expansion genotypes; TSV contains read-level spanning information; haplotagged CRAM for downstream use. **hg38 only.** Automatically enables `--snp`. |
| `--cnv` | [Spectre](https://github.com/nanoporetech/ont-spectre) | `<sample>.wf_cnv.vcf.gz`<br>`<sample>.wf-human-cnv-report.html` | Large CNVs (>100kb), annotated with SnpEff; includes sex chromosome karyotype prediction. Use `--cnv --use_qdnaseq` for shallow WGS or adaptive sampling. |
| `--mod` | [modkit](https://github.com/nanoporetech/modkit) | `<sample>.wf_mods.1.bedmethyl.gz`<br>`<sample>.wf_mods.2.bedmethyl.gz`<br>`<sample>.wf_mods.ungrouped.bedmethyl.gz` | Aggregated CpG methylation counts (bedMethyl format); per-haplotype methylation when used with `--phased`; bigWig for genome browser visualisation |
| `--phased` | [WhatsHap](https://github.com/whatshap/whatshap) | `<sample>.haplotagged.cram`<br>Phased SNP + SV VCFs | Haplotagged BAM/CRAM with HP and PS tags; phased VCFs for SNPs and SVs |

**Additional outputs (all runs):**

- `<sample>.wf-human-alignment-report.html` — Alignment QC summary
- `<sample>.mosdepth.summary.txt` — Mean coverage per chromosome
- `<sample>.regions.bed.gz` — Mean coverage per region (BED)
- `<sample>.stats.json` — Base statistics (reads, mappings, SNPs, SVs)

---

#Running nf-cavalier on ONT data

- Currently only tested on singletons (as samples are not joint called). Merging multiple SNV VCFs and SV VCFs from different samples is possible but not yet tested.
- CRAM input accepted in current dev branch
- Add the parameters below for ONT samples:

```bash
    --short_vcf = '/path/to/sample.wf_snp.vcf.gz'
    --struc_vcf = '/path/to/sample.wf_sv.vcf.gz'
    
    --struc_fill_tags = true
    --short_fill_tags = true 
```

---

# Common use cases

## 1) Variant phasing

Visualising phased variants in IGV:

1. Load haplotagged BAM/CRAM (i.e. `HP` and `PS` tags have been added by WhatsHap)
2. Navigate to gene/region of interest
3. Highlight variants of interest (via 'Define a region of interest')
4. "Group alignments by" → phase (i.e. `HP` haplotype tag)
5. "Color alignments by" → tag → `PS` (phase set — contiguous block over which variants are jointly phased)

## 2) Translocation detection

Visualising translocations (i.e. `BND` events called by Sniffles2) in IGV:

1. Load BAM/CRAM
2. Navigate to translocation breakpoint (e.g. via "Go to" → `chr1:12345`)
3. Click on a read supporting the translocation and view: a) 'Supplementary alignments' b) Clipping
4. "Group alignments by" → "Chimeric"
5. Right click on chimeric read → "Supplementary Reads Diagram" to view a schematic of the translocation event
6. Right click on chimeric read → "View chimeric alignments in split screen"

## 3) Methylation analysis

1. Load BAM/CRAM (must have `MM`/`ML` tags from Dorado basecalling with a modified-base model)
2. Navigate to gene/region of interest
3. "Color alignments by" → "Base modification 2 color (all)": blue = unmodified, red = modified (opacity = confidence of call)
4. Click on a histogram bar to see breakdown of methylation calls at that position

Load CpG island annotations:

**File → Load track from URL →** `https://data.broadinstitute.org/igvdata/annotations/hg38/cpgIslandExt.bed.gz`

Example methylation analysis scripts are available in the `methylation_scripts/` directory of this repo. These include:
- `methylation_compare.R`: Compare methylation profiles between two samples (used for comparing twin pairs), including pathway enrichment analysis of differentially methylated genes.
- Tested on modkit output with "traditional" preset

```bash
$ modkit pileup $input_bam $output_bed --ref $ref
  --modified-bases 5mC 5hmC \ 
  --combine-mods \ #Combine all modified bases into a single "modified" category for downstream analysis (e.g. DMR calling)
  --combine-strands --cpg
```

> **IGV tip:** Enable "Hide Small Indels" (e.g. <5bp) to reduce noise in screenshots.

---

## 🛠️ Common issues & troubleshooting

### Image pull fails in Nextflow

**Fix 1:** Increase Singularity/Apptainer pull timeout in `epi2me.config`:

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

Increase resources for the failing process in `epi2me.config`.
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

### Coverage below minimum threshold

epi2me wf-human-variation requires ≥20x coverage by default. If you have lower coverage and wish to proceed:

```bash
--bam_min_coverage 0   # disables the coverage check entirely
```

### Clair3 model not detected

If the basecall model cannot be detected from BAM headers, provide it explicitly:

```bash
--override_basecaller_cfg dna_r10.4.1_e8.2_400bps_sup@v5.0.0
```
---

## 📚 References

- epi2me wf-human-variation: https://github.com/epi2me-labs/wf-human-variation
- Dorado basecaller: https://github.com/nanoporetech/dorado
- Dorado models: https://software-docs.nanoporetech.com/dorado/latest/models/list/
- Clair3 (SNP): https://github.com/HKU-BAL/Clair3
- Sniffles2 (SV): https://github.com/fritzsedlazeck/Sniffles
- Straglr (STR): https://github.com/bcgsc/straglr
- Spectre (CNV): https://github.com/nanoporetech/ont-spectre
- modkit (methylation): https://github.com/nanoporetech/modkit
- WhatsHap (phasing): https://github.com/whatshap/whatshap
- pod5 tools: https://pod5-file-format.readthedocs.io/en/0.1.21/docs/tools.html
- GRCh38 no-alt reference guidance: https://lh3.github.io/2017/11/13/which-human-reference-genome-to-use