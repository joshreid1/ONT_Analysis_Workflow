#!/bin/bash
# Define the image pulls
declare -a pulls=(
    "apptainer pull /vast/scratch/users/reid.j/nextflow/singularity_cache/ontresearch-wf-human-variation-sv-shac591518dd32ecc3936666c95ff08f6d7474e9728.img docker://ontresearch/wf-human-variation-sv:shac591518dd32ecc3936666c95ff08f6d7474e9728"
    "apptainer pull /vast/scratch/users/reid.j/nextflow/singularity_cache/ontresearch-wf-common-sha338caea0a2532dc0ea8f46638ccc322bb8f9af48.img docker://ontresearch/wf-common:sha338caea0a2532dc0ea8f46638ccc322bb8f9af48"
    "apptainer pull /vast/scratch/users/reid.j/nextflow/singularity_cache/ontresearch-wf-human-variation-sha2b856c1f358ddf1576217a336bc0e9864b6dc0ed.img docker://ontresearch/wf-human-variation:sha2b856c1f358ddf1576217a336bc0e9864b6dc0ed"
    "apptainer pull /vast/scratch/users/reid.j/nextflow/singularity_cache/ontresearch-spectre-sha49a9fe474da9860f84f08f17f137b47a010b1834.img docker://ontresearch/spectre:sha49a9fe474da9860f84f08f17f137b47a010b1834"
    "apptainer pull /vast/scratch/users/reid.j/nextflow/singularity_cache/ontresearch-wf-human-variation-snp-sha17e686336bf6305f9c90b36bc52ff9dd1fa73ee9.img docker://ontresearch/wf-human-variation-snp:sha17e686336bf6305f9c90b36bc52ff9dd1fa73ee9"
)

# Submit each pull as a separate job
for pull in "${pulls[@]}"
do
    job_script=$(mktemp)
    cat <<EOT > $job_script
#!/bin/bash
#SBATCH --job-name=pull_apptainer_image
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --time=02:00:00
#SBATCH --mem=10G
#SBATCH --cpus-per-task=1

module load apptainer

echo $pull >&2

$pull
if [ $? -ne 0 ]; then
    echo "Failed to pull image: $pull" >&2
    exit 1
fi
EOT
    sbatch $job_script
    rm $job_script
done
