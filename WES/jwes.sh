#!/bin/bash

#SBATCH --job-name=WES
#SBATCH --output=logs/wes_%j.log
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --array=1-18%4  

# ==========================================
# Comprehensive WES Analysis Pipeline
# ==========================================

# Configuration
CONFIG_FILE="config.txt"
if [[ ! -f "$CONFIG_FILE" ]]; then
    echo "Error: Config file not found. Create config.txt first."
    exit 1
fi
source $CONFIG_FILE

# Get sample information from sample sheet
SAMPLE_ID=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${SAMPLE_SHEET} | cut -f1)
SAMPLE_TYPE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${SAMPLE_SHEET} | cut -f2)
FASTQ_R1=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${SAMPLE_SHEET} | cut -f3)
FASTQ_R2=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${SAMPLE_SHEET} | cut -f4)
NORMAL_SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${SAMPLE_SHEET} | cut -f5)

# Create sample output directory
SAMPLE_DIR="${OUTPUT_DIR}/${SAMPLE_ID}"
mkdir -p ${SAMPLE_DIR}/logs ${SAMPLE_DIR}/qc ${SAMPLE_DIR}/bam ${SAMPLE_DIR}/vcf

# Start logging
LOGFILE="${SAMPLE_DIR}/logs/${SAMPLE_ID}.log"
exec > >(tee -a ${LOGFILE}) 2>&1
echo "[$(date)] Starting WES analysis for sample ${SAMPLE_ID}"

# ====================
# QC and Preprocessing
# ====================
echo "[$(date)] Step 1: Quality control with FastQC..."
fastqc -o ${SAMPLE_DIR}/qc -t 8 ${FASTQ_R1} ${FASTQ_R2}

echo "[$(date)] Step 2: Adapter trimming with Trimmomatic..."
trimmomatic PE -threads 8 ${FASTQ_R1} ${FASTQ_R2} \
    ${SAMPLE_DIR}/qc/${SAMPLE_ID}_R1_paired.fastq.gz ${SAMPLE_DIR}/qc/${SAMPLE_ID}_R1_unpaired.fastq.gz \
    ${SAMPLE_DIR}/qc/${SAMPLE_ID}_R2_paired.fastq.gz ${SAMPLE_DIR}/qc/${SAMPLE_ID}_R2_unpaired.fastq.gz \
    ILLUMINACLIP:${ADAPTERS}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

# ====================
# Alignment & Processing
# ====================
echo "[$(date)] Step 3: Aligning reads with BWA MEM..."
bwa mem -M -t 16 -R "@RG\tID:${SAMPLE_ID}\tSM:${SAMPLE_ID}\tPL:ILLUMINA\tLB:lib1" \
    ${REFERENCE} \
    ${SAMPLE_DIR}/qc/${SAMPLE_ID}_R1_paired.fastq.gz \
    ${SAMPLE_DIR}/qc/${SAMPLE_ID}_R2_paired.fastq.gz \
    | samtools view -bS - > ${SAMPLE_DIR}/bam/${SAMPLE_ID}.bam

echo "[$(date)] Step 4: Sorting BAM file..."
samtools sort -@ 16 -m 4G -o ${SAMPLE_DIR}/bam/${SAMPLE_ID}.sorted.bam ${SAMPLE_DIR}/bam/${SAMPLE_ID}.bam
rm ${SAMPLE_DIR}/bam/${SAMPLE_ID}.bam

echo "[$(date)] Step 5: Marking duplicates..."
gatk MarkDuplicates \
    -I ${SAMPLE_DIR}/bam/${SAMPLE_ID}.sorted.bam \
    -O ${SAMPLE_DIR}/bam/${SAMPLE_ID}.sorted.markdup.bam \
    -M ${SAMPLE_DIR}/qc/${SAMPLE_ID}.metrics.txt \
    --CREATE_INDEX true

# ====================
# GATK Best Practices
# ====================
echo "[$(date)] Step 6: Base quality score recalibration..."
gatk BaseRecalibrator \
    -R ${REFERENCE} \
    -I ${SAMPLE_DIR}/bam/${SAMPLE_ID}.sorted.markdup.bam \
    --known-sites ${DBSNP} \
    --known-sites ${MILLS} \
    --known-sites ${INDELS} \
    -O ${SAMPLE_DIR}/bam/${SAMPLE_ID}.recal.table

gatk ApplyBQSR \
    -R ${REFERENCE} \
    -I ${SAMPLE_DIR}/bam/${SAMPLE_ID}.sorted.markdup.bam \
    --bqsr-recal-file ${SAMPLE_DIR}/bam/${SAMPLE_ID}.recal.table \
    -O ${SAMPLE_DIR}/bam/${SAMPLE_ID}.final.bam

echo "[$(date)] Step 7: Collecting alignment metrics..."
gatk CollectHsMetrics \
    -I ${SAMPLE_DIR}/bam/${SAMPLE_ID}.final.bam \
    -O ${SAMPLE_DIR}/qc/${SAMPLE_ID}.hsmetrics.txt \
    -R ${REFERENCE} \
    --BAIT_INTERVALS ${TARGETS_INTERVALS} \
    --TARGET_INTERVALS ${TARGETS_INTERVALS}

# ====================
# Variant Calling
# ====================
if [[ "${SAMPLE_TYPE}" == "tumor" && ! -z "${NORMAL_SAMPLE}" ]]; then
    # Somatic variant calling for tumor samples with matching normal
    NORMAL_BAM="${OUTPUT_DIR}/${NORMAL_SAMPLE}/bam/${NORMAL_SAMPLE}.final.bam"
    
    echo "[$(date)] Step 8: Somatic variant calling with Mutect2..."
    gatk Mutect2 \
        -R ${REFERENCE} \
        -I ${SAMPLE_DIR}/bam/${SAMPLE_ID}.final.bam \
        -I ${NORMAL_BAM} \
        -tumor ${SAMPLE_ID} \
        -normal ${NORMAL_SAMPLE} \
        --germline-resource ${GERMLINE_RESOURCE} \
        --panel-of-normals ${PON} \
        -L ${TARGETS_INTERVALS} \
        -O ${SAMPLE_DIR}/vcf/${SAMPLE_ID}.somatic.vcf.gz
    
    echo "[$(date)] Step 9: Filtering somatic variants..."
    gatk FilterMutectCalls \
        -R ${REFERENCE} \
        -V ${SAMPLE_DIR}/vcf/${SAMPLE_ID}.somatic.vcf.gz \
        -O ${SAMPLE_DIR}/vcf/${SAMPLE_ID}.somatic.filtered.vcf.gz
    
    echo "[$(date)] Step 10: Annotating somatic variants..."
    ${VEP_PATH}/vep \
        --input_file ${SAMPLE_DIR}/vcf/${SAMPLE_ID}.somatic.filtered.vcf.gz \
        --output_file ${SAMPLE_DIR}/vcf/${SAMPLE_ID}.somatic.annotated.vcf \
        --format vcf \
        --cache \
        --dir_cache ${VEP_CACHE} \
        --pick \
        --symbol \
        --hgvs \
        --af \
        --protein \
        --numbers \
        --regulatory \
        --canonical \
        --check_existing \
        --vcf \
        --fork 4

else
    # Germline variant calling for normal samples
    echo "[$(date)] Step 8: Germline variant calling with HaplotypeCaller..."
    gatk HaplotypeCaller \
        -R ${REFERENCE} \
        -I ${SAMPLE_DIR}/bam/${SAMPLE_ID}.final.bam \
        -L ${TARGETS_INTERVALS} \
        -O ${SAMPLE_DIR}/vcf/${SAMPLE_ID}.raw.vcf.gz
    
    echo "[$(date)] Step 9: Filtering variants with VQSR..."
    gatk VariantRecalibrator \
        -R ${REFERENCE} \
        -V ${SAMPLE_DIR}/vcf/${SAMPLE_ID}.raw.vcf.gz \
        --resource:hapmap,known=false,training=true,truth=true,prior=15.0 ${HAPMAP} \
        --resource:omni,known=false,training=true,truth=false,prior=12.0 ${OMNI} \
        --resource:1000G,known=false,training=true,truth=false,prior=10.0 ${PHASE1_SNP} \
        --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ${DBSNP} \
        -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
        -mode SNP \
        -O ${SAMPLE_DIR}/vcf/${SAMPLE_ID}.recal \
        --tranches-file ${SAMPLE_DIR}/vcf/${SAMPLE_ID}.tranches \
        --rscript-file ${SAMPLE_DIR}/vcf/${SAMPLE_ID}.plots.R
    
    gatk ApplyVQSR \
        -R ${REFERENCE} \
        -V ${SAMPLE_DIR}/vcf/${SAMPLE_ID}.raw.vcf.gz \
        -O ${SAMPLE_DIR}/vcf/${SAMPLE_ID}.filtered.vcf.gz \
        --truth-sensitivity-filter-level 99.5 \
        --tranches-file ${SAMPLE_DIR}/vcf/${SAMPLE_ID}.tranches \
        --recal-file ${SAMPLE_DIR}/vcf/${SAMPLE_ID}.recal \
        -mode SNP
    
    echo "[$(date)] Step 10: Annotating variants..."
    ${VEP_PATH}/vep \
        --input_file ${SAMPLE_DIR}/vcf/${SAMPLE_ID}.filtered.vcf.gz \
        --output_file ${SAMPLE_DIR}/vcf/${SAMPLE_ID}.annotated.vcf \
        --format vcf \
        --cache \
        --dir_cache ${VEP_CACHE} \
        --pick \
        --symbol \
        --hgvs \
        --af \
        --protein \
        --numbers \
        --regulatory \
        --canonical \
        --check_existing \
        --vcf \
        --fork 4
fi

# ====================
# Generate Reports
# ====================
echo "[$(date)] Step 11: Generating coverage reports..."
gatk DepthOfCoverage \
    -R ${REFERENCE} \
    -I ${SAMPLE_DIR}/bam/${SAMPLE_ID}.final.bam \
    -L ${TARGETS_INTERVALS} \
    -O ${SAMPLE_DIR}/qc/${SAMPLE_ID}.coverage \
    --summary-coverage-threshold 10 \
    --summary-coverage-threshold 20 \
    --summary-coverage-threshold 50 \
    --summary-coverage-threshold 100

echo "[$(date)] Step 12: Generating MultiQC report..."
multiqc -o ${SAMPLE_DIR}/qc ${SAMPLE_DIR}

echo "[$(date)] WES analysis for sample ${SAMPLE_ID} completed"