#!/bin/bash
#SBATCH --job-name=PLACER
#SBATCH --output=/mnt/home1/miska/hl725/scratch/projects/ciclid/urika/placer_out/%x_%j.out
#SBATCH --error=/mnt/home1/miska/hl725/scratch/projects/ciclid/urika/placer_out/%x_%j.err
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=256G
#SBATCH --partition=2204

# PLACER v0.2 - TE Insertion Detection
# Data: Cichlid D23
# BAM: final_merged_Yohann_D23.bam (23GB)
# Ref: Astatotilapia_calliptera.fAstCal1.3.dna.toplevel.fa
# TE Library: MWCichlidTE-3.2.splitted.namefixed_renamed_all_all.fa

set -e

# 输出目录
OUT_DIR=/mnt/home1/miska/hl725/scratch/projects/ciclid/urika/placer_out
mkdir -p $OUT_DIR

# 路径配置
BAM=/mnt/home1/miska/hl725/scratch/projects/ciclid/urika/data/final_merged_Yohann_D23.bam
REF=/mnt/home1/miska/hl725/scratch/projects/ciclid/urika/references/Astatotilapia_calliptera.fAstCal1.3.dna.toplevel.fa
TE=/mnt/home1/miska/hl725/scratch/projects/ciclid/urika/references/MWCichlidTE-3.2.splitted.namefixed_renamed_all_all.fa
PLACER=/mnt/home1/miska/hl725/projects/PLACER/build/placer

echo "========================================"
echo "PLACER Pipeline - Cichlid D23"
echo "========================================"
echo "BAM: $BAM"
echo "REF: $REF"
echo "TE: $TE"
echo "Output: $OUT_DIR"
echo "Date: $(date)"
echo "========================================"

# 生成索引（如果不存在）
if [ ! -f "${REF}.fai" ]; then
    echo "Creating reference index..."
    samtools faidx $REF
fi

if [ ! -f "${TE}.fai" ]; then
    echo "Creating TE index..."
    samtools faidx $TE
fi

# 运行 PLACER（使用 TE 库，Gate1 和 TE Reverse Index 都会启用）
cd $OUT_DIR
$PLACER $BAM $REF $TE

echo "========================================"
echo "PLACER completed at $(date)"
echo "========================================"
