#!/bin/bash

#SBATCH --time=24:00:00
#SBATCH -J cellranger_DB_AKC
#SBATCH -o process_%j.out
#SBATCH -e process_%j.err
#SBATCH -A XXXX
#SBATCH -N 1
#SBATCH --tasks-per-node=20

echo "Script used:"
echo "########################################"
cat $0
echo "########################################"
echo ""

module purge
module load CellRanger/7.0.0

cellranger count --id=DB_AKC --libraries library.csv --transcriptome /projects/fs1/razofz/DB_AKC/mouse_reference/refdata-gex-mm10-2020-A --chemistry SC3Pv3 --feature-ref feature-reference.csv --nosecondary --no-bam
