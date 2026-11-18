#!/usr/bin/env bash
# pySCENIC: run between 08_scenic_export.R and 10_scenic_import.R
# needs a python env with pyscenic; edit PROJECT_DIR when running elsewhere
set -euo pipefail

PROJECT_DIR="E:/my_project/70_intestine_code"
SCENIC_DIR="${PROJECT_DIR}/Single_cell/output/scenic"

LOOM="${SCENIC_DIR}/input.loom"
ADJ="${SCENIC_DIR}/adj.tsv"
REG="${SCENIC_DIR}/reg.csv"
OUT="${SCENIC_DIR}/scenic.loom"

TF_LIST="${PROJECT_DIR}/Sample_metdata/hs_hgnc_curated_tfs2pig.txt"

# aertslab resources, download once:
#   https://resources.aertslab.org/cistarget/motif2tf/motifs-v9-nr.hgnc-m0.001-o0.0.tbl
#   https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq-r80/10kb_up_and_down_tss/mc9nr/
MOTIF_RANKINGS="hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather"
MOTIF_ANNOT="motifs-v9-nr.hgnc-m0.001-o0.0.tbl"

WORKERS=10

# 1. gene regulatory network
pyscenic grn \
  --num_workers ${WORKERS} \
  --output "${ADJ}" \
  --method grnboost2 \
  "${LOOM}" \
  "${TF_LIST}"

# 2. regulon prediction from motif enrichment
pyscenic ctx \
  "${ADJ}" \
  "${MOTIF_RANKINGS}" \
  --annotations_fname "${MOTIF_ANNOT}" \
  --expression_mtx_fname "${LOOM}" \
  --mode "dask_multiprocessing" \
  --output "${REG}" \
  --num_workers ${WORKERS} \
  --mask_dropouts

# 3. regulon activity per cell
pyscenic aucell \
  "${LOOM}" \
  "${REG}" \
  --output "${OUT}" \
  --num_workers ${WORKERS}
