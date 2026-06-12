set -eu
fq1="/work/A_R1.fastq.gz"
fq2="/work/A_R2.fastq.gz"
barcodes="/work/barcodes.txt"

STAR --soloType CB_UMI_Complex --soloCBposition 0_0_0_9 0_10_0_19 --soloUMIposition 0_25_0_34 --soloCBmatchWLtype EditDist_2 --soloCBwhitelist $barcodes $barcodes \
    --runMode alignReads \
    --readFilesIn $fq2 $fq1 \
    --genomeDir /reference/human/refdata-gex-GRCh38-2024-A/fasta \
    --outFileNamePrefix star_outdir/ \
    --soloFeatures Gene GeneFull_Ex50pAS Velocyto \
    --quantMode GeneCounts \
    --soloCellFilter EmptyDrops_CR \
    --soloStrand Unstranded \
    --soloCellReadStats Standard \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM \
    --soloUMIdedup 1MM_CR --clipAdapterType CellRanger4 \
    --outFilterScoreMin 30 --soloUMIfiltering MultiGeneUMI_CR \
    --runThreadN 16 \
    --soloBarcodeReadLength 0
