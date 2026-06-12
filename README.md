# scTP-seq: A Single-cell Multi-omics Approach for Simultaneous Profiling of Transcriptome and Intracellular Proteins
# scTP-seq description
Single-cell RNA sequencing has revolutionized cell typing but cannot concurrently profile intracellular proteins, creating a critical blind spot that obscures post-transcriptional regulation and impedes functional biomarkers discovery. To address this, we developed scTP-seq (single-cell Transcriptome and Protein sequencing), a platform that enables simultaneous single-cell profiling of the transcriptome and intracellular proteins without prior cell manipulation, ensuring minimal technical bias. 

<img width="1788" height="1929" alt="图片1" src="https://github.com/user-attachments/assets/aa1c2584-1130-43d4-b897-0ea4ba4b0f33" />

## Step 1: Raw data processing and alignment
Use STARsolo with default parameters to process the raw FASTQ data, performing cell barcode identification and UMI counting.

```
sh 0_Preprocess/run_star.sh
```

## Step 2: Downstream analysis pipeline
The downstream analysis scripts from the scTP-seq paper have been uploaded to the code directories organized by figures.
