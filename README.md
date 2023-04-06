# rtHealthyTissues
readthrough in healthy tissues

# readthrough in healthy tissues
framework to quantify transcription readthrough in human samples from the GTEx project

### Accessing the Data
The primary entry point for accessing GTEx data is through the GTEx portal (https://gtexportal.org/). The GTEx Portal provides open access to data including gene expression, QTLs, and histology images. However, due to the nature of the donor consent agreement, raw data and attributes which might be used to identify the donors, such as raw sequencing data or variant calls, are not publicly available on the GTEx Portal. Accessing the raw data requires authorization from the NIH database of Genotypes and Phenotypes (dbGaP). One can apply for access to the data in the “Request Access” tab. In this project Authorization was granted to Ana Rita Grosso (dbGaP Accession phs000424.v8. p2), where NIH Genomic Data Sharing Policy policies are applied.

### Quantification of Transcription Readthrough
Our pipeline employs already existent tools, namely STAR (v2.7.8a) (Dobin et al., 2013), ARTDeco (Roth, Heinz and Benner, 2020) and bedtools (v2.30.0) (Quinlan and Hall, 2010). We first converted the bam files downloaded from the dbGaP back to fastq using samtools (v.1.10) (Danecek et al., 2021) and then re-aligned them to the reference genome using STAR (GRCh38.P13 assembly; annotation.gtf v33).

```
samtools collate -u -O filename.bam | samtools fastq -1 filename_1.fq -2 filename_2.fq -0 /dev/null -s /dev/null -n
```

```
STAR --runThreadN 8 –genomeDir ref_genome.fa.star.idx/ \
--outFileNamePrefix output_folder/filename --outSAMtype BAM Unsorted \
--readFilesIn filename_1.fastq filename_2.fastq
```

Then, to detect transcription readthrough (TRT) directly from the bam files generated from STAR, we used ARTDeco, a pipeline design to search for continuous coverage over a minimal length downstream of the 3’end of each gene locus using a rolling window approach. Transcription levels of the window must meet the thresholds to be considered as part of the readthrough tail: we used a rolling window of 500bp, a minimum length of 2000 bp and min coverage of 0.15 FPKM. 

```
ARTDeco -home-dir ARTDECO_DIR -bam-files-dir BAM_FILES_DIR -gtf-file GTF_FILE -cpu NUM_CPU -chrom-sizes-file CHROM_SIZES_FILE
```

ARTDeco returns a variety of metrics to measure readthrough. For the purpose of this work, the information of interest are contained inside the “quantification” and “dogs” folders (expression levels and novel transcripts created as a result of readthrough, respectively).

