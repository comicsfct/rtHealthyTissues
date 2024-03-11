# rtHealthyTissues
Repository dedicated to the code and data associated with the manuscript titled:
Transcription readthrough is prevalent in healthy human tissues and associated with inherent genomic features 
https://www.nature.com/articles/s42003-024-05779-5

# readthrough in healthy tissues
framework to quantify transcription readthrough in human samples from the GTEx project

### Accessing the Data
The primary entry point for accessing GTEx data is through the GTEx portal (https://gtexportal.org/). The GTEx Portal provides open access to data including gene expression, QTLs, and histology images. However, due to the nature of the donor consent agreement, raw data and attributes which might be used to identify the donors, such as raw sequencing data or variant calls, are not publicly available on the GTEx Portal. Accessing the raw data requires authorization from the NIH database of Genotypes and Phenotypes (dbGaP). One can apply for access to the data in the “Request Access” tab. In this project Authorization was granted to Ana Rita Grosso (dbGaP Accession phs000424.v8. p2), where NIH Genomic Data Sharing Policy policies are applied.

### Filtering the Data
We considered only paired-end samples with at least 60 million reads per sample and prepared with the Illumina TruSeq library construction protocol (non-strand specific polyA+ selected library). Cell culture samples and tissues containing fewer than 50 samples were excluded. Healthy subjects were selected by filtering samples for “violent and fast deaths" and "no terminal diseases". We obtained 2778 samples from 23 healthy human tissues that were used for downstream analyses. 

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

Then, to detect transcription readthrough (TRT) directly from the bam files generated from STAR, we used ARTDeco,
a pipeline design to search for continuous coverage over a minimal length downstream of the 3’end of each gene locus using 
a rolling window approach. Transcription levels of the window must meet the thresholds to be considered as part of the readthrough tail: 
we used a rolling window of 500bp, a minimum length of 2000 bp and min coverage of 0.15 FPKM. ARTDeco uses HOMER’s tools (Heinz et al., 2010) 
to select only uniquely mapped reads for downstream analysis and returns a variety of metrics to measure readthrough.

```
ARTDeco -home-dir ARTDECO_DIR -bam-files-dir BAM_FILES_DIR -gtf-file GTF_FILE -cpu NUM_CPU -chrom-sizes-file CHROM_SIZES_FILE
```

For the purpose of this work, the information of interest is contained inside the “quantification” and “dogs” folders
(expression levels and novel transcripts created as a result of readthrough, respectively).

### Filter ARTDeco output for non-stranded RNA-seq samples

GTEx samples were profiled using non-stranded RNAseq libraries. Since transcriptional signals can come from either direction, ARTDeco is ambiguous when inferring a true downstream transcript in some cases. Thus, a significant number of reads identified as downstream transcripts were in reality reads coming from genes being expressed in the opposite direction. To eliminate these false positives created by the lack of strandedness, we filtered the output from ARTDeco to report only entries that have no overlap with genes in the opposite strand using the intersect function from bedtools (v2.30.0) (Quinlan and Hall, 2010). 

for each .bed file inside the dog folder, we run bedtools as: 

```
bedtools intersect -a filename.bed -b genes_condensed.bed -v -S > filename.filtered.intersect
```

This approach discards RT transcripts with close downstream neighbors in the opposite strand but ensures that our list of readthrough genes is robust. 
To optimize this step, we created a script (`run-bedtools-intersect.sh`) that automatically runs this bedtool command to all files inside the `dogs` 
folder (artdeco output) and creates a new folder containing all the new filtered dog files (`dogs-filtered`). 

```
./filter-artdeco-output/run-bedtools-intersect.sh artdeco_output_folder ./filter-artdeco-output/dogs.filter.dubious.genes.py
```

The scripts takes two files as input: (1) standard artdeco output containing all subfolders; 
(2) custom python script (available inside the same folder) that filters all the remaining .txt files according to the genes filtered in the new dog-filtered bed file. 
A new folder is created containing these new files.

### Create a Master Table Containing all Data for Downstream Analysis

We created a script that takes as input a folder containing one or multiple ARTDeco output folders (one for each tissue in this scenario) 
and uses two specific subfolders ("quantification" and "dogs") to combine expression levels from genes and readthrough regions into a simple dataframe 
that can be used for all the downstream analysis.

```
~/python ./compact_artdeco_output.py artdeco-out-example/
```

The output is a table `master.table.GTEx.samples.RT.analysis.tab` contanining most of the information used for downstream analyis.
Here, we defined as expressed, genes with a FPKM > 1 in at least 25% of the samples within a given tissue.
Only RT transcripts coming from expressed genes in each given tissue were considered for downstream analysis.

Scripts to plot some of the figures in the manuscript can be found above.
