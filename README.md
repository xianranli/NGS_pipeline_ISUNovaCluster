# NGS_pipeline_ISUNovaCluster
The purpose of these two scripts is to automatically align publicly available NGS reads to the reference in Nova Cluster at Iowa State University.

One input file storing the SRR accessions is required.

The first perl script, run_bwa_1st.pl will read the target SRR list, download the SRA file, convert it to Fastq file and split the large file into chunk files, then the multiple chunk files will be submit by Slurm, to run the 2nd perl script, which is the core pipeline for QC and mapping through bwa for each chunk. After the last chunk is finished, the bam file from each chunk will be merged into a single bam file.
