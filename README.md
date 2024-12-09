# aparnaanandkumar.github.io
# Instructions to Create JBrowse

## I. Obtain Viral Genomes

1. **Obtain the H1N1 reference genome from NCBI:** [LINK](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_001343785.1/)

2. **Obtain genomes by segment from NCBI for multiple sequence alignment** (download as FASTA files):
   - **H1N1 Genome Assembly ViralMultiSegProj274766 by Segment**:
     - [Segment 1 (RefSeq ID: NC_026438.1)](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_001343785.1/)
     - [Segment 2 (RefSeq ID: NC_026435.1)](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_001343785.1/)
     - [Segment 3 (RefSeq ID: NC_026437.1)](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_001343785.1/)
     - [Segment 4 (RefSeq ID: NC_026433.1)](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_001343785.1/)
     - [Segment 5 (RefSeq ID: NC_026436.1)](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_001343785.1/)
     - [Segment 6 (RefSeq ID: NC_026434.1)](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_001343785.1/)
     - [Segment 7 (RefSeq ID: NC_026431.1)](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_001343785.1/)
     - [Segment 8 (RefSeq ID: NC_026432.1)](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_001343785.1/)

   - **H5N1 Genome Assembly ASM3868529v1 by Segment**:
     - [Segment 1 (GenBank ID: KJ907628.1)](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_038685295.1/)
     - [Segment 2 (GenBank ID: KJ907629.1)](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_038685295.1/)
     - [Segment 3 (GenBank ID: KJ907630.1)](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_038685295.1/)
     - [Segment 4 (GenBank ID: KJ907631.1)](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_038685295.1/)
     - [Segment 5 (GenBank ID: KJ907632.1)](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_038685295.1/)
     - [Segment 6 (GenBank ID: KJ907633.1)](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_038685295.1/)
     - [Segment 7 (GenBank ID: KJ907634.1)](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_038685295.1/)
     - [Segment 8 (GenBank ID: KJ907635.1)](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_038685295.1/)

   - **H1N2 Genome Assembly ASM3875498v1 by Segment**:
     - [Segment 1 (GenBank ID: KT225475.1)](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_038754985.1/)
     - [Segment 2 (GenBank ID: KT225474.1)](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_038754985.1/)
     - [Segment 3 (GenBank ID: KT225473.1)](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_038754985.1/)
     - [Segment 4 (GenBank ID: KT225468.1)](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_038754985.1/)
     - [Segment 5 (GenBank ID: KT225471.1)](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_038754985.1/)
     - [Segment 6 (GenBank ID: KT225470.1)](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_038754985.1/)
     - [Segment 7 (GenBank ID: KT225469.1)](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_038754985.1/)
     - [Segment 8 (GenBank ID: KT225472.1)](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_038754985.1/)

3. **Generate Multiple Sequence Alignments as BAM files**:
   - Install the following packages: `samtools`, `bwa`.
   - Run the following commands:
     ```bash
     cat h1n2_segment1.fa h5n1_segment1.fa > non_human.fa
     bwa index non_human.fasta
     samtools faidx non_human.fa
     samtools faidx h1n1_segment1.fa
     ```
   - Align merged FASTA to "reference" (H1N1 human version of segment):
     ```bash
     bwa mem h1n1_segment1.fa non_human.fasta > merged_seg#.sam
     ```

