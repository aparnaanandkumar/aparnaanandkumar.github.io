# README

## I. Obtain Viral Genomes
1. **Obtain the H1N1 reference genome from NCBI:** [LINK](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_001343785.1/)
2. **Obtain H1N1 annotation file:** [LINK](https://www.ncbi.nlm.nih.gov/)
3. **Obtain genomes by segment from NCBI for multiple sequence alignment** (download as FASTA files):
   - **H1N1 Genome Assembly ViralMultiSegProj274766 by Segment:**
     - [Segment 1 (RefSeq ID: NC_026438.1)](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_001343785.1/)
     - [Segment 2 (RefSeq ID: NC_026435.1)](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_001343785.1/)
     - [Segment 3 (RefSeq ID: NC_026437.1)](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_001343785.1/)
     - [Segment 4 (RefSeq ID: NC_026433.1)](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_001343785.1/)
     - [Segment 5 (RefSeq ID: NC_026436.1)](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_001343785.1/)
     - [Segment 6 (RefSeq ID: NC_026434.1)](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_001343785.1/)
     - [Segment 7 (RefSeq ID: NC_026431.1)](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_001343785.1/)
     - [Segment 8 (RefSeq ID: NC_026432.1)](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_001343785.1/)

   - **H5N1 Genome Assembly ASM3868529v1 by Segment:**
     - [Segment 1 (GenBank ID: KJ907628.1)](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_038685295.1/)
     - [Segment 2 (GenBank ID: KJ907629.1)](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_038685295.1/)
     - [Segment 3 (GenBank ID: KJ907630.1)](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_038685295.1/)
     - [Segment 4 (GenBank ID: KJ907631.1)](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_038685295.1/)
     - [Segment 5 (GenBank ID: KJ907632.1)](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_038685295.1/)
     - [Segment 6 (GenBank ID: KJ907633.1)](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_038685295.1/)
     - [Segment 7 (GenBank ID: KJ907634.1)](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_038685295.1/)
     - [Segment 8 (GenBank ID: KJ907635.1)](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_038685295.1/)

   - **H1N2 Genome Assembly ASM3875498v1 by Segment:**
     - [Segment 1 (GenBank ID: KT225475.1)](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_038754985.1/)
     - [Segment 2 (GenBank ID: KT225474.1)](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_038754985.1/)
     - [Segment 3 (GenBank ID: KT225473.1)](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_038754985.1/)
     - [Segment 4 (GenBank ID: KT225468.1)](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_038754985.1/)
     - [Segment 5 (GenBank ID: KT225471.1)](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_038754985.1/)
     - [Segment 6 (GenBank ID: KT225470.1)](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_038754985.1/)
     - [Segment 7 (GenBank ID: KT225469.1)](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_038754985.1/)
     - [Segment 8 (GenBank ID: KT225472.1)](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_038754985.1/)

4. **Generate Multiple Sequence Alignments as BAM files:**
   - Install the following packages: `samtools`, `bwa`.
   - Run the following commands to generate MSA BAM files (repeat the procedure below for all 8 segments).
   - The final outputs to upload to Jbrowse will be a sorted bam file and bam index file (e.g. sorted_seg1.bam and sorted_seg1.bam.bai):
     ```bash
     # Concatenate swine and duck segment files together
     cat h1n2_segment1.fa h5n1_segment1.fa > non_human.fa
     
     # Index files using bwa and samtools
     bwa index non_human.fasta
     bwa index h1n1_segment1.fa
     samtools faidx non_human.fa
     samtools faidx h1n1_segment1.fa
     
     # Align merged FASTA to the reference (H1N1 human segment)
     bwa mem h1n1_segment1.fa non_human.fa > aligned_segment1.sam
     
     # Convert the SAM file to BAM
     samtools view -S -b aligned_segment1.sam > seg1.bam
     
     # Sort and index the BAM file
     samtools sort seg1.bam -o sorted_seg1.bam
     samtools index sorted_seg1.bam
     ```
5. **Generate Phylogenetic Tree:**
---

## II. Launch an AWS EC2 Instance
1. **Log in to your AWS account and launch an EC2 instance:**
   - **Name:** Choose a name, such as bioe131_final_project.
   - **AMI:** Select Ubuntu 22.04 as the operating system.
   - **Key Pair:** Create a new key pair (e.g., bioe131) and save it securely.
   - **Network Settings:** Allow HTTP and HTTPS traffic.
   - **Storage Configuration:** Set storage to 30 GiB gp3.

2. **Once the instance is launched, click Connect to get connection instructions.**


## III. Set Up JBrowse2 in Terminal

1. **Connect to Your EC2 Instance via SSH**:
   - Follow the instructions from the AWS EC2 console to connect to your instance:
     ```bash
     ssh -i "your-key.pem" ubuntu@<public_IP>
     ```

2. **Switch to Root User**:
   - Run the following command to switch to the root user:
     ```bash
     sudo su -
     ```

3. **Set a Password for the Default User**:
   - Set a password for the default `ubuntu` user:
     ```bash
     passwd ubuntu
     ```
   - Enter a new password when prompted and save it securely.

4. **Exit Root User**:
   - Exit the root user session:
     ```bash
     exit
     ```

5. **Install Linuxbrew**:
   - Use the following command to install Linuxbrew:
     ```bash
     /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
     ```
   - Press Return/Enter to continue when prompted.
   - Enter the password you set earlier when prompted.

6. **Add Linuxbrew to Your Execution Path**:
   - Add Linuxbrew to your shell environment:
     ```bash
     echo 'eval "$(/home/linuxbrew/.linuxbrew/bin/brew shellenv)"' >> ~/.bashrc
     eval "$(/home/linuxbrew/.linuxbrew/bin/brew shellenv)"
     ```

7. **Install Node.js**:
   - Check if Node.js is installed:
     ```bash
     node -v
     ```
   - If Node.js is not installed, use the following commands to install it:
     ```bash
     sudo apt update
     sudo apt install unzip -y
     curl -fsSL https://fnm.vercel.app/install | bash
     source ~/.bashrc
     fnm use --install-if-missing 20
     ```
   - Verify the installation:
     ```bash
     node -v # Should print v20.x.x
     npm -v  # Should print the npm version
     ```

8. **Install npm**:
   - Install `npm` (Node.js package manager):
     ```bash
     sudo apt update
     sudo apt install npm -y
     ```

9. **Install JBrowse CLI**:
   - Use `npm` to install the JBrowse CLI globally:
     ```bash
     sudo npm install -g @jbrowse/cli
     ```
   - Verify the installation:
     ```bash
     jbrowse --version
     ```

10. **Install Apache2 and Dependencies**:
    - Install Apache2 and required tools:
      ```bash
      sudo apt install wget apache2 -y
      brew install samtools htslib
      ```

11. **Start the Apache Server**:
    - Start the Apache2 service:
      ```bash
      sudo service apache2 start
      ```

12. **Create and Configure JBrowse2**:
    - Create a temporary working directory:
      ```bash
      mkdir ~/tmp
      cd ~/tmp
      ```
    - Create a JBrowse2 app:
      ```bash
      jbrowse create output_folder
      ```
    - Move the JBrowse2 app to the Apache root directory:
      ```bash
      sudo mv output_folder /var/www/html/jbrowse2
      sudo chown -R $(whoami) /var/www/html/jbrowse2
      ```

13. **Verify the Installation**:
    - Open your browser and navigate to:
      ```
      http://<public_IP_of_your_instance>/jbrowse2/
      ```
    - You should see the message: **"JBrowse 2 is installed."**
   
---

## IV. Add Files to AWS Server

1. **Download FASTA Files**:
   - Add the FASTA files to the same directory as your `.pem` file for connecting to the AWS instance.

2. **Combine and Index Each Subtype (if needed)**:
   - Combine FASTA files into a single file:
     ```bash
     cat f1.fasta f2.fasta > combined.fasta
     ```

3. **Upload Files to the Server**:
   - Use `scp` to upload the files to your AWS instance:
     ```bash
     scp -i "your-key.pem" combined.fasta ubuntu@<public_IP>:/home/ubuntu/
     ```

---

## V. Upload Files to JBrowse

1. **Index the FASTA File**:
   - Use `samtools` to index the FASTA file on your AWS instance:
     ```bash
     samtools faidx combined.fasta
     ```

2. **Upload the FASTA File to JBrowse**:
   - Add the combined FASTA file as an assembly in JBrowse:
     ```bash
     sudo jbrowse add-assembly combined.fasta --out /var/www/html/jbrowse2 --load copy
     ```

3. **Upload Annotations**:
   - Add annotations to JBrowse using the following command:
     ```bash
     sudo jbrowse add-track annotations.gz --assemblyNames combined.fasta --out /var/www/html/jbrowse2 --load copy
     ```

---

## IV. Add Files to AWS Server

1. **Download FASTA Files**:
   - Add the FASTA files to the same directory as your `.pem` file for connecting to the AWS instance.

---

## V. Upload Files to JBrowse

1.  **Upload the FASTA File to JBrowse**:
   - Add the combined FASTA file as an assembly in JBrowse:
     ```bash
     sudo jbrowse add-assembly <combined.fasta> --out $APACHE_ROOT/jbrowse2 --load copy
     ```

3. **Upload Annotations**:
   - Add annotations to JBrowse using the following command:
     ```bash
     sudo jbrowse add-track annotations.gz --assemblyNames combined.fasta --out /var/www/html/jbrowse2 --load copy
     ```

---

## VI. Launch JBrowse

1. **Open the JBrowse2 Interface**:
   - Navigate to the following URL in your browser:
     ```
     http://<public_IP_of_your_instance>/jbrowse2/
     ```

2. **Verify Your Uploads**:
   - Ensure that the assembly and tracks are visible in the JBrowse2 interface.
   - You should see the genomic data and annotations displayed correctly.





