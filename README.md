# README

## I. Obtain Viral Genomes

1. **Obtain the H1N1 Reference Genome**:
   - [Download link](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/343/785/GCF_001343785.1_ViralMultiSegProj274766/)
   - Download as a FASTA file and index it using `samtools`:
     ```bash
     samtools faidx your_h1n1_genome.fa
     ```

2. **Obtain H1N1 Annotation File**:
   - [Download link](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/343/785/GCF_001343785.1_ViralMultiSegProj274766/)

3. **Obtain Genomes for Multiple Sequence Alignment**:
   - Download as FASTA files via FTP:
     - **H1N1 Genome Assembly**: [Download file and unzip]
     - RefSeq ID: GCF_001343785.1
(https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/343/785/GCF_001343785.1_ViralMultiSegProj274766/)
       - Segment 1: NC_026438.1
       - Segment 2: NC_026435.1
       - Segment 3: NC_026437.1
       - Segment 4: NC_026433.1
       - Segment 5: NC_026436.1
       - Segment 6: NC_026434.1
       - Segment 7: NC_026431.1
       - Segment 8: NC_026432.1

     - **H5N1 Genome Assembly**: [Download file and unzip](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/038/685/295/GCA_038685295.1_ASM3868529v1/)
     - GenBank ID: GCA_038685295.1
       - Segment 1: KJ907628.1
       - Segment 2: KJ907629.1
       - Segment 3: KJ907630.1
       - Segment 4: KJ907631.1
       - Segment 5: KJ907632.1
       - Segment 6: KJ907633.1
       - Segment 7: KJ907634.1
       - Segment 8: KJ907635.1

     - **H1N2 Genome Assembly**: [Download file and unzip](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/038/754/985/GCA_038754985.1_ASM3875498v1/)
     - GenBank ID: GCA_038754985.1
       - Segment 1: KT225475.1
       - Segment 2: KT225474.1
       - Segment 3: KT225473.1
       - Segment 4: KT225468.1
       - Segment 5: KT225471.1
       - Segment 6: KT225470.1
       - Segment 7: KT225469.1
       - Segment 8: KT225472.1

4. **Generate Multiple Sequence Alignments as BAM Files**:
   - Install `samtools` and `bwa`:
     ```bash
     sudo apt install samtools bwa
     ```
   - Run the following commands:
     ```bash
     # Concatenate swine and duck segment files
     cat h1n2.fna h5n1.fna > non_human.fa

     # Index files
     bwa index non_human.fa
     bwa index h1n1.fna
     samtools faidx non_human.fa
     samtools faidx h1n1.fna

     # Align merged FASTA to reference (H1N1 human segment)
     bwa mem h1n1.fna non_human.fa > aligned_genomes.sam

     # Convert SAM to BAM
     samtools view -S -b aligned_genomes.sam > aligned_genomes.bam

     # Sort and index BAM
     samtools sort aligned_genomes.bam -o sorted_aligned.bam
     samtools index sorted_aligned.bam
     ```

## II. Generate Phylogenetic Tree

### 1. Download ~20 Genomes from NCBI
- Visit the [NCBI Genomes page](https://www.ncbi.nlm.nih.gov/datasets/genome/?taxon=11320&typical_only=true&search_text=human) to download approximately 20 genomes.

### 2. Align and Form the Tree in Terminal

1. **Enter the Directory Where Genomes are Stored**:
   - Navigate to the directory containing the downloaded genome FASTA files:
     ```bash
     cd /path/to/genomes/
     ```

2. **Combine All FASTA Files**:
   - Concatenate all `.fasta` files into one:
     ```bash
     cat *.fasta > combined_fasta.fasta
     ```

3. **Open a Linux Environment**:
   - Ensure you are in a Linux terminal or SSH into your AWS instance.

4. **Install `Mafft` and `FastTree`**:
   - Use the following command to install the necessary tools:
     ```bash
     sudo apt install mafft fasttree
     ```

5. **Perform Alignment**:
   - Align the combined FASTA file:
     ```bash
     mafft --auto combined_fasta.fasta > aligned_sequences.fasta
     ```

6. **Create the Phylogenetic Tree**:
   - Generate the tree using `FastTree`:
     ```bash
     fasttree -nt aligned_sequences.fasta > phylo_tree.nwk
     ```

---

### 3. Generate Tree Visualization on Interactive Tree of Life (iTOL)

1. **Upload the `.nwk` File**:
   - Go to [Interactive Tree of Life (iTOL)](https://itol.embl.de/) and upload the `phylo_tree.nwk` file.

2. **Add Annotations and Labels**:
   - Customize the tree by adding annotations and labels as needed.

3. **Download as a PDF**:
   - Save the visualized tree as a PDF file.

---

### 4. Create an HTML File for the Phylogenetic Tree

1. **Upload the PDF to Google Drive**:
   - Upload the PDF file to your Google Drive.
   - Share it with everyone and copy the shareable link.

2. **Create an HTML File**:
   - Open a text editor and create a file named `phylogenetic_tree.html` with the following content:
     ```html
     <!DOCTYPE html>
     <html>
     <head>
         <title>Phylogenetic Tree Viewer</title>
     </head>
     <body>
         <h1>Human Influenza Phylogenetic Tree</h1>
         <p>Click on the link below to view the detailed phylogenetic tree for human influenza:</p>
         <a href="https://drive.google.com/file/d/1Qvz1DAiVgQsVuPqRZI3vUBPJf03JdLBh/view?usp=sharing" target="_blank">View Phylogenetic Tree</a>
     </body>
     </html>
     ```

3. **Upload the HTML File to GitHub or Another Secure Server**:
   - Copy the HTML file to your GitHub Pages repository or another secure server.

---

### 5. Add the HTML to JBrowse

1. **Update the JSON File for the Track**:
   - Add the following entry to the JSON configuration file for JBrowse:
     ```json
     {
       "label": "CustomHTML",
       "key": "Custom HTML Track",
       "type": "JBrowse/View/Track/HTMLFeatures",
       "urlTemplate": "http://yourserver.com/path/to/phylogenetic_tree.html"
     }
     ```

2. **Reload JBrowse**:
   - Navigate to your JBrowse instance and verify that the custom HTML track is displayed correctly.

---

## III. Launch an AWS EC2 Instance

1. **Log in to your AWS account and launch an EC2 instance with the following configuration:**
   - **Name**: bioe131_final_project
   - **AMI**: Ubuntu 22.04
   - **Key Pair**: Create a new key pair (e.g., `bioe131.pem`).
   - **Network Settings**: Allow HTTP and HTTPS traffic.
   - **Storage Configuration**: Set storage to 30 GiB gp3.

2. **Connect to your instance. Press Connect.**

## IV. Set Up JBrowse2 in Terminal

1. **Switch to Root User**:
   - Run the following command to switch to the root user:
     ```bash
     sudo su -
     ```

2. **Set a Password for the Default User**:
   - Set a password for the default `ubuntu` user:
     ```bash
     passwd ubuntu
     ```
   - Enter a new password when prompted and save it securely.

3. **Exit Root User**:
   - Exit the root user session:
     ```bash
     exit
     ```

4. **Install Linuxbrew**:
   - Use the following command to install Linuxbrew:
     ```bash
     /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
     ```
   - Press Return/Enter to continue when prompted.
   - Enter the password you set earlier when prompted.

5. **Add Linuxbrew to Your Execution Path**:
   - Add Linuxbrew to your shell environment:
     ```bash
     echo 'eval "$(/home/linuxbrew/.linuxbrew/bin/brew shellenv)"' >> ~/.bashrc
     eval "$(/home/linuxbrew/.linuxbrew/bin/brew shellenv)"
     ```

6. **Install Node.js**:
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

7. **Install npm**:
   - Install `npm` (Node.js package manager):
     ```bash
     sudo apt update
     sudo apt install npm -y
     ```

8. **Install JBrowse CLI**:
   - Use `npm` to install the JBrowse CLI globally:
     ```bash
     sudo npm install -g @jbrowse/cli
     ```
   - Verify the installation:
     ```bash
     jbrowse --version
     ```

9. **Install Apache2 and Dependencies**:
    - Install Apache2 and required tools:
      ```bash
      sudo apt install wget apache2 -y
      brew install samtools htslib
      ```

10. **Start the Apache Server**:
    - Start the Apache2 service:
      ```bash
      sudo service apache2 start
      ```

11. **Create and Configure JBrowse2**:
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

12. **Verify the Installation**:
    - Open your browser and navigate to:
      ```
      http://<public_IP_of_your_instance>/jbrowse2/
      ```
    - You should see the message: **"JBrowse 2 is installed."**
   
---

## V. Add Files to AWS Server

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

## VI. Upload Files to JBrowse

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





