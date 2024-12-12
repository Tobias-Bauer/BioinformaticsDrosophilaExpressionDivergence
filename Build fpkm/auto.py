import os
import subprocess

# Define paths and file lists
output_dir = "rna_seq_output"
genome_dir = os.path.join(output_dir, "genomes")
annotation_dir = os.path.join(output_dir, "annotations")
fastq_dir = os.path.join(output_dir, "fastq")
genome_index_dir = os.path.join(output_dir, "star_genome_indices")
threads = 8

# URLs for genome and annotation files
genome_files = {
    "dmel": {
        "fasta": "ftp://ftp.flybase.net/releases/FB2012_06/dmel_r5.48/fasta/dmel-all-chromosome-r5.48.fasta.gz",
        "gff": "ftp://ftp.flybase.net/releases/FB2012_06/dmel_r5.48/gff/dmel-all-r5.48.gff.gz"
    },
    "dsec": {
        "fasta": "ftp://ftp.flybase.net/releases/FB2012_06/dsec_r1.3/fasta/dsec-all-chromosome-r1.3.fasta.gz",
        "gff": "ftp://ftp.flybase.net/releases/FB2012_06/dsec_r1.3/gff/dsec-all-r1.3.gff.gz"
    },
    "dsim": {
        "fasta": "https://ftp.flybase.net/releases/FB2012_06/dsim_r1.4/fasta/dsim-all-chromosome-r1.4.fasta.gz",
        "gff": "https://ftp.flybase.net/releases/FB2012_06/dsim_r1.4/gff/dsim-all-r1.4.gff.gz"
    }
}

# Sample SRA URLs
sample_ids = {
    # Dmel male TW rep1, Dmel female TW rep1
    "dmel": ["SRR5638358", "SRR5638359"],
    # Dsec male TW rep1, Dsec male TW rep2, Dsec male TW rep3, Dsec female TW rep1, Dsec female TW rep2, Dsec female TW rep3, Dsec male JP rep1, Dsec male JP rep2, Dsec male JP rep3, Dsec female JP rep1, Dsec female JP rep2, Dsec female JP rep3
    "dsec": ["SRR1952772", "SRR1952773", "SRR1952774", "SRR1952775", "SRR1952776", "SRR1952777", "SRR1973486", "SRR1973487", "SRR1973488", "SRR1973489", "SRR1973490", "SRR1973491"],
    # Dsim male JP rep1, Dsim male JP rep2, Dsim male JP rep3, Dsim female JP rep1, Dsim female JP rep2, Dsim female JP rep3
    "dsim": ["SRR1973492", "SRR1973493", "SRR1973494", "SRR1973495", "SRR1973496", "SRR1973497"]
}

# Ensure output directories exist
os.makedirs(genome_dir, exist_ok=True)
os.makedirs(annotation_dir, exist_ok=True)
os.makedirs(fastq_dir, exist_ok=True)
os.makedirs(genome_index_dir, exist_ok=True)

def run_command(command, error_message):
    """Utility function to run a shell command."""
    try:
        subprocess.run(command, shell=True, check=True)
    except subprocess.CalledProcessError:
        print(f"Error: {error_message}")
        exit(1)

# Download genome and annotation files
print("Downloading genome and annotation files...")
file_paths = {}
for species, urls in genome_files.items():
    file_paths[species] = {}
    for file_type, url in urls.items():
        target_dir = genome_dir if file_type == "fasta" else annotation_dir
        file_name = url.split('/')[-1]
        file_path = os.path.join(target_dir, file_name)
        file_paths[species][file_type] = file_path.replace('.gz', '')

        if not os.path.exists(file_path):
            run_command(f"wget -P {target_dir} {url}", f"Failed to download {file_type} for {species}.")
        else:
            print(f"{file_name} already exists, skipping download.")

        if not os.path.exists(file_paths[species][file_type]):
            run_command(f"gunzip -c {file_path} > {file_paths[species][file_type]}", f"Failed to unzip {file_type} for {species}.")
        else:
            print(f"{file_paths[species][file_type]} already exists, skipping unzip.")

# Generate genome indices
for species, paths in file_paths.items():
    index_dir = os.path.join(genome_index_dir, species)
    os.makedirs(index_dir, exist_ok=True)
    if not os.listdir(index_dir):
        print(f"Indexing genome for {species}...")
        run_command(
            f"STAR --runThreadN {threads} --runMode genomeGenerate "
            f"--genomeDir {index_dir} "
            f"--genomeFastaFiles {paths['fasta']} "
            f"--sjdbGTFfile {paths['gff']} "
            f"--sjdbOverhang 100 "
            f"--genomeSAindexNbases 12 "
            f"--limitGenomeGenerateRAM 16000000000 ",
            f"Failed to index genome for {species}."
        )
    else:
        print(f"Genome index for {species} already exists, skipping indexing.")

# Process each sample
for species, sra_ids in sample_ids.items():
    index_dir = os.path.join(genome_index_dir, species)
    for sample_id in sra_ids:
        print(f"Processing sample {sample_id} for species {species}...")

        # Download FASTQ files using prefetch and fastq-dump
        fastq_file_1 = os.path.join(fastq_dir, f"{sample_id}_1.fastq.gz")
        fastq_file_2 = os.path.join(fastq_dir, f"{sample_id}_2.fastq.gz")

        # Fetch the sample using prefetch
        if not os.path.exists(fastq_file_1) or not os.path.exists(fastq_file_2):
            run_command(f"prefetch {sample_id}", f"Failed to fetch SRA file {sample_id}.")

            # Dump the FASTQ files using fastq-dump
            run_command(
                f"fastq-dump --split-files --gzip --outdir {fastq_dir} {sample_id}",
                f"Failed to convert SRA file {sample_id} to FASTQ."
            )
        else:
            print(f"FASTQ files for {sample_id} already exist, skipping download.")

        # Mapping with STAR
        print(f"Running STAR for {sample_id}...")
        bam_file = os.path.join(output_dir, f"{sample_id}_Aligned.sortedByCoord.out.bam")
        run_command(
            f"STAR --runThreadN {threads} --genomeDir {index_dir} --readFilesIn {fastq_file_1} {fastq_file_2} "
            f"--readFilesCommand zcat --outFileNamePrefix {os.path.join(output_dir, sample_id + '_')} "
            f"--outSAMtype BAM SortedByCoordinate "
            f"--outSAMattributes NH HI AS nM MD XS ",
            f"Failed to process {sample_id}."
        )

        # Calculate FPKM with Cufflinks
        gff_file = file_paths[species]["gff"]
        run_command(
            f"cufflinks -p {threads} -G {gff_file} -o {os.path.join(output_dir, sample_id + '_cufflinks')} {bam_file}",
            f"Failed to calculate FPKM for {sample_id}."
        )

        # Delete FASTQ files to save space
        if os.path.exists(fastq_file_1):
            os.remove(fastq_file_1)
            print(f"Deleted {fastq_file_1} to save space.")

        if os.path.exists(fastq_file_2):
            os.remove(fastq_file_2)
            print(f"Deleted {fastq_file_2} to save space.")

print("Pipeline completed successfully.")