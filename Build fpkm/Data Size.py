import subprocess

# SRA Sample IDs with Descriptive Names
sample_ids = {
    "dmel": {
        "Dmel male TW rep1": "SRR5638358",
        "Dmel female TW rep1": "SRR5638359"
    },
    "dsec": {
        "Dsec male TW rep1": "SRR1952772",
        "Dsec male TW rep2": "SRR1952773",
        "Dsec male TW rep3": "SRR1952774",
        "Dsec female TW rep1": "SRR1952775",
        "Dsec female TW rep2": "SRR1952776",
        "Dsec female TW rep3": "SRR1952777",
        "Dsec male JP rep1": "SRR1973486",
        "Dsec male JP rep2": "SRR1973487",
        "Dsec male JP rep3": "SRR1973488",
        "Dsec female JP rep1": "SRR1973489",
        "Dsec female JP rep2": "SRR1973490",
        "Dsec female JP rep3": "SRR1973491"
    },
    "dsim": {
        "Dsim male JP rep1": "SRR1973492",
        "Dsim male JP rep2": "SRR1973493",
        "Dsim male JP rep3": "SRR1973494",
        "Dsim female JP rep1": "SRR1973495",
        "Dsim female JP rep2": "SRR1973496",
        "Dsim female JP rep3": "SRR1973497"
    }
}

# Function to get sizes for each SRA sample
def get_sample_sizes(sample_ids):
    results = {}
    for species, samples in sample_ids.items():
        results[species] = {}
        for name, sra_id in samples.items():
            try:
                # Run vdb-dump to fetch size
                command = f"vdb-dump {sra_id} --info"
                result = subprocess.run(command, shell=True, capture_output=True, text=True)
                if result.returncode == 0:
                    output = result.stdout
                    # Parse size in bytes
                    for line in output.splitlines():
                        if line.startswith("size   :"):
                            size_bytes = int(line.split(":")[1].strip().replace(",", ""))
                            size_gb = size_bytes / (1024 ** 3)  # Convert bytes to GB
                            results[species][name] = f"{round(size_gb, 2)} GB"
                else:
                    results[species][name] = "Error"
            except Exception as e:
                results[species][name] = f"Error: {e}"
    return results

# Get and print sizes
sizes = get_sample_sizes(sample_ids)
for species, samples in sizes.items():
    print(f"--- {species.upper()} ---")
    for name, size in samples.items():
        print(f"{name}: {size}")

# Calculate and display total size for all datasets
def calculate_total_size(sizes):
    total_size_gb = 0
    for species, samples in sizes.items():
        for size in samples.values():
            try:
                # Only add valid numeric sizes
                total_size_gb += float(size.split()[0])
            except ValueError:
                continue
    print(f"\nTotal Size: {round(total_size_gb, 2)} GB")

# Calculate total size
calculate_total_size(sizes)
