import os
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import glob
from scipy.cluster.hierarchy import linkage, dendrogram

# Step 1: Load the orthologous genes file
orthologs_file = "../../all_3spe-exist+noDup.csv"
orthologs = pd.read_csv(orthologs_file)

# Extract gene families
orthologs['#gene'] = orthologs['#gene'].fillna('')
ir_genes = orthologs[orthologs['#gene'].str.match(r'^Ir\d.*')]
obp_genes = orthologs[orthologs['#gene'].str.startswith('Obp')]
or_genes = orthologs[orthologs['#gene'].str.startswith('Or')]
gr_genes = orthologs[orthologs['#gene'].str.startswith('Gr')]

# Function to prepare gene family data for any species
def get_family_ids_and_names(gene_family_df, species='Dsec'):
    gene_id_name_map = {}
    for _, row in gene_family_df.iterrows():
        if pd.notnull(row[species]):  # Dynamically use the species passed as argument
            gene_id_name_map[row[species]] = row['#gene']
    return gene_id_name_map

# Create ID-to-name mappings for each family for all species
species_list = ['Dsec', 'Dsim', 'Dmel']
ir_map = {species: get_family_ids_and_names(ir_genes, species) for species in species_list}
obp_map = {species: get_family_ids_and_names(obp_genes, species) for species in species_list}
or_map = {species: get_family_ids_and_names(or_genes, species) for species in species_list}
gr_map = {species: get_family_ids_and_names(gr_genes, species) for species in species_list}

# Step 2: Load FPKM tracking files
fpkm_files = glob.glob("GSM*.fpkm_tracking")

# Group files by species, gender, and country based on filenames
def group_files_by_species_gender_and_country(files):
    grouped_files = {}
    for file in files:
        parts = file.split("_")
        species = parts[2]  # Species is always in position 2
        gender = "Male" if "Mal" in parts else "Female"
        country = parts[4]  # Country is at position 4
        key = f"{species}-{gender}-{country}"
        if key not in grouped_files:
            grouped_files[key] = []
        grouped_files[key].append(file)
    return grouped_files

grouped_files = group_files_by_species_gender_and_country(fpkm_files)

# Load all FPKM data
def load_fpkm_file(file, gene_map):
    try:
        df = pd.read_csv(file, sep="\t")
        print(df)
        # Filter for relevant genes
        df = df[df['gene_id'].isin(gene_map.keys())]
        df['gene_name'] = df['gene_id'].map(gene_map)
        # Extract sample name from filename
        sample_name = "_".join(file.split("_")[2:4])  # Adjust based on filename structure
        df.rename(columns={'FPKM': sample_name}, inplace=True)
        df = df[['gene_name', sample_name]].set_index('gene_name')
        return df
    except Exception as e:
        print(f"Error processing file {file}: {e}")
        return pd.DataFrame()

# Step 3: Combine data for each group (species-gender-country)
def prepare_combined_data(gene_maps, grouped_files):
    combined_data = {}
    for group, files in grouped_files.items():
        group_data = []
        species = group.split('-')[0]  # Extract species from the group key
        gene_map = gene_maps[species]  # Get the gene map for the species
        for file in files:
            data = load_fpkm_file(file, gene_map)
            if not data.empty:
                group_data.append(data)
        # Average FPKM values across replicates for each group
        if group_data:
            combined_data[group] = pd.concat(group_data, axis=1).mean(axis=1).rename(group)
    # Combine all groups into a single DataFrame
    if combined_data:
        final_data = pd.concat(combined_data.values(), axis=1)
        return final_data
    else:
        return pd.DataFrame()

# Get data for each family
ir_data = prepare_combined_data(ir_map, grouped_files)
obp_data = prepare_combined_data(obp_map, grouped_files)
or_data = prepare_combined_data(or_map, grouped_files)
gr_data = prepare_combined_data(gr_map, grouped_files)


# Step 4: Filter genes with FPKM < 0.05 in all species
def filter_genes_by_expression_threshold(*datasets, threshold=0.05):
    """
    Filter out genes that have FPKM < threshold in all species.
    """
    combined_data = pd.concat(datasets, axis=1)
    # Keep genes that have at least one FPKM >= threshold in any species
    filtered_data = combined_data.loc[combined_data.ge(threshold).any(axis=1)]
    return filtered_data

# Apply the filter to keep only genes with FPKM >= 0.05 in at least one species
filtered_ir_data = filter_genes_by_expression_threshold(ir_data)
filtered_obp_data = filter_genes_by_expression_threshold(obp_data)
filtered_or_data = filter_genes_by_expression_threshold(or_data)
print(len(filtered_ir_data))
#print("\n".join(filtered_ir_data.index))
print(len(filtered_obp_data))
#print("\n".join(filtered_obp_data.index))
print(len(filtered_or_data.index))
print("\n".join(filtered_or_data.index))

# Step 5: Rank All Genes Across All Data
def rank_all_genes_across_data(*datasets):
    """
    Rank genes across all gene families globally, not just within each family.
    """
    # Concatenate all data for global ranking
    combined_data = pd.concat(datasets, axis=0)
    
    # Rank all genes across all species globally
    ranked_data = combined_data.rank(axis=0, method='min', ascending=True).fillna(1)  # FPKM=0 as rank 1
    return ranked_data

print(filtered_ir_data)
# Rank all genes (combine ir_data, obp_data, or_data for global ranking)
all_ranked_data = rank_all_genes_across_data(filtered_ir_data, filtered_obp_data, filtered_or_data, gr_data)

# Step 6: Apply Global Ranking to Family Data
def get_family_ranked_data(family_data, all_ranked_data):
    """
    Apply global rankings to each family (ir_data, obp_data, or_data).
    """
    return all_ranked_data.loc[family_data.index]

# Apply the global ranking to each family
ir_family_ranked = get_family_ranked_data(filtered_ir_data, all_ranked_data)
obp_family_ranked = get_family_ranked_data(filtered_obp_data, all_ranked_data)
or_family_ranked = get_family_ranked_data(filtered_or_data, all_ranked_data)

# Step 7: Heatmap and Dendrogram
def plot_heatmap_and_dendrogram(data, title, filename):
    if data.empty:
        print(f"No data to plot for {title}")
        return

    data_clean = data.fillna(0)

    # Heatmap
    plt.figure(figsize=(12, 8))
    sns.heatmap(data_clean, cmap="YlOrRd", xticklabels=True, yticklabels=True, cbar_kws={'label': 'Rank'})
    plt.title(f"{title} Heatmap")
    plt.xlabel('Groups')
    plt.ylabel('Genes')
    plt.tight_layout()
    plt.savefig(f"{filename}_heatmap.png", dpi=300)
    plt.show()

    # Dendrogram
    plt.figure(figsize=(10, 8))
    linkage_matrix = linkage(data_clean.T, method='complete', metric='euclidean')
    dendrogram(linkage_matrix, labels=data.columns, leaf_rotation=90)
    plt.title(f"{title} Dendrogram")
    plt.tight_layout()
    plt.savefig(f"{filename}_dendrogram.png", dpi=300)
    plt.show()

# Step 8: Plot for each gene family
plot_heatmap_and_dendrogram(ir_family_ranked, "IR Gene Family", "IR_gene_family")
plot_heatmap_and_dendrogram(obp_family_ranked, "OBP Gene Family", "OBP_gene_family")
plot_heatmap_and_dendrogram(or_family_ranked, "OR Gene Family", "OR_gene_family")
# Export ranked data to CSV for each gene family
ir_family_ranked.to_csv("ir_gene_family_ranked.csv")
obp_family_ranked.to_csv("obp_gene_family_ranked.csv")
or_family_ranked.to_csv("or_gene_family_ranked.csv")
