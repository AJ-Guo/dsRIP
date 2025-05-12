import os
import pandas as pd
import re
from Bio import SeqIO

constant_files = "/home/dcedden/PycharmProjects/dsRIP_web/main_site/constant_files/"
dsRIP_orthology_path = os.path.join(constant_files, "dsRIP_orthology/")

def extract_tc_ids(tc_id_string):
    # This function extracts 'TC' followed by 6 numbers from a given string
    return re.findall(r'TC\d{6}', tc_id_string)

def clean_gene_id(gene_id):
    # Removes the '.pX' part at the end of gene IDs
    return re.sub(r'\.p\d+$', '', gene_id)

def filter_species_genes(species_name, excel_file=(constant_files + "lethal_gene_threshold.xlsx"), list_criteria=None, go_criteria=None, kegg_criteria=None, without_paralog=False):
    # Load gene threshold data from Excel
    gene_thresholds = pd.read_excel(excel_file, engine='openpyxl')
    gene_thresholds['Tribolium ID'] = gene_thresholds['Tribolium ID'].apply(lambda x: x.split(' ')[0])  # Clean up IDs

    # Convert 'list' to integer for comparison and filter based on the criteria
    gene_thresholds['list'] = pd.to_numeric(gene_thresholds['list'], errors='coerce')
    if list_criteria is not None:
        gene_thresholds = gene_thresholds[gene_thresholds['list'] >= list_criteria]

    # Filter by KEGG if criteria are provided
    if kegg_criteria is not None and isinstance(kegg_criteria, list):
        kegg_pattern = '|'.join(kegg_criteria)
        gene_thresholds = gene_thresholds[gene_thresholds['KEGG'].str.contains(kegg_pattern, na=False)]

    tribolium_ids = set(gene_thresholds['Tribolium ID'])

    # Find the correct .tsv file in the dsRIP_orthology directory
    species_file = next((os.path.join(dsRIP_orthology_path, file)
                         for file in os.listdir(dsRIP_orthology_path)
                         if species_name in file and file.endswith('.tsv')), None)

    if not species_file:
        print("No file found for the given species.")
        return None

    # Read the species .tsv file
    species_data = pd.read_csv(species_file, sep="\t", header=0)

    # Collect rows for final output to avoid repeated DataFrame concatenations
    output_rows = []

    # Process each row in the .tsv file
    for _, row in species_data.iterrows():
        tribolium_column = row.iloc[-2]
        species_gene_id = clean_gene_id(row.iloc[-1])  # Clean the species gene ID

        # Extract TC IDs and match with threshold data
        tc_ids = extract_tc_ids(tribolium_column)
        matched_ids = [tc_id for tc_id in tc_ids if tc_id in tribolium_ids]

        for tc_id in matched_ids:
            # Fetch annotation for the current Tribolium ID
            annotation_data = gene_thresholds[gene_thresholds['Tribolium ID'] == tc_id].iloc[0]

            gene_ids = species_gene_id.split(', ')
            for gene_id in gene_ids:
                if without_paralog and len(gene_ids) > 1:
                    continue

                output_rows.append({
                    'Annotation': annotation_data['Annotation'],
                    'Tribolium Gene ID': tc_id,
                    'Pest Gene ID': clean_gene_id(gene_id),
                    'Number of transfers': annotation_data['list'],
                    'GO': annotation_data['GO'],
                    'KEGG': annotation_data['KEGG'],
                    'Have_Paralogs': len(gene_ids) > 1
                })

    # Create DataFrame from collected rows
    output_data = pd.DataFrame(output_rows)

    output_data = output_data.sort_values(by='Number of transfers', ascending=False)

    # Remove duplicate rows based on 'Pest Gene ID'
    output_data = output_data.drop_duplicates(subset='Pest Gene ID')

    return output_data

def filter_fasta_by_ids(species_name, output_dir, list_criteria, kegg_criteria, paralog, pest_corr_dir=constant_files + "/pest_corr"):
    """
    Find the .fna file for a given species, filter it to keep only specified gene IDs,
    and save the result as a .fasta file in the specified output directory.

    Args:
    species_name (str): Name of the species to filter.
    list_criteria (int or None): Criteria for list threshold.
    kegg_criteria (list or None): KEGG criteria for filtering.
    paralog (bool): Whether to include paralogs.
    pest_corr_dir (str): Directory where the .fna files are stored.
    output_dir (str): Directory where the filtered .fasta files should be saved.
    """
    # Clear the output directory before proceeding
    for filename in os.listdir(output_dir):
        file_path = os.path.join(output_dir, filename)  # Get the full path of the file
        if os.path.isfile(file_path) or os.path.islink(file_path):  # Check if it's a file or a symlink
            os.remove(file_path)  # Remove the file

    # Get filtered Tribolium IDs
    filtered_ids = filter_species_genes(species_name, constant_files + "lethal_gene_threshold.xlsx", list_criteria=list_criteria, kegg_criteria=kegg_criteria, without_paralog=paralog)

    if filtered_ids is None or filtered_ids.empty:
        print("No matching IDs found.")
        return

    # Extract the relevant IDs from the Pest Gene ID column
    ids_to_filter = set(filtered_ids['Pest Gene ID'])

    # Locate the .fna file for the species
    fna_file = next((os.path.join(pest_corr_dir, file)
                     for file in os.listdir(pest_corr_dir)
                     if file.startswith(species_name) and file.endswith('.fna')), None)

    if not fna_file:
        print("No FNA file found for the given species.")
        return

    # Prepare output file path
    output_fasta = os.path.join(output_dir, f"{species_name}_target_genes.fasta")

    # Filter sequences without using the entire FASTA index in memory
    with open(fna_file, "r") as infile, open(output_fasta, "w") as outfile:
        write_flag = False
        for line in infile:
            if line.startswith(">"):
                header = line[1:].strip().split()[0]  # Extract the header ID
                write_flag = header in ids_to_filter  # Set flag if header matches filtered IDs
            if write_flag:
                outfile.write(line)

    return filtered_ids


def extract_sclero_ids(tc_id_string):
    # This function extracts gene names starting with 'XM_' followed by exactly 11 characters (total length 14)
    return re.findall(r'XM_\S{11}', tc_id_string)

def filter_fasta_by_ids_pathogen(species_name, output_dir, list_criteria=0, kegg_criteria=[], paralog=False, pest_corr_dir=constant_files + "/pest_corr"):
    """
    Find the .fna file for a given species, filter it to keep only specified gene IDs,
    and save the result as a .fasta file in the specified output directory.

    Args:
    species_name (str): Name of the species to filter.
    list_criteria (int or None): Criteria for list threshold.
    kegg_criteria (list or None): KEGG criteria for filtering.
    paralog (bool): Whether to include paralogs.
    pest_corr_dir (str): Directory where the .fna files are stored.
    output_dir (str): Directory where the filtered .fasta files should be saved.
    """
    # Clear the output directory before proceeding
    for filename in os.listdir(output_dir):
        file_path = os.path.join(output_dir, filename)  # Get the full path of the file
        if os.path.isfile(file_path) or os.path.islink(file_path):  # Check if it's a file or a symlink
            os.remove(file_path)  # Remove the file

    # Get filtered Tribolium IDs
    filtered_ids = filter_species_genes_pathogen(species_name, constant_files + "lethal_gene_threshold_pathogen.xlsx")

    if filtered_ids is None or filtered_ids.empty:
        print("No matching IDs found.")
        return

    # Extract the relevant IDs from the Pest Gene ID column
    ids_to_filter = set(filtered_ids['Pest Gene ID'])

    # Locate the .fna file for the species
    fna_file = next((os.path.join(pest_corr_dir, file)
                     for file in os.listdir(pest_corr_dir)
                     if file.startswith(species_name) and file.endswith('.fna')), None)

    if not fna_file:
        print("No FNA file found for the given species.")
        return

    # Prepare output file path
    output_fasta = os.path.join(output_dir, f"{species_name}_target_genes.fasta")

    # Filter sequences without using the entire FASTA index in memory
    with open(fna_file, "r") as infile, open(output_fasta, "w") as outfile:
        write_flag = False
        for line in infile:
            if line.startswith(">"):
                header = line[1:].strip().split()[0]  # Extract the header ID
                write_flag = header in ids_to_filter  # Set flag if header matches filtered IDs
            if write_flag:
                outfile.write(line)

    return filtered_ids



def filter_species_genes_pathogen(species_name, excel_file=(constant_files + "lethal_gene_threshold_pathogen.xlsx")):
    # Load gene threshold data from Excel
    gene_thresholds = pd.read_excel(excel_file, engine='openpyxl')
    gene_thresholds['Sclerotinia'] = gene_thresholds['Sclerotinia'].apply(lambda x: x.split(' ')[0])  # Clean up IDs

    # Convert 'list' to integer for comparison and filter based on the criteria
    gene_thresholds['list'] = pd.to_numeric(gene_thresholds['list'], errors='coerce')

    gene_thresholds = gene_thresholds[gene_thresholds['list'] >= 0]



    tribolium_ids = set(gene_thresholds['Sclerotinia'])

    # Find the correct .tsv file in the dsRIP_orthology directory
    species_file = next((os.path.join((dsRIP_orthology_path + "fungi"), file)
                         for file in os.listdir((dsRIP_orthology_path + "fungi"))
                         if species_name in file and file.endswith('.tsv')), None)

    if not species_file:
        print("No file found for the given species.")
        return None

    # Read the species .tsv file
    species_data = pd.read_csv(species_file, sep="\t", header=0)

    # Collect rows for final output to avoid repeated DataFrame concatenations
    output_rows = []

    # Process each row in the .tsv file
    for _, row in species_data.iterrows():
        tribolium_column = row.iloc[-2]
        species_gene_id = clean_gene_id(row.iloc[-1])  # Clean the species gene ID

        # Extract TC IDs and match with threshold data
        tc_ids = extract_sclero_ids(tribolium_column)
        #print(tc_ids)
        matched_ids = [tc_id for tc_id in tc_ids if tc_id in tribolium_ids]
        #print(matched_ids)
        for tc_id in matched_ids:

            # Fetch annotation for the current Tribolium ID
            annotation_data = gene_thresholds[gene_thresholds['Sclerotinia'] == tc_id].iloc[0]


            gene_ids = species_gene_id.split(', ')
            for gene_id in gene_ids:
                print(gene_id)

                output_rows.append({
                    'Annotation': annotation_data['Annotation'],
                    'Sclerotinia ID': tc_id,
                    'Pest Gene ID': clean_gene_id(gene_id),
                    'Have_Paralogs': len(gene_ids) > 1
                })

    # Create DataFrame from collected rows
    output_data = pd.DataFrame(output_rows)

    # Remove duplicate rows based on 'Pest Gene ID'
    output_data = output_data.drop_duplicates(subset='Pest Gene ID')

    return output_data







