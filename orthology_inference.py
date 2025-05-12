import re
from openpyxl import Workbook
import csv
import gc
import os
import shutil
import os
import pandas as pd
import subprocess
import mmap
from dsRNA_efficiency import check_valid_dna
##### ADD to the filter to remove NNNNN

absolute_path_main_site = "/home/dcedden/PycharmProjects/dsRIP_web/main_site/"
constant_files = "/home/dcedden/PycharmProjects/dsRIP_web/main_site/constant_files/"

def extract_tc_ids(tc_id_string):
    # This function extracts 'TC' followed by 6 numbers from a given string
    return re.findall(r'TC\d{6}', tc_id_string)

def clean_gene_id(gene_id):
    # Removes the '.pX' part at the end of gene IDs
    return re.sub(r'\.p\d+$', '', gene_id)

def find_latest_results_dir(user_orthofinder_dir):
    # List all directories in the user_orthofinder_dir
    all_dirs = [d for d in os.listdir(user_orthofinder_dir) if os.path.isdir(os.path.join(user_orthofinder_dir, d))]

    # Use regex to match directories with the pattern "Results_XXYY"
    results_dirs = []
    for d in all_dirs:
        match = re.match(r'Results_\D*(\d{1,2})$', d)
        if match:
            yy = int(match.group(1))
            results_dirs.append((d, yy))

    # Get the directory with the highest YY value
    if not results_dirs:
        raise ValueError("No valid Results_XXYY directory found.")

    latest_results_dir = max(results_dirs, key=lambda x: x[1])[0]
    return os.path.join(user_orthofinder_dir, latest_results_dir)


def filter_fasta_by_length(user_dir, fasta="user_species.fasta", min_length=800):
    """
    Read a FASTA file, filter out sequences shorter than a specified length, and remove non-DNA/RNA
    nucleotides ('N', 'n', and any other characters except A, T, C, G, U). Save the filtered sequences
    to a new file called 'fasta_filtered.fasta'.

    Parameters:
    - user_dir: str, the directory path where the FASTA file is stored.
    - fasta: str, the name of the input FASTA file (default: "user_species.fasta").
    - min_length: int, the minimum length of nucleotide sequences to keep (default: 800).

    Returns:
    - filtered_sequences: dict, a dictionary of filtered sequences with headers as keys.
    """
    valid_nucleotides = set("ATCGUatcgu")  # Set of valid DNA/RNA nucleotides
    filtered_sequences = {}
    current_header = None
    current_sequence = []

    # Define the output file path
    output_file_path = user_dir + "/user_orthofinder/fasta_filtered.fasta"

    with open(user_dir + "/user_orthofinder" + "/" + fasta, 'r') as fasta_file:
        for line in fasta_file:
            line = line.strip()
            if line.startswith('>'):  # Header line
                if current_header:  # Save the previous sequence
                    full_sequence = ''.join(current_sequence)
                    if len(full_sequence) >= min_length:
                        filtered_sequences[current_header] = full_sequence
                # Start a new sequence
                current_header = line
                current_sequence = []
            else:  # Sequence line
                # Filter out invalid characters (e.g., N, n, or any non-DNA/RNA characters)
                cleaned_sequence = ''.join([char for char in line if char in valid_nucleotides])
                current_sequence.append(cleaned_sequence)

        # Don't forget the last sequence in the file
        if current_header:
            full_sequence = ''.join(current_sequence)
            if len(full_sequence) >= min_length:
                filtered_sequences[current_header] = full_sequence

    # Write filtered sequences to output file
    with open(output_file_path, 'w') as output_file:
        for header, sequence in filtered_sequences.items():
            output_file.write(f"{header}\n")
            # Write sequence in 80-character wide lines (FASTA format)
            for i in range(0, len(sequence), 80):
                output_file.write(f"{sequence[i:i + 80]}\n")

    return None


### cd-hit transcriptome - 0.9 similarity elimination ####



def cd_hit_transcriptome(user_dir,similarity_threshold):



    cd_hit_command = [
        'cd-hit',
        '-i', (user_dir+"/fasta_filtered.fasta.transdecoder.pep"),
        '-o', (user_dir+"/fasta_filtered_transdecoder_cdhit.pep"),
        '-c', '0.99',# 'str(similarity_threshold),
        '-M', '3000',
        '-T', '2'
    ]


    subprocess.run(cd_hit_command, text=True)







####  generate protein.fasta - longest isoform ####


def transdecoder_all(fasta_file_path, dir_path,transdecoder_path = "/home/dcedden/PycharmProjects/dsRIP_web/main_site/ext_packages/TransDecoder-TransDecoder-v5.7.1/"):

    # Prepare the commands
    transdecoder_longest_command = [
        transdecoder_path + 'TransDecoder.LongOrfs',
        '-t', fasta_file_path,
        '-m', '200',
        '--complete_orfs_only',
        '--output_dir', dir_path
    ]

    transdecoder_predict_command = [
        transdecoder_path + 'TransDecoder.Predict',
        '-t', fasta_file_path,
        '--single_best_only',
        '--no_refine_starts',
        '-T', '100',
        '--output_dir', dir_path
    ]

    # Run the commands
    subprocess.run(transdecoder_longest_command, text=True, check=True)
    subprocess.run(transdecoder_predict_command, text=True, check=True)



def orthofinder_after_run(user_dir, user_species_name):


    orthology_result_dir_find = find_latest_results_dir(user_dir+"/orthofinder_input/OrthoFinder")


    orthology_result = orthology_result_dir_find +"/Orthologues/Orthologues_Tribolium/Tribolium__v__" + user_species_name + ".tsv"


    return orthology_result

def orthology_inference_preprocess(user_dir, user_species_name="new_species", Tc_protein_file =  (absolute_path_main_site + "constant_files/Tc_proteins/Tc_protein_dsRIP_91_C1.fasta"), fasta="user_species.fasta", short_threshold= 800, similar_threshold= 0.95):

    try:
        Tc_protein = Tc_protein_file

        os.makedirs(user_dir+ "/user_orthofinder", exist_ok=True)

        os.makedirs(user_dir + "/user_orthofinder/orthofinder_input", exist_ok=True)



        shutil.copy(Tc_protein,  user_dir + "/user_orthofinder/orthofinder_input/Tribolium.fasta")


        filter_fasta_by_length(user_dir, min_length=int(short_threshold))


        #cd_hit_transcriptome(user_dir)


        transdecoder_all(user_dir+"/user_orthofinder/"+ "fasta_filtered.fasta",user_dir)

        cd_hit_transcriptome(user_dir,similar_threshold)

        shutil.copy(user_dir+"/fasta_filtered_transdecoder_cdhit.pep", user_dir+ "/user_orthofinder/orthofinder_input/" + user_species_name +".fasta")


        user_orthofinder_dir = (user_dir+ "/user_orthofinder/orthofinder_input")

        orthofinder_command = [
            'orthofinder',
            '-f', user_orthofinder_dir,
            '-t',"3"
        ]

        subprocess.run(orthofinder_command, text=True)



        return None
    except:
        print("Orthology inference failed")
        return None


def filter_species_genes_new_orthology(species_file, species_name, excel_file=(constant_files + "lethal_gene_threshold.xlsx"), list_criteria=0, go_criteria=None, kegg_criteria=None, without_paralog=False):
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

def filter_fasta_by_ids_new_orthology(user_dir_path,orthology_tsv_path, user_species_name):


    target_genes_dir = os.path.join(user_dir_path, "target_genes")
    os.makedirs(target_genes_dir, exist_ok=True)
    # Clear the output directory before proceeding
    for filename in os.listdir(target_genes_dir):
        file_path = os.path.join(target_genes_dir, filename)  # Get the full path of the file
        if os.path.isfile(file_path) or os.path.islink(file_path):  # Check if it's a file or a symlink
            os.remove(file_path)  # Remove the file

    # Get filtered Tribolium IDs
    filtered_ids = filter_species_genes_new_orthology(orthology_tsv_path,user_species_name)

    if filtered_ids is None or filtered_ids.empty:
        print("No matching IDs found.")
        return

    # Extract the relevant IDs from the Pest Gene ID column
    ids_to_filter = set(filtered_ids['Pest Gene ID'])


    fna_file = user_dir_path + "/user_orthofinder/fasta_filtered.fasta"

    if not fna_file:
        print("No FNA file found for the given species.")
        return

    # Prepare output file path
    output_fasta = os.path.join(target_genes_dir, f"{user_species_name}_target_genes.fasta")

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


