import os
import subprocess
import tempfile
import pandas as pd
import pysam
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import csv
from brokenaxes import brokenaxes
from random import randrange



def generate_siRNAs(sequence, si_length):
    siRNA_sequences = []
    for i in range(0, len(sequence) - si_length + 1):
        kmer = sequence[i:i+si_length]
        siRNA_sequences.append(kmer)
    return siRNA_sequences


def run_bowtie1(siRNA_names, siRNA_sequences, index_prefix, siRNA_length, off_target_dir, mis_input):
    with tempfile.NamedTemporaryFile(mode='w+', delete=False, suffix='.fna') as tmp:
        for siRNA_name, siRNA_sequence in zip(siRNA_names, siRNA_sequences):
            tmp.write(f">{siRNA_name}\n{siRNA_sequence}\n")
        tmp_fasta_name = tmp.name

    bowtie1_command = [
        "/home/dcedden/PycharmProjects/dsRIP_web/main_site/ext_packages/bowtie-1.3.1-linux-x86_64/bowtie-align-s",
        "-n", str(mis_input),
        "-a",
        "-x", index_prefix,
        "-f", tmp_fasta_name,
        "-S", off_target_dir + "/siRNAs.sam"
    ]

    try:
        subprocess.run(bowtie1_command, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error occurred: {str(e)}")
    finally:
        os.unlink(tmp_fasta_name)

    return off_target_dir + "/siRNAs.sam"


def read_lethal_genes(species_list):
    all_lethals = {}
    for species in species_list:
        file_path = os.path.join("all_lethals", species + "_all_lethals")
        #print(f"Processing: {file_path}")
        try:
            with open(file_path, 'r') as file:
                reader = csv.reader(file, delimiter='\t')
                all_lethals[species] = [row[2] for row in reader if
                                        len(row) > 2 and row[2]]  # Access the third column (index 2)
                #print(f"Lethal genes for {species}: {all_lethals[species]}")

        except FileNotFoundError:
            print(f"No lethal genes file found for {species}.")
        except IndexError:
            print(f"Error: The file for {species} does not contain enough columns.")
        except Exception as e:
            print(f"An error occurred while processing {species}: {e}")

    return all_lethals

def predict_off_target(input_siRNAS, species_list, kmer_length, all_lethals, off_target_dir, mis_input="1", custom_index_prefix=None):
    alignment_data = {}
    siRNA_sequences = input_siRNAS

    filtered_siRNA_sequences = []
    siRNA_names = []
    for idx, kmer in enumerate(siRNA_sequences):
        #polyAT_count = kmer.count('A') + kmer.count('T')
        #if polyAT_count / kmer_length < 2:
        filtered_siRNA_sequences.append(kmer)
        siRNA_names.append(f"siRNA_{idx + 1}")

    siRNA_stats = {name: {"total_maps": 0, "total_lethal_maps": 0} for name in siRNA_names}

    for species in species_list:
        if custom_index_prefix:
            index_prefix = custom_index_prefix  # Use the custom index prefix for user transcriptome
        else:
            index_prefix = os.path.join("constant_files/indexes", species, species)

        sam_file_path = run_bowtie1(siRNA_names, filtered_siRNA_sequences, index_prefix, kmer_length, off_target_dir, mis_input=mis_input)

        try:
            with pysam.AlignmentFile(sam_file_path, "r") as samfile:
                for read in samfile.fetch():
                    if not read.is_unmapped:
                        siRNA_name = read.query_name
                        mismatches = read.get_tag("NM")
                        ref_name = read.reference_name
                        query_sequence = read.query_sequence
                        perfect_matches = len(query_sequence) - mismatches
                        strand = "reverse" if read.is_reverse else "forward"
                        gene_name = ref_name

                        is_lethal = gene_name in all_lethals.get(species, [])
                        if is_lethal:
                            gene_name += "_LETHAL"
                            siRNA_stats[siRNA_name]["total_lethal_maps"] += 1

                        siRNA_stats[siRNA_name]["total_maps"] += 1

                        alignment_result = {
                            "query_sequence": query_sequence,
                            "fasta_file": os.path.basename(index_prefix).replace("_index", ".fna"),
                            "gene_name": gene_name,
                            "perfect_matches": perfect_matches,
                            "mismatches": mismatches,
                            "strand": strand
                        }

                        if species not in alignment_data:
                            alignment_data[species] = {}
                        if siRNA_name not in alignment_data[species]:
                            alignment_data[species][siRNA_name] = []

                        alignment_data[species][siRNA_name].append(alignment_result)
        finally:
            #os.remove(sam_file_path)
            pass

    return alignment_data, [(name, stats["total_maps"], stats["total_lethal_maps"]) for name, stats in siRNA_stats.items()]



def main(input_siRNAs, input_length, name, species_input, off_target_dir, user_dir, siRNA_length, mismatch_input, new_off_target_species_input_name):
    seq_name = name
    species_list = species_input.split(',') if species_input else []
    kmer_length = siRNA_length
    mismatch = str(mismatch_input)
    lethal_plot_data = {}
    all_plot_data = {}
    all_siRNAs_report = []
    all_lethals = {}

    # Process predefined species
    alignment_data = {}
    if species_list:
        try:
            all_lethals = read_lethal_genes(species_list)
            alignment_data, all_siRNAs_report = predict_off_target(
                input_siRNAs,
                species_list,
                kmer_length,
                all_lethals,
                off_target_dir,
                mis_input=mismatch
            )
        except Exception as e:
            print(f"Error processing predefined species: {e}")

    # Process additional species from user_transcriptome
    transcriptome_dir = os.path.join(user_dir, ("efficiency/" + new_off_target_species_input_name))
    additional_all_siRNAs_report = []

    if os.path.exists(transcriptome_dir):
        for species_file in os.listdir(transcriptome_dir):
            if species_file.endswith(".ebwt"):  # Only process Bowtie index files with .ebwt extensions
                species_basename = species_file.split('.')[0]  # Get the prefix without extension
                species_index_prefix = os.path.join(transcriptome_dir, species_basename)

                print(f"Processing species from user_transcriptome: {species_basename}")

                # Skip lethal gene reading for this species, only run Bowtie alignment
                additional_alignment_data, additional_all_siRNAs_report = predict_off_target(
                    input_siRNAs,
                    [species_basename],
                    kmer_length,
                    all_lethals,
                    off_target_dir,
                    mis_input=mismatch,
                    custom_index_prefix=species_index_prefix
                )

                # Merge additional alignment data
                alignment_data.update(additional_alignment_data)
                break  # Exit the loop after processing the first species

    # Combine all_siRNAs_report and additional_all_siRNAs_report
    # Combine all_siRNAs_report and additional_all_siRNAs_report intelligently
    if additional_all_siRNAs_report:
        if all_siRNAs_report:
            # Convert all_siRNAs_report and additional_all_siRNAs_report to dictionaries for easy merging
            combined_report_dict = {name: [total_maps, total_lethal_maps] for name, total_maps, total_lethal_maps in
                                    all_siRNAs_report}

            for name, total_maps, total_lethal_maps in additional_all_siRNAs_report:
                if name in combined_report_dict:
                    # Sum values for the same siRNA
                    combined_report_dict[name][0] += total_maps
                    combined_report_dict[name][1] += total_lethal_maps
                else:
                    # Add new siRNA entries from the additional report
                    combined_report_dict[name] = [total_maps, total_lethal_maps]

            # Convert the dictionary back to a list
            all_siRNAs_report = [(name, values[0], values[1]) for name, values in combined_report_dict.items()]
        else:
            # If no predefined species report exists, use the additional report directly
            all_siRNAs_report = additional_all_siRNAs_report



    # Process and plot the alignment data
    for species, siRNA_results in alignment_data.items():
        print(species)
        if species not in lethal_plot_data:
            lethal_plot_data[species] = {}
        if species not in all_plot_data:
            all_plot_data[species] = {}

        species_df = pd.DataFrame()

        for siRNA_name, data in siRNA_results.items():
            if data:
                df = pd.DataFrame(data)
                df['siRNA_name'] = siRNA_name
                species_df = pd.concat([species_df, df])

                for entry in data:
                    position = int(siRNA_name.replace('siRNA_', ''))
                    all_plot_data[species].setdefault(position, 0)
                    all_plot_data[species][position] += 1

                    if '_LETHAL' in entry['gene_name']:
                        lethal_plot_data[species].setdefault(position, 0)
                        lethal_plot_data[species][position] += 1

        if not species_df.empty:
            species_df.to_csv(off_target_dir + f"/{seq_name}_{species}.txt", sep='\t', index=False)

    # Create a color map for species
    all_species = set(list(alignment_data.keys()) + list(lethal_plot_data.keys()))
    color_map = {species: cm.tab10(i % 10) for i, species in enumerate(all_species)}

    plt.figure(figsize=(20, 5))

    # First subplot
    plt.subplot(1, 2, 1)
    for species, siRNA_results in alignment_data.items():
        plot_data = {}
        for siRNA_name, data in siRNA_results.items():
            for entry in data:
                position = int(siRNA_name.replace('siRNA_', ''))
                plot_data.setdefault(position, 0)
                plot_data[position] += 1
        if plot_data:
            plt.scatter(plot_data.keys(), plot_data.values(), label=species.replace("_", " "), s=3,
                        color=color_map[species])

    plt.xlabel('siRNA Position')
    plt.ylabel('Number of genes')
    plt.title(f'All off-targets for {seq_name}')
    plt.legend(loc='lower left', bbox_to_anchor=(1, 1))

    plt.tight_layout()

    # Second subplot
    plt.subplot(1, 2, 2)
    for species in lethal_plot_data:
        species_data = lethal_plot_data[species]
        if species_data:
            positions = list(species_data.keys())
            counts = list(species_data.values())
            plt.scatter(positions, counts, label=species, s=3, color=color_map[species])

    plt.xlabel('siRNA Position')
    plt.ylabel('Number of genes')
    plt.title(f'Essential gene off-targets for {seq_name}')
    plt.tight_layout()
    plt.subplots_adjust(wspace=0.1)

    return lethal_plot_data, int(input_length), all_plot_data, plt, all_siRNAs_report




def clean_off_target_per_species(lethal_targets, all_targets):
    summary_data = []

    all_species = set(lethal_targets.keys()) | set(all_targets.keys())

    for species in all_species:
        lethal_positions_counts = lethal_targets.get(species, {})
        all_positions_counts = all_targets.get(species, {})

        lethal_off_targets_sum = sum(lethal_positions_counts.values())
        all_off_targets_sum = sum(all_positions_counts.values())

        summary_data.append({
            "species": species,
            "all_off_targets": all_off_targets_sum,
            "lethal_off_targets": lethal_off_targets_sum
        })
    return summary_data


def aggregate_lethal_counts_from_dict(lethal_plot_data):
    aggregated_counts = {}
    for species, positions in lethal_plot_data.items():
        for position, count in positions.items():
            aggregated_counts[position] = aggregated_counts.get(position, 0) + count
    #print(aggregated_counts)
    return aggregated_counts




def sequence_within_range(range_tuple, sequence):
    """
    Extracts a subsequence defined by the given range from the sequence.

    Parameters:
    - range_tuple: A tuple (start, end) defining the range.
    - sequence: The sequence from which to extract the subsequence.

    Returns:
    The subsequence of `sequence` defined by `range_tuple`.
    """
    start, end = range_tuple
    return sequence[start:end]


# Usage Example
def off_target_siRNA_all(siRNAs_DNA_off_target, seq_length, name, off_target_species, user_dir,siRNA_length, user_mismatch,new_off_target_species_input_name):
    off_target_id = "id_" + name
    off_target_dir = os.path.join(user_dir, "efficiency/off_target", off_target_id)
    os.makedirs(off_target_dir, exist_ok=True)

    lethal_plot_data, input_seq_length, all_plot_data, off_targets_plot, all_siRNAs_report = main(
        siRNAs_DNA_off_target, seq_length, name, off_target_species, off_target_dir, user_dir, siRNA_length= siRNA_length, mismatch_input= user_mismatch,new_off_target_species_input_name=new_off_target_species_input_name)

    per_species_summary = clean_off_target_per_species(lethal_plot_data, all_plot_data)

    all_off_targets_sum = sum(item['all_off_targets'] for item in per_species_summary)
    print(all_off_targets_sum)
    lethal_off_targets_sum = sum(item['lethal_off_targets'] for item in per_species_summary)

    return per_species_summary, all_off_targets_sum, lethal_off_targets_sum, all_siRNAs_report, off_targets_plot



