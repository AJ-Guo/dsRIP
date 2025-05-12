import primer3
from Bio import SeqIO
import pandas as pd

# Constants for the T7 promoter sequence
T7_PROMOTER = 'TAATACGACTCACTATAGGGAGA'
SP6_PROMOTER = 'ATTTAGGTGACACTATAGAAGAG'
T3_PROMOTER = 'AATTAACCCTCACTAAAGGGAGA'

EcoRI = "GAATTC"
NotI = "GCGGCCGC"

EcoRI_UPS = "TAAGCA"
NotI_UPS = "TGCTTA"



def generate_well_positions():
    """Generate well positions for a 96-well plate."""
    rows = 'ABCDEFGH'
    cols = range(1, 13)
    for row in rows:
        for col in cols:
            yield f'{row}{col}'

def primer_design(sequence, min_length, max_length, opt_length_pri, min_length_pri, max_length_pri, min_tm, max_tm, target_tm, min_gc, max_gc):
    product_size_range = [min_length, max_length]
    sequence_args = {
        'SEQUENCE_ID': 'example',
        'SEQUENCE_TEMPLATE': sequence,
        'SEQUENCE_INCLUDED_REGION': [0, len(sequence)]
    }
    global_args = {
        'PRIMER_OPT_SIZE': opt_length_pri,
        'PRIMER_MIN_SIZE': min_length_pri,
        'PRIMER_MAX_SIZE': max_length_pri,
        'PRIMER_OPT_TM': target_tm,
        'PRIMER_MIN_TM': min_tm,
        'PRIMER_MAX_TM': max_tm,
        'PRIMER_MIN_GC': min_gc,
        'PRIMER_MAX_GC': max_gc,
        'PRIMER_MAX_POLY_X': 100,
        'PRIMER_SALT_MONOVALENT': 50.0,
        'PRIMER_DNA_CONC': 50.0,
        'PRIMER_MAX_NS_ACCEPTED': 0,
        'PRIMER_MAX_SELF_ANY': 7,
        'PRIMER_MAX_SELF_END': 7,
        'PRIMER_PAIR_MAX_COMPL_ANY': 7,
        'PRIMER_PAIR_MAX_COMPL_END': 7,
        'PRIMER_PRODUCT_SIZE_RANGE': product_size_range
    }
    data = primer3.design_primers(sequence_args, global_args)
    return select_best_primer(data, target_tm, sequence)


def select_best_primer(data, target_tm, sequence):
    best_primer = None
    best_penalty = float('inf')
    best_tm_diff = float('inf')
    input_sequence_length = len(sequence)  # Length of the input sequence
    for i in range(data['PRIMER_PAIR_NUM_RETURNED']):
        left_primer_tm = data[f'PRIMER_LEFT_{i}_TM']
        right_primer_tm = data[f'PRIMER_RIGHT_{i}_TM']
        pair_penalty = data[f'PRIMER_PAIR_{i}_PENALTY']
        avg_tm = (left_primer_tm + right_primer_tm) / 2
        tm_diff = abs(avg_tm - target_tm)
        if pair_penalty < best_penalty and tm_diff < best_tm_diff:
            best_penalty = pair_penalty
            best_tm_diff = tm_diff
            best_primer = i
    if best_primer is None:
        return None

    start = data['PRIMER_LEFT_'+str(best_primer)][0]
    end = data['PRIMER_RIGHT_'+str(best_primer)][0]
    amplicon_length = end - start + 1
    amplicon_sequence = sequence[start:end+1]

    # Calculate non-overlapping sequences and their lengths
    non_overlapping_1 = sequence[:start]
    non_overlapping_2 = sequence[end+1:]

    percentage_length = (amplicon_length / input_sequence_length) * 100
    return {
        'left_primer': data[f'PRIMER_LEFT_{best_primer}_SEQUENCE'],
        'right_primer': data[f'PRIMER_RIGHT_{best_primer}_SEQUENCE'],
        'avg_TM' :avg_tm,
        'amplicon_length': amplicon_length,
        'amplicon_sequence': amplicon_sequence,
        'input_sequence_length': input_sequence_length,
        'percentage_length': percentage_length,
        'amplicon_start': start,
        'amplicon_end': end,
        'non_overlapping_1': non_overlapping_1,
        'non_overlapping_2': non_overlapping_2,
    }


def design_primers_from_fasta(fasta_file, output_file, primer_params):
    well_position_generator = generate_well_positions()
    results = []

    for record in SeqIO.parse(fasta_file, "fasta"):
        gene_name = record.id
        sequence = str(record.seq)
        input_sequence_length = len(sequence)

        # Check if the sequence is valid (basic validation for nucleotide sequences)
        if not set(sequence.upper()).issubset({"A", "T", "G", "C", "N"}):
            print(f"Skipping {gene_name}: Sequence contains invalid characters.")
            continue

        # Calculate dynamic min and max amplicon size based on fraction
        min_length = int(int(primer_params['dsRNA_length_min_perc']) * input_sequence_length) / 100
        max_length = input_sequence_length

        # Use the user-defined primer parameters
        result = primer_design(
            sequence,
            min_length,
            max_length,
            primer_params['primer_size_opt'],
            primer_params['primer_size_min'],
            primer_params['primer_size_max'],
            primer_params['annealing_temp_min'],
            primer_params['annealing_temp_max'],
            primer_params['annealing_temp_opt'],
            primer_params['primer_GC_min'],
            primer_params['primer_GC_max']
        )

        if result:
            primer_types = []
            sequences = []

            # Check if two primer pairs are required
            if primer_params['two_primer_pairs']:
                # Add primers without overhangs
                primer_types.extend(['F', 'R'])
                sequences.extend([result['left_primer'], result['right_primer']])

                # Add primers with overhangs
                primer_types.extend(['F_ovr', 'R_ovr'])
                sequences.extend([
                    str(primer_params['forward_overhang']) + result['left_primer'],
                    str(primer_params['reverse_overhang']) + result['right_primer']
                ])
            else:
                # Only add primers with overhangs
                primer_types.extend(['F_ovr', 'R_ovr'])
                sequences.extend([
                    str(primer_params['forward_overhang']) + result['left_primer'],
                    str(primer_params['reverse_overhang']) + result['right_primer']
                ])

            for primer_type, seq in zip(primer_types, sequences):
                well_position = next(well_position_generator)
                results.append({
                    'Well Position': well_position,
                    'Name': f'{gene_name}_{primer_type}',
                    'Sequence': seq,
                    'Average Primer TM (gene complementary region)': result["avg_TM"],
                    'Amplicon Length': result['amplicon_length'],
                    'Amplicon Sequence': result['amplicon_sequence'] if primer_type in ['F', 'R'] else '',
                    'Coverage of optimized dsRNA region (%)': f'{result["percentage_length"]:.2f}%'
                })

    # Create a DataFrame and remove duplicate rows based on the 'Sequence' column
    df = pd.DataFrame(results).drop_duplicates(subset=['Sequence'])

    # Write the output to Excel
    df.to_excel(output_file, index=False)




# Usage example
#fasta_file = "/Users/dogacedden/Desktop/PhD_drive/dsRIP test files/dsRIP_outputs/best_dsRNA_regions.fasta"  # Update with the actual path to your FASTA file
#output_file = "/Users/dogacedden/Desktop/PhD_drive/dsRIP test files/dsRIP_outputs/Pbra_best_dsRNA_regions_primers.xlsx"  # Update with the desired path for the output Excel file

#design_primers_from_fasta(fasta_file, output_file)
