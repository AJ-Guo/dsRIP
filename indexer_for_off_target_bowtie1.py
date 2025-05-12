import os
import subprocess


def filter_fasta_by_length(user_dir, min_length=200):
    """
    Read a FASTA file, filter out sequences shorter than a specified length, and remove non-DNA/RNA
    nucleotides ('N', 'n', and any other characters except A, T, C, G, U). Optionally, save the
    filtered sequences to a new file.

    Parameters:
    - user_dir: str, the directory containing the input FASTA file.
    - fasta: str, name of the input FASTA file.
    - min_length: int, the minimum length of nucleotide sequences to keep.

    Returns:
    - filtered_fasta_path: str, path to the filtered fasta file.
    """
    valid_nucleotides = set("ATCGUatcgu")  # Set of valid DNA/RNA nucleotides
    filtered_sequences = {}
    current_header = ""




    with open(user_dir + "/user_index.fasta", 'r') as fasta_file:
        for line in fasta_file:
            line = line.strip()
            if line.startswith('>'):  # Header line
                current_header = line
            else:  # Sequence line
                # Filter out invalid characters (e.g., N, n, or any non-DNA/RNA characters)
                cleaned_sequence = ''.join([char for char in line if char in valid_nucleotides])

                if current_header in filtered_sequences:
                    filtered_sequences[current_header] += cleaned_sequence
                else:
                    filtered_sequences[current_header] = cleaned_sequence

    # Filter sequences by length
    filtered_sequences = {header: seq for header, seq in filtered_sequences.items() if len(seq) >= min_length}

    # Write filtered sequences to a new FASTA file
    with open(user_dir + "/index_filtered.fna", 'w') as output_file:
        for header, sequence in filtered_sequences.items():
            output_file.write(f"{header}\n")
            output_file.write(f"{sequence}\n")





def index_fna_file(user_dir, user_species_name, bowtie_build_path="/home/dcedden/PycharmProjects/dsRIP_web/main_site/ext_packages/bowtie-1.3.1-linux-x86_64/bowtie-build"):

    fna_file = user_dir + "/index_cdhit.fna"
    # Create a directory for the .fna file if it doesn't exist
    index_dir = os.path.join(os.path.dirname(fna_file), user_species_name)
    if not os.path.exists(index_dir):
        os.makedirs(index_dir)

    # Run bowtie-build to create the index
    try:
        subprocess.run([bowtie_build_path, fna_file, (index_dir + "/" + user_species_name)], check=True)

        return "indexed"
    except subprocess.CalledProcessError:

        return "error"



def cd_hit_transcriptome(user_dir,  c_param):
    """
    Run CD-HIT on a given fasta file.

    Parameters:
    - user_dir (str): Directory path where the fasta file is located.
    - fasta (str): Fasta filename to be processed.
    - output_file (str): Output file name for CD-HIT.
    - c_param (float): The sequence identity threshold for CD-HIT.

    Returns:
    - str: Path to the output fasta file.
    """
    input_fasta = user_dir + "/index_filtered.fna"

    cd_hit_command = [
        'cd-hit-est',
        '-i', input_fasta,
        '-o', user_dir + "/index_cdhit.fna",
        '-c', str(c_param),
        '-n', '8',
        '-M', '4000',
        '-T', '2'
    ]

    try:
        subprocess.run(cd_hit_command, check=True)

    except subprocess.CalledProcessError:

        return None




def run_cd_hit_and_indexer(user_dir_path, user_species_name= "user_transcriptome", c_param=0.99, min_length=800,
                           bowtie_build_path='bowtie-build'):
    """
    Runs filtering, CD-HIT, and then Bowtie-Build indexer on the output file.

    Parameters:
    - user_dir_path (str): Path to the user's directory.
    - fasta (str): The fasta file to be processed. Defaults to 'user_species.fasta'.
    - c_param (float): CD-HIT sequence identity threshold. Defaults to 0.95.
    - min_length (int): Minimum sequence length for filtering. Defaults to 800.
    - bowtie_build_path (str): Path to bowtie-build executable. Defaults to 'bowtie-build'.

    Returns:
    - None
    """
    # Step 1: Filter sequences by length
    filtered_fasta = filter_fasta_by_length(user_dir_path, min_length)



    # Step 3: Run CD-HIT on the filtered sequences
    cd_hit_transcriptome(user_dir_path, c_param)



    # Step 4: Rename CD-HIT output to have .fna extension for bowtie-build compatibility



    index_fna_file(user_dir_path,user_species_name)