from Bio import SeqIO
from Bio.Seq import Seq
from orffinder import orffinder
from Bio.SeqRecord import SeqRecord
import pandas as pd


def get_dict_with_highest_length(list_of_dicts, gene):
    if not list_of_dicts:
        return None

    dict_with_highest_length = max(list_of_dicts, key=lambda d: d.get('length', 0))

    if dict_with_highest_length["sense"] == "+":
        dict_with_highest_length["gene_sense_direction"] = gene.seq  # gene.seq is the correct access point
        dict_with_highest_length["id"] = gene.id
        return dict_with_highest_length
    elif dict_with_highest_length["sense"] == "-":
        gene_seq_sense = gene.seq.reverse_complement()  # no need to use Seq() around it
        dict_with_highest_length["gene_sense_direction"] = gene_seq_sense
        dict_with_highest_length["id"] = gene.id
        return dict_with_highest_length

def process_ORF_output(gene):
    # Pass the gene to orffinder and get the ORF details
    abc = orffinder.getORFNucleotides(gene, remove_nested=True, return_loci=True)
    return get_dict_with_highest_length(abc, gene)

def process_sequence(sequence, gene_name):
    gene = SeqRecord(Seq(sequence.upper()), id=gene_name)
    processed_orf = process_ORF_output(gene)

    if processed_orf is None:
        return sequence  # Return the original sequence if no ORF is found

    return (processed_orf["gene_sense_direction"], processed_orf["start"], processed_orf["end"],
            processed_orf["sense"], processed_orf["nucleotide"])

def get_ORF(sequence, gene_name):
    gene = SeqRecord(Seq(sequence.upper()), id=gene_name)
    processed_orf = process_ORF_output(gene)

    if processed_orf is None:
        return sequence  # Return the original sequence if no ORF is found

    if int(processed_orf["start"]) < int(processed_orf["end"]):
        return (processed_orf["gene_sense_direction"], processed_orf["start"], processed_orf["end"], processed_orf["sense"], processed_orf["nucleotide"])
    else:
        return process_sequence(str(processed_orf["gene_sense_direction"]), gene_name)


def correct_sequence_based_on_ORF(sequence, gene_name):
    """
    This function corrects the input sequence based on the direction of the longest ORF.
    If no valid ORF is found, or something fails, it returns the original sequence.

    Args:
    sequence (str): The input sequence to be corrected.
    gene_name (str): The name of the gene.

    Returns:
    str: The corrected sequence or the original sequence if no ORF is found or an error occurs.
    """
    try:
        # Process the ORF and get the best match
        gene = SeqRecord(Seq(sequence.upper()), id=gene_name)
        processed_orf = process_ORF_output(gene)

        # If no ORF found, return original sequence
        if processed_orf is None:
            return sequence

        # If ORF has valid start and end, check the direction and return the corrected sequence
        if int(processed_orf["start"]) < int(processed_orf["end"]):
            return str(processed_orf["gene_sense_direction"])
        else:
            # Reverse the sequence if needed and return
            corrected_sequence = process_sequence(str(processed_orf["gene_sense_direction"]), gene_name)
            return str(corrected_sequence[0])  # Return only the corrected gene sequence
    except Exception as e:
        print(f"An error occurred: {e}")
        return sequence  # Return the original sequence if anything goes wrong

# Example usage
#print(correct_sequence_based_on_ORF("CTGTGCTTGTAATGTGCTTGACATCTCCAGTTTTGACTTAAAATTTTTCGAAAGTCATGAAATGAATTGTAATGTTTCATGTTTTAATTGTTTCAATTCTAACGAAATCGTGGTTTGTGGTTAAATATTGTAATTTCTAAGTCGCGGTAAAACATTAAGATATAACAATGGAGCTTACAGCAGAAGGTAAAACTCCAGAGTACATGGCATTAGCCGGAATAAAGTTTAAGTTGACAATCTCGGAATTGAAAAACGACCCGGTTCTAAGGGATCAGCTACTTCAAGGTATTAAAGCTGGTAACATGGCTCCATATTACAAAGAAGTTTGTACGGACTTAGGGTGGACCTTTGACCAAAAGTTGTTTGATGAAATGGCTAAAGAAAACCAAGACAGACTTTCAAAGTTCCAGGAAGACGACTCTGAAACTCCTGTATGGCAAGATAGATTAGATTATCTGTGTTCCATTGGTGACAAAGATGCTGCAACTGCGTTAGCCCAATCTAAATATGAAGACTCCACATTGACAACCAACAGGAGACTTGATGCTATATTTGCTCTATTCAGAATTTCATACTTCCATGGATGTAATGTCAAAGATATGGGCAAATACATCAATAAGGCCCATGAGTTAGTAGACAAAGGCGGTGATTGGAGGTCCCGCAACAAACTTAAGGCATATGAGGCTATCTACTGTTTAGCAGTGAGGGATTACAGTCGTTCAGCAGATCTTTTTATTGACTGTGTTTCAACATTTGAGTCATATGAATTGGTTGATTTTGGGACCATAATTCAATACTGCGTGTTGGCTTGTGCGTTGGCTCTCGAGCGCCACGCGTTGCAAGCGGCTTTGCGTCGCCAAGGCGCGGCGGTGCAAGCTCTACGCTCGCGATTCCCAGAACTAAGGGAATTGGTTGAGTCACTTCATGAGTGTCGGTACGCTGATTTTATGAAGAGTCTTGCTTGGGTGGAAACTCAGATTTGCGTGGACCCCGTGTTCCGTCCGCACTACCAGCACTACGTGCGCGAAGCGCGCATTAAGGCGTACGTTCAACTGTTACGCGCGTATCGCTCGCTCAGTCTTGACAACATTGCCGACACCTTTGGAGTAACCAGGGAATTCATTGAAGATGAGATTTCTACTTTTACAGCGGCCGGCCGTCTCCAATGTCGTATAGATGCAGTGGCGGGTTGTGTGGTCACAGGCGCTGGTCGCGGCGCTGATGCCGATCGTTCGCATCTCTACCAAGCCACCATTCGCGAGGGTGACTTACTGCTTAATAGAGTCAAGAAACTCGCTAGTGTCATTAACTTCTAGTTTTATATGGATGTCTAATTATCATACTTGAATACTTCTTACATAAATGAAATAGAGTCTTTAATAAAACTTTTTTTTATTTTGTAGAATGGTTTGTCTATTGTATGAACTAATACTAACTCTGTGGTCTTCACAGATCACAGTGTATATTTAATGTAATGGTTAAGATGATTTCCTTGCGAATACAACAAATACATTACTTTTTTGTACTAAAG","dsa"))
