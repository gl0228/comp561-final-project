import random
from Bio import SeqIO

def fasta_to_array(fasta_file):
    # Read the FASTA file and extract the sequence
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequence = str(record.seq)
        # Convert the sequence into a list (array) of nucleotides
        sequence_array = list(sequence)
        return sequence_array  # Return the array of nucleotides
    
    # If the FASTA file is empty, return an empty list
    return []

def pick_random_segment(sequence_array, segment_length, fasta_filename="query_sequence.fasta"):
    """
    Randomly pick a segment of the specified length from the given sequence array 
    and write it to a FASTA file. Ensures the FASTA file is cleared before writing.
    
    Parameters:
        sequence_array (list): The long sequence represented as a list of characters (e.g., ['A', 'C', 'T', 'G']).
        segment_length (int): The length of the segment to extract.
        fasta_filename (str): The name of the FASTA file to write the segment to.
    
    Returns:
        tuple: A tuple containing the segment (str), start index (int), and end index (int).
    """
    if segment_length > len(sequence_array):
        raise ValueError("Segment length cannot be longer than the input sequence array.")
    
    # Randomly pick a start index such that the segment fits within the sequence array
    start_index = random.randint(0, len(sequence_array) - segment_length)
    end_index = start_index + segment_length  # Calculate the end index (exclusive)
    
    # Extract the segment
    segment = ''.join(sequence_array[start_index:end_index])  # Convert the array slice back to a string

    with open("query_sequence.fasta", "w") as fasta_file:
        fasta_file.write(">random_query_sequence\n")
        fasta_file.write(segment + "\n")
    
    return segment, start_index, end_index

database_array = fasta_to_array("database.fasta")

segment, start, end = pick_random_segment(database_array, 1000)

print("The segment randomly picked is", segment)
print("The start index is ", start)
# The end index is exclusive
print("The end index is ", end - 1)