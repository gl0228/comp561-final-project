from Bio.Blast.Applications import NcbiblastnCommandline
from Bio import SeqIO
import csv
import math

# Define BLAST command with tabular output
blastn_cline = NcbiblastnCommandline(
    query="target_sequence.fasta",
    db="my_database",
    evalue=10,
    word_size=4,
    perc_identity=10,  # 70% match
    outfmt=6,  # Tabular format
    out="blast_results.tsv",
    strand="plus"
)

def input_length(fasta_file):
    for record in SeqIO.parse(fasta_file, "fasta"):
        return len(str(record.seq))
    
input_sequence_length = input_length("target_sequence.fasta")

# Execute BLAST
stdout, stderr = blastn_cline()

# Define header for tabular format
headers = [
    "qseqid", "sseqid", "pident", "length",
    "mismatch", "gapopen", "qstart", "qend",
    "sstart", "send", "evalue", "bitscore"
]

def fasta_to_array(fasta_file):
    # Read the FASTA file and extract the sequence
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequence = str(record.seq)
        # Convert the sequence into a list (array) of nucleotides
        sequence_array = list(sequence)
        return sequence_array  # Return the array of nucleotides
    
    # If the FASTA file is empty, return an empty list
    return []

def probability_file_to_array(file_path):
    with open(file_path, 'r') as file:
        content = file.read()
        # Split the content by any whitespace character
        tokens = content.split()
        float_tokens = [float(token) for token in tokens]
        return float_tokens
    
def get_positional_probabilities(database_array, probability_array): 
    # Initialize the outer dictionary
    index_dict = {} 
    nucleotides = ['A', 'T', 'G', 'C']
    for i in range(len(database_array)):
        nucleotide_probabilities = {}
        other_probability = round((1 - probability_array[i]) / 3, 3)
        for nucleotide in nucleotides:
            if nucleotide == database_array[i]:
                nucleotide_probabilities[database_array[i]] = probability_array[i]
            else:
                nucleotide_probabilities[nucleotide] = other_probability
        index_dict[i] = nucleotide_probabilities
    return index_dict

database_array = fasta_to_array("database.fasta")
probability_array = probability_file_to_array("database_probs.txt")
positional_probabilities = get_positional_probabilities(database_array, probability_array)

print(positional_probabilities[0])

# Read the tabular BLAST output and write to CSV
with open("blast_results.tsv") as tsv_file, open("blast_results.csv", mode='w', newline='') as csv_file:
    tsv_reader = csv.reader(tsv_file, delimiter='\t')
    csv_writer = csv.writer(csv_file)

    # Write headers
    csv_writer.writerow(headers)
    epsilon = 1e-10

    max_prob = math.log(epsilon)
    print("epsilon prob is", max_prob)
    max_index = -1
    
    match_probability = {}
    max_match = (0,0)
    # Write all rows from TSV to CSV
    for row in tsv_reader:
        csv_writer.writerow(row)
        probability = math.log(1)
        # Extract qstart and qend from the row
        query_start = int(row[6]) - 1
        print("query start is ", query_start)
        query_end = int(row[7]) - 1
        print("query end is ", query_end)
        database_start = int(row[8]) - 1 # Index 6 corresponds to 'qstart'
        print("database start is ", database_start)
        database_end = int(row[9]) - 1   # Index 7 corresponds to 'qend'
        print("database end is ", database_end)
        print("input_sequence_length is", input_sequence_length)
        start_end_tuple = (database_start + 1, database_end + 1)
        left_extension = query_start
        right_extension = input_sequence_length - query_end + 1
        print("I am before for loop, the range is ", start_end_tuple)

        for index in range(database_start, database_end):
            print("I am in for loop")
            nucleotide = database_array[index]
            nucleotide_prob = positional_probabilities[index][nucleotide]
            print("The index is ", index)
            print("The nucleotide is ", nucleotide)
            print("The nucleotide prob is ", nucleotide_prob)
            probability += math.log(nucleotide_prob) if nucleotide_prob > 0 else math.log(epsilon)
            # print("The positional probabilities are ", positional_probabilities[index])
        
        print("I am after for loop")
        print("Match robability for", start_end_tuple,  "is ", probability)
        if probability > max_prob:
            print("I am here")
            max_prob = probability
            max_match = start_end_tuple
        
        match_probability[start_end_tuple] = probability

        # Print the indices
        # print(f"Match in target sequence from nucleotide {input_start} to {input_end}")

    print("The max prob is", math.exp(max_prob))
    print("The max match is ", max_match)
    print(match_probability)



