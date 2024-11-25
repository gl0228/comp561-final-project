from Bio.Blast.Applications import NcbiblastnCommandline
from Bio import SeqIO
import csv
import math
import random

def fasta_to_array(fasta_file):
    # Read the FASTA file and extract the sequence
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequence = str(record.seq)
        # Convert the sequence into a list (array) of nucleotides
        sequence_array = list(sequence)
        return sequence_array  # Return the array of nucleotides
    
    # If the FASTA file is empty, return an empty list
    return []

database_array = fasta_to_array("database.fasta")
query_sequence = fasta_to_array("query_sequence.fasta")
database_length = len(database_array)
query_length = len(query_sequence)

HEURISTIC = True
# The number of consecutive low probability nucleotides before the algorithm discards the sequence
HEURISTIC_NUM = 10
# The probability that is considered low for the heuristic
HEURISTIC_PROB = 0.1

# print("The segment randomly picked is", segment)
# print("The start index is ", start)
# print("The end index is ", end)


# Define BLAST command with tabular output
blastn_cline = NcbiblastnCommandline(
    query="query_sequence.fasta",
    db="database/my_database",
    evalue=1000,
    word_size=4,
    perc_identity=0,
    outfmt=6,  # Tabular format
    out="blast_results.tsv",
    ungapped=True,
    strand="plus"
)

# Execute BLAST
stdout, stderr = blastn_cline()

# Define header for tabular format
headers = [
    "qseqid", "sseqid", "pident", "length",
    "mismatch", "gap_open", "query_start", "query_end",
    "database_start", "database_end", "evalue", "bitscore"
]

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

    max_prob = float('-inf')
    # print("epsilon prob is", max_prob)
    max_index = -1
    
    match_probability = {}
    max_match = (0,0)
    unlikely_match = False
    # Write all rows from TSV to CSV
    for row in tsv_reader:
        csv_writer.writerow(row)
        probability = math.log(1)
        query_start = int(row[6]) - 1
        # print("query_start", query_start)
        query_end = int(row[7]) - 1
        # print("query_end", query_end)
        database_start = int(row[8]) - 1
        # print("database_start", database_start)
        diff = database_start - query_start
        database_end = int(row[9]) - 1
        # print("database_end", database_end)
        # print("database diff", database_end - database_start)
       #  print("query diff", query_end - query_start)
        left_extension = query_start if query_start <= database_start else database_start
        right_extension = query_length - query_end if query_length - query_end <= database_length - database_end else database_length - database_end

        extended_database_start = database_start - left_extension 
        extended_database_end = database_end + right_extension - 1

        # print("extended database diff", extended_database_end - extended_database_start)
        extended_query_start = query_start - left_extension
        extended_query_end = query_end + right_extension -1 
        # print("extended query diff", extended_query_end - extended_query_start)
        # print("query length", len(query_sequence))
        # print("extended_query_start", extended_query_start)
        # print("extended_query_end", extended_query_end)
        start_end_tuple = (extended_database_start, extended_database_end )

        consecutive_low_prob_count = 0 
        for index in range(extended_database_start, extended_database_end):
            # print("The length is ", extended_database_end - extended_database_start)
            # print(len(query_sequence))
            # print("the diff is", diff)
            # print("the index is ", index)
            # print("index is", index)
            # print("extended_query_start is", extended_query_start)
            nucleotide = query_sequence[extended_query_start]
            # print("The nucleotide is ", nucleotide)
            nucleotide_prob = positional_probabilities[index][nucleotide]
            if HEURISTIC: 
                if nucleotide_prob < HEURISTIC_PROB:
                    consecutive_low_prob_count += 1
                else:
                    consecutive_low_prob_count = 0
                if consecutive_low_prob_count >= HEURISTIC_NUM :
                    # print("Unlikely match, skip!!!")
                    unlikely_match = True
                    consecutive_low_prob_count = 0
                    continue
            # print("The nucleotide prob is ", nucleotide_prob)
            probability += math.log(nucleotide_prob) if nucleotide_prob > 0 else math.log(epsilon)
            extended_query_start+=1
            # print("The positional probabilities are ", positional_probabilities[index])
        # print("The sequence is ", start_end_tuple)
        # print("The match probability is ", probability)
        if unlikely_match:
            unlikely_match = False
            continue

        if probability > max_prob:
            max_prob = probability
            max_match = start_end_tuple
        
        match_probability[start_end_tuple] = probability

        # Print the indices
        # print(f"Match in target sequence from nucleotide {input_start} to {input_end}")

def write_probabilities_to_csv(log_probabilities, output_file):
    """
    Writes a dictionary of sequence log probabilities to a CSV file.

    Parameters:
        log_probabilities (dict): A dictionary where keys are tuples (start, end) and values are log probabilities.
        output_file (str): The name of the output CSV file.
    """
    counter = 1
    # Define the headers
    headers = ["probability_ranking","qseqid", "sseqid", "database_start", "database_end", "log_probability"]

    # Open the file for writing
    with open(output_file, mode='w', newline='') as csv_file:
        csv_writer = csv.writer(csv_file)

        # Write the headers
        csv_writer.writerow(headers)

        # Write each dictionary entry to the CSV
        for (start, end), probability in log_probabilities.items():
            csv_writer.writerow([counter,"random_query_sequence", "database", start, end, probability])
            counter+=1

    print("The max prob is", math.exp(max_prob))
    print("The max match is ", max_match)
    # print(match_probability)

write_probabilities_to_csv(match_probability, "final_result.csv")





