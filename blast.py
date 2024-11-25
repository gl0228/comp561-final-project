from Bio.Blast.Applications import NcbiblastnCommandline
from Bio import SeqIO
import csv
import math

# Takes a fasta file and return an array
def fasta_to_array(fasta_file):
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequence = str(record.seq)
        sequence_array = list(sequence)
        return sequence_array  # Return the array of nucleotides
    
    return []

# Get array for database sequence
database_array = fasta_to_array("database.fasta")
# Get array for query sequence
query_sequence = fasta_to_array("query_sequence.fasta")
# Get length for database sequence
database_length = len(database_array)
# Get length for query sequence
query_length = len(query_sequence)

# Variable to turn heuristic on/off
# The heuristic eliminates a match if it exceeds a certain number of low probability nucleotides in a row
HEURISTIC = False
# The number of consecutive low probability nucleotides before the algorithm discards the sequence
HEURISTIC_NUM = 10
# The probability that is considered low for the heuristic
HEURISTIC_PROB = 0.1

# BLAST command
# To get more matches, increase evalue or word_size 
blastn_cline = NcbiblastnCommandline(
    query="query_sequence.fasta",
    db="database/my_database",
    evalue=1000,
    word_size=4,
    # Percentage of mismatches
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
    
# Returns a dictionary with keys as indexes in the database array,
# and values as dictionaries containing the probability of each nucleotide
def get_positional_probabilities(database_array, probability_array): 
    index_dict = {} 
    nucleotides = ['A', 'T', 'G', 'C']

    # For each position in the database, get the corresponding probabilities of all the nucleotides
    for i in range(len(database_array)):
        nucleotide_probabilities = {}
        # The probability of the other 3 nucleotides
        other_probability = round((1 - probability_array[i]) / 3, 3)
        for nucleotide in nucleotides:
            # The nucleotide written in the database file gets the probability written in the probability file
            if nucleotide == database_array[i]:
                nucleotide_probabilities[database_array[i]] = probability_array[i]
            # The unwritten nucleotides get (leftover probability / 3)
            else:
                nucleotide_probabilities[nucleotide] = other_probability
        index_dict[i] = nucleotide_probabilities
    return index_dict

probability_array = probability_file_to_array("database_probs.txt")
# Get the dictionary with the probabilities of each nucleotide at each position
positional_probabilities = get_positional_probabilities(database_array, probability_array)

# Read the tabular BLAST output and write to CSV
with open("blast_results.tsv") as tsv_file, open("blast_results.csv", mode='w', newline='') as csv_file:
    tsv_reader = csv.reader(tsv_file, delimiter='\t')
    csv_writer = csv.writer(csv_file)

    # Write headers
    csv_writer.writerow(headers)
    # Is used instead of 0 for log probabilities
    epsilon = 1e-10

    max_prob = float('-inf')
    max_index = -1
    
    match_probability = {}
    max_match = (0,0)
    # Becomes true when the limit of consecutive low probability nucleotides is reached for the heuristic
    unlikely_match = False

    # Write all rows from TSV to CSV
    for row in tsv_reader:
        csv_writer.writerow(row)
        probability = math.log(1)

        # Start and end indices of the BLAST match in the query sequence
        query_start = int(row[6]) - 1
        query_end = int(row[7]) - 1

        # Start and end indices of the BLAST match in the database sequence
        database_start = int(row[8]) - 1
        database_end = int(row[9]) - 1
        diff = database_start - query_start

        # The number of nucleotides to extend on the left of the seed, to match the length of the query sequence
        left_extension = query_start if query_start <= database_start else database_start
        # The number of nucleotides to extend on the right of the seed
        right_extension = query_length - query_end if query_length - query_end <= database_length - database_end else database_length - database_end

        # The start and end indices of the database after extending, to match the length of the query sequence
        extended_database_start = database_start - left_extension 
        extended_database_end = database_end + right_extension - 1

        # The start and end indices of the query after extending, to match the length of the query sequence
        extended_query_start = query_start - left_extension
        extended_query_end = query_end + right_extension - 1


        start_end_tuple = (extended_database_start, extended_database_end )

        consecutive_low_prob_count = 0 
        for index in range(extended_database_start, extended_database_end):
            nucleotide = query_sequence[extended_query_start]
            nucleotide_prob = positional_probabilities[index][nucleotide]

            if HEURISTIC: 
                if nucleotide_prob < HEURISTIC_PROB:
                    consecutive_low_prob_count += 1
                else:
                    consecutive_low_prob_count = 0
                if consecutive_low_prob_count >= HEURISTIC_NUM :
                    unlikely_match = True
                    consecutive_low_prob_count = 0
                    continue

            probability += math.log(nucleotide_prob) if nucleotide_prob > 0 else math.log(epsilon)
            extended_query_start+=1

        if unlikely_match:
            unlikely_match = False
            continue

        if probability > max_prob:
            max_prob = probability
            max_match = start_end_tuple
        
        match_probability[start_end_tuple] = probability

def write_probabilities_to_csv(log_probabilities, output_file):
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

    print("The max match is ", max_match)

    found_string = ""
    for index in range(max_match[0], max_match[1] + 1):
        found_string += database_array[index]
    
    print("The found string is ", found_string)

# Sort the final matches from highest to lowest probability
sorted_match_probabilities = dict(sorted(match_probability.items(), key=lambda item: item[1], reverse=True))
write_probabilities_to_csv(sorted_match_probabilities, "final_result.csv")





