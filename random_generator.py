import random
from Bio import SeqIO
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

def probability_file_to_array(file_path):
    with open(file_path, 'r') as file:
        content = file.read()
        # Split the content by any whitespace character
        tokens = content.split()
        float_tokens = [float(token) for token in tokens]
        return float_tokens
# Get array for database sequence
database_array = fasta_to_array("database.fasta")
probability_array = probability_file_to_array("database_probs.txt")

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

positional_probabilities = get_positional_probabilities(database_array, probability_array)
max_index = 604465
min_length = 10
max_length = 1000

def generate_random_sequence(probabilities, max_position=604465):
    # Randomly choose a starting position between 0 and max_position
    start_position = random.randint(0, max_position)

    # Randomly choose a length between 10 and the maximum possible length
    length = random.randint(min_length, max_length)

    sequence = []
    end_position = start_position + length - 1

    for position in range(start_position, end_position + 1):
        if position not in probabilities:
            raise ValueError(f"Position {position} is out of range of the probability dictionary.")
        
        nucleotide_probs = probabilities[position]
        nucleotides = list(nucleotide_probs.keys())
        probs = list(nucleotide_probs.values())

        # Randomly choose a nucleotide based on the probabilities
        chosen_nucleotide = random.choices(nucleotides, weights=probs, k=1)[0]
        sequence.append(chosen_nucleotide)

    return ''.join(sequence), start_position, end_position, length

start_position = 289385
length = 100

generated_sequence, start_position, end_position, length = generate_random_sequence(positional_probabilities)

with open("query_sequence.fasta", "w") as fasta_file:
    fasta_file.write(">random_query_sequence\n")
    fasta_file.write(generated_sequence + "\n")

print("Generated sequence:", generated_sequence)
print("The start position is ", start_position)
print("The end position is ", end_position)
print("The query length is ", length)


