from Bio.Blast.Applications import NcbiblastnCommandline
import csv

# Define BLAST command with tabular output
blastn_cline = NcbiblastnCommandline(
    query="target_sequence.fasta",
    db="my_database",
    evalue=10,
    word_size=11,
    outfmt=6,  # Tabular format
    out="blast_results.tsv"
)

# Execute BLAST
stdout, stderr = blastn_cline()

# Define header for tabular format
headers = [
    "qseqid", "sseqid", "pident", "length",
    "mismatch", "gapopen", "qstart", "qend",
    "sstart", "send", "evalue", "bitscore"
]

# Read the tabular BLAST output and write to CSV
with open("blast_results.tsv") as tsv_file, open("blast_results.csv", mode='w', newline='') as csv_file:
    tsv_reader = csv.reader(tsv_file, delimiter='\t')
    csv_writer = csv.writer(csv_file)
    
    # Write headers
    csv_writer.writerow(headers)
    
    # Write all rows from TSV to CSV
    for row in tsv_reader:
        csv_writer.writerow(row)

print("BLAST results have been successfully converted to 'blast_results.csv'.")


